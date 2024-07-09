#ifndef VNTR_H_
#define VNTR_H_


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <vector>
#include <string>
#include <iostream>
#include <cassert>
#include <algorithm>
#include <map>
#include <stdint.h>
#include <inttypes.h>
#include "read.h"
#include "option.h"

using namespace std;




/**
 * @brief class for one vntr motif
 */
class MOTIF
{
public:
    string seq;     // the sequence of the motif
    int len;        // the lenght of the motif

    /// @brief Construct a new MOTIF object
    MOTIF() {};
    /**
     * @brief Construct a new MOTIF object
     * @param Seq input motif sequence
     */
    MOTIF(string &Seq) : seq(Seq) { len = Seq.length();};
    /// @brief Destroy the MOTIF object
    ~MOTIF() {};
};




/**
 * @brief class for the stringDecomposer dynamic programming tables (i.e., book)
 */
class SDTables
{
public:
    vector<vector<vector<int> > > pathMat;      // DP path book (outer: motif, middle: motif len, inner: vntr len)
    vector<vector<vector<int> > > nMatchMat;    // # of match book (outer: motif, middle: motif len, inner: vntr len)
    vector<vector<vector<int> > > nDelMat;      // # of del book (outer: motif, middle: motif len, inner: vntr len)
    vector<vector<vector<int> > > scoreMat;     // DP score book (outer: motif, middle: motif len, inner: vntr len)
    vector<int> pathRow0;                       // DP path book general first row (length: vntr len)
    vector<int> scoreRow0;                      // DP score book general first row (length: vntr len)
    vector<vector<int> > nMatchOnOptPath;       //
    vector<vector<int> > nDelOnOptPath;         //

    /**
     * @brief Initialize and size all DP records
     * @param motifs input vector of all motifs
     * @param seq input vntr sequence
     */
    void Init(vector<MOTIF> &motifs, string &seq)
    {
        try
        {
            // each motif has one table (outer layer indexed by motif)
            pathMat.resize(motifs.size());
            nMatchMat.resize(motifs.size());
            nDelMat.resize(motifs.size());
            scoreMat.resize(motifs.size());
            int totalSeqLen=0; // maintain a total motif length

            for (auto m=0; m < motifs.size(); m++)
            {
                // middle layer (col index) by vntr seq length (1 addition)
                pathMat[m].resize(seq.size()+1);
                nMatchMat[m].resize(seq.size()+1);
                nDelMat[m].resize(seq.size()+1);
                scoreMat[m].resize(seq.size()+1);

                for (auto s=0; s < seq.size() + 1; s++) 
                {
                    // inner layer (row index) indexced by motif length
                    totalSeqLen += motifs[m].len;
                    pathMat[m][s].resize(motifs[m].len);
                    nMatchMat[m][s].resize(motifs[m].len);
                    nDelMat[m][s].resize(motifs[m].len);
                    scoreMat[m][s].resize(motifs[m].len);

                    // fill all cells with 0
                    fill(pathMat[m][s].begin(), pathMat[m][s].end(), 0);
                    fill(nMatchMat[m][s].begin(), nMatchMat[m][s].end(), 0);
                    fill(nDelMat[m][s].begin(), nDelMat[m][s].end(), 0);
                    fill(scoreMat[m][s].begin(), scoreMat[m][s].end(), 0);
                }
            }

            // DP book general first row records
            pathRow0.resize(seq.size() + 1);
            scoreRow0.resize(seq.size() + 1);
        }
        catch (std::exception const & e)
        {
            cerr << "Exception allocating tables " << e.what() << endl;
            assert(0);
        }
    }
};




/**
 * @brief class for one vntr sequence
 * @note annotations are string of motif indeces
 */
class VNTR
{
public:

    string chr;                             // The VNTR reference chromosome
    uint32_t ref_start;                     // The VNTR reference start coordinate
    uint32_t ref_end;                       // The VNTR reference end coordinate
    int len;                                // The VNTR reference length
    string region;                          // The VNTR reference region in format: "chr:ref_start:ref_end"
    vector<MOTIF> motifs;                   // The pre-configured motifs for this VNTR
    int mappedContigLength;                 // The length of the contig that mapped this sequence.
                                            // When mapping from contigs, keep annotation from the longer contig, if multiple overlap.
    vector<READ *> reads;                   // all reads overlapping this VNTR
    int nreads;                             // number of reads
    int cur_len;                            // the sample sequence length (averaged over all reads)
    vector<uint8_t> haps;
    vector<vector<uint8_t> > annos;         // the motif annotations for each read
    vector<vector<uint8_t> > consensus;     // the final consensus annotation

    vector<int> h1_reads;                   // index of all haplotype1 reads (manipulated in consensusReadForHapByABpoa)
    vector<int> h2_reads;                   // index of all haplotype2 reads (manipulated in consensusReadForHapByABpoa)
    vector<READ *> Hap_seqs;
    bool het;                               // if the VNTR is heterozygous
    bool skip;                              // if skip the VNTR
    bool nullAnno;                          // if the annotation can't be generated
    vector<bool> nullAnnos;
    int len_h1;                             // 
    int len_h2;                             // 
    bool readsArePhased;                    // 
  int index;

    /// @brief Construct a new VNTR object
    VNTR()
    {
      index = 0;
        het = false;
        nreads = 0;
        skip = false;
        nullAnno = false;
        readsArePhased = false;
	mappedContigLength = 0;
    };


    /**
     * @brief Construct a new VNTR object
     * @param Chr VNTR reference chromosome
     * @param Start VNTR reference start coordinate
     * @param End VNTR reference end coordinate
     * @param Len VNTR reference length
     */
    VNTR(string Chr, uint32_t Start, uint32_t End, uint32_t Len)
    :chr(Chr), ref_start(Start), ref_end(End), len(Len) 
    {
        string s = ":" + to_string(ref_start);
        string e = "-" + to_string(ref_end);
        region = Chr + s + e;
        het = false;
        nreads = 0;
        skip = false;
        nullAnno = false;
        len_h1 = 0;
        len_h2 = 0;
	mappedContigLength = 0;
    };


    /// @brief Destroy the VNTR object
    ~VNTR() {};


    /**
     * @brief Clear reads and haplotype sequences to free up memory
     */
    void clearReads()
    {
        for (size_t i = 0; i < reads.size(); ++i) 
        { 
            delete reads[i];
        }
        reads.clear();

        if (!het and nreads == 1) return;

	Hap_seqs.clear();

        return;
    }


    /**
     * @brief Annotate each sequence by motifs
     * @param opt general options
     * @param sdTables stringDecomposer tables
     * @param mismatchCI mismatch confidence intervals
     */
    void motifAnnoForOneVNTR(const OPTION &opt, SDTables &sdTables, vector<int> &mismatchCI);


    /**
     * @brief Get the reads consensus (nt sequence) of one/two haplotypes by abpoa msa
     * @param opt general options
     */
    void consensusReadForHapByABpoa (const OPTION &opt);


    /**
     * @brief Get the consensus of all annotated strings by abpoa msa
     * @param opt general options
     */
    void consensusMotifAnnoForOneVNTRByABpoa (const OPTION &opt);


    /**
     * @brief Get the comma seperated motif annotation string for the vcf output
     * @param h1 if this is haplotype 1 consensus
     * @param motif_rep the output comma seperated motif annotation string
     */
    void commaSeparatedMotifAnnoForConsensus (bool h1, string &motif_rep);
};




/**
 * @brief debug function to print the consensus
 * @param consensus 
 */
void outputConsensus (vector<uint8_t> &consensus);


class CompareVNTRPos {
public:
  vector<VNTR*>* vntrs;
  bool operator()(const int &a, const int &b) const {
    return (*vntrs)[a]->ref_start < (*vntrs)[b]->ref_start;
  }
};

class UpperBoundSearchVNTRPos {

public:
  vector<VNTR*>* vntrs;
  bool operator()(const int &a, const int &b) const {
    return a < (*vntrs)[b]->ref_start;
  }
};

class LowerBoundSearchVNTRPos {
public:
  vector<VNTR*>* vntrs;
  bool operator()(const int &a, const int &b) const {
    return  (*vntrs)[a]->ref_start < b;
  }
};

void SortVNTRIndexByPos(vector<VNTR*> &vntrs, vector<int> &idx);
void ChromToVNTRMap(vector<VNTR*> &vntrs, map<string, vector<int> > &vntrMap);

#endif
