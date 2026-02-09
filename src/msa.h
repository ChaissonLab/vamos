#ifndef MSA_H_
#define MSA_H_

#include "read.h"
#include "abpoa.h"
#include "option.h"
#include "spoa/spoa.hpp"

class MSA
{
public:
    vector<string> seqs;
    string consensus;

    MSA() 
    {
    };

  ~MSA() {
    cleanMSA();
  }
    /**
     * @brief Construct a new MSA object for msa of anno sequences - initialize all anno seqs as vector of motif indices
     * @param annos annotations (by motif index vector) of all reads for this VNTR
     * @param motifs_size number of motifs for this VNTR
     */

    /**
     * @brief Construct a new MSA object for msa of raw sequences - initialize all raw nt seqs as 0-4 nt encoding
     * @param gp vector of indices of all reads in a haplotype group
     * @param reads vector of all reads
     */
    MSA(const vector<int> &gp, const vector<READ *> &reads)
    {
	for (int i = 0; i < gp.size(); i++) {
	  if (reads[gp[i]]->seq.size() > 0) {
	    seqs.push_back(reads[gp[i]]->seq);
	  }
	}
    };


    /**
     * @brief Run the abpoa msa
     */
    void runConsensus();


    /**
     * @brief Convert the consensus (by motif indices + 1) back to the original motif indices string
     * @param consensus converted final MSA consensus (by motif indices)
     * @param motifs_size number of motifs
     */
    void extractConsensusAnno(vector<uint8_t> &consensus, int motifs_size);


    /**
     * @brief Convert the consensus (by 0-4 nt encoding) back to the original nt sequence
     * @param Hap_seq converted final MSA consensus (by nt AaCcGgTtNn)
     */
    void extractConsensusRawSeq(READ * Hap_seq);


    /**
     * @brief Generate consensus sequence for a group of motif anno sequences through abpoa
     * @param motifs_size number of motifs
     * @param annos annotations (by motif index vector) of all reads for this VNTR
     * @param consensus converted final MSA consensus (by motif indices)
     */
    void MSA_anno_group(int motifs_size, const vector<vector<uint8_t>> &annos, vector<uint8_t> &consensus);


    /**
     * @brief Generate consensus sequence for a group of read sequences through abpoa
     * @param gp indeces of all input reads for MSA
     * @param reads all input reads
     * @param Hap_seq converted final MSA consensus (by nt AaCcGgTtNn)
     */
    void MSA_seq_group(const vector<int> &gp, const vector<READ *> &reads, READ * Hap_seq);


    /**
     * @brief Clear memory
     */
    void cleanMSA();
};


#endif
