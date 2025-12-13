#ifndef IO_H_
#define IO_H_
#include <iostream>
#include <stdlib.h>
#include <istream>
#include <fstream>
#include <sstream>
#include <assert.h>
#include <string>
#include <tuple> 
#include <vector>
#include "read.h"
#include "vntr.h"
#include "vcf.h"
#include "option.h"
#include <mutex>
#include <map>
#include <list>
#include <queue>
#include <htslib/faidx.h>

using namespace std;

void GetItOfOverlappingVNTRs(vector<VNTR*> &vntrs,
		   map<string, vector<int> > &vntrMap, string &chrom, int regionStart, int regionEnd,
		   vector<int>::iterator &startIt,
		   vector<int>::iterator &endIt);

class PiledRead {
public:
  int refStartPos;
  int refEndPos;
  char flag;
  string readToRefSeq, readSeq;
  string readName;
  int phased;
  int hap;
  PiledRead(int initRefPos, string &initReadSeq ) {
    refStartPos = initRefPos;
    readToRefSeq = initReadSeq;
    refEndPos = refStartPos + initReadSeq.size();
  }
  
  int GetCharAt(int pos, char &res) {
    if (pos >= refStartPos and pos < refEndPos) {
      int offset = pos - refStartPos;
      res = readToRefSeq[offset];
      return 1;
    }
    else {
      res=255;
      return 0;
    }
  }
  int GetSubstr(int start, int end, string &res) {
    if (start >= refStartPos and end <= refEndPos) {
      start-=refStartPos;
      end-=refStartPos;
      res=readToRefSeq.substr(start, end-start);
      return 1;
    }
    else {
      return 0;
    }
  }
};

class Consensus {
public:
  int counts[6];
  int operator[](int i) { return counts[i];}
  char a,b;
  bool isCluster;
  int pos;
  int cov;
  Consensus() {
    counts[0] = counts[1] =counts[2]=counts[3]=counts[4]=counts[5]=0;
    a=b='N';
    isCluster=false;
    cov=0;
  }
  int Inc(char c, char nucIndex[]) {
    cov++;

    assert(c <= 255 && c >= 0);
    assert(nucIndex[c] < 6);
    counts[nucIndex[c]]++;
    return counts[nucIndex[c]];
  }
  void Dec() {
    cov--;
  }
  int Calc(int minSNVCov, int minAltCov ) {
    int colSum=0;
    int snvSum=0;
    for (auto i=0; i < 5; i++) {
	colSum += counts[i];
    }
    snvSum=colSum - counts[4];
    vector<float> frac(4,0);
    vector<int> snvIdx;
    vector<char> snvNuc;
    if (snvSum < minSNVCov) {
      a = 'N';
      b = 'N';
      return 0;
    }
    
    for (auto i=0; i < 4; i++) {
	frac[i] = ((float)counts[i]/colSum);
	if (frac[i] > 0.25) {
	    snvIdx.push_back(i);
	    snvNuc.push_back("ACGT"[i]);
	}
    }
    
    if (snvIdx.size() == 2) {
      if (counts[snvIdx[0]] < minAltCov or counts[snvIdx[1]]< minAltCov) {
	a='N';
	b='N';
	return 0;
      }
      else {
	a=snvNuc[0];
	b=snvNuc[1];
	return 1;
      }
    }
    return 0;
  }
};

void InitNucIndex(char nucIndex[]);
class HetSNV {
public:
  int pos;
  char a, b;
  HetSNV() {}
  HetSNV(int initPos) { pos=initPos;}
};


class CompPosWithSNV {
public:
  bool operator()(const HetSNV &lhs, const HetSNV &rhs ) {
    return lhs.pos < rhs.pos;
  }
};

class EndPosPair {
public:
  int pos;
  list<PiledRead>::iterator it;
  EndPosPair(int initPos, list<PiledRead>::iterator &initIt) {
    pos=initPos;
    it = initIt;
  }
};

class SortEndPosPair {
public:
  bool operator()(const EndPosPair &lhs, const EndPosPair &rhs) const {
    return lhs.pos > rhs.pos;
  }
};
    
class Pileup {
public:
  list<PiledRead> reads;
  map<int, Consensus > consensus;
  int consensusStart;
  int consensusEnd;
  char nucIndex[256];
  vector<HetSNV> hetSNVs;
  priority_queue< EndPosPair, vector<EndPosPair>, SortEndPosPair>  endPosPQueue;
  Pileup() {
    InitNucIndex(nucIndex);
    consensusStart=-1;
  }
  
  void AddRead(PiledRead &read ) {

    reads.push_back(read);
    //
    // No consensus window defined.
    //
    if (consensusStart == -1) {
      consensusStart = read.refStartPos;
      consensusEnd = consensusStart;
    }
    //
    // By now, a consensus window has been defined. It must have been
    // initialized, and should have a position before the current read.
    //
    assert(consensusStart <= read.refStartPos );
    //
    // Expand the consensus to fit the new read.
    if (consensusEnd < read.refEndPos ) {
      //      cerr << "Expanding consensus from " << consensusStart << "-" << consensusEnd << " to " << read.refEndPos << endl;
      while (consensusEnd < read.refEndPos) {
	consensus[consensusEnd] = Consensus();
	consensus[consensusEnd].pos = consensusEnd;
	consensusEnd++;
      }
    }
    //    cerr << "Adding read to range " << read.refStartPos << " - " << read.refEndPos << endl;
    for (int i=0; i < read.readToRefSeq.size(); i++) {
      int alignedPos = i + read.refStartPos;
      //      assert(consensus.find(alignedPos) != consensus.end());
      consensus[alignedPos].Inc(read.readToRefSeq[i], nucIndex);
    }
    auto it=reads.end();
    it--;
    endPosPQueue.push(EndPosPair(read.refEndPos, it));
    /*
    for (auto it : consensus) {
      if (it.first % 100 == 0 && it.first < 44838100 ) {
	cerr << it.first << "\t" << it.second.cov << endl;
      }
    }
    */
      
  }
  void RemoveRead(list<PiledRead>::iterator it) {
    for (int i=it->refStartPos; i < it->refEndPos; i++ ){
      //      cerr << "dec before " << consensus[i].cov << " " << consensus[i].pos << endl;
      consensus[i].Dec();
      //      cerr << "dec after " << consensus[i].cov << " " << consensus[i].pos << endl;
    }
    reads.erase(it);
  }
  
  void ProcessUntil(int pos, OPTION &opts) {
    while (consensusStart < pos) {
      Consensus& cons=consensus[consensusStart];
      if (consensus[consensusStart].Calc(opts.minSNVCoverage, opts.minAltCoverage)) {
	HetSNV snv;
	snv.pos = consensusStart;
	snv.a = consensus[consensusStart].a;
	snv.b = consensus[consensusStart].b;
	/*
	cerr << "Found het snv: " << snv.pos << "\t" << snv.a << "\t" << snv.b
	     << "\t" << consensus[consensusStart].counts[0]
	     << "\t" << consensus[consensusStart].counts[1]
	     << "\t" << consensus[consensusStart].counts[2]
	     << "\t" << consensus[consensusStart].counts[3]
	     << "\t" << consensus[consensusStart].counts[4]
	     << "\t" << consensus[consensusStart].counts[5] 	   << endl;
	*/
	hetSNVs.push_back(snv);
      }
      consensus.erase(consensusStart);
      consensusStart++;
    }
  }

  void PurgeBefore(int pos) {    
    auto readIt = reads.begin();
    while (readIt != reads.end()) {
      if (readIt->refEndPos <= pos ) {
	readIt = reads.erase(readIt);
      }
      else {
	readIt++;
      }
    }
    if (reads.size() == 0) {
      consensus.clear();
    }
    else {
      cerr << "Removing consensus from " << consensusStart << " to " << reads.front().refStartPos << endl;
      while (consensusStart < reads.front().refStartPos ) {
	assert(consensus.size() > 0);
	consensus.erase(consensusStart);
	consensusStart++;
      }
    }
  }
};

class IOAlignBuff {
public:
  int buffSize;
  vector<bam1_t*> alignments;
  int curAln;
  hts_itr_t *itr;
  samFile * fp_in;
  mutex *ioLock;
  int nRead;
  void Free() {
    for (; curAln < alignments.size(); curAln++) {
      bam_destroy1(alignments[curAln]);
    }
  }
  IOAlignBuff(samFile *initFp_in, hts_itr_t *initIter, mutex * initIoLock, int initBuffSize) {
    itr = initIter;
    fp_in = initFp_in;
    ioLock = initIoLock;
    buffSize=initBuffSize;
    nRead=0;
    curAln=0;
  }
  int GetNext(bam1_t* &next) {
    if (curAln >= 0 and curAln < nRead ) {
      next=alignments[curAln];
      ++curAln;
      return curAln;
    }
    else {
      curAln = 0;
      nRead=0;
      if (ioLock != NULL) { ioLock->lock();}
      bam1_t * aln;// = bam_init1(); //initialize an alignment
      if (alignments.size() == 0 ){
	alignments.resize(buffSize);
      }
      //      cerr << "Refreshing IO Buffer " << endl;      
      while (nRead < buffSize) {
	aln = bam_init1(); //initialize an alignment		  
	if (bam_itr_next(fp_in, itr, aln) >= 0) {
	  alignments[nRead] = aln;
	  ++nRead;
	}
	else {
	  bam_destroy1(aln);
	  break;
	}
      }
      // Last alignment not needed.

      if (nRead > 0) {
	//
	// Allocated one extra
	next=alignments[curAln];
	curAln++;
	if (ioLock != NULL) { ioLock->unlock();}
	return curAln;
      }
      if (ioLock != NULL) { ioLock->unlock();}
      return nRead;
    }
  }
};

class IO
{
public:
    string region_and_motifs;
    string input_fasta;
    string input_bam;
    string reference;
    string vntr_bed;
    string motif_csv;
    string out_vcf;
    string sampleName;
    string version;
    OutWriter outWriter;
    string bai;
    samFile * fp_in;
    bam_hdr_t * bamHdr;
    hts_idx_t * idx;
    mutex *ioLock;
    int *numProcessed;
    int thread;
    int minMapQV;
  int minChrY;
    uint64_t nChrY;
    string curChromosome;
    vector<string> chromosomeNames;
    vector<int> chromosomeLengths;
    faidx_t *fai;
    int maxLength;
    std::map<string, vector<int> > *vntrMap;

    // Not the best place to put this, but since the IO is batched and we
    // do not want to store the flanking sequences at each locus for
    // the duration of the run of the program, it's passed through
    // to here for now. Eventually if sites are processed immediately
    // after overlapping reads are read, this can move.
    int phaseFlank;
    IO() 
    {
      version = "2.1.8";
        region_and_motifs = "";
        input_bam = "";
        reference = "";
        vntr_bed = "";
        motif_csv = "";
        out_vcf = "";
        sampleName = "";
        phaseFlank = 0;
        ioLock = NULL;
	thread = 0;
	curChromosome="";
	minMapQV=3;
	maxLength=10000;
	minChrY = nChrY = 0;
    };

  void clear() {
    cerr << "Clearing IO" << endl;
      if (bamHdr != NULL) {
	//
	// The structure was read using sam_hdr_read, so destroy it with this.
	//
	bam_hdr_destroy(bamHdr);
      }
      if (idx != NULL) {
	hts_idx_destroy(idx);
      }
  }    
    ~IO() 
    {
    };


    int readRegionAndMotifs (vector<VNTR*> &vntrs);


    int readMotifsFromCsv (vector<VNTR *> &vntrs);


    int read_tsv(vector<vector<string>> &items);


    /* get the sequences from input_bam_file that overlapping with chr:start-end */
    //void readSeqFromBam (vector<VNTR *> &vntrs, int nproc, int cur_thread, int sz);
    void initializeBam();
  void initializeRefFasta();

    void closeBam();


  //    void readSeqFromBam(vector<VNTR *>&vntrs, int pos, OPTION &opts);

    // void readSeqFromBam (vector<READ*> &reads, string &chr, const uint32_t &ref_VNTR_start, 
    //    const uint32_t &ref_VNTR_end, const uint32_t &VNTR_len, string &region);
    
    // int outputVCF (vector<VNTR *> &vntrs);


    int writeVCFHeader_locuswise(ofstream& out);


    int writeVCFBody_locuswise(ofstream& out, vector<VNTR *> &vntrs, int tid, int nproc);


    int writeBEDHeader_readwise(ofstream& out);


    int writeBEDBody_readwise(ofstream& out, vector<VNTR *> &vntrs, int tid, int nproc);


    void writeFa(ofstream& out, vector<VNTR *> &vntrs);


    void readSeqFromFasta(vector<VNTR *> &vntrs);

  //  void StoreVNTRLoci(vector<VNTR*> &vntrs, vector<int> &vntrIndex, Pileup &pileup, int &refAlignPos, int &curVNTR);

  int CallSNVs(string &chrom, int regionStart, int regionEnd, vector<VNTR*> &vntrs, map<string, vector<int> > &vntrMap, Pileup &pileup, bool& readsArePhased, OPTION &opts);  
  void StoreReadsOnChrom(string &chrom, int regionStart, int regionEnd, vector<VNTR*> &vntrs, map<string, vector<int> > &vntrMap, Pileup &pileup, int thread, bool readsArePhased);

  void StoreReadSeqAtRefCoord(bam1_t *aln, string &seq, string &toRef, vector<int> &map);
  void ProcessOneContig(bam1_t *aln, vector<VNTR*> &vntrs, map<string, vector<int> > &vntrIndex);
  void StoreAllContigs(vector<VNTR*> &vntrs, map<string, vector<int> > &vntrIndex);
  void StoreSeq(bam1_t *aln, string &seq);
  string MakeRegion(string chrom, int start, int end) {
    stringstream sstrm;
    sstrm << chrom << ":" << start << "-" << end;
    return sstrm.str();
  }
};


#endif
