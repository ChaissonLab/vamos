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

using namespace std;

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
  vector<int> counts;
  int operator[](int i) { return counts[i];}
  static int minCov;
  static char*nucIndex;
  char a,b;
  bool isCluster;
  int pos;
  int cov;
  Consensus() {
    counts.resize(6,0);
    a=b='N';
    isCluster=false;
    cov=0;
  }
  int Inc(char c) {
    cov++;
    counts[nucIndex[c]]++;
    return counts[nucIndex[c]];
  }
  void Dec() {
    cov--;
  }
  int Calc() {
    int colSum=0;    
    for (auto i=0; i < 5; i++) {
	colSum += counts[i];
    }
    vector<float> frac(4,0);
    vector<int> suffIdx;
    vector<char> snvNuc;
    if (colSum < minCov) {
      a = 'N';
      b = 'N';
      return 0;
    }
    
    for (auto i=0; i < 4; i++) {
	frac[i] = ((float)counts[i]/colSum);
	if (frac[i] > 0.25) {
	    suffIdx.push_back(i);
	    snvNuc.push_back("ACGT"[i]);
	}
    }
    if (suffIdx.size() == 2) {
      a=snvNuc[0];
      b=snvNuc[1];
      return 1;
    }
    return 0;
  }
};

void InitNucIndex(char nucIndex[]);
class HetSNV {
public:
  int pos;
  char a, b;
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
    Consensus::nucIndex = new char[256];
    InitNucIndex(Consensus::nucIndex);
    Consensus::minCov = 10;
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
      consensus[alignedPos].Inc(read.readToRefSeq[i]);
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
  
  void ProcessUntil(int pos) {
    while (consensusStart < pos) {
      Consensus& cons=consensus[consensusStart];
      if (consensus[consensusStart].Calc()) {
	HetSNV snv;
	snv.pos = consensusStart;
	snv.a = consensus[consensusStart].a;
	snv.b = consensus[consensusStart].b;	
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
  IOAlignBuff(samFile *initFp_in, hts_itr_t *initIter, mutex * initIoLock, int initBuffSize) {
    itr = initIter;
    fp_in = initFp_in;
    ioLock = initIoLock;
    buffSize=initBuffSize;
    nRead=0;
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
      bam1_t * aln = bam_init1(); //initialize an alignment
      if (alignments.size() == 0 ){
	alignments.resize(buffSize);
      }
      //      cerr << "Refreshing IO Buffer " << endl;      
      while (nRead < buffSize) {
	if (bam_itr_next(fp_in, itr, aln) >= 0) {
	  alignments[nRead] = aln;
	  aln = bam_init1(); //initialize an alignment		  
	  ++nRead;
	}
	else {
	  break;
	}
      }
      if (nRead > 0) {
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
    char * region_and_motifs;
    char * input_fasta;
    char * input_bam;
    char * vntr_bed;
    char * motif_csv;
    char * out_vcf;
    char * sampleName;
    char * version;
    OutWriter outWriter;
    char *bai;
    samFile * fp_in;
    bam_hdr_t * bamHdr;
    hts_idx_t * idx;
    mutex *ioLock;
    int *numProcessed;
    int thread;
    int curChromosome;
    vector<string> chromosomeNames;
    std::map<string, vector<int> > *vntrMap;

    // Not the best place to put this, but since the IO is batched and we
    // do not want to store the flanking sequences at each locus for
    // the duration of the run of the program, it's passed through
    // to here for now. Eventually if sites are processed immediately
    // after overlapping reads are read, this can move.
    int phaseFlank;
    IO() 
    {
        version = (char *) malloc(7);
        strcpy(version, "1.2.8");
        region_and_motifs = NULL;
        input_bam = NULL;
        vntr_bed = NULL;
        motif_csv = NULL;
        out_vcf = NULL;
        sampleName = NULL;
        phaseFlank = 0;
        ioLock = NULL;
	thread = 0;
	curChromosome=0;
    };


    ~IO() 
    {
        free(region_and_motifs);
        free(version);
        free(input_bam);
        free(vntr_bed);
        free(motif_csv);
        free(out_vcf);
        free(sampleName);
    };


    int readRegionAndMotifs (vector<VNTR*> &vntrs);


    int readMotifsFromCsv (vector<VNTR *> &vntrs);


    int read_tsv(vector<vector<string>> &items);


    void readVNTRFromBed (vector<VNTR *> &vntrs);


    /* get the sequences from input_bam_file that overlapping with chr:start-end */
    //void readSeqFromBam (vector<VNTR *> &vntrs, int nproc, int cur_thread, int sz);
    void initializeBam();


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

  void StoreVNTRLoci(vector<VNTR*> &vntrs, vector<int> &vntrIndex, Pileup &pileup, int &refAlignPos, int &curVNTR);

  void ProcessReadsOnChrom(string &chrom, vector<VNTR*> &vntrs, map<string, vector<int> > &vntrMap);
  void CallSNVs(string &chrom, vector<VNTR*> &vntrs, map<string, vector<int> > &vntrMap, Pileup &pileup);  

  void StoreReadSeqAtRefCoord(bam1_t *aln, string &seq, string &toRef, vector<int> &map);
  void ProcessOneContig(bam1_t *aln, vector<VNTR*> &vntrs, map<string, vector<int> > &vntrIndex);
  void StoreAllContigs(vector<VNTR*> &vntrs, map<string, vector<int> > &vntrIndex);
  void StoreSeq(bam1_t *aln, string &seq);
};

#endif
