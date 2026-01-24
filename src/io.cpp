#include <iostream>
#include <stdlib.h>
#include <istream>
#include <fstream>
#include <sstream>
#include <assert.h>
#include <string>
#include <cctype>
#include <algorithm>
#include <tuple> 
#include <vector>
#include <regex>
#include "io.h"
#include "read.h"
#include "vntr.h"
#include "vcf.h"
#include "phase.h"
#include "abpoa.h"
// #include <seqan/align.h>
// #include <seqan/graph_msa.h>
#include "htslib/hts.h"
#include "htslib/sam.h"
#include <zlib.h>  
#include <queue>
#include "htslib/kseq.h"  
#include <parasail.h>
#include <parasail/matrix_lookup.h>
#include <parasail/io.h>
extern int naive_flag;
extern int debug_flag;
extern int hclust_flag;
extern int seqan_flag;
extern int liftover_flag;
extern int conseq_anno_flag;
extern int locuswise_prephase_flag;
extern int locuswise_flag;

//0123456789ABCDEF
//=ACMGRSVTWYHKDBN  aka seq_nt16_str[]
//=TGKCYSBAWRDMHVN  comp1ement of seq_nt16_str
//084C2A6E195D3B7F
// static int seqi_rc[] = { 0,8,4,12,2,10,6,14,1,9,5,13,3,11,7,15 };

const char rcseq_nt16_str[] = "!TGKCYSBAWRDMHVN";

// STEP 1: declare the type of file handler and the read() function  
KSEQ_INIT(gzFile, gzread)



int IO::readMotifsFromCsv (vector<VNTR *> &vntrs) 
{
    ifstream ifs(motif_csv);
    if (ifs.fail()) 
    {
        cerr << "Unable to open file " << motif_csv << endl;
        return 1;
    }    

    string line;
    size_t numOfLine = 0;
    while (getline(ifs, line)) 
    {
        stringstream ss(line);
        string tmp;
        while(getline(ss, tmp, ',')) 
        {
            vntrs[numOfLine]->motifs.push_back(MOTIF(tmp));
            // cerr << vntrs[numOfLine]->motifs.back().seq << endl;
        }
        numOfLine += 1; // 0-indexed
    }
    assert(vntrs.size() == numOfLine);
    return 0;
}

int IO::read_tsv(vector<vector<string>> &items) 
{
    items.clear();
    ifstream ifs(vntr_bed);
    if (ifs.fail()) 
    {
        cerr << "Unable to open file " << vntr_bed << endl;
        return 1;
    }

    string line;
    while (getline(ifs, line)) 
    {
        stringstream ss(line);
        vector<string> item;
        string tmp;
        while(getline(ss, tmp, '\t')) 
            item.push_back(tmp);

        // for (auto &i : item)
        //     cerr << i << "\t";
        // cerr << endl;
        items.push_back(item);
    }
    return 0;
}



int IO::readRegionAndMotifs (vector<VNTR*> &vntrs) 
{
    ifstream ifs(region_and_motifs);
    if (ifs.fail()) 
    {
        cerr << "Unable to open file " << motif_csv << endl;
        return 1;
    }    

    string line;
    size_t numOfLine = 0;
    string prevChrom="";
    int lineNumber=1;
    int prevEnd =-1;
    while (getline(ifs, line)) 
    {
        stringstream ss(line);
        string motifs;
        string tmp;
        // string tmp_r;
        string chrom;
        int start;
        int end;
	string version, svtype;
        ss >> chrom >> start >> end >> motifs >> version >> svtype;
	if (end - start > maxLength) {
	  cerr << "WARNING, locus " << chrom << ":" << start << "-" << end << " has length greater than " << maxLength << " and will be ignored." << endl;
	  continue;
	}
	if (start > end) {
	  cout << "ERROR on line " << lineNumber << " of motif file: start > end: " << start << "\t" << end << endl;
	  exit(1);
	}
	if (chrom == prevChrom and prevEnd > start) {
	  cout << "ERROR on line " << lineNumber << " of motif file: start is before previous end: " << start << "\t" << prevEnd << endl;
	  exit(1);
	}
	prevChrom=chrom;
	prevEnd = end;
	lineNumber++;
        VNTR * vntr = new VNTR(chrom, start, end, end-start,svtype);
        vntrs.push_back(vntr);
	vntr->index = vntrs.size()-1;
        stringstream mm(motifs);
        while(getline(mm, tmp, ',')) 
        {
            // tmp_r = regex_replace(tmp, std::regex("^\\t+"), std::string(""));
            vntrs[numOfLine]->motifs.push_back(MOTIF(tmp));
            // cerr << vntrs[numOfLine]->motifs.back().seq << endl;
        }
        numOfLine += 1; // 0-indexed
    }
    assert(vntrs.size() == numOfLine);
    return 0;
}



int FRONT=0;
int BACK=1;

void IO::StoreReadSeqAtRefCoord(bam1_t *aln, string &readSeq, string &readSeqAtRefCoord, vector<int> &refToReadMap) {
    uint32_t readPos=0;
    uint32_t refEnd=bam_endpos(aln);
    int refPos=aln->core.pos;
    uint32_t refStart=refPos;    
    uint32_t *cigar = bam_get_cigar(aln);
    //
    // Check for unmapped sequence.
    //
    uint32_t readLen = aln->core.l_qseq;    
    if (refPos < 0 or readLen == 0) {
      return;
    }
    refToReadMap.resize(refEnd-refStart);
    
    int op, type, len;
    uint8_t *read;

    // If there is a prefix gap, add those.
    read=bam_get_seq(aln);
    char *name = bam_get_qname(aln);
    readSeqAtRefCoord.resize(refEnd - refPos);
    int refOffset=0;
    uint32_t k;
    for (k = 0; k < aln->core.n_cigar && refPos < refEnd && readPos < readLen; k++) 
    {
        assert(readPos < readLen);
        op = bam_cigar_op(cigar[k]);
        type = bam_cigar_type(op);
        len = bam_cigar_oplen(cigar[k]);

        if (!(type & 1) and !(type & 2)) // Hard clip
            continue;
        // A match or mismatch type
        if ((type & 1) && (type & 2)) 
        {
            for (int i=0; i < len; i++ ) 
            {
	      readSeqAtRefCoord[refOffset] = readSeq[readPos];
	      refToReadMap[refPos-refStart] = readPos;
	      refOffset++;
	      readPos++;
	      refPos++;
            }
        }
        else if ( (type & 1) == 0 && (type & 2) != 0 ) {
            for (int i=0; i < len; i++) 
            {
	      readSeqAtRefCoord[refOffset] = '-';
	      refToReadMap[refPos-refStart] = readPos;
	      refPos +=1;
	      refOffset++;
            }
        }
        else if ((type & 1) != 0 && (type & 2) == 0) {
            readPos+=len;
        }
    }
    assert(refOffset == refEnd-refStart);
}


int MappedStartPosInRead(vector<int> &refIndex, int refStart, int refPos) {
  //
  // For handling flank regions, the ref pos may have been set beyond the end of the read.
  //
  if (refPos < refStart) {
    refPos = refStart;
  }
  
  assert(refPos - refStart < refIndex.size());
  int i=refPos-refStart;
  while (i+ 1 < refIndex.size() and refIndex[i] == refIndex[i+1]) {
    //    cerr << "skipping a gap character at " << i << endl;
    i++;
  }
  //
  // Check to see if there was an insertion right before the matched pos. If so, return the
  // start of the insertion.
  // The case for an insertion is the two adjacent mapped positions have a larger gap than a single base.
  if (i > 0 and refIndex[i-1] + 1  < refIndex[i]) {
    return(refIndex[i-1]+1);
  }
  else {
    return(refIndex[i]);
  }
}

int MappedEndPosInRead(vector<int> &refIndex, int refStart, int refPos) {
  //
  // For handling flank regions, the ref pos may have been set beyond the end of the read.
  //
  if (refPos >= refStart + refIndex.size()) {
    refPos = refStart + refIndex.size() - 1;
  }
  assert(refPos - refStart < refIndex.size());
  int i=refPos-refStart;
  while (i - 1 > 0 and refIndex[i-1] == refIndex[i]) {
    //    cerr << "skipping a gap character at end " << i << endl;    
    i--;
  }
  //
  // Check for an inserted sequence at the end of the match
  if (i + 1 < refIndex.size() and refIndex[i] + 1 < refIndex[i+1]) {
    return refIndex[i+1]-1;
  }
  else {
    return(refIndex[i]);
  }
}

void IO::StoreAllContigs(vector<VNTR*> &vntrs, map<string, vector<int> > &vntrMap) {
  initializeBam();
  initializeRefFasta();
  bam1_t * aln = bam_init1(); //initialize an alignment
  while (sam_read1(fp_in, bamHdr, aln) >= 0) {
    ProcessOneContig(aln, vntrs, vntrMap);
    bam_destroy1(aln);
    aln=bam_init1();    
  }
  bam_destroy1(aln);
}


void IO::ProcessOneContig(bam1_t *aln, vector<VNTR*> &vntrs, map<string, vector<int> > &vntrMap) {

  vector<int> readMap;
  string readSeq, readToRef;
  StoreSeq(aln, readSeq);
  StoreReadSeqAtRefCoord(aln, readSeq, readToRef, readMap);
  if (readSeq.size() == 0) {
    return;
  }
  int tid=aln->core.tid;
  if (tid == -1) {
    return;
  }
  string refName(bamHdr->target_name[tid]);
  //
  // Check if any vntrs found in this chromosome (possibly a problem for an assembly)
  //
  if (vntrMap.find(refName) == vntrMap.end()) {
    return;
  }
  LowerBoundSearchVNTRPos lbComp;
  lbComp.vntrs=&vntrs;
  int refAlnStart = aln->core.pos;
  int refAlnEnd = bam_endpos(aln);
  vector<int> &vntrIndex=vntrMap[refName];
  vector<int>::iterator it = std::lower_bound(vntrIndex.begin(), vntrIndex.end(), refAlnStart, lbComp);

  //
  // Second check if any mapped by this contig.
  //
  if (it == vntrMap[refName].end()) {
    return;
  }
  int i=it-vntrMap[refName].begin();
  //
  while (i > 0 and vntrs[vntrMap[refName][i-1]]->ref_start == refAlnStart) {
    i--;
  }
  
  while (i < vntrMap[refName].size() and
	 vntrs[vntrMap[refName][i]]->ref_start >= refAlnStart and
	 vntrs[vntrMap[refName][i]]->ref_end < refAlnEnd) {

    int idx= vntrMap[refName][i];
    if (vntrs[idx]->mappedContigLength < refAlnEnd - refAlnStart) {
      int readStart = MappedStartPosInRead(readMap, refAlnStart, vntrs[idx]->ref_start-1);
      int readEnd   = MappedEndPosInRead(readMap, refAlnStart, vntrs[idx]->ref_end-1);
      string vntrSeq;
      /*
      int preadStart   = MappedStartPosInRead(readMap, refAlnStart, vntrs[idx]->ref_start-12);
      int preadEnd = MappedEndPosInRead(readMap, refAlnStart, vntrs[idx]->ref_start-2);
      int sreadStart = MappedEndPosInRead(readMap, refAlnStart, vntrs[idx]->ref_end);
      int sreadEnd   = MappedStartPosInRead(readMap, refAlnStart, vntrs[idx]->ref_end+10);
      */
      if (readEnd == readStart) {
        vntrSeq = readSeq.substr(readStart, readEnd-readStart);
      }
      else {
        vntrSeq = readSeq.substr(readStart, readEnd-readStart+1);
      }
      /*
      string pSeq, sSeq;
      pSeq=readSeq.substr(preadStart, preadEnd-preadStart);
      sSeq=readSeq.substr(sreadStart, sreadEnd - sreadStart);
      cout << refName << "\t" << vntrs[idx]->ref_start-1 << "\t" << vntrs[idx]->ref_end-1 << "\t" << pSeq << "\t" << vntrSeq << "\t" << sSeq << endl;
      */
      
      READ *read = new READ;
      read->qname = refName; //new char[refName.size()+1];
      //      memcpy(read->qname, refName.c_str(), refName.size());
      //      read->qname[refName.size()]='\0';
      read->l_qname = read->qname.size();
      read->seq = vntrSeq;
      read->len = vntrSeq.size();
      read->flag = aln->core.flag;
      if (vntrs[idx]->reads.size() == 0) {      
	vntrs[idx]->reads.push_back(read);
	READ *hsRead = new READ;
	*hsRead = *read;
	vntrs[idx]->Hap_seqs.push_back(hsRead);	
      }
      else {
	delete vntrs[idx]->reads[0];	
	vntrs[idx]->reads[0] = read;
	READ *hsRead = new READ;
	*hsRead = *read;
	vntrs[idx]->Hap_seqs[0] = hsRead;	
      }
      assert(vntrs[idx]->reads.size() == 1);
    }
    i++;
  }
}
int GetHap(char *s, int l) {
    int si=0, nc=0;
    for (; si < l; si++) 
    {
        if (s[si] == ':') { nc+=1;}
        if (nc == 2) { si++; break; }
    }
    if (nc == 2) 
    {
        return atoi(&s[si]);
    }
    return -1;
}

int QueryAndSetPhase(bam1_t* aln, int &setHap) {
  //    kstring_t auxStr = { 0, 0, NULL };
  char *auxStr;
    uint8_t* auxRes= bam_aux_get(aln, "HP");
    bool foundAux=false;
    if (auxRes != NULL) {
      int hap = bam_aux2i(auxRes);
      setHap=-1;    
      bool isPhased = false;
      int si=0, nc=0;

      setHap=hap;
      foundAux = true;	

    }
    //    free(auxRes);
    return foundAux;
}


void IO::StoreSeq(bam1_t *aln, string &readSeq) {
  uint8_t *read = bam_get_seq(aln);  
  uint32_t readLen = aln->core.l_qseq;
  readSeq.resize(readLen);
  for (auto i=0;i<readLen; i++) {
    readSeq[i] = seq_nt16_str[bam_seqi(read,i)];
  }
}

void GetItOfOverlappingVNTRs(vector<VNTR*> &vntrs,
		   map<string, vector<int> > &vntrMap, string &chrom, int regionStart, int regionEnd,
		   vector<int>::iterator &startIt,
		   vector<int>::iterator &endIt) {

  vector<int> &vntrIndex=vntrMap[chrom];
  LowerBoundSearchVNTRPos lbComp;
  UpperBoundSearchVNTRPos ubComp;   
  lbComp.vntrs = &vntrs;
  ubComp.vntrs = &vntrs;
  startIt = std::lower_bound(vntrIndex.begin(), vntrIndex.end(), regionStart, lbComp);
  endIt   = std::upper_bound(vntrIndex.begin(), vntrIndex.end(), regionEnd, ubComp);
}

		   

void IO::StoreReadsOnChrom(string &chrom, int regionStart, int regionEnd, vector<VNTR*> &vntrs,
			   map<string, vector<int> > &vntrMap,
			   Pileup &pileup, int thread, bool readsArePhased) {

   bam1_t * aln;
   hts_itr_t * itr;

   uint32_t isize; // observed template size
   uint32_t * cigar;
   uint8_t * s; // pointer to the read sequence
   bool rev;
   char * name;

   uint32_t ref_aln_start, ref_aln_end;
   uint32_t read_aln_start, read_len;
   uint32_t ref_len;
   uint32_t VNTR_s, VNTR_e;
   unsigned char base;
   uint32_t total_len;

   int i;
   VNTR * vntr = NULL;
   uint32_t origPhaseFlank = phaseFlank;
   cerr << "Storing TR sequences on " << chrom << ":" << regionStart << "-" << regionEnd << endl;
   if (vntrMap.find(chrom) == vntrMap.end()) {
     return;
   }
   // For now process all reads on chrom

   string region=MakeRegion(chrom, regionStart, regionEnd);
   itr = bam_itr_querys(idx, bamHdr, region.c_str());
   
   IOAlignBuff alignBuff(fp_in, itr, ioLock, 1000);
   int curVNTR=0;

   vector<int> &vntrIndex=vntrMap[chrom];

   int match = 3;
   int mismatch = -1;
   int gap_open = 3;
   int gap_extend = 3;
   const parasail_matrix_t* matrix = parasail_matrix_create("ACGT", match, mismatch);

   
   int nRead=0;
   CompPosWithSNV posSnvComp;
   int chromLen = faidx_seq_len(fai, chrom.c_str());
   stringstream regionStream;
   regionStream << chrom << ":1-" << chromLen;
   string chromRegion=regionStream.str();
   int chromSeqLen=0;
   char *chromSeq = fai_fetch(fai, chromRegion.c_str(), &chromSeqLen);
   string chromSeqStr(chromSeq);
   
   while (alignBuff.GetNext(aln) > 0) {
     string readSeq, readToRef;
     vector<int> readMap;
     if ((aln->core.flag & 256) or (aln->core.flag & 2048) ) {
       bam_destroy1(aln);
       continue;
     }
     if (aln->core.qual < minMapQV) {
       bam_destroy1(aln);
       continue;
     }
     
     StoreSeq(aln, readSeq);     
     StoreReadSeqAtRefCoord(aln, readSeq, readToRef, readMap);
     PiledRead read(aln->core.pos, readToRef);
     read.flag = aln->core.flag;
     int hap;
     bool readIsPhased = QueryAndSetPhase(aln, hap);

     if (chrom == "chrX" or chrom=="X" and minChrY > 0 and nChrY > minChrY) {
       readIsPhased = true;
       hap = 0;
     }

     int refAlnStart = aln->core.pos;
     int refAlnEnd = bam_endpos(aln);
     vector<SNV> snvs;
     if (readsArePhased == false ) {
       auto snvIt = lower_bound(pileup.hetSNVs.begin(), pileup.hetSNVs.end(), HetSNV(refAlnStart), posSnvComp);
       auto snvEnd = upper_bound(pileup.hetSNVs.begin(), pileup.hetSNVs.end(), HetSNV(refAlnEnd), posSnvComp);

       // Collect a list of snvs for this read.
       
       for (; snvIt != snvEnd; snvIt++) {
	 if (readToRef[snvIt->pos - refAlnStart] == snvIt->a) {
	   snvs.push_back(SNV(snvIt->pos, snvIt->a));
	 }
	 else if (readToRef[snvIt->pos - refAlnStart] == snvIt->b) {
	   snvs.push_back(SNV(snvIt->pos, snvIt->b));
	 }
       }
     }
     //
     // For each vntr that overlaps this read, add the read to the VNTR.
     //
     vector<int>::iterator it, endIt;
     GetItOfOverlappingVNTRs(vntrs, vntrMap, chrom, refAlnStart, refAlnEnd, it, endIt);
     
     while ( it != vntrIndex.end() and vntrs[*it]->ref_end <= refAlnEnd ) {
       //
       // Need to skip this vntr since it is outside the region. Here the loop exits because
       // all additional VNTRs will not overlap this read.
       //
       if (vntrs[*it]->ref_end < regionStart or vntrs[*it]->ref_start > regionEnd ) {
	 break;
       }

       //
       // Handle the corner case where a read starts at the same position as the tandem repeat.
       // Since this is not informative for the full length of the TR, do not process it.
       //
       if (vntrs[*it]->ref_start == refAlnStart) {
	 ++it;
	 continue;
       }

       int FLANK_REF_LEN=64;

       int refFlankStart = max((int)(0), (int)(vntrs[*it]->ref_start-1 -FLANK_REF_LEN));
       int refFlankEnd   = min((int)chromSeqLen, (int)( vntrs[*it]->ref_end-1 + FLANK_REF_LEN));

       int readStart = MappedStartPosInRead(readMap, refAlnStart, vntrs[*it]->ref_start-1);
       int readEnd   = MappedEndPosInRead(readMap, refAlnStart, vntrs[*it]->ref_end-1);

       int readFlankStart = MappedStartPosInRead(readMap, refAlnStart, refFlankStart);
       int readFlankEnd   = MappedEndPosInRead(readMap, refAlnStart, refFlankEnd);
       
       string vntrSeq;


       string refStartFlankSeq= chromSeqStr.substr(refFlankStart, vntrs[*it]->ref_start-1 - refFlankStart);
       string refEndFlankSeq= chromSeqStr.substr( vntrs[*it]->ref_end-1, refFlankEnd - vntrs[*it]->ref_end-1);

       for (auto ch=0; ch< refStartFlankSeq.size(); ch++) { refStartFlankSeq[ch] = toupper(refStartFlankSeq[ch]);}
       for (auto ch=0; ch< refEndFlankSeq.size(); ch++) { refEndFlankSeq[ch] = toupper(refEndFlankSeq[ch]);}       
       
       string readStartFlankSeq = readSeq.substr(readFlankStart, readStart - readFlankStart);
       string readEndFlankSeq   = readSeq.substr(readEnd, readFlankEnd - readEnd);

       parasail_result_t* startResult = NULL;
       parasail_cigar_t *startCigar = NULL;       
       int refinedReadStart = readStart;
       int refinedReadEnd = readEnd;
       if (readStartFlankSeq.size() >0 and refStartFlankSeq.size() >0){ 
	   startResult = parasail_sg_dx_trace_striped_sat( refStartFlankSeq.c_str(), refStartFlankSeq.length(),
						      readStartFlankSeq.c_str(), readStartFlankSeq.length(),						      
						      gap_open, gap_extend,
						      matrix);
	   startCigar = parasail_result_get_cigar(startResult, refStartFlankSeq.c_str(), refStartFlankSeq.length(),
								    readStartFlankSeq.c_str(), readStartFlankSeq.length(), matrix);

	   refinedReadStart=readFlankStart + startCigar->beg_ref;
	   if (startCigar) {
	     for (auto c=0; c < startCigar->len; c++) {
	       int cl = parasail_cigar_decode_len(startCigar->seq[c]);
	       char op = parasail_cigar_decode_op(startCigar->seq[c]);
	       if (op == '=' or op == 'X' or op=='M' or op=='D') {
		 refinedReadStart += cl;
	       }
	     }
	   }
       }
       parasail_result_t* endResult=NULL;
       parasail_cigar_t *endCigar = NULL;       
       if (readEndFlankSeq.size() >0 and refEndFlankSeq.size() >0){        
	 endResult = parasail_sg_dx_trace_striped_sat( refEndFlankSeq.c_str(), refEndFlankSeq.length(),
								       readEndFlankSeq.c_str(), readEndFlankSeq.length(),
								       gap_open, gap_extend,
								       matrix);
	 
	 endCigar=  parasail_result_get_cigar(endResult, refEndFlankSeq.c_str(), refEndFlankSeq.length(),
								readEndFlankSeq.c_str(), readEndFlankSeq.length(), matrix);
	 
	 refinedReadEnd = readEnd;
	 if (endCigar) {
	   for (auto c=0; c < endCigar->len; c++) {
	     int cl = parasail_cigar_decode_len(endCigar->seq[c]);
	     char op = parasail_cigar_decode_op(endCigar->seq[c]);
	     if (op == 'M' or op == '=' or op == 'X' or op == 'I') {
	       break;
	     }
	     if (op == 'D') {
	       refinedReadEnd += cl;
	       break;
	     }
	   }
	 }
       }


       if (startResult and endResult) {
	 readStart=refinedReadStart;
	 readEnd = refinedReadEnd;
	 parasail_result_free(startResult);
	 parasail_result_free(endResult);
	 parasail_cigar_free(startCigar);
	 parasail_cigar_free(endCigar);	 
       }
       else {

       }
       vntrSeq = readSeq.substr(readStart, readEnd-readStart);

       READ *read = new READ;
       read->qname = bam_get_qname(aln);
       read->l_qname = read->qname.size();
       read->seq = vntrSeq;
       read->len = vntrSeq.size();
       read->flag = aln->core.flag;

       //
       // Copy over all snvs
       //
       read->snvs = snvs;
       vntrs[*it]->reads.push_back(read);
       if (readIsPhased) {
	 read->haplotype=hap;
       }
       ++it;
     } // End looping over vntr loci
     bam_destroy1(aln);
   } // End looping over reads.
   bam_itr_destroy(itr);   
   if ( readsArePhased == false) {
     cerr << "Phasing reads in " << chrom << ":" << regionStart << "-" << regionEnd << endl;
     vector<int>::iterator it,end;
     GetItOfOverlappingVNTRs(vntrs, vntrMap, chrom, regionStart, regionEnd, it, end);
     
     int nProc=0;
     int total=end-it;
     int itPos = *it;

     int sz = vntrIndex.size();
     if (total < 0 ){
       int endSearch=itPos;
       while (endSearch < vntrs.size() and vntrs[endSearch]->ref_start < regionEnd) {
	 endSearch++;
       }
       //       cerr << "Ended search from " << itPos << " at " << endSearch << endl;
       return;
     }
     assert(total >=0);
     
     while (it != end ) {
       MaxCutPhase(vntrs[*it], phaseFlank);
       ++it;
       ++nProc;
     }
   }
   free(chromSeq);
}

int IO::CallSNVs(string &chrom, int regionStart, int regionEnd,  vector<VNTR*> &vntrs,
		  map<string, vector<int> > &vntrMap,
		 Pileup &pileup, bool & readsArePhased, OPTION &options) {
   bam1_t * aln;
   hts_itr_t * itr;

   uint32_t isize; // observed template size
   uint32_t * cigar;
   uint8_t * s; // pointer to the read sequence
   bool rev;
   char * name;

   uint32_t ref_aln_start, ref_aln_end;
   uint32_t read_aln_start, read_len;
   uint32_t ref_len;
   uint32_t VNTR_s, VNTR_e;
   unsigned char base;
   uint32_t total_len;

   int i;
   VNTR * vntr = NULL;
   uint32_t origPhaseFlank = phaseFlank;

   if (vntrMap.find(chrom) == vntrMap.end()) {
     return 0;
   }
   // For now process all reads on chrom
   cerr << "Calling snvs: " << chrom << ":" << regionStart << "-" << regionEnd << endl;
   string region=MakeRegion(chrom, regionStart, regionEnd);
   itr = bam_itr_querys(idx, bamHdr, region.c_str());
   IOAlignBuff alignBuff(fp_in, itr, ioLock, 1000);
   int curVNTR=0;

   vector<int> &vntrIndex=vntrMap[chrom];
   int nRead=0;
   while (alignBuff.GetNext(aln) > 0) {
     int tid=aln->core.tid;
     if (tid != -1) {
       string refName(bamHdr->target_name[tid]);
     }
     
     ++nRead;
     int refStartPos = aln->core.pos;
     int refAlnEnd = bam_endpos(aln);
     // If not overlapping anything, continue reading.
     if ((aln->core.flag & 256) or (aln->core.flag & 2048) ) {
       bam_destroy1(aln);
       //       aln=bam_init1();
       continue;
     }
     if (aln->core.qual < minMapQV) {
       bam_destroy1(aln);
       continue;
     }
     //
     // Add a read to the pileup. This increments all counts of snps
     //
     string readSeq, readToRef;
     vector<int> readMap;
     StoreSeq(aln, readSeq);     
     StoreReadSeqAtRefCoord(aln, readSeq, readToRef, readMap);
     PiledRead read(aln->core.pos, readToRef);
     read.flag = aln->core.flag;
     int hap;
     bool readIsPhased = QueryAndSetPhase(aln, hap);
     read.phased = readIsPhased;
     read.hap = hap;
     if (readIsPhased) {
       readsArePhased= true;
       alignBuff.Free();
       cerr << "Reads are phased. Skipping SNV calling." << endl;
       return -1;
     }

     //
     // Short circuit phasing for chrX and passing y chrom limit.
     //
     if (chrom == "chrX" or chrom=="X" and minChrY > 0 and nChrY > minChrY) {
       readIsPhased = true;
       hap = 0;
       read.phased = true;
       read.hap = 0;
       return -1;
     }
     
     read.readName = bam_get_qname(aln);
     pileup.AddRead(read);
     while(pileup.endPosPQueue.top().pos <= refStartPos) {
       //
       // Found a read that ends before the current read starts. That means
       // no additional phasing information will be added by future alignments
       // and positions up to this read.
       int nextEndPos = pileup.endPosPQueue.top().pos;
       pileup.RemoveRead(pileup.endPosPQueue.top().it);
       pileup.endPosPQueue.pop();
     }
     
     int finishedCons=pileup.consensusStart;
     int tmpCov = pileup.consensus[finishedCons].cov;
     while (pileup.consensus[finishedCons].cov == 0) {
       finishedCons++;
     }
     int before=pileup.consensus.size();
     if (readIsPhased == false) {
       pileup.ProcessUntil(finishedCons, options);
     }
     bam_destroy1(aln);
   }
   bam_itr_destroy(itr);
   pileup.ProcessUntil(pileup.consensusEnd, options);
   return 1;
}

/*
 read alignment from bam
 liftover ref_VNTR_start, ref_VNTR_end of every vntr
 get the subsequence 
*/
string &FlankSeq(READ* read, int side) {
    if (side == 0) {
        return read->upstream;
    }
    else {
        return read->downstream;
    }
}

void IO::initializeRefFasta() {
  return;
  if (reference =="") {
    cerr << "ERROR. Vamos now requires a reference to be specified with -R ref.fasta ." << endl;
  }
  fai = fai_load(reference.c_str());
  if (fai == NULL) {
    cerr << "Could not load reference " << reference << endl;
  }
}

void IO::initializeBam() {
    fp_in = hts_open(input_bam.c_str(), "r"); //open bam file
    if (!reference.empty()) {
        hts_set_opt(fp_in, CRAM_OPT_REFERENCE, reference.c_str());
    }
    bamHdr = sam_hdr_read(fp_in); //read header
    idx = sam_index_load(fp_in, input_bam.c_str()); //samtools will implicitly search for the appropriate header if given the filename
    for (int i = 0; i < bamHdr->n_targets; ++i) {
      chromosomeNames.push_back(bamHdr->target_name[i]);
      chromosomeLengths.push_back(bamHdr->target_len[i]);      
    }

    if (minChrY > 0) {
      for (int tid = 0; tid < bamHdr->n_targets; ++tid) {
        if (strcmp(bamHdr->target_name[tid], "chrY") == 0) {
	  uint64_t unmapped;
	  hts_idx_get_stat(idx, tid, &nChrY, &unmapped);
	  break;
        }
      }
    }
}

void InitNucIndex(char nucIndex[]) {
    memset(nucIndex, 254, 256);
    nucIndex['a']=0;  nucIndex['A']=0;
    nucIndex['c']=1;  nucIndex['C']=1;
    nucIndex['g']=2;  nucIndex['G']=2;
    nucIndex['t']=3;  nucIndex['T']=3;
    nucIndex['-']=4;
    // Special character for flanking sequence. 
    nucIndex['F']=5;
}


void IO::closeBam() {
    bam_hdr_destroy(bamHdr);
    hts_idx_destroy(idx);
    sam_close(fp_in); 
}
/* read consensus read in */
void IO::readSeqFromFasta(vector<VNTR *> &vntrs)
{
    int i;
    assert(vntrs.size() == 1);
    gzFile fp;  
    kseq_t *seq;  
    int l; 
    fp = gzopen(input_bam.c_str(), "r"); // STEP 2: open the file handler  
    seq = kseq_init(fp); // STEP 3: initialize seq  
    while ((l = kseq_read(seq)) >= 0) { // STEP 4: read sequence  
        READ * read = new READ();
        read->l_qname = seq->name.l;
        read->qname = seq->name.s; //bam_get_qnamenew char[read->l_qname + 1];
	read->seq = seq->seq.s;
	//        read->seq = new char[seq->seq.l + 1]; // read sequence array
	//        strcpy(read->seq, seq->seq.s);
        read->len = seq->seq.l;
	//	read->seq[seq->seq.l] = '\0';

        vntrs[0]->reads.push_back(read); 
        if (debug_flag)
        {
             cerr << "read_name: "; 
             cout << read->qname;
             cerr << endl;
             cerr << "read length: " << read->len << endl;
             cout << read->seq;
             cout << endl; 
        }
    }  
    vntrs[0]->nreads = vntrs[0]->reads.size();
    assert(vntrs[0]->nreads == 1);
    vntrs[0]->cur_len = vntrs[0]->reads[0]->len;
    kseq_destroy(seq); // STEP 5: destroy seq  
    gzclose(fp); // STEP 6: close the file handler  
    return;
}

int IO::writeVCFHeader_locuswise(ofstream &out)
{
    outWriter.init(input_bam, version, sampleName);
    outWriter.writeHeader_locuswise(out);
    return 0;
}

int IO::writeVCFBody_locuswise(ofstream& out, vector<VNTR *> &vntrs, int tid, int nproc)
{
    outWriter.writeBody_locuswise(vntrs, out, tid, nproc);
    return 0;
}


int IO::writeBEDHeader_readwise(ofstream &out)
{
    outWriter.init(input_bam, version, sampleName);
    outWriter.writeHeader_readwise(out);
    return 0;
}

int IO::writeBEDBody_readwise(ofstream& out, vector<VNTR *> &vntrs, int tid, int nproc)
{
    outWriter.writeBody_readwise(vntrs, out, tid, nproc);
    return 0;
}

void IO::writeFa(ofstream& out, vector<VNTR *> &vntrs)
{
    for (auto &vntr : vntrs)
    {
        for (auto &read : vntr->reads)
        {
	  out << ">" << read->qname << endl;
	  out << read->seq << endl;
        }
    }
    return;
}

