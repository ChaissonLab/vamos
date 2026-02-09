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

void IO::UnrefinedStoreReadsOnChrom(string &chrom, int regionStartIndex, int regionEndIndex, vector<VNTR*> &vntrs,
				    map<string, vector<int> > &vntrMap,
				    Pileup &pileup, int thread, bool readsArePhased, int oneOffset) {

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
   int regionStart =    vntrs[regionStartIndex]->ref_start;
   int regionEnd   =vntrs[regionEndIndex]->ref_end;
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

   
   int nRead=0;
   CompPosWithSNV posSnvComp;


   while (alignBuff.GetNext(aln) > 0) {
     string readSeq, readToRef, readName;
     readName = bam_get_qname(aln);
     
     vector<int> readMap;
     if ((aln->core.flag & 256) or (aln->core.flag & 2048) ) {
       bam_destroy1(aln);
       //       cout << "Skipping secondary alignment " << readName << endl;
       continue;
     }
     if (aln->core.qual < minMapQV) {
       bam_destroy1(aln);
       //       cout << "Skipping minmapqv " << readName << endl;
       continue;
     }
     
     StoreSeq(aln, readSeq);     
     StoreReadSeqAtRefCoord(aln, readSeq, readToRef, readMap);
     PiledRead read(aln->core.pos, readToRef);
     read.flag = aln->core.flag;
     int hap=0;
     bool readIsPhased = QueryAndSetPhase(aln, hap);
     if ((chrom == "chrX" or chrom=="X") and minChrY > 0 and nChrY > minChrY) {
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
     //     cout << readName << "\t" << nRead << " read " << refAlnStart << "\t" << refAlnEnd << "\t" << hap << endl;
     ++nRead;
     while ( it != vntrIndex.end() and vntrs[*it]->ref_end <= refAlnEnd ) {
       //
       // Need to skip this vntr since it is outside the region. Here the loop exits because
       // all additional VNTRs will not overlap this read.
       //
       //       cout << "Trying for vntr  " << vntrs[*it]->region << endl;
       if (vntrs[*it]->ref_end < regionStart or vntrs[*it]->ref_start > regionEnd ) {
	 //	 cout << "skipping end before" << endl;
	 break;
       }

       //
       // Handle the corner case where a read starts at the same position as the tandem repeat.
       // Since this is not informative for the full length of the TR, do not process it.
       //
       if (vntrs[*it]->ref_start == refAlnStart) {
	 //	 cout << "Skipping starts at aln start - possible that entire seq is not captured. " << endl;
	 ++it;
	 continue;
       }

       string vntrSeq;
       int startOffset=1;
       if (vntrs[*it]->ref_start-1 - startOffset <= 0) {
	 startOffset = 0;
       }
       int endOffset = 1;
       
       int readStart = MappedStartPosInRead(readMap, refAlnStart, vntrs[*it]->ref_start-oneOffset-startOffset);
       int readEnd = MappedEndPosInRead(readMap, refAlnStart, vntrs[*it]->ref_end-1+endOffset);
       if (readEnd - readStart <= 0) {
	 //	 cout << "Skipping due to read end less than or equal to start" << endl;
	 ++it;	 
	 continue;
       }

       int normalStart = MappedStartPosInRead(readMap, refAlnStart, vntrs[*it]->ref_start-oneOffset);
       int normalEnd = MappedEndPosInRead(readMap, refAlnStart, vntrs[*it]->ref_end-1);
       
       while (readEnd < 0 and endOffset >= -1) {
	 endOffset--;
	 readEnd = MappedEndPosInRead(readMap, refAlnStart, vntrs[*it]->ref_end-1+endOffset);
       }

       if (readEnd < 0) {
	 // No suitable alignment found
	 ++it;	 
	 continue;
       }
       vntrSeq = readSeq.substr(readStart+startOffset, readEnd-readStart-endOffset);
       if (vntrSeq.size() == 0) {
	 //	 cout << "Empty vntr" << endl;
	 ++it;	 
	 continue;
       }
       if (vntrSeq.size() > 40000) {
	 ++it;
	 continue;
       }
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
       //       cout << "Added TR " << vntrSeq << "\t" << hap << "\t" << (int) readIsPhased << endl;
       vntrs[*it]->reads.push_back(read);
       if (readIsPhased) {
	 read->haplotype=hap;
       }
       ++it;
     } // End looping over vntr loci for this alignment.
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
}
