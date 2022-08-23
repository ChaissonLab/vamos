#include <iostream>
#include <stdlib.h>
#include <istream>
#include <fstream>
#include <sstream>
#include <assert.h>
#include <string>
#include <algorithm>
#include <tuple> 
#include <vector>
#include <regex>
#include "io.h"
#include "read.h"
#include "vntr.h"
#include "vcf.h"
#include "abpoa.h"
// #include <seqan/align.h>
// #include <seqan/graph_msa.h>
#include "htslib/hts.h"
#include "htslib/sam.h"
#include <zlib.h>  
#include "htslib/kseq.h"  
#include "phase.h"

extern int naive_flag;
extern int debug_flag;
extern int hclust_flag;
extern int seqan_flag;
extern int liftover_flag;
extern int conseq_anno_flag;
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
    while (getline(ifs, line)) 
    {
        stringstream ss(line);
        string tmp;
        string tmp_r;
        string chrom;
        int start;
        int end;
        ss >> chrom >> start >> end;
        VNTR * vntr = new VNTR(chrom, start, end, end-start);
        vntrs.push_back(vntr);
        while(getline(ss, tmp, ',')) 
        {
            tmp_r = regex_replace(tmp, std::regex("^\\t+"), std::string(""));
            vntrs[numOfLine]->motifs.push_back(MOTIF(tmp_r));
            // cerr << vntrs[numOfLine]->motifs.back().seq << endl;
        }
        numOfLine += 1; // 0-indexed
    }
    assert(vntrs.size() == numOfLine);
    return 0;
}



/* read vntrs coordinates from file `vntr_bed`*/
void IO::readVNTRFromBed (vector<VNTR*> &vntrs)
{
    vector<vector<string>> items;
    read_tsv(items);
    uint32_t start, end, len;
    for (auto &it : items)
    {
        start = stoi(it[1]);
        end = stoi(it[2]);
        len = start < end ? end - start : 0;
        VNTR * vntr = new VNTR(it[0], start, end, len);
        vntrs.push_back(vntr);
    }
    return;
}

int FRONT=0;
int BACK=1;

void StoreReadSeqAtRefCoord(bam1_t *aln, uint32_t refRangeStart, uint32_t refRangeEnd, string &seq) {
  uint32_t readPos=0;
  uint32_t refEnd=bam_endpos(aln);
  uint32_t refPos=aln->core.pos;
  uint32_t refStartPos=refPos;    
  uint32_t *cigar = bam_get_cigar(aln);
  int op, type, len;
  uint8_t *read;
  uint32_t readLen = aln->core.l_qseq;
  // If there is a prefix gap, add those.
  for (int i=refRangeStart; i < refPos; i++) {
    seq.push_back('F');
  }
  read=bam_get_seq(aln);
  char *name = bam_get_qname(aln);
  string nucSeq;
  nucSeq.resize(readLen);
  for (auto i=0;i<readLen; i++) {
    nucSeq[i] = seq_nt16_str[bam_seqi(read,i)];
  }
  for (uint32_t k = 0; k < aln->core.n_cigar && refPos < refRangeEnd && readPos < readLen; k++) {
    assert(readPos < readLen);
    op = bam_cigar_op(cigar[k]);
    type = bam_cigar_type(op);
    len = bam_cigar_oplen(cigar[k]);

    if (!(type & 1) and !(type & 2)) // Hard clip
      continue;
    // A match or mismatch type
    if ((type & 1) && (type & 2)) {
      for (int i=0; i < len; i++ ) {
	if (refPos >= refRangeStart and refPos < refRangeEnd) {	  
	  seq.push_back(nucSeq[readPos]);
	}
	readPos++;
	refPos++;
      }
    }
    else if ( (type & 1) == 0 && (type & 2) != 0 ) {
      for (int i=0; i < len; i++) {
	if (refPos >= refRangeStart and refPos < refRangeEnd) {	  	
	  seq.push_back('F');
	}
	refPos +=1;
      }
    }
    else if ((type & 1) != 0 && (type & 2) == 0) {
      readPos+=len;
    }
  }
  if (refPos < refRangeStart) { refPos = refRangeStart;}
  while (refPos < refRangeEnd) {
    seq.push_back('F');
    refPos++;
  }
  assert(seq.size() == refRangeEnd-refRangeStart);
}

pair<uint32_t, bool> processCigar(bam1_t * aln, uint32_t * cigar, uint32_t &CIGAR_start, uint32_t target_crd, uint32_t &ref_aln_start, uint32_t &read_aln_start)
{
    assert(read_aln_start <= (uint32_t) aln->core.l_qseq);
    /* trivial case: when ref_aln_start equals target_crd*/
    if (target_crd == ref_aln_start) return make_pair(read_aln_start, 1);

    uint32_t cigar_start = CIGAR_start;
    int op, type, len;

    for (uint32_t k = cigar_start; k < aln->core.n_cigar; k++) 
    {
        op = bam_cigar_op(cigar[k]);
        type = bam_cigar_type(op);
        len = bam_cigar_oplen(cigar[k]);

        if (ref_aln_start > target_crd) 
            return make_pair(read_aln_start, 0); // skip the out-of-range alignment

        else if (ref_aln_start == target_crd)
            return make_pair(read_aln_start, 1);

        else if (!(type & 1) and !(type & 2)) // Hard clip
            continue;

        if (type & 2 or op == BAM_CMATCH or op == BAM_CEQUAL or op == BAM_CDIFF) 
        {
            if (target_crd < ref_aln_start + len) 
            {
                 if (op == BAM_CMATCH or op == BAM_CEQUAL) 
                    return make_pair(read_aln_start + target_crd - ref_aln_start, 1);
                else
                    return make_pair(read_aln_start, 1);                
            }
              
        }

        CIGAR_start = k + 1;
        if (type & 1) read_aln_start += len;        
        if (type & 2) ref_aln_start += len;

        // // for debug 
        // uint32_t callen = bam_cigar2qlen(k + 1, cigar);
        // assert(read_aln_start == callen);
    }

    assert(read_aln_start <= (uint32_t) aln->core.l_qseq);
    //
    // Reached the end of the cigar but did not find a match overlapping
    // the target_crd. Return a pair that flags invalid search.
    //
    return make_pair(read_aln_start, 0);
}

int GetHap(char *s, int l) {
  int si=0, nc=0;
  for (; si < l; si++) {
    if (s[si] == ':') { nc+=1;}
    if (nc == 2) { si++; break; }
  }
  if (nc == 2) {
    return atoi(&s[si]);
  }
  return -1;
}

/*
 read alignment from bam
 liftover ref_VNTR_start, ref_VNTR_end of every vntr
 get the subsequence 
*/
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

string &FlankSeq(READ* read, int side) {
  if (side == 0) {
    return read->upstream;
  }
  else {
    return read->downstream;
  }
}
void SimpleSNV(VNTR *vntr, uint32_t start, uint32_t end, int side=0) {
/* vamos -in read.bam -vntr vntrs.bed -motif motifs.csv -o out.vcf */
  char nucIndex[256];
  InitNucIndex(nucIndex);
  vector<vector<int> > counts;
  int len=end-start;
  int nChars=6;
  counts.resize(nChars);
  for (auto i=0; i <nChars; i++) {
    counts[i].resize(len, 0);
  }
  for (auto i=0; i < vntr->nreads; i++) {
    for (auto p=0; p < FlankSeq(vntr->reads[i], side).size(); p++) {
      counts[nucIndex[FlankSeq(vntr->reads[i], side)[p]]][p] += 1;
    }
  }
  //  cout << "simpsnv ";
  vector<int> snvIdx;
  vector<bool> cluster(counts[0].size(), false);
  int prevIndex=-1;
  for (auto j=0; j < counts[0].size(); j++) {
    int colSum=0;
    // Exclude flanking sequences in total count.
    for (auto i=0; i < 5; i++) {
      colSum += counts[i][j];
    }
    if (colSum > 0) {
      vector<float> frac(4,0);
      vector<int> suffIdx;
      vector<char> snvNuc;
      for (auto i=0; i < 4; i++) {
	frac[i] = ((float)counts[i][j]/colSum);
	if (frac[i] > 0.3) {
	  suffIdx.push_back(i);
	  snvNuc.push_back("ACGT"[i]);
	}
      }
      if (suffIdx.size() == 2) {
	if (prevIndex >= 0 and j - prevIndex < 50) {
	  cluster[prevIndex] = true;
	  cluster[j] = true;
	}
	prevIndex=j;
      }
    }
  }
  
  for (auto j=0; j < counts[0].size(); j++) {
    int colSum=0;
    // Exclude flanking sequences in total count.
    for (auto i=0; i < 5; i++) {
      colSum += counts[i][j];
    }
    if (colSum > 0) {
      vector<float> frac(4,0);
      vector<int> suffIdx;
      vector<char> snvNuc;
      for (auto i=0; i < 4; i++) {
	frac[i] = ((float)counts[i][j]/colSum);
	if (frac[i] > 0.3) {
	  suffIdx.push_back(i);
	  snvNuc.push_back("ACGT"[i]);
	}
      }
      if (suffIdx.size() == 2 and cluster[j] == false) {
	/*	for (int k=0; k < vntr->nreads; k++) {
	  cout << vntr->reads[k]->upstream[j];
	}
	cout << endl;*/
	snvIdx.push_back(j);
	for (auto ri=0; ri < vntr->nreads; ri++) {
	  SNV snv;
	  snv.nuc='-';
	  snv.pos=start + j;
	  if (FlankSeq(vntr->reads[ri], side)[j] == snvNuc[0]) {
	    snv.nuc=snvNuc[0];
	  }
	  if (FlankSeq(vntr->reads[ri], side)[j] == snvNuc[1]) {
	    snv.nuc=snvNuc[1];
	  }
	  if (snv.nuc != '-') {
	    vntr->reads[ri]->snvs.push_back(snv);
	  }
	}      
	//	cout << "*";
	//
	// Found a het site
	//	cout << "HET: " << vntr->chr << ":" << start + j << "-" << start + j + 1 << " " << "ACGT"[suffIdx[0]] << " " << counts[suffIdx[0]][j] << " " << "ACGT"[suffIdx[1]] << " " << counts[suffIdx[1]][j] << endl;
      }
      else {
	//	cout << " ";
      }
    }
  }
  //  cout << endl;
  for (auto si=0; si < snvIdx.size(); si++ ) {
    vector<float> frac(4,0);
    vector<int> suffIdx;
    int colSum=0;
    for (auto i=0; i < 4; i++) {
      colSum+= counts[i][snvIdx[si]];
    }
    for (auto i=0; i < 4; i++) {
      frac[i] = ((float)counts[i][snvIdx[si]]/colSum);
      if (frac[i] > 0.3) {
	suffIdx.push_back(i);
      }
    }
    if (suffIdx.size() == 2) {
      //      cout << "HET: " << vntr->chr << ":" << start + snvIdx[si] << "-" << start + snvIdx[si] + 1 << " " << "ACGT"[suffIdx[0]] << " " << counts[suffIdx[0]][snvIdx[si]] << " " << "ACGT"[suffIdx[1]] << " " << counts[suffIdx[1]][snvIdx[si]] << endl;
    }
  }
}


void IO::readSeqFromBam (vector<VNTR *> &vntrs, int nproc, int cur_thread, int sz) 
{
    char * bai = (char *) malloc(strlen(input_bam) + 4 + 1); // input_bam.bai
    strcpy(bai, input_bam);
    strcat(bai, ".bai");

    samFile * fp_in = hts_open(input_bam, "r"); //open bam file
    bam_hdr_t * bamHdr = sam_hdr_read(fp_in); //read header
    hts_idx_t * idx = sam_index_load(fp_in, bai);
    bam1_t * aln = bam_init1(); //initialize an alignment
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
    // vector<READ *> initial_reads;
    uint32_t total_len;

    int i;
    VNTR * vntr = NULL;

    for (int j = cur_thread; j < sz; j += nproc) 
    {
        vntr = vntrs[j];
        total_len = 0;
        itr = bam_itr_querys(idx, bamHdr, vntr->region.c_str());
        while(bam_itr_next(fp_in, itr, aln) >= 0)
        {
            rev = bam_is_rev(aln);
            if (aln->core.flag & BAM_FSECONDARY or aln->core.flag & BAM_FUNMAP) 
                continue; // skip secondary alignment / unmapped reads

            if (aln->core.qual < 20) 
                continue; // skip alignment with mapping qual < 20

            isize = aln->core.isize;
            if (isize == 0) // 9th col: Tlen is unavailable in bam
                isize = bam_endpos(aln);

            if (isize > 0 and isize < (uint32_t) vntr->len)
                continue; // skip alignment with length < vntr->len

            cigar = bam_get_cigar(aln);
            ref_aln_start = aln->core.pos;
            ref_aln_end = ref_aln_start + isize;
            ref_len = bamHdr->target_len[aln->core.tid]; 

            read_aln_start = 0;
            read_len = aln->core.l_qseq;
	    s=bam_get_seq(aln);
            uint32_t tmp;
            VNTR_s = vntr->ref_start;
            VNTR_e = vntr->ref_end;
	    //	    kstring_t s;
	    kstring_t auxStr = KS_INITIALIZE;
	    int auxStat =  bam_aux_get_str(aln, "HP", &auxStr);
	    int hap = GetHap(auxStr.s, auxStr.l);
            if (VNTR_s < ref_aln_start or VNTR_e > ref_aln_end) // the alignment doesn't fully cover the VNTR locus
                continue;

            /*
            reference: VNTR_s, VNTR_e
            read alignment: ref_aln_start, ref_aln_end; read_aln_start, read_aln_end;
            */
            uint32_t cigar_start = 0;
            auto [liftover_read_s, iflift_s] = processCigar(aln, cigar, cigar_start, VNTR_s, ref_aln_start, read_aln_start);
            auto [liftover_read_e, iflift_e] = processCigar(aln, cigar, cigar_start, VNTR_e, ref_aln_start, read_aln_start);
	    string upstream, downstream;
            // [liftover_read_s, liftover_read_e]
            if (iflift_s and iflift_e and liftover_read_e > liftover_read_s and liftover_read_e <= read_len)
            {
                liftover_read_s = (liftover_read_s == 0) ? 0 : liftover_read_s - 1;
                liftover_read_e = (liftover_read_e == 0) ? 0 : liftover_read_e - 1;

                assert(liftover_read_e < read_len);

                READ * read = new READ();
                read->chr = bamHdr->target_name[aln->core.tid]; 
                read->qname = (char *) malloc(aln->core.l_qname + 1);
                name = bam_get_qname(aln);
                strcpy(read->qname, name);
                string tmp_name(read->qname);
                read->l_qname = tmp_name.length();
                read->len = liftover_read_e - liftover_read_s + 1; // read length
                read->seq = (char *) malloc(read->len + 1); // read sequence array
                read->rev = rev;

		StoreReadSeqAtRefCoord(aln, vntr->ref_start - phaseFlank, vntr->ref_start, read->upstream);
		StoreReadSeqAtRefCoord(aln, vntr->ref_end, vntr->ref_end+phaseFlank, read->downstream);
		//		cout << "prefix: " << read->upstream << endl;
		if (hap >= 0) {
		  read->haplotype=hap;
		}

		kstring_t auxStr = KS_INITIALIZE;
		int auxStat =  bam_aux_get_str(aln, "HP", &auxStr);
		if (auxStat ) {
		  
		  int si=0, nc=0;
		  int hap=GetHap(auxStr.s, auxStr.l);
		  if (hap >= 0 ) {
		    read->haplotype=hap;
		  }
		}
                for(i = 0; i < read->len; i++)
                {
                    assert(i + liftover_read_s < read_len);
                    base = bam_seqi(s, i + liftover_read_s);
                    assert(0 < base < 16);
                    read->seq[i] = seq_nt16_str[base]; //gets nucleotide id and converts them into IUPAC id.
                } 
                vntr->reads.push_back(read); 

                total_len += read->len;

                if (debug_flag)
                {
                     cerr << "read_name: " << bam_get_qname(aln) << endl; 
                     // cerr << "vntr->ref_start: " << vntr->ref_start << " vntr->ref_end: " << vntr->ref_end << endl;
                     // cerr << "liftover_read_s: " << liftover_read_s << " liftover_read_e: " << liftover_read_e << endl;
                     cerr << "read length: " << read->len << endl;
                     cout.write(read->seq, read->len);
                     cout << endl; 
                }
            }
        }
        vntr->nreads = vntr->reads.size();
	SimpleSNV(vntr, vntr->ref_start- phaseFlank, vntr->ref_start, 0);
	SimpleSNV(vntr, vntr->ref_end, vntr->ref_end + phaseFlank, 1);
	MaxCutPhase(vntr);
	
        vntr->cur_len = (vntr->nreads == 0) ? 0 : (total_len / vntr->nreads);
    }
    free(bai);
    bam_destroy1(aln);
    bam_hdr_destroy(bamHdr);
    hts_itr_destroy(itr);
    hts_idx_destroy(idx);
    sam_close(fp_in); 
    return;  
}

/* read consensus read in */
void IO::readSeqFromFasta(vector<VNTR *> &vntrs)
{
    int i;
    assert(vntrs.size() == 1);
    gzFile fp;  
    kseq_t *seq;  
    int l; 
    fp = gzopen(input_bam, "r"); // STEP 2: open the file handler  
    seq = kseq_init(fp); // STEP 3: initialize seq  
    while ((l = kseq_read(seq)) >= 0) { // STEP 4: read sequence  
        READ * read = new READ();
        read->l_qname = seq->name.l;
        read->qname = (char *) malloc(read->l_qname + 1);
        strcpy(read->qname, seq->name.s);

        read->seq = (char *) malloc(seq->seq.l + 1); // read sequence array
        strcpy(read->seq, seq->seq.s);
        read->len = seq->seq.l;

        vntrs[0]->reads.push_back(read); 
        if (debug_flag)
        {
             cerr << "read_name: "; 
             cout.write(read->qname, read->l_qname);
             cerr << endl;
             cerr << "read length: " << read->len << endl;
             cout.write(read->seq, read->len);
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
            out << ">"; 
            out.write(read->qname, read->l_qname);
            out << "\n";
            out.write(read->seq, read->len);
            out << "\n";
        }
    }
    return;
}
