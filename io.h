#ifndef IO_H_
#define IO_H_

#include <iostream>
#include <stdlib.h>
#include <istream>
#include <fstream>
#include <assert.h>
#include <string>
#include <tuple> 
#include "read.h"
#include "htslib/htslib/sam.h"
using namespace std;

void readMotifsFromCsv (const string &input_motif_csv, vector<VNTR> &vntrs) 
{
    int vntr_size = vntrs.size();

    ifstream ifs(input_motif_csv);
    if (ifs.fail()) 
    {
        cerr << "error" << endl;
        return;
    }    

    string line;
    int numOfLine = 0;
    while (getline(ifs, line)) {
        stringstream ss(line);
        string tmp;
        while(getline(ss, tmp, ',')) {
            vntrs[numOfLine].motifs.push_back(tmp);
        }
        numOfLine += 1; // 0-indexed
    }
    assert(vntrs.size() == numOfLine.size() + 1);
    return;
}

void read_tsv(const string &input_vntr_bed, vector<vector<string>> &items) 
{
    items.clear();
    ifstream ifs(input_vntr_bed);
    if (ifs.fail()) 
    {
        cerr << "error" << endl;
        return;
    }

    string line;
    while (getline(ifs, line)) {
        stringstream ss(line);
        vector<string> item;
        string tmp;
        while(getline(ss, tmp, '\t')) {
            item.push_back(tmp);
        }
        items.push_back(item);
    }
    return;
}

pair<uint32_t, bool> processCigar(bam1_t * aln, uint32_t * cigar, uint32_t &cigar_start, uint32_t target_crd, uint32_t &ref_aln_start, uint32_t &read_aln_start)
{
    uint32_t k;
    for (k = cigar_start; k < aln->core.n_cigar; ++k) 
    {
        int op = bam_cigar_op(cigar[k]);
        int l = bam_cigar_oplen(cigar[k]);

        if (ref_aln_start > target_crd) return make_pair(read_aln_start, 0); // skip the out-of-range alignment
        else if (ref_aln_start <= target_crd and target_crd < ref_aln_start + l)
        {
            if (op == BAM_CMATCH) read_aln_start += target_crd - ref_aln_start;
            break;
        }

        switch (op) 
        {
        case BAM_CDEL:
            read_aln_start += l;  
            break;    

         case BAM_CINS:
            ref_aln_start += l;
            break;            

         case BAM_CMATCH:
            read_aln_start += l;
            ref_aln_start += l;
            break;

         case BAM_CDIFF:
            read_aln_start += l;
            ref_aln_start += l;
            break;

         case BAM_CHARD_CLIP:
            break;

         case BAM_CSOFT_CLIP:
            read_aln_start += l;
            break;
        }
    }
    return make_pair(read_aln_start, 1);
}

/*
 read alignment from bam overlapping a VNTR locus
 liftover ref_VNTR_start, ref_VNTR_end
 get the subsequence 
*/
void readSeqFromBam (const string &input_bam_file, string &chr, uint32_t ref_VNTR_start, uint32_t ref_VNTR_end, uint32_t VNTR_len) 
{
    bamFile * fp_in = bam_open(input_bam_file, "r"); //open bam file
    bam_hdr_t * bamHdr = bam_hdr_read(fp_in); //read header
    hts_idx_t * idx = bam_index_load(fp, input_bam_file + ".bai");
    bam1_t * aln = bam_init1(); //initialize an alignment
    hts_itr_t * itr = bam_itr_querys(idx, bamHdr, chr + "." + to_string(ref_VNTR_start) + '-' + to_string(ref_VNTR_end));

    uint16_t flag;
    uint32_t mapq;
    uint32_t isize; // observed template size
    uint32_t * cigar;
    uint8_t * s; // pointer to the read sequence
    uint32_t * cigar;

    uint32_t ref_aln_start, ref_aln_end;
    uint32_t read_aln_start, read_aln_end;
    uint32_t ref_len;
    uint32_t VNTR_s, VNTR_e;
    uint32_t k;
    uint32_t liftover_read_s, liftover_read_e;

    while(bam_itr_next(fp_in, itr, aln) >= 0)
    {
        flag = aln->core.flag;
        if (flag & BAM_FSECONDARY or flag & BAM_FUNMAP) continue; // skip secondary alignment / unmapped reads

        mapq = aln->core.qual;
        if (mapq < 20) continue; // skip alignment with mapping qual < 20

        isize = ali->core.isize;
        if (isize < VNTR_len) continue; // skip alignment with length < VNTR_len

        cigar = bam_get_cigar(b);
        ref_aln_start = aln->core.pos;
        ref_aln_end = ref_aln_start + isize;

        read_aln_start = 0;
        read_aln_end = aln->core.l_qseq;

        ref_len = *(bamHdr->target_len[aln->core.tid]); 

        uint32_t tmp;
        VNTR_s = ref_VNTR_start;
        VNTR_e = ref_VNTR_end;
        if (bam_is_rev(aln)) 
        {   
            tmp = VNTR_s;
            VNTR_s = ref_len - VNTR_e;
            VNTR_e = ref_len - tmp - 1;

            tmp = ref_aln_start;
            ref_aln_start = ref_len - ref_aln_end;
            ref_aln_end = ref_len - ref_aln_start - 1;
        }

        /*
        reference: VNTR_s, VNTR_e
        read alignment: ref_aln_start, ref_aln_end; read_aln_start, read_aln_end;
        */
        k = 0;
        auto [liftover_read_s, iflift_s] = processCigar(aln, cigar, k, VNTR_s, ref_aln_start, read_aln_start);
        auto [liftover_read_e, iflift_e] = processCigar(aln, cigar, k, VNTR_e, ref_aln_start, read_aln_start);

        if (iflift_s and iflift_e and liftover_read_e > liftover_read_s)
        {

            READ * read = new READ();
            read->chr = bamHdr->target_name[aln->core.tid]; 
            read->len = liftover_read_e - liftover_read_s; // read length
            read->seq = (char *) malloc(read->len); // read sequence array
            s = bam_get_seq(aln); 

            for(uint32_t i = liftover_read_s; i < liftover_read_e; i++)
            {
                read->seq[i] = bam_nt16_rev_table[bam_seqi(s, i)] //gets nucleotide id and converts them into IUPAC id.
            }
            reads.push_back(read);           
        }
    }
    
    bam_destroy1(aln);
    hts_itr_destroy(itr);
    bam_hdr_destroy(bamHdr);
    bam_close(fp_in);   

    return &read;
}

void VcfWriteHeader(ostream& out, const VcfWriter & vcfWriter)
{
    vcfWriter.writeHeader(out);
    return;
}


void VCFWriteBody(const vector<VNTR> &vntrs, const VcfWriter & vcfWriter, ostream& out)
{
    vcfWriter.writeBody(vntrs, out);
    return;
}

#endif

