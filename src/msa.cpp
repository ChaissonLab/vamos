#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <algorithm>
#include <utility>
#include <numeric>
#include <iostream>
#include <fstream>
#include "msa.h"
#include "read.h"
#include "abpoa.h"
#include "option.h"


// 65,97=>A, 67,99=>C, 71,103=>G, 84,85,116,117=>T, else=>N
const char nt256_table[256] = {
       'A', 'C', 'G', 'T',  'N', '-', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', '-',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'A', 'N', 'C',  'N', 'N', 'N', 'G',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'T', 'T', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'A', 'N', 'C',  'N', 'N', 'N', 'G',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'T', 'T', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N'
};


void MSA::runConsensus()
{
    int i, j;

    // alignment parameters

    abpt->align_mode = 2;               // 0:global, 1:local, 2:extension
    abpt->match = 1;                    // match score
    abpt->mismatch = 1;                 // mismatch penalty
    abpt->gap_mode = ABPOA_AFFINE_GAP;  // gap penalty mode
    abpt->gap_open1 = 1;                // gap open penalty #1
    abpt->gap_ext1 = 1;                 // gap extension penalty #1
    abpt->gap_open2 = 1;                // gap open penalty #2
    abpt->gap_ext2 = 1;                 // gap extension penalty #2
    // gap_penalty = min{gap_open1 + gap_len * gap_ext1, gap_open2 + gap_len * gap_ext2}
    // abpt->is_diploid = 0;

    abpt->out_msa = 0;                  // generate Row-Column multiple sequence alignment(RC-MSA), set 0 to disable
    abpt->out_cons = 1;                 // generate consensus sequence, set 0 to disable
    abpt->use_qv = 0;
    abpt->progressive_poa = 1;
    abpt->ret_cigar = 0;
    abpt->rev_cigar = 0;

    abpt->out_pog = 0;
    abpt->out_gfa = 0;
    // abpt->out_msa_header = 0;
    abpt->use_read_ids = 0;

    abpoa_post_set_para(abpt);
    abpoa_msa(ab, abpt, n_seqs, NULL, seq_lens, bseqs, NULL, NULL);
    
    cons_seq = ab->abc->cons_base;
    cons_cov = ab->abc->cons_cov;
    cons_l   = ab->abc->cons_len;
    cons_n   = ab->abc->n_cons;
    msa_seq  = ab->abc->msa_base;
    msa_l    = ab->abc->msa_len;    
    // abpoa_msa(ab, abpt, n_seqs, NULL, seq_lens, bseqs, NULL, &cons_seq, &cons_cov, NULL);

    return;
}


void MSA::extractConsensusAnno(vector<uint8_t> &consensus, int motifs_size)
{
    assert((cons_n > 0) && "No consensus sequence exists!"); // cons_n == 0 means no consensus sequence exists
    for (int j = 0; j < cons_l[0]; ++j)
    {
        assert(cons_seq[0][j] <= motifs_size);
        consensus.push_back((cons_seq[0][j] - 1)); // -1 to cancel out +1
    }

    assert(consensus.size() > 0);

    return;
}


void MSA::extractConsensusRawSeq(READ * Hap_seq)
{
    assert((cons_n > 0) && "No consensus sequence exists!"); // cons_n == 0 means no consensus sequence exists
    Hap_seq->len = cons_l[0]; // read length
    assert(Hap_seq->len < 1000000);
    //    Hap_seq->seq = (char *) malloc(Hap_seq->len + 1); // read sequence array
    Hap_seq->seq.resize(Hap_seq->len);
    for(int i = 0; i < Hap_seq->len; i++)
        Hap_seq->seq[i] = nt256_table[cons_seq[0][i]]; // get nucleotide id and convert them into IUPAC id.

    //    Hap_seq->seq[Hap_seq->len] = '\0';
}


void MSA::MSA_anno_group(int motifs_size, const vector<vector<uint8_t> > &annos, vector<uint8_t> &consensus)
{
    if (annos.empty()) return;

    if (annos.size() == 1)
    {
        consensus = annos[0];
        return;
    }
    
    // encodeSeqsFromAnno(gp, annos, n_seqs, motifs_size);
    runConsensus ();
    extractConsensusAnno (consensus, motifs_size);
    return;
}


void MSA::MSA_seq_group(const vector<int> &gp, const vector<READ *> &reads, READ * Hap_seq)
{
    if (gp.empty()) return;

    if (n_seqs == 1)
    {
        Hap_seq->seq = reads[gp[0]]->seq;
        Hap_seq->len = reads[gp[0]]->len;	
        return;
    }

    // encodeSeqsFromRaw(gp, reads);
    if (n_seqs > 0) {
      runConsensus ();
      extractConsensusRawSeq (Hap_seq);
    }
    else if (n_empty_seqs  > 0) {
      Hap_seq->seq = "";
      Hap_seq->len = 0;
    }
    return;
}


void MSA::cleanMSA()
{
    int i;
    if (n_seqs > 0 and bseqs != NULL) {
      for(i = 0; i < n_seqs; ++i) delete[] bseqs[i];    
      delete [] bseqs;
      delete [] seq_lens;
    }
    bseqs = NULL;
    seq_lens = NULL;
    if (ab) {
      abpoa_free(ab);
      ab=NULL;
    }
    if (abpt) {
      abpoa_free_para(abpt);
      abpt=NULL;
    }

    return;
}

