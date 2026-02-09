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

#include <stdio.h>
void MSA::runConsensus()
{

  auto alignment_engine = spoa::AlignmentEngine::Create(
    spoa::AlignmentType::kSW,  // global alignment
    3,                         // match score
    -3,                        // mismatch penalty
    -3,                        // gap open
    -3                         // gap extension
);
  spoa::Graph graph {};
  if (seqs.size()  == 0) {
    return;
  }
  vector<int> seqLengths;

  for (const auto &seq : seqs) {
    int s=seq.size();
    seqLengths.push_back(s);
  }
  sort(seqLengths.begin(), seqLengths.end());
  int medianLength=seqLengths[seqLengths.size()/2];
  int minLength=seqLengths[0];
  int maxLength=seqLengths[seqLengths.size()-1];
  int l=seqLengths.size()/2;
  int h=l;
  int m=seqLengths.size()/2;
  /*
  while (l> 0 and seqLengths[m] == seqLengths[l-1]) { l--;}
  while (h+1 < seqLengths.size() and seqLengths[m] == seqLengths[h+1]) { h++;}
  if (h - l < 6) {
    bool expanded=true;
    while(expanded) {
      if (
  }
  */
  for (const auto &seq : seqs) {
    // Align seq to the graph
    if (seq.size() == medianLength or (seq.size() != minLength and seq.size() != maxLength)) {
      auto alignment = alignment_engine->Align(seq, graph);
      
      // Add alignment to the graph
      //    cout << "MSA adding " << seq << endl;
      graph.AddAlignment(alignment, seq);
    }
  }
  //  auto msa = graph.generate_msa();
  int minCoverage=min(4,(int)seqs.size());
  vector<uint32_t> summary;
  consensus = graph.GenerateConsensus(minCoverage);
  /*
  cout << "Summary" << endl;
  for (auto c: consensus) {
    cout << " " << c;
  }
  cout << endl;
  auto msa = graph.GenerateMultipleSequenceAlignment(true);
  for (std::uint32_t i = 0; i < msa.size(); ++i) {
    std::cout << ">" << i << std::endl
	      << msa[i] << std::endl;
  }
    
  cout << "cons" << endl << consensus << endl;
  */
  /*
    int i, j;

    // alignment parameters

    abpt->align_mode = 0;               // 0:global, 1:local, 2:extension
    abpt->match = 1;                    // match score
    abpt->mismatch = -2;                 // mismatch penalty
    abpt->gap_mode = ABPOA_AFFINE_GAP;  // gap penalty mode
    abpt->gap_open1 = 2;                // gap open penalty #1
    abpt->gap_ext1 = 1;                 // gap extension penalty #1
    abpt->gap_open2 = 2;                // gap open penalty #2
    abpt->gap_ext2 = 1;                 // gap extension penalty #2
    // gap_penalty = min{gap_open1 + gap_len * gap_ext1, gap_open2 + gap_len * gap_ext2}
    // abpt->is_diploid = 0;

    abpt->out_msa = 1;                  // generate Row-Column multiple sequence alignment(RC-MSA), set 0 to disable
    abpt->out_cons = 1;                 // generate consensus sequence, set 0 to disable
    abpt->use_qv = 0;
    abpt->progressive_poa = 0;
    abpt->ret_cigar = 0;
    abpt->rev_cigar = 0;

    abpt->out_pog = 0;
    abpt->out_gfa = 0;
    // abpt->out_msa_header = 0;
    abpt->use_read_ids = 0;

    abpoa_post_set_para(abpt);
    abpoa_msa(ab, abpt, n_seqs, NULL, seq_lens, bseqs, NULL, NULL);

    abpoa_output(ab, abpt, stdout);    

    

    cons_seq = ab->abc->cons_base;
    cons_cov = ab->abc->cons_cov;
    cons_l   = ab->abc->cons_len;
    cons_n   = ab->abc->n_cons;
    msa_seq  = ab->abc->msa_base;
    msa_l    = ab->abc->msa_len;    
    // abpoa_msa(ab, abpt, n_seqs, NULL, seq_lens, bseqs, NULL, &cons_seq, &cons_cov, NULL);
    */
    return;
}


void MSA::extractConsensusAnno(vector<uint8_t> &consensus, int motifs_size)
{
    return;
}


void MSA::extractConsensusRawSeq(READ * Hap_seq)
{
    Hap_seq->len=consensus.size();
    Hap_seq->seq = consensus;
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

    if (seqs.size() == 1)
    {
        Hap_seq->seq = reads[gp[0]]->seq;
        Hap_seq->len = reads[gp[0]]->len;	
        return;
    }

    // encodeSeqsFromRaw(gp, reads);
    if (seqs.size() > 0) {
      runConsensus ();
      extractConsensusRawSeq (Hap_seq);
    }
    else {
      Hap_seq->seq = "";
      Hap_seq->len = 0;
    }
    return;
}


void MSA::cleanMSA()
{
    return;
}

