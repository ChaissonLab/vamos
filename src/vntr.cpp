#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <algorithm>
#include <utility>
#include "vntr.h"
#include "bounded_anno.cpp"
#include "naive_anno.cpp"
#include "edlib.h"
#include "dataanalysis.h"
#include "abpoa.h"
#include "option.h"

// #include "seqan/sequence.h"
// #include "seqan/align.h"
// #include "seqan/score.h"
void VNTR::motifAnnoForOneVNTR (const OPTION &opt) 
{
	if (nreads == 0) return;

	annos.resize(nreads);
	for (int i = 0; i < nreads; ++i)
	{
		/* 
			apply Bida's code here
			input:   string : reads[i].seq; vector<MOTIF> : motifs
			output:  vector<vector<int>> annos[i]
		*/
		if (opt.fasterAnnoAlg)
			bounded_anno(annos[i], motifs, reads[i]->seq, reads[i]->len);
		else
			naive_anno(annos[i], motifs, reads[i]->seq, reads[i]->len);
		// cerr << "finish for reads: " << i << endl;
	}
	return;
}

// skip "-" (45) for encoding
static void encode (int motifs_size, vector<uint8_t> &annoNum, string &annoStr)
{
	string tmp_s;
	uint8_t num;
	for (size_t i = 0; i < annoNum.size(); ++i)
	{
		num = annoNum[i];
		assert(num < motifs_size);

		// tmp_s.assign(1, (char) num);
		// annoStr.replace(i, 1, tmp_s);  

		// for outputing string to stdout
		if (num < 45)
		{
			tmp_s.assign(1, (char) (num + 33));
			annoStr.replace(i, 1, tmp_s);  
		}
		else if (num >= 45 and num <= 221)
		{
			tmp_s.assign(1, (char) (num + 34));
			annoStr.replace(i, 1, tmp_s);  		
		}		
	}
	return;
}

void annoTostring_helper (int motifs_size, string &annoStr, vector<uint8_t> &annoNum)
{
	annoStr.clear();
	annoStr.resize(annoNum.size(), '+');
	encode(motifs_size, annoNum, annoStr);
	return;
}

void VNTR::annoTostring (const OPTION &opt)
{
	if (nreads == 0) return;
	annoStrs.resize(nreads);
	for (int r = 0; r < nreads; ++r)
	{
		annoTostring_helper(motifs.size(), annoStrs[r], annos[r]);
		if (opt.debug) cerr << "string " << r << ":  " << annoStrs[r] << endl;
	}
	return;
}

size_t VNTR::getAnnoStrLen (int i)
{
	assert(i < nreads);
	return annoStrs[i].length();
}

int VNTR::hClust (vector<int> &gp1, vector<int> &gp2, double edist [])
{
	/*
	 compute the edit distance matrix for annotation strings 
	 input: vector<string> annoStrs
	 output: double edist[nreads * nreads] = {...}
	*/
    gp1.clear(); gp2.clear();
	if (nreads == 1)
	{
	   gp1.push_back(0);
	   return 0;
	}

	int i, j;
	EdlibAlignResult result;
	for (i = 0; i < nreads; ++i)
	{
		for (j = i + 1; j < nreads; ++j) 
		{
			result = edlibAlign(annoStrs[i].c_str(), getAnnoStrLen(i), annoStrs[j].c_str(), getAnnoStrLen(j), edlibDefaultAlignConfig());
            if (result.status == EDLIB_STATUS_OK) 
            {
                edist[i * nreads + j] = (double) result.editDistance; // [i][j]
                edist[j * nreads + i] = (double) result.editDistance; // [j][i]
            }
            else 
            {
                printf("edlib goes wrong!\n");
                return 1;
            }
		}
	}

	/* do hierarchical clustering based on the edit ditance matrix */
    alglib::clusterizerstate s;
    alglib::ahcreport rep;
    alglib::clusterizercreate(s);

    alglib::real_2d_array d;
    d.setcontent(nreads, nreads, edist);
    alglib::clusterizersetdistances(s, d, true);
    alglib::clusterizerrunahc(s, rep);

	/* get the cluster results when nclusts == 2 */
    alglib::integer_1d_array cidx;
    alglib::integer_1d_array cz;    
    clusterizergetkclusters(rep, 2, cidx, cz);

    gp1.clear(); gp2.clear();
    for (i = 0; i < nreads; ++i)
    {
    	if (cidx[i] == 0)
    		gp1.push_back(i);
    	else
    		gp2.push_back(i);
    }
    return 0;
}

void outputConsensus (vector<uint8_t> &consensus, const OPTION &opt)
{
	size_t i = 0;
    for (auto &it : consensus)
    {
    	if (opt.debug)
    	{
			if (i < consensus.size() - 1)
				cerr << to_string(it) << "->";
			else
				cerr << to_string(it);    		
    	}
		i += 1;
    }
	if (opt.debug) cerr << endl;
	return;
}

/* use abPOA to do MSA */
void MSA_helper (int motifs_size, int n_seqs, vector<uint8_t> &consensus, int *seq_lens, uint8_t **bseqs)
{
	int i, j;
    // initialize variables
    abpoa_t *ab = abpoa_init();
    abpoa_para_t *abpt = abpoa_init_para();

    // alignment parameters
    abpt->align_mode = 0; // 0:global 1:local, 2:extension
    abpt->match = 1;      // match score
    abpt->mismatch = 1;   // mismatch penalty
    abpt->gap_mode = ABPOA_AFFINE_GAP; // gap penalty mode
    abpt->gap_open1 = 1;  // gap open penalty #1
    abpt->gap_ext1 = 1;   // gap extension penalty #1
    abpt->gap_open2 = 1; // gap open penalty #2
    abpt->gap_ext2 = 1;   // gap extension penalty #2
                    	  // gap_penalty = min{gap_open1 + gap_len * gap_ext1, gap_open2 + gap_len * gap_ext2}
    // abpt->bw = 10;        // extra band used in adaptive banded DP
    // abpt->bf = 0.01; 

    abpt->is_diploid = 0;
	// abpt->min_freq = 0.8; 
    abpt->out_msa = 1; // generate Row-Column multiple sequence alignment(RC-MSA), set 0 to disable
    abpt->out_cons = 1; // generate consensus sequence, set 0 to disable
    // abpt->w = 6, abpt->k = 2; abpt->min_w = 10; // minimizer-based seeding and partition
    abpt->progressive_poa = 1;

    // variables to store result
    uint8_t **cons_seq; int **cons_cov, *cons_l, cons_n = 0;
    uint8_t **msa_seq; int msa_l = 0;

    abpoa_post_set_para(abpt);

    // 1. output to stdout
    // fprintf(stdout, "=== output to stdout ===\n");

    // perform abpoa-msa
    // abpoa_msa(ab, abpt, n_seqs, NULL, seq_lens, bseqs, stdout, NULL, NULL, NULL, NULL, NULL, NULL);

    // ab->abs->n_seq = 0; // To re-use ab, n_seq needs to be set as 0
    abpoa_msa(ab, abpt, n_seqs, NULL, seq_lens, bseqs, NULL, &cons_seq, &cons_cov, &cons_l, &cons_n, &msa_seq, &msa_l);

    // assert(cons_n > 0); // cons_n == 0 means no consensus sequence exists
    for (j = 0; j < cons_l[0]; ++j)
    {
    	assert(cons_seq[0][j] < motifs_size);
    	consensus.push_back((uint8_t)cons_seq[0][j]);
    }
    assert(consensus.size() > 0);

    // fprintf(stdout, "=== output to variables ===\n");
    // for (i = 0; i < cons_n; ++i) {
    //     fprintf(stdout, ">Consensus_sequence\n");
    //     for (j = 0; j < cons_l[i]; ++j)
    //         fprintf(stdout, "%i", cons_seq[i][j]);
    //     fprintf(stdout, "\n");
    // }
    // fprintf(stdout, ">Multiple_sequence_alignment\n");
    // for (i = 0; i < n_seqs; ++i) {
    //     for (j = 0; j < msa_l; ++j) {
    //         fprintf(stdout, "%c", "ACGTN-"[msa_seq[i][j]]);
    //     }
    //     fprintf(stdout, "\n");
    // }

    if (cons_n) {
        for (i = 0; i < cons_n; ++i) 
        {
            free(cons_seq[i]); 
            free(cons_cov[i]);
        } 
        free(cons_seq); 
        free(cons_cov); 
        free(cons_l);
    }

    if (msa_l) 
    {
        for (i = 0; i < n_seqs; ++i) 
        {
        	free(msa_seq[i]); 
        	free(msa_seq);
        }
    }

    abpoa_free(ab); 
    abpoa_free_para(abpt); 

    return;
}

void MSA (int motifs_size, const vector<int> &gp, const vector<vector<uint8_t>> &annos, vector<uint8_t> &consensus)
{
	if (gp.empty()) return;

    int n_seqs = gp.size();
	if (n_seqs == 1)
	{
		consensus = annos[gp[0]];
		return;
	}

    // collect sequence length
    int *seq_lens = (int*)malloc(sizeof(int) * n_seqs);
    uint8_t **bseqs = (uint8_t**)malloc(sizeof(uint8_t*) * n_seqs);
    int i, j;
    for (i = 0; i < n_seqs; ++i) {
        seq_lens[i] = annos[gp[i]].size();
        bseqs[i] = (uint8_t*)malloc(sizeof(uint8_t) * seq_lens[i]);
        for (j = 0; j < seq_lens[i]; ++j)
        {
        	assert(annos[gp[i]][j] < motifs_size);
        	bseqs[i][j] = annos[gp[i]][j];

        	// if (annos[gp[i]][j] >= 0 and annos[gp[i]][j] < 45) 
        	// {
        	// 	bseqs[i][j] = annos[gp[i]][j] + 33;
        	// }
        	// else 
        	// {
        	// 	bseqs[i][j] = annos[gp[i]][j] + 34;
        	// }
            
        }
    }

	MSA_helper (motifs_size, n_seqs, consensus, seq_lens, bseqs);

    for (i = 0; i < n_seqs; ++i) free(bseqs[i]); 
   	free(bseqs); 
    free(seq_lens);
    return;
}

void MSA (int motifs_size, const vector<vector<uint8_t>> &annos, vector<uint8_t> &consensus)
{
    int n_seqs = annos.size();
	if (n_seqs == 1)
	{
		consensus = annos[0];
		return;
	}

    // collect sequence length
    int *seq_lens = (int*)malloc(sizeof(int) * n_seqs);
    uint8_t **bseqs = (uint8_t**)malloc(sizeof(uint8_t*) * n_seqs);
    int i, j;
    for (i = 0; i < n_seqs; ++i) {
        seq_lens[i] = annos[i].size();
        bseqs[i] = (uint8_t*)malloc(sizeof(uint8_t) * seq_lens[i]);
        for (j = 0; j < seq_lens[i]; ++j) 
        {
        	assert(annos[i][j] < motifs_size); 
        	bseqs[i][j] = annos[i][j];
        }
    }

	MSA_helper (motifs_size, n_seqs, consensus, seq_lens, bseqs);
    for (i = 0; i < n_seqs; ++i) free(bseqs[i]); 
   	free(bseqs); 
    free(seq_lens);
    return;
}

double avgSilhouetteCoeff (int nreads, double edist [], vector<int> &gp1, vector<int> &gp2)
{
	double a1 = 0.0;
	double b1 = 10000000.0;
	double sc1 = 0.0;

	double a2 = 0.0;
	double b2 = 10000000.0;
	double sc2 = 0.0;

	for (auto &i : gp1)
	{
		a1 = 0.0; b1 = 10000000.0;
		for (auto &j : gp1) 
		{
			if (i == j) continue;
			a1 += edist[i * nreads + j];
		}

		if (gp1.size() > 1)
			a1 /= gp1.size() - 1;

		for (auto &j : gp2)
			b1 = min(edist[i * nreads + j], b1);

		sc1 += (b1 - a1) / max(a1, b1);
	}

	for (auto &i : gp2)
	{
		a2 = 0.0; b2 = 10000000.0;
		for (auto &j : gp2) 
		{
			if (i == j) continue;
			a2 += edist[i * nreads + j];
		}

		if (gp2.size() > 1)
			a2 /= gp2.size() - 1;

		for (auto &j : gp1)
			b2 = min(edist[i * nreads + j], b2);

		sc2 += (b2 - a2) / max(a2, b2);
	}

	return (sc1 + sc2) / (gp1.size() + gp2.size());
}

/* based on vector<vector<uint8_t>> annos, get a consensus representation of the current locus 
   input:   vector<vector<uint8_t>> annos
   output:  vector<vector<uint8_t>> consensus
*/
void VNTR::concensusMotifAnnoForOneVNTR (const OPTION &opt)
{
	/* compute the MSA for each cluster when nclusts == 2 */
	vector<int> gp1, gp2;
	double edist [nreads * nreads];
	hClust(gp1, gp2, edist); 

	/* compare which is better: nclusts == 2 or nclusts == 1 */
	double sc_2clusts = avgSilhouetteCoeff(nreads, edist, gp1, gp2);
	if (sc_2clusts > 0.0 and gp1.size() > 0 and gp2.size() > 0)
	{
		het = true;
		consensus.resize(2);
		MSA (motifs.size(), gp1, annos, consensus[0]);
		MSA (motifs.size(), gp2, annos, consensus[1]);
		assert(consensus[0].size() > 0 and consensus[1].size() > 0);
	}
	else 
	{
		het = false;
		consensus.resize(1);
		MSA (motifs.size(), annos, consensus[0]);
		assert(consensus[0].size() > 0);
	}

	if (het)
	{
		if (opt.debug) {cerr << "het! " << endl; cerr << "h1 anno: " << endl;}
		outputConsensus(consensus[0], opt);
		if (opt.debug) cerr << "h2 anno: " << endl;
		outputConsensus(consensus[1], opt);
	}
	else
	{
		if (opt.debug) {cerr << "hom! " << endl; cerr << "h1&h2 anno: " << endl;}
		outputConsensus(consensus[0], opt);		
	}
	return;
}

void VNTR::concensusMotifAnnoForOneVNTRUsingABpoa (const OPTION &opt)
{
    int n_seqs = annos.size();
    int motifs_size = motifs.size();
	if (n_seqs == 1)
	{
		consensus.resize(1);
		consensus[0] = annos[0];
		return;
	}

    // collect sequence length
    int *seq_lens = (int*)malloc(sizeof(int) * n_seqs);
    uint8_t **bseqs = (uint8_t**)malloc(sizeof(uint8_t*) * n_seqs);
    int i, j;
    for (i = 0; i < n_seqs; ++i) {
        seq_lens[i] = annos[i].size();
        bseqs[i] = (uint8_t*)malloc(sizeof(uint8_t) * seq_lens[i]);
        for (j = 0; j < seq_lens[i]; ++j) 
        {
        	assert(annos[i][j] < motifs_size); 
        	bseqs[i][j] = annos[i][j];
        }
    }

    // initialize variables
    abpoa_t *ab = abpoa_init();
    abpoa_para_t *abpt = abpoa_init_para();

    // alignment parameters
    abpt->align_mode = 0; // 0:global 1:local, 2:extension
    abpt->match = 1;      // match score
    abpt->mismatch = 1;   // mismatch penalty
    abpt->gap_mode = ABPOA_AFFINE_GAP; // gap penalty mode
    abpt->gap_open1 = 1;  // gap open penalty #1
    abpt->gap_ext1 = 1;   // gap extension penalty #1
    abpt->gap_open2 = 1; // gap open penalty #2
    abpt->gap_ext2 = 1;   // gap extension penalty #2
                    	  // gap_penalty = min{gap_open1 + gap_len * gap_ext1, gap_open2 + gap_len * gap_ext2}
    // abpt->bw = 10;        // extra band used in adaptive banded DP
    // abpt->bf = 0.01; 

    abpt->is_diploid = 0;
	// abpt->min_freq = 0.8; 
    abpt->out_msa = 1; // generate Row-Column multiple sequence alignment(RC-MSA), set 0 to disable
    abpt->out_cons = 1; // generate consensus sequence, set 0 to disable
    // abpt->w = 6, abpt->k = 2; abpt->min_w = 10; // minimizer-based seeding and partition
    abpt->progressive_poa = 1;

    // variables to store result
    uint8_t **cons_seq; int **cons_cov, *cons_l, cons_n = 0;
    uint8_t **msa_seq; int msa_l = 0;

    abpoa_post_set_para(abpt);

    abpoa_msa(ab, abpt, n_seqs, NULL, seq_lens, bseqs, NULL, &cons_seq, &cons_cov, &cons_l, &cons_n, &msa_seq, &msa_l);

    assert(cons_n > 0); // cons_n == 0 means no consensus sequence exists

    int numHap = cons_n > 1 ? 2 : 1;
	consensus.resize(numHap);
	for (i = 0; i < numHap; ++i)
	{
		for (j = 0; j < cons_l[i]; ++j)
		{
			assert(cons_seq[i][j] < motifs_size);
			consensus[i].push_back((uint8_t)cons_seq[i][j]);
		}
		assert(consensus[i].size() > 0); 		
	}

    if (cons_n) {
        for (i = 0; i < cons_n; ++i) 
        {
            free(cons_seq[i]); 
            free(cons_cov[i]);
        } 
        free(cons_seq); 
        free(cons_cov); 
        free(cons_l);
    }

    if (msa_l) 
    {
        for (i = 0; i < n_seqs; ++i) 
        {
        	free(msa_seq[i]); 
        	free(msa_seq);
        }
    }

    abpoa_free(ab); 
    abpoa_free_para(abpt); 

    for (i = 0; i < n_seqs; ++i) free(bseqs[i]); 
   	free(bseqs); 
    free(seq_lens);

	if (het)
	{
		if (opt.debug) {cerr << "het! " << endl; cerr << "h1 anno: " << endl;}
		outputConsensus(consensus[0], opt);
		if (opt.debug) cerr << "h2 anno: " << endl;
		outputConsensus(consensus[1], opt);
	}
	else
	{
		if (opt.debug) {cerr << "hom! " << endl; cerr << "h1&h2 anno: " << endl;}
		outputConsensus(consensus[0], opt);		
	}

	return;
}

void VNTR::commaSeparatedMotifAnnoForConsensus (bool h1, string &motif_anno)
{
	if (h1)
		for (auto &it : consensus[0]) { motif_anno += "MOTIF_" + to_string(it) + ",";}
	else if (het)
		for (auto &it : consensus[1]) { motif_anno += "MOTIF_" + to_string(it) + ",";}	
	else 
		for (auto &it : consensus[0]) { motif_anno += "MOTIF_" + to_string(it) + ",";}

	if (!motif_anno.empty()) motif_anno.pop_back();
	// assert(!motif_anno.empty());
	return;
}
