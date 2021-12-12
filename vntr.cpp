#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <algorithm>
#include "htslib/sam.h"
#include "vntr.h"
#include "edlib/include/edlib.h"

// function to compute the S_i scores (the naive occurrence)
void anno(vector<int> &optMotifs, const vector<string> &motifs, const string &vntr) {

    int vntr_len = vntr.size(), motif_len = motifs.size();
    double max;
    vector<double> score(vntr_len + 1, 0); /* score[i]: the score of the best annotation of [vntr[0], vntr[i - 1]], score[0] means no sequnece, thus 0 */
    vector<int> traceI(vntr_len + 1, 0); /* traceI[i]: index of best j */
    vector<int> traceM(vntr_len + 1, 0); /* traceM[i]: index of best motif */

    traceI[0] = -1;
    traceM[0] = -1;

    // propagate
    int i, j, m, k;
    for (i = 1; i < vntr_len; i++) 
    {
        double best_score = 0;
        double cur_score;

        for (j = 0; j < i; j++) 
        {
            for (m = 0; m < motif_len; m++) 
            {
                char * motif = new char[motifs[m].length() + 1];
                char * vntr_seq = new char[i - j + 1];
                strcpy(motif, motifs[m].c_str());
                strcpy(vntr_seq, vntr.substr(j, i - j).c_str());
                /*the score of [vntr[0], vntr[j - 1]] + the alignment score of substring [j, i - 1]*/
                EdlibAlignResult result = edlibAlign(motif, motifs[m].length(), vntr_seq, i - j,
                                     edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_DISTANCE, NULL, 0));
                cur_score = score[j] + (double) result.editDistance;
                if (cur_score > best_score) 
                {
                    best_score = cur_score;
                    traceI[i] = j;
                    traceM[i] = m;
                }
                delete [] motif;
                delete [] vntr_seq;
            }
        }
        score[i] = best_score;
    }

    // traceback
    int cur_idx = vntr_len;
    while (cur_idx > 0)
    {
        optMotifs.push_back(traceM[cur_idx]);
        cur_idx = traceI[cur_idx];
    }
    reverse(optMotifs.begin(), optMotifs.end());
    return;
}

void VNTR::motifRepresentationForOneVNTR () 
{
	/* modify vector<vector<int>> reps
	   for now, just using dummy code
	*/
	reps.resize(reads.size());
	for (int i = 0; i < reads.size(); ++i)
	{
		/* 
			apply Bida's code here
			input: string : reads[i].seq
				  vector<MOTIF> : motifs
			output: vector<vector<int>> reps[i]
		*/
		anno(reps[i], motifs, reads[i]->seq);
	}
	return;
}

void VNTR::concensusMotifRepForOneVNTR ()
{
	/* based on vector<int> rep, get a concensus representation of the current locus 
	   for now, just using some dummy code 
	   input: vector<vector<int>> reps
	   output: vector<int> concensus_h1
	           vector<int> concensus_h2
	*/
	concensus_h1.resize(10, 1);
	concensus_h2.resize(10, 1);
	return;
}

void VNTR::commaSeparatedMotifRepForConsensus (bool h1, string &motif_rep)
{
	if (h1)
	{
		for (auto &it : concensus_h1) { motif_rep += "VNTR_" + to_string(it) + ",";}
	}
	else
	{
		for (auto &it : concensus_h2) { motif_rep += "VNTR_" + to_string(it) + ",";}
	}
	if (!motif_rep.empty()) motif_rep.pop_back();
}
