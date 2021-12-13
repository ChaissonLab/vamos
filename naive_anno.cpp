#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <algorithm>
#include <iostream>
#include "edlib/include/edlib.h"

// function to compute the S_i scores (the naive occurrence)
int anno(vector<int> &optMotifs, vector<MOTIF> &motifs, char * vntr) {

    /* get the maximum length of motifs */
    int max_len = 0;
    for (auto & motif : motifs)
    {
        max_len = max(max_len, motif.len);
    }

    int vntr_len = strlen(vntr), motif_len = motifs.size();
    vector<double> score(vntr_len + 1, 0); /* score[i]: the score of the best annotation of [vntr[0], vntr[i - 1]], score[0] means no sequnece, thus 0 */
    vector<int> traceI(vntr_len + 1, 0); /* traceI[i]: index of best j */
    vector<int> traceM(vntr_len + 1, 0); /* traceM[i]: index of best motif */

    traceI[0] = -1;
    traceM[0] = -1;

    // propagate
    int i, j, m, k;
    double best_score, cur_score;
    for (i = 1; i <= vntr_len; i++) 
    {
        best_score = 10000000.0;

        for (j = 0; j < i; j++) 
        {
            // char * vntr_subseq = (char *) malloc(i - j + 1); 
            // memcpy(vntr_subseq, vntr + j, i - j);
            // vntr_subseq[i - j] = 0;
            for (m = 0; m < motif_len; m++) 
            {
                /*the score of [vntr[0], vntr[j - 1]] + the alignment score of substring [j, i - 1]*/
                EdlibAlignResult result = edlibAlign(motifs[m].seq.c_str(), motifs[m].len, &vntr[j], i - j, edlibDefaultAlignConfig());
                // EdlibAlignResult result = edlibAlign(motifs[m].seq.c_str(), motifs[m].len, vntr_subseq, i - j, edlibDefaultAlignConfig());

                if (result.status == EDLIB_STATUS_OK) 
                {
                    cur_score = score[j] + (double) result.editDistance;

                    if (cur_score < best_score) 
                    {
                        best_score = cur_score;
                        traceI[i] = j;
                        traceM[i] = m;
                    }
                }
                else 
                {
                    cerr << "edlib goes wrong!" << endl;
                    return 1;
                }
            }
            // free(vntr_subseq);

        }
        // if (i % 1000 == 0)
        //     cerr << "i: " << i << " score: " << best_score << " traceI[i]: " << traceI[i] << " traceM[i]: " << traceM[i] <<  endl;
        score[i] = best_score;
    }

    // traceback
    int cur_idx = vntr_len;
    while (cur_idx >= 0)
    {
        if (traceM[cur_idx] >= 0)
            optMotifs.push_back(traceM[cur_idx]);
        cur_idx = traceI[cur_idx];
    }

    reverse(optMotifs.begin(), optMotifs.end());
    return 0;
}