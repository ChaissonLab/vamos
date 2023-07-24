#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "option.h"
#include <math.h>
using namespace std;


/* deprecated
static int scoring(vector<MOTIF> &motifs, vector<int> &matches, int top, double globalError, vector<int> &mismatchCI)
{

    if (matches[top] == 0) return 0;

    int M = motifs.size(), M1 = motifs[top].len;
    int mismatches[M];
    double p = 0, numer = 0, denom = 0, score = 0;

    for (auto m=0; m<M; m++)
    {
        mismatches[m] = motifs[m].len - matches[m];
    }

    if (mismatches[top] > mismatchCI[motifs[top].len])
        return 0.5;

    p = matches[top] / double(M1);
    if (p > 1 - globalError)
        p = 1 - globalError;

    numer = pow(p, matches[top]) * pow((1-p), (mismatches[top]));

    for (int m=0; m<M; m++)
    {
        denom += pow(p, matches[m]) * pow((1-p), (mismatches[m]));
    }
    score = numer / denom;

    return -10*log10(1-score);
}
*/


void string_decomposer(vector<uint8_t> &optMotifs, vector<int> &motifQV, vector<int> &starts, vector<int> &ends,
    vector<vector<int> > &motifNMatches, vector<MOTIF> &motifs, char *vntr, int N, const OPTION &opt,
    SDTables &sdTables, vector<int > &mismatchCI)
{

    string vntrS(vntr, N);
    sdTables.Init(motifs, vntrS);
    int left=0;
    int down=1;
    int diag=2;
    int nInf=-99999999;

    // First column is a gap. First row is handled by a separate vector
    for (auto m=0; m < motifs.size(); m++)
    {
        int curGap=-opt.penalty_indel;
        for (auto mi=0; mi < motifs[m].len; mi++)
        {
            sdTables.scoreMat[m][0][mi] = curGap;
            sdTables.pathMat[m][0][mi]  = down;
            curGap -= opt.penalty_indel;
        }
    }


    for (auto s=0; s < N; s++)
    {
        for (auto m=0; m < motifs.size(); m++)
        {

            // Current column in matrices is s+1, string s.
            // Current row in block m is mi, string mi
            for (auto mi=0; mi < motifs[m].len; mi++)
            {
                int diagScore, insScore;
                if (mi == 0)
                {
                    int ms=0;
                    if (motifs[m].seq[0] == vntr[s]) { ms=0;} else { ms =-opt.penalty_mismatch;}
                    diagScore= sdTables.scoreRow0[s] + ms;
                    insScore = sdTables.scoreRow0[s] - opt.penalty_indel;
                }
                else
                {
                    int ms=0;
                    if (motifs[m].seq[mi] == vntr[s]) { ms=0;} else { ms =-opt.penalty_mismatch;}    
                    diagScore= sdTables.scoreMat[m][s][mi-1] + ms;
                    insScore = sdTables.scoreMat[m][s+1][mi-1] - opt.penalty_indel;
                }

                int delScore = sdTables.scoreMat[m][s][mi] - opt.penalty_indel; //prev index =s, cur=s+1;
                int optScore=0;
                int optPath=0;
                if (diagScore >= insScore and diagScore >= delScore)
                {
                    optScore=diagScore;
                    optPath=diag;
                }
                else if (insScore >= diagScore and insScore >= delScore)
                {
                    optScore=insScore;
                    optPath=down;
                }
                else
                {
                    optScore=delScore;
                    optPath=left;
                }

                sdTables.scoreMat[m][s+1][mi] = optScore;
                sdTables.pathMat[m][s+1][mi]  = optPath;

                if (optScore == diagScore and motifs[m].seq[mi] == vntr[s] )
                {
                    if (mi > 0) {
                        sdTables.nMatchMat[m][s+1][mi] = sdTables.nMatchMat[m][s][mi-1] + 1;
                    }
                    else {
                        sdTables.nMatchMat[m][s+1][mi] = 1;
                    }
                    sdTables.nDelMat[m][s+1][mi] = sdTables.nDelMat[m][s][mi];
                    // cout << "nm " << m << " " << s+1<< " " << mi << " " <<  sdTables.nMatchMat[m][s+1][mi] << endl;
                }
                else if (optScore == insScore)
                {
                    if (mi > 0)
                    {
                        sdTables.nMatchMat[m][s+1][mi] = sdTables.nMatchMat[m][s+1][mi-1];
                        sdTables.nDelMat[m][s+1][mi] = sdTables.nDelMat[m][s+1][mi-1];
                    }
                }
                else
                {
                    sdTables.nMatchMat[m][s+1][mi] = sdTables.nMatchMat[m][s][mi];
                    sdTables.nDelMat[m][s+1][mi] = sdTables.nDelMat[m][s][mi]+1;
                }
            }
        }

        // score loop back.

        int maxScore = sdTables.scoreMat[0][s+1][motifs[0].len-1];
        int maxIndex = 0;
        for (auto m=1; m < motifs.size(); m++)
        {
            if (sdTables.scoreMat[m][s+1][motifs[m].len-1] > maxScore )
            {
                maxScore = sdTables.scoreMat[m][s+1][motifs[m].len-1];
                maxIndex = m;
            }
        }
        sdTables.scoreRow0[s+1] = maxScore;
        sdTables.pathRow0[s+1]  = maxIndex;
    }

    // Trace back;
    int maxScore;
    int maxIndex;
    for (auto m=0; m < motifs.size(); m++)
    {
        if (m == 0 or sdTables.scoreMat[m][N][motifs[m].len-1] > maxScore )
        {
            maxScore = sdTables.scoreMat[m][N][motifs[m].len-1];
            maxIndex = m;
        }
    }

    int cm, cs, cmi;
    cmi=motifs[maxIndex].len-1;
    cs=N;
    cm=maxIndex;
    ends.push_back(cs);
    optMotifs.push_back((uint8_t)cm);

    /* remove scoring
    motifNMatches.push_back(vector<int>());
    for (auto tm=0; tm < motifs.size(); tm++)
    {
        motifNMatches[motifNMatches.size()-1].push_back(
            max(0, sdTables.nMatchMat[tm][cs][motifs[tm].len-1] - sdTables.nDelMat[tm][cs][motifs[tm].len-1]));
    }*/

    while (cs > 0)
    {
        assert(optMotifs.size() < 100000);
        int arrow=sdTables.pathMat[cm][cs][cmi];

        // At the beginning of a motif, may need to wrap around
        if (cmi == 0)
        {
            if (arrow == left)
            {
                assert(cs > 0);
                cs--;
            }
            else
            {
                cm = sdTables.pathRow0[cs-1];
                cmi= motifs[cm].len-1;
                if (cs > 1)
                {
                    optMotifs.push_back((uint8_t)cm);    
                    ends.push_back(cs-1);
                    starts.push_back(cs-1);

                    /* remove scoring
                    motifNMatches.push_back(vector<int>());
                    for (auto tm=0; tm < motifs.size(); tm++)
                    {
                        motifNMatches[motifNMatches.size()-1].push_back(max(0, 
                            sdTables.nMatchMat[tm][cs][motifs[tm].len-1] - sdTables.nDelMat[tm][cs][motifs[tm].len-1]));
                    }*/
                }
                cs--;
            }
        }
        else
        {
            if (arrow == diag)
            {
                cs--;
                cmi--;
            }
            else if (arrow == down)
            {
                cmi--;
            }
            else
            {
                cs--;
            }
        }
        assert(cmi >= 0);
        assert(cm >=0);
    }

    starts.push_back(cs);
    assert(starts.size() == ends.size());
    reverse(optMotifs.begin(), optMotifs.end());
    reverse(starts.begin(), starts.end());
    reverse(ends.begin(), ends.end());
    reverse(motifNMatches.begin(), motifNMatches.end());

    // for (auto m=0; m < optMotifs.size(); m++)
    // {
    //     int qv=scoring(motifs, motifNMatches[m], optMotifs[m], 1-opt.accuracy, mismatchCI);
    //     motifQV.push_back(qv);
    // }
}

