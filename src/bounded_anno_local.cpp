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

static int scoring(vector<MOTIF> &motifs, vector<int> &matches, int top, double globalError, vector<int> &mismatchCI) {

    if (matches[top] == 0) return 0;

    int M = motifs.size(), M1 = motifs[top].len;
    int mismatches[M];
    double p = 0, numer = 0, denom = 0, score = 0;

    for (auto m=0; m<M; m++) {
      mismatches[m] = motifs[m].len - matches[m];
    }
    if (mismatches[top] > mismatchCI[motifs[top].len]) {
      return 0.5;
    }
    p = matches[top] / double(M1);
    if (p > 1 - globalError)
        p = 1 - globalError;

    numer = pow(p, matches[top]) * pow((1-p), (mismatches[top]));
//    numer=pow((1-p), (mismatches[top]));
    for (int m=0; m<M; m++) {
        denom += pow(p, matches[m]) * pow((1-p), (mismatches[m]));
//      denom += pow((1-p), mismatches[m]);
//      cout << numer << "\t" << mismatches[m] << "\t" << pow((1-p), mismatches[m]) << "\t" << denom;
//      if (m == top) { cout << " TOP";}
//      cout << endl;
    }
    score = numer / denom;
//    score = numer / denom;

    return -10*log10(1-score);
}


// function to compute the distance of given pair of characters
static double distance(char a, char b, const OPTION &opt) {

    double distance;
    if (a == b) {
        distance = 0; // distance for match
    } else if (a == '-' || b == '-') {
        distance = opt.penalty_indel; // distance for indel
    } else {
        distance = opt.penalty_mismatch; // distance for mismatch
    }
    return distance;
}


// function to find min edit distance in dynamic programming
static uint8_t compare(vector<double> &update, uint8_t l) {

    double min = update[0];
    uint8_t k = 0;

    for (uint8_t i = 1; i < l; i++) {
        if (update[i] <= min) {
            min = update[i];
            k = i;
        }
    }
    return k;
}

// function to compute the variant N-W alignment (linear space, only computes opt distance and traces starting position)
// (leading gaps cost 0)
void global(const string &motif, char *vntr, vector<vector<double>> &dists, vector<vector<int>> &starts, uint8_t m, int N, const OPTION &opt) { 

    int i, j, oldStart, tempStart; // i as row (string A) index, j as column (string B) index, k as distance index for "compare" function
    int M = motif.length();
    double oldDist;
    vector<double> update(3, 0.0);
    uint8_t k;

    // initialization for the first row (dist=0, start=j for first row, 0 penalty for leading gaps)
    for (j=0; j<=N; j++) {
        dists[m][j] = 0;
        starts[m][j] = j;
    }

    // propagate to later rows (dist=i, start=0 for first column)
    for (i = 1; i <= M; i++) {
        oldDist = dists[m][0];
        oldStart = starts[m][0];
        dists[m][0] = i;
        starts[m][0] = 0;

        for (j=1; j<=N; j++) {

            // dynamic programming
            update[0] = dists[m][j-1] + distance(vntr[j-1], '-', opt); // use current value: dist[j-1] after update
            update[1] = dists[m][j] + distance('-', motif[i-1], opt); // use old value: dist[j] before update
            update[2] = oldDist + distance(vntr[j-1], motif[i-1], opt); // use old value: dist[j-1] before update

            k = compare(update, 3);

            if (k == 0) {
                tempStart = starts[m][j-1];
            } else if (k == 2) {
                tempStart = oldStart;
            } else {
                tempStart = starts[m][j];
            }

            oldDist = dists[m][j]; // save old value: dist[j-1] before update
            oldStart = starts[m][j]; // save old value: start[j-1] before update

            // update dist and start
            dists[m][j] = update[k];
            starts[m][j] = tempStart;
        }
    }
}

void string_decomposer(vector<uint8_t> &optMotifs, vector<int> &motifQV, vector<int> &starts, vector<int> &ends, vector<vector<int> > &motifNMatches, vector<MOTIF> &motifs, char *vntr, int N, const OPTION &opt, SDTables &sdTables, vector<int > &mismatchCI) {
    string vntrS(vntr);
    sdTables.Init(motifs, vntrS);
    int left = 0;
    int down = 1;
    int diag = 2;
    int nInf = -99999999;
    //
    // First column is a gap. First row is handled by a separate vector
    //
    for (auto m = 0; m < motifs.size(); m++) 
    {
        int curGap = -opt.penalty_indel;
        for (auto mi = 0; mi < motifs[m].len; mi++) 
        {
          sdTables.scoreMat[m][0][mi] = curGap;
          sdTables.pathMat[m][0][mi] = down;
          curGap -= opt.penalty_indel;
        }
    }
      
    for (auto s = 0; s < N; s++) 
    {
        for (auto m = 0; m < motifs.size(); m++) 
        {
            //
            // Current column in matrices is s+1, string s.
            // Current row in block m is mi, string mi
            // 
            for (auto mi = 0; mi < motifs[m].len; mi++) 
            {
                int diagScore, insScore;
                if (mi == 0) 
                {
                    int ms = 0;
                    if (motifs[m].seq[0] == vntr[s]) {ms = 0;} else { ms = -opt.penalty_mismatch;}
                    diagScore = sdTables.scoreRow0[s] + ms;
                    insScore = sdTables.scoreRow0[s] - opt.penalty_indel;
                }
                else 
                {
                    int ms = 0;
                    if (motifs[m].seq[mi] == vntr[s]) { ms = 0;} else { ms = -opt.penalty_mismatch;}	  
                    diagScore = sdTables.scoreMat[m][s][mi-1] + ms;
                    insScore = sdTables.scoreMat[m][s+1][mi-1] - opt.penalty_mismatch;
                }

                int delScore = sdTables.scoreMat[m][s][mi] - opt.penalty_indel; //prev index =s, cur=s+1;
                int optScore =0;
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
                    if (mi > 0) 
                        sdTables.nMatchMat[m][s+1][mi] = sdTables.nMatchMat[m][s][mi-1] + 1;
                    else 
                        sdTables.nMatchMat[m][s+1][mi] = 1;
                    

                    assert(0 <= m and  m < sdTables.nDelMat.size());
                    assert(0 <= s + 1 and s + 1 < sdTables.nDelMat[m].size());
                    assert(0 <= mi and mi < sdTables.nDelMat[m][s+1].size());

                    assert(0 <= s and s < sdTables.nDelMat[m].size());
                    assert(0 <= mi and mi < sdTables.nDelMat[m][s].size());
                    sdTables.nDelMat[m][s+1][mi] = sdTables.nDelMat[m][s][mi];
                    //	  cout << "nm " << m << " " << s+1<< " " << mi << " " <<  sdTables.nMatchMat[m][s+1][mi] << endl;
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
        //
        // score loop back.
        //
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
    //
    // Trace back;
    //
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
    motifNMatches.push_back(vector<int>());
    for (auto tm = 0; tm < motifs.size(); tm++) 
    {
        motifNMatches[motifNMatches.size()-1].push_back(max(0, sdTables.nMatchMat[tm][cs][motifs[tm].len-1] - sdTables.nDelMat[tm][cs][motifs[tm].len-1] ) );
    }
    while (cs > 0) 
    {
        int arrow = sdTables.pathMat[cm][cs][cmi];
        //
        // At the beginning of a motif, may need to wrap around
        //
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
                    motifNMatches.push_back(vector<int>());
                    for (auto tm=0; tm < motifs.size(); tm++) 
                    {
                        motifNMatches[motifNMatches.size()-1].push_back(max(0, sdTables.nMatchMat[tm][cs][motifs[tm].len-1] - sdTables.nDelMat[tm][cs][motifs[tm].len-1] ) );
                    }
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
                cmi--;
            else  
                cs--;
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
    for (auto m=0; m < optMotifs.size(); m++) 
    {
        int qv=scoring(motifs, motifNMatches[m], optMotifs[m], 1-opt.accuracy, mismatchCI);
        motifQV.push_back(qv);
    }
}

// function to compute the S_i distances (the naive occurrence)
void bounded_anno(vector<uint8_t> &optMotifs, vector<MOTIF> &motifs, char *vntr, int N, const OPTION &opt) 
{
    vector<uint8_t> optAnno;
    optAnno.clear();
    int i; // i as position index (of vntr), m as motif index, k as distance index for "compare" function
    uint8_t M = motifs.size();
    uint8_t k, m;

    vector<double> dist(N + 1);
    vector<double> update(M);
    vector<bool> gap(M);
    vector<int> traceI(N + 1);
    vector<uint8_t> traceM(N + 1);
    double ratio = 0.80;

    // generate dp matrix for each motif
    vector<vector<double>> dists(M);
    vector<vector<int>> starts(M);
    for (i = 0; i < M; ++i) {dists[i].resize(N + 1, 0);}
    for (i = 0; i < M; ++i) {starts[i].resize(N + 1, 0);}

    for (m = 0; m < M; m++) {
        global(motifs[m].seq, vntr, dists, starts, m, N, opt);
    }

    // calculate S_i distances
    // initialization for i=0
    dist[0] = 0;
    traceI[0] = 0;
    traceM[0] = 0;

    // propagate
    for (i = 1; i <= N; i++) {
        for (m = 0; m < M; m++) {
            if (dists[m][i] > motifs[m].len * ratio){
                update[m] = dist[starts[m][i]] + motifs[m].len * ratio;
                gap[m] = 1;
            } else {
                update[m] = dist[starts[m][i]] + dists[m][i];
                gap[m] = 0;
            }
        }
        k = compare(update, M); // k gives the index of the picked motif
        dist[i] = update[k];
        traceI[i] = starts[k][i];
        if (gap[k] == 1) {
            traceM[i] = M;
        } else {
            traceM[i] = k;
        }
    }

    // traceback
    i = N;
    while (i > 0) {
        optAnno.push_back(traceM[i]);
        // assert(traceI[i] != i);
        if (traceI[i] == i) break;
        i = traceI[i];
    }
    reverse(optAnno.begin(), optAnno.end());

    // skip the gap 
    optMotifs.clear();
    for (auto &it : optAnno)
    {
        if (it < M)
            optMotifs.push_back(it);
    }

    // assert(optMotifs.size() > 0);
    return;
}
