#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

using namespace std;

// function to compute the distance of given pair of characters
static double distance(char a, char b) {

    double distance;
    if (a == b) {
        distance = 0; // distance for match
    } else if (a == '-' || b == '-') {
        distance = 1; // distance for indel
    } else {
        distance = 1; // distance for mismatch
    }
    return distance;
}


// function to find min edit distance in dynamic programming
static int compare(vector<double> &update, int l) {

    double min = update[0];
    int k = 0;

    for (int i=1; i<l; i++) {
        if (update[i] <= min) {
            min = update[i];
            k = i;
        }
    }
    return k;
}

// function to compute the variant N-W alignment (linear space, only computes opt distance and traces starting position)
// (leading gaps cost 0)
static void
global(const string &motif, char *vntr, vector<vector<double>> &dists, vector<vector<int>> &starts, int m, int N) { 

    uint8_t i, k;
    int j, oldStart, tempStart; // i as row (string A) index, j as column (string B) index, k as distance index for "compare" function
    uint8_t M = motif.size();
    double oldDist;
    vector<double> update(3);

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
            update[0] = dists[m][j-1] + distance(vntr[j-1], '-'); // use current value: dist[j-1] after update
            update[1] = dists[m][j] + distance('-', motif[i-1]); // use old value: dist[j] before update
            update[2] = oldDist + distance(vntr[j-1], motif[i-1]); // use old value: dist[j-1] before update

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

// function to compute the S_i distances (the naive occurrence)
void bounded_anno(vector<uint8_t> &optMotifs, vector<MOTIF> &motifs, char *vntr, int N) 
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
        global(motifs[m].seq, vntr, dists, starts, m, N);
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

    assert(optMotifs.size() > 0);
    return;
}
