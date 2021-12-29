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
static int compare(double update[], int l) {

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
global(const string &motif, char *vntr, double dist[], uint8_t start[]) { 

    int i, j, k;
    uint8_t oldStart, tempStart; // i as row (string A) index, j as column (string B) index, k as distance index for "compare" function
    const int m = motif.size(), n = strlen(vntr);
    double oldDist, update[3];

    // initialization for the first row (dist=0, start=j for first row, 0 penalty for leading gaps)
    for (j=0; j<=n; j++) {
        dist[j] = 0;
        start[j] = j;
    }

    // propagate to later rows (dist=i, start=0 for first column)
    for (i=1; i<=m; i++) {
        oldDist = dist[0];
        oldStart = start[0];
        dist[0] = i;
        start[0] = 0;

        for (j=1; j<=n; j++) {

            // dynamic programming
            update[0] = dist[j-1] + distance(vntr[j-1], '-'); // use current value: dist[j-1] after update
            update[1] = dist[j] + distance('-', motif[i-1]); // use old value: dist[j] before update
            update[2] = oldDist + distance(vntr[j-1], motif[i-1]); // use old value: dist[j-1] before update

            k = compare(update, 3);

            if (k == 0) {
                tempStart = start[j-1];
            } else if (k == 2) {
                tempStart = oldStart;
            } else {
                tempStart = start[j];
            }

            oldDist = dist[j]; // save old value: dist[j-1] before update
            oldStart = start[j]; // save old value: start[j-1] before update

            // update dist and start
            dist[j] = update[k];
            start[j] = tempStart;
        }
    }
}


// function to compute the S_i distances (the naive occurrence)
int bounded_anno(vector<uint8_t> &optMotifs, vector<MOTIF> &motifs, char * vntr) 
{
    optMotifs.clear();
    int i, k, m; // i as position index (of vntr), m as motif index, k as distance index for "compare" function
    int n = strlen(vntr), M = motifs.size();
    double dist[n+1], update[M];
    uint8_t traceI[n+1], traceM[n+1];

    // generate dp matrix for each motif
    double dists[M][n+1];
    uint8_t starts[M][n+1];
    for (int m=0; m<M; m++) 
        global(motifs[m].seq, vntr, dists[m], starts[m]);

    // calculate S_i distances
    // initialization for i=0
    dist[0] = 0;
    traceI[0] = 0;
    traceM[0] = 0;

    // propagate
    for (i=1; i<=n; i++) {
        for (m=0; m<M; m++) {
			update[m] = dist[starts[m][i]] + dists[m][i];
        }
        k = compare(update, M); // k gives the index of the picked motif
        dist[i] = update[k];
        traceI[i] = starts[k][i];
        traceM[i] = k;
    }

    // traceback
    i = n;
    while (i > 0) {
        optMotifs.push_back(traceM[i]);
        i = traceI[i];
    }
    reverse(optMotifs.begin(), optMotifs.end());
    return 0;
}



