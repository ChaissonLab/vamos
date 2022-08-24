#include "phase.h"
#include <stdlib.h>
#include <assert.h>
#include <iostream>
#include <set>

extern int debug_flag;


int  CountDiff(READ* r1, READ* r2) {
  int r1i=0, r2i=0;
  int nDiff=0;
  while (r1i < r1->snvs.size() and r2i < r2->snvs.size()) {
      while (r1i < r1->snvs.size() and r2i < r2->snvs.size() and
	     r1->snvs[r1i].pos < r2->snvs[r2i].pos) {
	r1i+=1;
      }
      while (r1i < r1->snvs.size() and r2i < r2->snvs.size() and
	     r2->snvs[r2i].pos < r1->snvs[r1i].pos) {
	r2i+=1;
      }
      if (r1i < r1->snvs.size() and r2i < r2->snvs.size()) {
	if (r1->snvs[r1i].pos == r2->snvs[r2i].pos) {
	  if (r1->snvs[r1i].nuc != r2->snvs[r2i].nuc) { nDiff++;}
	}
	r1i++;
	r2i++;
      }
  }
  return nDiff;
}

void MaxCutPhase(VNTR *vntr) {
  vector<READ*> &reads = vntr->reads;
  vector<vector<int> > diffs;
  int maxDiff=-1;
  int max1i=0,max2i=0;
  int nSNVs=0;
  for (auto i=0; i < reads.size(); i++) {
    nSNVs += reads[i]->snvs.size();
  }
  if (debug_flag and reads.size() > 0) {
    cerr << "Average " << ((float)nSNVs) / reads.size() << endl;
  }
  if (debug_flag and reads.size() == 0) {
    cerr << vntr->region << " NoReads"<< endl << endl << endl << endl;
    return;
  }
  vector<int> totDiffs(reads.size(), 0);
    
  diffs.resize(reads.size());
  for (auto i=0; i < reads.size(); i++) {
    diffs[i].resize(reads.size());
  }
  for (auto i=0; i < reads.size()-1; i++) {
    diffs[i][i] = 0;
    for (auto j=i+1; j< reads.size(); j++ ) {
      diffs[i][j] = CountDiff(reads[i], reads[j]);
      totDiffs[i] += diffs[i][j];
      totDiffs[j] += diffs[i][j];
      diffs[j][i] = diffs[i][j];
      if (diffs[i][j] > maxDiff) {
      	maxDiff = diffs[i][j];
      	max1i=i;
      	max2i=j;
      }
    }
  }
  vector<int> cut1, cut2, cut0;
  /*  int maxTotDiff=0, maxTotDiffIndex=0;
  
  for (auto i = 0; i < reads.size(); i++ ) {
    if (totDiffs[i] > maxTotDiff) {
      maxTotDiff = totDiffs[i];
      maxTotDiffIndex=i;
    }
  }
  int maxDiffCut=0;
  int maxDiffCutIndex=0;
  for (auto i = 0; i< reads.size(); i++ ) {
    if (diffs[maxTotDiffIndex][i] > maxDiffCut) {
      maxDiffCutIndex = i;
      maxDiffCut = diffs[maxTotDiffIndex][i];
    }
    }
  max1i=maxTotDiffIndex;
  max2i=maxDiffCutIndex;
  */
  cut1.push_back(max1i);
  cut2.push_back(max2i);
  
  for (auto j = 0; j < reads.size(); j++) {
    /*
    int diff1=0, diff2=0;
    if (j == max1i or j == max2i) {
      continue;
    }
    for (auto &i: cut1) {
      diff1+=diffs[i][j];
    }
    for (auto &i: cut2) {
      diff2 += diffs[i][j];
    }
    if (diff1 == diff2 ) {
      cut0.push_back(j);
    }
    else if (diff1 > diff2) {
      // i has the largest difference to members in 1, so add to 2
      cut2.push_back(j);
    }
    else {
      cut1.push_back(j);
    }
    */

    if (maxDiff <= 1 or diffs[j][max1i] == diffs[j][max2i]) {
      cut0.push_back(j);
    }
    else if (diffs[j][max1i] < diffs[j][max2i]) {
      // J is closer to max1i than max2i
      cut1.push_back(j);
    }
    else {
      cut2.push_back(j);
    }
  }

  if (((float)cut0.size() )/(cut0.size() + cut1.size() + cut2.size()) > 0.25) {
    for (auto &v: cut1) {
      cut0.push_back(v);
    }
    for (auto &v: cut2) {
      cut0.push_back(v);
    }
    cut1.resize(0);
    cut2.resize(0);
  }
  
  if (debug_flag) {
    cerr << vntr->region << " cut 1 ";
    for (auto &idx: cut1 ) { cerr << idx << "/" << reads[idx]->haplotype << " "; }
    cerr << endl;
    cerr << vntr->region << " cut 2 ";
    for (auto &idx: cut2 ) { cerr << idx << "/" << reads[idx]->haplotype << " "; }
    cerr << endl;
    cerr << vntr->region << " cut 0 ";  
    for (auto &idx: cut0 ) { cerr << idx << "/" << reads[idx]->haplotype << " "; }
    cerr << endl;    
  }

  
}
