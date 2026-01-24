#include "phase.h"
#include <stdlib.h>
#include <assert.h>
#include <iostream>
#include <set>

extern int debug_flag;


int  CountDiff(READ* r1, READ* r2, int phaseFlank, int start, int end) {
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
	if (r1->snvs[r1i].pos == r2->snvs[r2i].pos and
	    ((r1->snvs[r1i].pos < start and start - r1->snvs[r1i].pos < phaseFlank) or
	     (r1->snvs[r1i].pos > end and r1->snvs[r1i].pos - end < phaseFlank))) {
	  if (r1->snvs[r1i].nuc != r2->snvs[r2i].nuc) { nDiff++;}
	}
	r1i++;
	r2i++;
      }
  }
  return nDiff;
}

void MaxCutPhase(VNTR *vntr, int phaseFlank) {
  if (vntr->reads.size() == 0) {
    return;
  }
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
  int nReads = reads.size();
  diffs.resize(nReads);
  for (auto i=0; i < nReads; i++) {
    diffs[i].resize(nReads);
    assert(diffs[i].size() == nReads);
  }
  for (auto i=0; i + 1 < nReads; i++) {
    diffs[i][i] = 0;
    for (auto j=i+1; j< nReads; j++ ) {
      assert(i < nReads);
      assert(j < nReads);			  
      diffs[i][j] = CountDiff(reads[i], reads[j], phaseFlank, vntr->ref_start, vntr->ref_end);
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

  cut1.push_back(max1i);
  cut2.push_back(max2i);
  
  for (auto j = 0; j < reads.size(); j++) {
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
  for (auto &v: cut0) {
    reads[v]->haplotype = 0;
  }
  for (auto &v: cut1) {
    reads[v]->haplotype = 1;
  }
  for (auto &v : cut2) {
    reads[v]->haplotype = 2;
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
  if (debug_flag) {
    cerr << vntr->region << " cut 1 ";
    for (auto &idx: cut1 ) { cerr << reads[idx]->qname << "   " << reads[idx]->haplotype << "\n"; }
    cerr << endl;
    cerr << vntr->region << " cut 2 ";
    for (auto &idx: cut2 ) { cerr << reads[idx]->qname << "   " << reads[idx]->haplotype << "\n"; }
    cerr << endl;
    cerr << vntr->region << " cut 0 ";  
    for (auto &idx: cut0 ) { cerr << reads[idx]->qname << "   " << reads[idx]->haplotype << "\n"; }
    cerr << endl;    
  }
  
}
