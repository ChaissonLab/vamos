#include "acc_lookup_table.h"

using boost::math::poisson;
using boost::math::cdf;


void CreateAccLookupTable(vector<VNTR*> &vntrs, double accuracy, vector<int > &mismatchCI, double ci=0.99) {
  int maxMotifLength=0;
  for (auto i = 0; i < vntrs.size(); i++) {
    VNTR *v=vntrs[i];
    for (auto m=0; m < v->motifs.size(); m++ ) {
      if (v->motifs[m].len > maxMotifLength) {
	maxMotifLength = v->motifs[m].len;
      }
    }
  }
  mismatchCI.resize(maxMotifLength+1);
  mismatchCI[0] = 0;
  for (auto i=1; i < mismatchCI.size(); i++ ) {
    double lambda=(1-accuracy)*i;
    const poisson dist(lambda);
    double cdf=0;
    for (auto j=0; j <= i; j++ ) {      
      double c=pdf(dist,j);
      cdf+=c;
      if (cdf > ci) {
	mismatchCI[i] = j;
	break;
      }
    }
  }
}

