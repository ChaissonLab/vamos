#ifndef ACC_LOOKUP_TABLE
#define ACC_LOOKUP_TABLE

#include <boost/math/distributions/poisson.hpp>
#include <boost/math/distributions/negative_binomial.hpp>
#include <boost/math/distributions/binomial.hpp>

#include "vntr.h"

void CreateAccLookupTable(vector<VNTR*> &vntrs, double accuracy, vector<int > &mismatchCI, double ci);

#endif
