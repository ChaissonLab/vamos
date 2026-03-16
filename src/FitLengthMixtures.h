#include <cmath>
#include <vector>
#include <numeric>
#include <stdexcept>
#include <limits>

// ---------- helpers ----------

// Log-likelihood of a Poisson(lambda) for one observation x
double poissonLogPMF(int x, double lambda);

// ---------- single-Poisson fit ----------

struct SingleFit {
    double lambda;
    double logLikelihood;
};

SingleFit fitSingle(const std::vector<int>& data);

// ---------- two-component Poisson mixture fit (EM) ----------

struct MixtureFit {
    double lambda1, lambda2;
    double pi1;           // mixing weight of component 1
    double logLikelihood;
};

MixtureFit fitMixture(const std::vector<int>& data,
                      int    maxIter = 500,
                      double tol     = 1e-8);
// ---------- BIC-based model selection ----------

enum class PoissonModel { Single, Mixture };

struct DetectionResult {
    PoissonModel model;
    SingleFit    single;
    MixtureFit   mixture;
    double       bicSingle;
    double       bicMixture;
};

DetectionResult detectPoissonModel(const std::vector<int>& data);
