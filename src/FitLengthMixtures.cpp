#include <cmath>
#include <vector>
#include <numeric>
#include <stdexcept>
#include <limits>
#include <algorithm>
#include "FitLengthMixtures.h"
// ---------- helpers ----------

// Log-likelihood of a Poisson(lambda) for one observation x
double poissonLogPMF(int x, double lambda) {
    if (lambda <= 0.0) return -std::numeric_limits<double>::infinity();
    // log P(X=x) = x*log(lambda) - lambda - log(x!)
    double logFactX = 0.0;
    for (int i = 2; i <= x; ++i) logFactX += std::log(static_cast<double>(i));
    return x * std::log(lambda) - lambda - logFactX;
}

// ---------- single-Poisson fit ----------

SingleFit fitSingle(const std::vector<int>& data) {
    double mean = std::accumulate(data.begin(), data.end(), 0.0) / data.size();
    double ll = 0.0;
    for (int x : data) ll += poissonLogPMF(x, mean);
    return {mean, ll};
}



DetectionResult detectPoissonModel(const std::vector<int>& data) {
    if (data.size() < 3)
        throw std::invalid_argument("Need at least 3 data points.");
    double n = static_cast<double>(data.size());


    auto sf = fitSingle(data);
    auto mf = fitMixture(data);

    // BIC = -2 * logL + k * ln(n)
    // single Poisson: 1 free parameter (lambda)
    // mixture of two Poissons: 3 free parameters (lambda1, lambda2, pi)
    double bicSingle  = -2.0 * sf.logLikelihood + 1.0 * std::log(n);
    double bicMixture = -2.0 * mf.logLikelihood + 3.0 * std::log(n);

    PoissonModel chosen = (bicMixture < bicSingle)
                          ? PoissonModel::Mixture
                          : PoissonModel::Single;

    return {chosen, sf, mf, bicSingle, bicMixture};
}

DetectionResult detectPoissonModel(const std::vector<int>& guess1,
                                    const std::vector<int>& guess2,
                                    int    maxIter,
                                    double tol ) {
    if (guess1.empty() || guess2.empty())
        throw std::invalid_argument("Both guess vectors must be non-empty.");

    // Combine for single-model fit
    std::vector<int> data;
    data.insert(data.end(), guess1.begin(), guess1.end());
    data.insert(data.end(), guess2.begin(), guess2.end());

    double n = static_cast<double>(data.size());
    if (n < 3)
        throw std::invalid_argument("Need at least 3 data points total.");

    auto sf = fitSingle(data);
    auto mf = fitMixture(guess1, guess2, maxIter, tol);

    double bicSingle  = -2.0 * sf.logLikelihood + 1.0 * std::log(n);
    double bicMixture = -2.0 * mf.logLikelihood + 3.0 * std::log(n);

    PoissonModel chosen = (bicMixture < bicSingle)
                          ? PoissonModel::Mixture
                          : PoissonModel::Single;

    return {chosen, sf, mf, bicSingle, bicMixture};
}

// --- internal core: runs EM given initial lambda estimates ---
MixtureFit fitMixtureFromInit(const std::vector<int>& data,
                               double initLam1, double initLam2,
                               int    maxIter,
                               double tol ) {
    const int n = static_cast<int>(data.size());

    double lam1 = initLam1;
    double lam2 = initLam2;
    double pi1  = 0.5;

    std::vector<double> r(n);
    double prevLL = -std::numeric_limits<double>::infinity();

    for (int iter = 0; iter < maxIter; ++iter) {
        // E-step
        double curLL = 0.0;
        for (int i = 0; i < n; ++i) {
            double p1    = pi1         * std::exp(poissonLogPMF(data[i], lam1));
            double p2    = (1.0 - pi1) * std::exp(poissonLogPMF(data[i], lam2));
            double total = p1 + p2;
            r[i]   = (total > 0.0) ? p1 / total : 0.5;
            curLL += std::log(total > 0.0 ? total : 1e-300);
        }

        // M-step
        double R1 = 0.0, R2 = 0.0, S1 = 0.0, S2 = 0.0;
        for (int i = 0; i < n; ++i) {
            R1 += r[i];
            R2 += (1.0 - r[i]);
            S1 += r[i]         * data[i];
            S2 += (1.0 - r[i]) * data[i];
        }
        pi1  = R1 / n;
        lam1 = (R1 > 1e-12) ? S1 / R1 : 1e-6;
        lam2 = (R2 > 1e-12) ? S2 / R2 : 1e-6;

        if (lam1 > lam2) {
            std::swap(lam1, lam2);
            pi1 = 1.0 - pi1;
            for (auto& ri : r) ri = 1.0 - ri;
        }

        if (std::abs(curLL - prevLL) < tol) break;
        prevLL = curLL;
    }

    return {lam1, lam2, pi1, prevLL};
}

// --- two-vector version: user supplies initial guess partitions ---
MixtureFit fitMixture(const std::vector<int>& guess1,
                      const std::vector<int>& guess2,
		      int maxIter, double tol) {
    if (guess1.empty() || guess2.empty())
        throw std::invalid_argument("Both guess vectors must be non-empty.");

    double initLam1 = std::accumulate(guess1.begin(), guess1.end(), 0.0) / guess1.size();
    double initLam2 = std::accumulate(guess2.begin(), guess2.end(), 0.0) / guess2.size();

    // Combine into one dataset for EM
    std::vector<int> data;
    data.insert(data.end(), guess1.begin(), guess1.end());
    data.insert(data.end(), guess2.begin(), guess2.end());

    return fitMixtureFromInit(data, initLam1, initLam2, maxIter, tol);
}

// --- one-vector version: no guess, uses sorted split heuristic ---
MixtureFit fitMixture(const std::vector<int>& data, int maxIter, double tol) {
    const int n = static_cast<int>(data.size());

    std::vector<int> sorted(data);
    std::sort(sorted.begin(), sorted.end());

    int    mid      = n / 2;
    double initLam1 = std::accumulate(sorted.begin(),        sorted.begin() + mid, 0.0) / mid;
    double initLam2 = std::accumulate(sorted.begin() + mid,  sorted.end(),         0.0) / (n - mid);

    return fitMixtureFromInit(data, initLam1, initLam2, maxIter, tol);
}
