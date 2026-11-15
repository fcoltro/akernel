// akernel — convergence-rate estimators for the diagnostics loop
// (docs/BLUEPRINT.md §2.2 d).
//
// A smooth problem converges at its theoretical order. An observed order
// below that is evidence of a nearby singularity, available before any
// hard failure threshold is reached.
#pragma once

#include <cmath>
#include <cstddef>
#include <vector>

namespace ak::diag {

struct RateFit {
    double order = 0.0;    ///< p in  error ~ C h^p
    double constant = 0.0; ///< C
    bool ok = false;
};

/// Least-squares fit of log(error) = log C + p log h.
inline RateFit fit_rate(const std::vector<double>& h, const std::vector<double>& error) {
    RateFit r;
    const std::size_t n = h.size();
    if (n < 2 || error.size() != n) return r;
    double sx = 0, sy = 0, sxx = 0, sxy = 0;
    for (std::size_t i = 0; i < n; ++i) {
        if (!(h[i] > 0.0) || !(error[i] > 0.0)) return r;
        const double x = std::log(h[i]), y = std::log(error[i]);
        sx += x; sy += y; sxx += x * x; sxy += x * y;
    }
    const double den = n * sxx - sx * sx;
    if (!(std::fabs(den) > 0.0)) return r;
    r.order = (n * sxy - sx * sy) / den;
    r.constant = std::exp((sy - r.order * sx) / n);
    r.ok = true;
    return r;
}

/// Observed order q of an iteration, from e_{k+1} ~ C e_k^q, using the last
/// three errors that are still above the rounding floor. Newton gives q = 2
/// at a regular root and q = 1 at a singular (tangential) one.
inline double iteration_order(const std::vector<double>& errors, double floor = 1e-13) {
    std::size_t last = errors.size();
    while (last > 0 && !(errors[last - 1] > floor)) --last;
    if (last < 3) return 0.0;
    const double e0 = errors[last - 3], e1 = errors[last - 2], e2 = errors[last - 1];
    const double den = std::log(e1 / e0);
    return den != 0.0 ? std::log(e2 / e1) / den : 0.0;
}

} // namespace ak::diag
