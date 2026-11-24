// akernel — B-spline basis functions, knot utilities and quadrature.
// Algorithms A2.1 and A2.3 of Piegl & Tiller, "The NURBS Book".
#pragma once

#include <algorithm>
#include <cmath>
#include <vector>

#include <akernel/num/detmath.hpp>

namespace ak::geom {

/// Open (clamped) uniform knot vector on [0,1] for `ncp` control points of degree p.
inline std::vector<double> open_uniform_knots(int ncp, int p) {
    std::vector<double> U(static_cast<std::size_t>(ncp + p + 1));
    const int spans = ncp - p;
    for (int i = 0; i <= ncp + p; ++i) {
        if (i <= p) U[i] = 0.0;
        else if (i >= ncp) U[i] = 1.0;
        else U[i] = static_cast<double>(i - p) / spans;
    }
    return U;
}

/// Knot span index i with U[i] <= u < U[i+1]; n is the last control-point index.
inline int find_span(int n, int p, double u, const std::vector<double>& U) {
    if (u >= U[n + 1]) return n;
    if (u <= U[p]) return p;
    int low = p, high = n + 1, mid = (low + high) / 2;
    while (u < U[mid] || u >= U[mid + 1]) {
        if (u < U[mid]) high = mid; else low = mid;
        mid = (low + high) / 2;
    }
    return mid;
}

/// ders[k][j] = k-th derivative of basis function N_{span-p+j,p} at u,
/// for k = 0..nders and j = 0..p. Derivatives of order > p are zero.
inline std::vector<std::vector<double>> ders_basis_funs(int span, double u, int p, int nders,
                                                        const std::vector<double>& U) {
    std::vector<std::vector<double>> ders(nders + 1, std::vector<double>(p + 1, 0.0));
    std::vector<std::vector<double>> ndu(p + 1, std::vector<double>(p + 1, 0.0));
    std::vector<double> left(p + 1), right(p + 1);

    ndu[0][0] = 1.0;
    for (int j = 1; j <= p; ++j) {
        left[j] = u - U[span + 1 - j];
        right[j] = U[span + j] - u;
        double saved = 0.0;
        for (int r = 0; r < j; ++r) {
            ndu[j][r] = right[r + 1] + left[j - r];
            const double temp = ndu[r][j - 1] / ndu[j][r];
            ndu[r][j] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        ndu[j][j] = saved;
    }
    for (int j = 0; j <= p; ++j) ders[0][j] = ndu[j][p];

    const int nd = std::min(nders, p);
    std::vector<std::vector<double>> a(2, std::vector<double>(p + 1, 0.0));
    for (int r = 0; r <= p; ++r) {
        int s1 = 0, s2 = 1;
        a[0][0] = 1.0;
        for (int k = 1; k <= nd; ++k) {
            double d = 0.0;
            const int rk = r - k, pk = p - k;
            if (r >= k) {
                a[s2][0] = a[s1][0] / ndu[pk + 1][rk];
                d = a[s2][0] * ndu[rk][pk];
            }
            const int j1 = rk >= -1 ? 1 : -rk;
            const int j2 = r - 1 <= pk ? k - 1 : p - r;
            for (int j = j1; j <= j2; ++j) {
                a[s2][j] = (a[s1][j] - a[s1][j - 1]) / ndu[pk + 1][rk + j];
                d += a[s2][j] * ndu[rk + j][pk];
            }
            if (r <= pk) {
                a[s2][k] = -a[s1][k - 1] / ndu[pk + 1][r];
                d += a[s2][k] * ndu[r][pk];
            }
            ders[k][r] = d;
            std::swap(s1, s2);
        }
    }
    double fac = p;
    for (int k = 1; k <= nd; ++k) {
        for (int j = 0; j <= p; ++j) ders[k][j] *= fac;
        fac *= (p - k);
    }
    return ders;
}

/// Greville abscissae: the natural collocation sites of a B-spline space.
inline std::vector<double> greville(int ncp, int p, const std::vector<double>& U) {
    std::vector<double> g(ncp);
    for (int i = 0; i < ncp; ++i) {
        double s = 0.0;
        for (int k = 1; k <= p; ++k) s += U[i + k];
        g[i] = s / p;
    }
    return g;
}

/// Gauss-Legendre nodes and weights on [-1, 1], exact for degree 2n-1.
inline void gauss_legendre(int n, std::vector<double>& x, std::vector<double>& w) {
    x.assign(n, 0.0);
    w.assign(n, 0.0);
    const double pi = 3.14159265358979323846;
    for (int i = 0; i < (n + 1) / 2; ++i) {
        double z = num::det::cos(pi * (i + 0.75) / (n + 0.5));
        double pp = 1.0;
        for (int it = 0; it < 100; ++it) {
            double p1 = 1.0, p2 = 0.0;
            for (int j = 1; j <= n; ++j) {
                const double p3 = p2;
                p2 = p1;
                p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3) / j;
            }
            pp = n * (z * p1 - p2) / (z * z - 1.0);
            const double dz = p1 / pp;
            z -= dz;
            if (std::fabs(dz) < 1e-16) break;
        }
        x[i] = -z;
        x[n - 1 - i] = z;
        w[i] = w[n - 1 - i] = 2.0 / ((1.0 - z * z) * pp * pp);
    }
}

} // namespace ak::geom
