// akernel — variational four-sided fill (docs/BLUEPRINT.md §1.7 and §4.4).
//
// Unknown: a tensor-product B-spline patch on [0,1]^2. On every side the
// cross-boundary jet up to order k is prescribed, which fixes the k+1 rows
// of control points next to that side. The remaining control points minimise
//     E = sum_m w_m * integral |D^m S|^2 ,
// a quadratic form assembled exactly by Gauss quadrature. With the jets taken
// from the neighbouring surfaces in a common parameterisation this gives
// C^k, hence G^k, joins.
//
// Order rule: prescribing jets of order 0..k is well posed for the energy of
// order m = k+1 (thin plate for G1, m = 3 for G2, m = 4 for G3).
//
// Not covered here: gauge functions other than the identity (true G^k with a
// reparameterisation), N-sided regions, inequality constraints.
#pragma once

#include <array>
#include <cmath>
#include <functional>
#include <vector>

#include <akernel/geom/bspline.hpp>
#include <akernel/geom/nurbs_surface.hpp>
#include <akernel/num/dense.hpp>

namespace ak::var {

using geom::Vec3d;

enum Side : int { U0 = 0, U1 = 1, V0 = 2, V1 = 3 };

/// jet(side, j, s): j-th derivative across the side, with respect to the
/// patch's own parameter (u for U0/U1, v for V0/V1), at position s along it.
using JetFn = std::function<Vec3d(int side, int j, double s)>;

struct FillSpec {
    int degree = 7;
    int ncp = 14;       ///< control points per direction
    int continuity = 3; ///< k: jets of order 0..k are prescribed on all four sides
    std::array<double, 5> energy_weights = {0, 0, 0, 0, 1}; ///< w_m for m = 0..4
};

struct FillResult {
    bool ok = false;
    geom::NurbsSurface surface;
    int free_dofs = 0;
    double corner_mismatch = 0.0;      ///< largest disagreement between sides at shared control points
    std::array<double, 5> energy{};    ///< integral |D^m S|^2 of the solution, m = 0..4
};

namespace detail {

/// Gram matrix A[i][j] = integral over [0,1] of N_i^(a) N_j^(a).
inline std::vector<double> gram(int n, int p, int a, const std::vector<double>& U) {
    std::vector<double> A(static_cast<std::size_t>(n * n), 0.0);
    if (a > p) return A;
    std::vector<double> gx, gw;
    geom::gauss_legendre(p + 1, gx, gw);
    for (int span = p; span < n; ++span) {
        const double a0 = U[span], a1 = U[span + 1];
        if (!(a1 > a0)) continue;
        for (std::size_t g = 0; g < gx.size(); ++g) {
            const double u = 0.5 * (a0 + a1) + 0.5 * (a1 - a0) * gx[g];
            const double w = 0.5 * (a1 - a0) * gw[g];
            const auto d = geom::ders_basis_funs(span, u, p, a, U);
            for (int i = 0; i <= p; ++i)
                for (int j = 0; j <= p; ++j)
                    A[static_cast<std::size_t>((span - p + i) * n + (span - p + j))] += w * d[a][i] * d[a][j];
        }
    }
    return A;
}

inline double binomial(int n, int k) {
    double c = 1.0;
    for (int i = 1; i <= k; ++i) c = c * (n - k + i) / i;
    return c;
}

/// integral |D^m S|^2 for m = 0..4, by Gauss quadrature of the derivatives
/// themselves. A sum of squares, so it cannot go negative through rounding
/// (the quadratic form x^T K x can, when the true energy is near zero).
inline std::array<double, 5> energies(const std::vector<Vec3d>& cp, int n, int p,
                                      const std::vector<double>& U) {
    struct Node {
        int span;
        double w;
        std::vector<std::vector<double>> d;
    };
    std::vector<double> gx, gw;
    geom::gauss_legendre(p + 1, gx, gw);
    std::vector<Node> nodes;
    for (int span = p; span < n; ++span) {
        const double a0 = U[span], a1 = U[span + 1];
        if (!(a1 > a0)) continue;
        for (std::size_t g = 0; g < gx.size(); ++g) {
            const double u = 0.5 * (a0 + a1) + 0.5 * (a1 - a0) * gx[g];
            nodes.push_back({span, 0.5 * (a1 - a0) * gw[g], geom::ders_basis_funs(span, u, p, 4, U)});
        }
    }
    std::array<double, 5> E{};
    for (const Node& qu : nodes)
        for (const Node& qv : nodes) {
            Vec3d D[5][5]{};
            for (int i = 0; i <= p; ++i)
                for (int j = 0; j <= p; ++j) {
                    const Vec3d& P = cp[static_cast<std::size_t>((qu.span - p + i) * n + (qv.span - p + j))];
                    for (int a = 0; a <= 4; ++a)
                        for (int b = 0; b <= 4 - a; ++b) D[a][b] += P * (qu.d[a][i] * qv.d[b][j]);
                }
            for (int m = 0; m <= 4; ++m) {
                double s = 0.0;
                for (int a = 0; a <= m; ++a) s += binomial(m, a) * geom::dot(D[a][m - a], D[a][m - a]);
                E[m] += qu.w * qv.w * s;
            }
        }
    return E;
}

} // namespace detail

inline FillResult solve_fill(const FillSpec& spec, const JetFn& jet) {
    FillResult out;
    const int p = spec.degree, n = spec.ncp, k = spec.continuity;
    if (p < 1 || k < 0 || k > 3 || p < k || n < 2 * (k + 1) || n < p + 1) return out;

    const std::vector<double> U = geom::open_uniform_knots(n, p);
    const std::vector<double> gv = geom::greville(n, p, U);
    auto idx = [n](int i, int j) { return static_cast<std::size_t>(i * n + j); };

    // End derivatives of the basis: D0[j][i] = N_i^(j)(0), D1[j][i] = N_{n-1-i}^(j)(1).
    const auto d0 = geom::ders_basis_funs(p, 0.0, p, k, U);
    const auto d1 = geom::ders_basis_funs(n - 1, 1.0, p, k, U);

    // Collocation matrix at the Greville sites, factored once per use.
    std::vector<double> colloc(static_cast<std::size_t>(n * n), 0.0);
    for (int r = 0; r < n; ++r) {
        const int span = geom::find_span(n - 1, p, gv[r], U);
        const auto b = geom::ders_basis_funs(span, gv[r], p, 0, U);
        for (int i = 0; i <= p; ++i) colloc[idx(r, span - p + i)] = b[0][i];
    }

    // Control points fixed by the jets. Corner blocks are set by two sides.
    std::vector<Vec3d> cp(static_cast<std::size_t>(n * n));
    std::vector<int> hits(static_cast<std::size_t>(n * n), 0);
    for (int side = 0; side < 4; ++side) {
        const bool at_end = side == U1 || side == V1;
        // rows[i][r]: i-th control row from the boundary, sampled at Greville site r
        std::vector<double> rhs(static_cast<std::size_t>(n * (k + 1) * 3));
        for (int r = 0; r < n; ++r) {
            std::vector<Vec3d> row(k + 1);
            for (int j = 0; j <= k; ++j) {
                Vec3d v = jet(side, j, gv[r]);
                for (int i = 0; i < j; ++i) v -= row[i] * (at_end ? d1[j][p - i] : d0[j][i]);
                row[j] = v / (at_end ? d1[j][p - j] : d0[j][j]);
            }
            for (int i = 0; i <= k; ++i) {
                rhs[static_cast<std::size_t>(r * (k + 1) * 3 + i * 3 + 0)] = row[i].x;
                rhs[static_cast<std::size_t>(r * (k + 1) * 3 + i * 3 + 1)] = row[i].y;
                rhs[static_cast<std::size_t>(r * (k + 1) * 3 + i * 3 + 2)] = row[i].z;
            }
        }
        std::vector<double> M = colloc;
        if (!num::solve_linear(M, rhs, n, (k + 1) * 3)) return out;
        for (int l = 0; l < n; ++l)
            for (int i = 0; i <= k; ++i) {
                const std::size_t o = static_cast<std::size_t>(l * (k + 1) * 3 + i * 3);
                const Vec3d v{rhs[o], rhs[o + 1], rhs[o + 2]};
                const int across = at_end ? n - 1 - i : i;
                const std::size_t c = (side == U0 || side == U1) ? idx(across, l) : idx(l, across);
                if (hits[c] > 0)
                    out.corner_mismatch = std::max(out.corner_mismatch, geom::norm(cp[c] / hits[c] - v));
                cp[c] += v;
                ++hits[c];
            }
    }
    for (std::size_t c = 0; c < cp.size(); ++c)
        if (hits[c] > 1) cp[c] = cp[c] / static_cast<double>(hits[c]);

    // Energy matrices K_m = sum_a C(m,a) A^(a) (x) A^(m-a).
    std::array<std::vector<double>, 5> G;
    for (int a = 0; a <= 4; ++a) G[a] = detail::gram(n, p, a, U);
    const int N = n * n;
    auto energy_entry = [&](int m, int i, int j, int i2, int j2) {
        double s = 0.0;
        for (int a = 0; a <= m; ++a)
            s += detail::binomial(m, a) * G[a][idx(i, i2)] * G[m - a][idx(j, j2)];
        return s;
    };

    std::vector<int> free_ids;
    for (int c = 0; c < N; ++c)
        if (hits[static_cast<std::size_t>(c)] == 0) free_ids.push_back(c);
    const int nf = static_cast<int>(free_ids.size());
    out.free_dofs = nf;

    if (nf > 0) {
        std::vector<double> Kff(static_cast<std::size_t>(nf * nf), 0.0);
        std::vector<double> rhs(static_cast<std::size_t>(nf * 3), 0.0);
        for (int a = 0; a < nf; ++a) {
            const int i = free_ids[a] / n, j = free_ids[a] % n;
            int col = 0;
            for (int c = 0; c < N; ++c) {
                const int i2 = c / n, j2 = c % n;
                // B-spline supports: entries vanish unless both index gaps are <= p.
                if (std::abs(i - i2) > p || std::abs(j - j2) > p) {
                    if (hits[static_cast<std::size_t>(c)] == 0) ++col;
                    continue;
                }
                double kij = 0.0;
                for (int m = 0; m <= 4; ++m)
                    if (spec.energy_weights[m] != 0.0)
                        kij += spec.energy_weights[m] * energy_entry(m, i, j, i2, j2);
                if (hits[static_cast<std::size_t>(c)] == 0) {
                    Kff[static_cast<std::size_t>(a * nf + col)] = kij;
                    ++col;
                } else {
                    const Vec3d& x = cp[static_cast<std::size_t>(c)];
                    rhs[static_cast<std::size_t>(a * 3 + 0)] -= kij * x.x;
                    rhs[static_cast<std::size_t>(a * 3 + 1)] -= kij * x.y;
                    rhs[static_cast<std::size_t>(a * 3 + 2)] -= kij * x.z;
                }
            }
        }
        if (!num::cholesky(Kff, nf)) return out;
        num::cholesky_solve(Kff, rhs, nf, 3);
        for (int a = 0; a < nf; ++a)
            cp[static_cast<std::size_t>(free_ids[a])] = {rhs[static_cast<std::size_t>(a * 3)],
                                                         rhs[static_cast<std::size_t>(a * 3 + 1)],
                                                         rhs[static_cast<std::size_t>(a * 3 + 2)]};
    }

    out.surface = geom::NurbsSurface(p, p, n, n, U, U);
    for (int c = 0; c < N; ++c) {
        const Vec3d& x = cp[static_cast<std::size_t>(c)];
        out.surface.cp[static_cast<std::size_t>(c)] = {x.x, x.y, x.z, 1.0};
    }

    out.energy = detail::energies(cp, n, p, U);
    out.ok = true;
    return out;
}

} // namespace ak::var
