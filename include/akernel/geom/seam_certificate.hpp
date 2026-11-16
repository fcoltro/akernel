// akernel — certified bounds on the geometric-continuity jumps across a seam
// between two polynomial Bezier patches (docs/BLUEPRINT.md §1.7; RESULTS.md
// §4, "Certified seam bounds").
//
// geom/continuity.hpp measures the jumps of the normal, of the second
// fundamental form b and of its covariant derivative C at sample points.
// Here they are bounded over the whole seam, rigorously.
//
// Setting. A and B are polynomial patches (all weights 1) in one common
// parameterisation: the seam is the same side of [0,1]^2 for both (for
// example v = 0), and B's jets along it are compared with A's at the same
// parameter. A neighbour that meets A along its own opposite side is
// brought to this form by restricting its polynomial beyond [0,1]^2 (the
// restriction is exact; see PolyPatch::restricted). This certifies jumps for
// joins that are C^k up to a residual in that parameterisation, which is
// how the variational fill (var/fill.hpp) is built; a G^k join through a
// nontrivial reparameterisation is not covered.
//
// Method, on each of `pieces` equal parts of the seam:
//  1. Residuals. Every jet entry D^(k,l) along the seam, k + l <= 3, is a
//     polynomial in the seam parameter. Its Bernstein coefficients on the
//     part are computed for A and for B by exact-structure de Casteljau
//     steps in interval arithmetic, raised to a common degree, and
//     subtracted. By the convex-hull property the largest coefficient
//     bounds the residual Delta = J_B - J_A on the part, per coordinate.
//  2. Mean value. Each jump is a smooth function H(J_A, J_B) with
//     H(J_A, J_A) = 0. So H(J_A, J_B) = G . Delta, with G = dH/dJ_B taken
//     at some point between J_A and J_B. G is enclosed over the part by
//     forward differentiation in interval arithmetic, with J_A ranging over
//     the hull of A's coefficients and J_B over the hull of both.
//  3. |jump| <= sum_i |G_i| eps_i, rounded upward.
// The bounds therefore shrink with the residual, however small it is,
// where enclosing each surface's invariants separately could not resolve
// a jump below the variation of the invariants along a part.
//
// Normal: the bound is on |n_B - n_A|, turned into an angle by
// angle = 2 asin(d / 2) <= d / sqrt(1 - d^2 / 4). b and C: their jumps in the
// (tau, nu) frame of A, tau the unit seam tangent and nu = n_A x tau, as in
// seam_jump. Given G1, G2 <=> db_nn = 0; given G2, G3 <=> dC_nnn = 0.
#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

#include <akernel/geom/bernstein.hpp>
#include <akernel/num/interval.hpp>

namespace ak::geom {

enum SeamSide : int { kSeamU0 = 0, kSeamU1 = 1, kSeamV0 = 2, kSeamV1 = 3 };

namespace seam_detail {

using num::Interval;

/// Interval value with an interval gradient over N variables (forward mode).
template <int N>
struct Dual {
    Interval v{};
    std::array<Interval, N> d{};
    Dual() = default;
    Dual(double x) : v(x) {}
    Dual(const Interval& x) : v(x) {}
};

template <int N> Dual<N> operator+(const Dual<N>& a, const Dual<N>& b) {
    Dual<N> r(a.v + b.v);
    for (int i = 0; i < N; ++i) r.d[i] = a.d[i] + b.d[i];
    return r;
}
template <int N> Dual<N> operator-(const Dual<N>& a, const Dual<N>& b) {
    Dual<N> r(a.v - b.v);
    for (int i = 0; i < N; ++i) r.d[i] = a.d[i] - b.d[i];
    return r;
}
template <int N> Dual<N> operator-(const Dual<N>& a) {
    Dual<N> r(-a.v);
    for (int i = 0; i < N; ++i) r.d[i] = -a.d[i];
    return r;
}
template <int N> Dual<N> operator*(const Dual<N>& a, const Dual<N>& b) {
    Dual<N> r(a.v * b.v);
    for (int i = 0; i < N; ++i) r.d[i] = a.d[i] * b.v + a.v * b.d[i];
    return r;
}
template <int N> Dual<N> operator/(const Dual<N>& a, const Dual<N>& b) {
    Dual<N> r(a.v / b.v);
    for (int i = 0; i < N; ++i) r.d[i] = (a.d[i] - r.v * b.d[i]) / b.v;
    return r;
}
template <int N> Dual<N>& operator+=(Dual<N>& a, const Dual<N>& b) { return a = a + b; }
template <int N> Dual<N>& operator-=(Dual<N>& a, const Dual<N>& b) { return a = a - b; }
template <int N> Dual<N> sqrt(const Dual<N>& a) {
    Dual<N> r(num::sqrt(a.v));
    const Interval two_r = Interval(2.0) * r.v;
    for (int i = 0; i < N; ++i) r.d[i] = a.d[i] / two_r;
    return r;
}

/// The invariants of geom/continuity.hpp, generic over the scalar type.
template <class T>
struct Frame {
    Vec3<T> Su, Sv, n;
    T ginv[2][2];
    T b[2][2];
    T C[2][2][2];
};

template <class T>
Frame<T> frame(const std::array<Vec3<T>, kJetSize>& J) {
    Frame<T> f{};
    f.Su = J[1];
    f.Sv = J[2];
    const Vec3<T> N = cross(f.Su, f.Sv);
    f.n = N / sqrt(dot(N, N));
    const Vec3<T> S1[2] = {f.Su, f.Sv};
    auto S2 = [&](int a, int c) -> const Vec3<T>& { return J[static_cast<std::size_t>(3 + a + c)]; };
    auto S3 = [&](int a, int c, int e) -> const Vec3<T>& { return J[static_cast<std::size_t>(6 + a + c + e)]; };
    T g[2][2];
    for (int a = 0; a < 2; ++a)
        for (int c = 0; c < 2; ++c) {
            g[a][c] = dot(S1[a], S1[c]);
            f.b[a][c] = dot(f.n, S2(a, c));
        }
    const T det = g[0][0] * g[1][1] - g[0][1] * g[1][0];
    f.ginv[0][0] = g[1][1] / det;
    f.ginv[1][1] = g[0][0] / det;
    f.ginv[0][1] = f.ginv[1][0] = -g[0][1] / det;
    T G[2][2][2];
    for (int d = 0; d < 2; ++d)
        for (int a = 0; a < 2; ++a)
            for (int c = 0; c < 2; ++c) {
                T s(0.0);
                for (int e = 0; e < 2; ++e) s += f.ginv[d][e] * dot(S1[e], S2(a, c));
                G[d][a][c] = s;
            }
    for (int a = 0; a < 2; ++a)
        for (int b2 = 0; b2 < 2; ++b2)
            for (int c = 0; c < 2; ++c) {
                T s = dot(f.n, S3(a, b2, c));
                for (int d = 0; d < 2; ++d)
                    s -= f.b[c][d] * G[d][a][b2] + f.b[d][b2] * G[d][c][a] + f.b[a][d] * G[d][c][b2];
                f.C[a][b2][c] = s;
            }
    return f;
}

template <class T, class F>
std::array<T, 2> to_param(const Frame<F>& f, const Vec3<Interval>& d) {
    const T r0 = dot(f.Su, Vec3<F>{F(d.x), F(d.y), F(d.z)}), r1 = dot(f.Sv, Vec3<F>{F(d.x), F(d.y), F(d.z)});
    return {f.ginv[0][0] * r0 + f.ginv[0][1] * r1, f.ginv[1][0] * r0 + f.ginv[1][1] * r1};
}

template <class T>
T form2(const Frame<T>& f, const std::array<T, 2>& x, const std::array<T, 2>& y) {
    T s(0.0);
    for (int a = 0; a < 2; ++a)
        for (int c = 0; c < 2; ++c) s += f.b[a][c] * x[a] * y[c];
    return s;
}
template <class T>
T form3(const Frame<T>& f, const std::array<T, 2>& x, const std::array<T, 2>& y, const std::array<T, 2>& z) {
    T s(0.0);
    for (int a = 0; a < 2; ++a)
        for (int c = 0; c < 2; ++c)
            for (int e = 0; e < 2; ++e) s += f.C[a][c][e] * x[a] * y[c] * z[e];
    return s;
}

using Curve = std::vector<Vec3<Interval>>; ///< Bernstein coefficients in the seam parameter

/// Jet entry (k, l) of P along `side`, as a curve in the seam parameter.
inline Curve jet_curve(const PolyPatchI& P, int k, int l, int side) {
    PolyPatchI d = P;
    for (int i = 0; i < k; ++i) d = d.deriv_u();
    for (int i = 0; i < l; ++i) d = d.deriv_v();
    Curve c;
    if (side == kSeamV0 || side == kSeamV1) {
        const int j = side == kSeamV0 ? 0 : d.q;
        for (int i = 0; i <= d.p; ++i) c.push_back(xyz(d.at(i, j)));
    } else {
        const int i = side == kSeamU0 ? 0 : d.p;
        for (int j = 0; j <= d.q; ++j) c.push_back(xyz(d.at(i, j)));
    }
    return c;
}

/// Degree elevation to `n` (interval coefficients).
inline Curve elevate(Curve c, int n) {
    while (static_cast<int>(c.size()) - 1 < n) {
        const int m = static_cast<int>(c.size()); // new degree
        Curve e(static_cast<std::size_t>(m + 1));
        e[0] = c[0];
        e[static_cast<std::size_t>(m)] = c[static_cast<std::size_t>(m - 1)];
        for (int i = 1; i < m; ++i) {
            const Interval a = Interval(static_cast<double>(i)) / Interval(static_cast<double>(m));
            e[static_cast<std::size_t>(i)] = c[static_cast<std::size_t>(i - 1)] * a + c[static_cast<std::size_t>(i)] * (Interval(1.0) - a);
        }
        c = std::move(e);
    }
    return c;
}

/// The curve restricted to [a, b] of its parameter.
inline Curve restrict_curve(Curve c, double a, double b) {
    if (c.size() <= 1 || (a == 0.0 && b == 1.0)) return c;
    detail::restrict1d<Vec3<Interval>, Interval>(c.data(), static_cast<int>(c.size()) - 1, a, b);
    return c;
}

inline Interval hull_of(const Curve& c, int coord) {
    auto get = [coord](const Vec3<Interval>& p) { return coord == 0 ? p.x : coord == 1 ? p.y : p.z; };
    Interval h = get(c[0]);
    for (const auto& p : c) h = num::hull(h, get(p));
    return h;
}

/// Upper bound of the square root of a non-negative interval's upper end.
inline double sqrt_up(const Interval& a) { return num::sqrt(Interval(0.0, std::max(0.0, a.hi))).hi; }

} // namespace seam_detail

/// Certified upper bounds over the whole seam (see the file comment).
struct SeamCertificate {
    bool ok = false;
    std::string why;           ///< why it failed, if it did
    double gap = 0.0;          ///< |S_B - S_A|
    double normal_angle = 0.0; ///< radians
    double db_tt = 0.0, db_tn = 0.0, db_nn = 0.0;
    double dC_ttt = 0.0, dC_ttn = 0.0, dC_tnn = 0.0, dC_nnn = 0.0;
    int pieces = 0;
};

/// Bounds for the seam `side` of A and B, which share the parameterisation.
/// Both must be polynomial (every weight exactly 1).
inline SeamCertificate certify_seam(const PolyPatchI& A, const PolyPatchI& B, int side, int pieces = 8) {
    using namespace seam_detail;
    constexpr int N = 3 * kJetSize;
    SeamCertificate out;
    for (const PolyPatchI* P : {&A, &B})
        for (const auto& cp : P->c)
            if (!(cp.w.lo == 1.0 && cp.w.hi == 1.0)) {
                out.why = "a patch is rational";
                return out;
            }
    // Jet entries along the seam, raised to a common degree per entry.
    std::array<Curve, kJetSize> ca, cb;
    for (int k = 0; k <= 3; ++k)
        for (int l = 0; l <= 3 - k; ++l) {
            const int e = jet_index(k, l);
            ca[static_cast<std::size_t>(e)] = jet_curve(A, k, l, side);
            cb[static_cast<std::size_t>(e)] = jet_curve(B, k, l, side);
            const int n = static_cast<int>(std::max(ca[static_cast<std::size_t>(e)].size(), cb[static_cast<std::size_t>(e)].size())) - 1;
            ca[static_cast<std::size_t>(e)] = elevate(ca[static_cast<std::size_t>(e)], n);
            cb[static_cast<std::size_t>(e)] = elevate(cb[static_cast<std::size_t>(e)], n);
        }
    const bool along_u = side == kSeamV0 || side == kSeamV1;
    const int tangent_entry = along_u ? jet_index(1, 0) : jet_index(0, 1);
    for (int m = 0; m < pieces; ++m) {
        const double t0 = static_cast<double>(m) / pieces, t1 = static_cast<double>(m + 1) / pieces;
        std::array<Vec3<Interval>, kJetSize> JA, HB;
        std::array<double, N> eps{};
        for (int e = 0; e < kJetSize; ++e) {
            const Curve a = restrict_curve(ca[static_cast<std::size_t>(e)], t0, t1);
            const Curve b = restrict_curve(cb[static_cast<std::size_t>(e)], t0, t1);
            Interval box_a[3], box_b[3];
            for (int c = 0; c < 3; ++c) {
                box_a[c] = hull_of(a, c);
                box_b[c] = hull_of(b, c);
                double r = 0.0;
                for (std::size_t i = 0; i < a.size(); ++i) {
                    const Vec3<Interval> d = b[i] - a[i];
                    r = std::max(r, num::mag(c == 0 ? d.x : c == 1 ? d.y : d.z));
                }
                eps[static_cast<std::size_t>(3 * e + c)] = r;
            }
            JA[static_cast<std::size_t>(e)] = {box_a[0], box_a[1], box_a[2]};
            HB[static_cast<std::size_t>(e)] = {num::hull(box_a[0], box_b[0]), num::hull(box_a[1], box_b[1]),
                                               num::hull(box_a[2], box_b[2])};
        }
        // A's frame over the part, and the seam frame (tau, nu) from it.
        const Frame<Interval> fa = frame(JA);
        const Vec3<Interval> T = JA[static_cast<std::size_t>(tangent_entry)];
        Vec3<Interval> tau = T - fa.n * dot(fa.n, T);
        tau = tau / sqrt(dot(tau, tau));
        const Vec3<Interval> nu = cross(fa.n, tau);
        // B's jet as the variables, over the hull of both.
        std::array<Vec3<Dual<N>>, kJetSize> JB;
        for (int e = 0; e < kJetSize; ++e) {
            Vec3<Dual<N>>& v = JB[static_cast<std::size_t>(e)];
            v = {Dual<N>(HB[static_cast<std::size_t>(e)].x), Dual<N>(HB[static_cast<std::size_t>(e)].y),
                 Dual<N>(HB[static_cast<std::size_t>(e)].z)};
            v.x.d[static_cast<std::size_t>(3 * e)] = Interval(1.0);
            v.y.d[static_cast<std::size_t>(3 * e + 1)] = Interval(1.0);
            v.z.d[static_cast<std::size_t>(3 * e + 2)] = Interval(1.0);
        }
        const Frame<Dual<N>> fb = frame(JB);
        // Orientation: the normals must agree (no flipped parameterisation).
        const Interval c_nn = dot(fa.n, Vec3<Interval>{fb.n.x.v, fb.n.y.v, fb.n.z.v});
        if (!(c_nn.lo > 0.0)) {
            out.why = "the normals are not proven to agree in orientation on part " + std::to_string(m);
            return out;
        }
        const auto tB = to_param<Dual<N>>(fb, tau), nB = to_param<Dual<N>>(fb, nu);
        // Each jump is (B's value as a Dual) - (A's value, a constant here).
        auto bound = [&](const Dual<N>& hb) {
            Interval s(0.0);
            for (int i = 0; i < N; ++i) s = s + Interval(num::mag(hb.d[static_cast<std::size_t>(i)])) * Interval(eps[static_cast<std::size_t>(i)]);
            return std::isfinite(s.hi) ? s.hi : std::numeric_limits<double>::infinity();
        };
        // Normal: |n_B - n_A| per coordinate, then the angle.
        const double dnx = bound(fb.n.x), dny = bound(fb.n.y), dnz = bound(fb.n.z);
        const Interval d2 = Interval(dnx) * Interval(dnx) + Interval(dny) * Interval(dny) + Interval(dnz) * Interval(dnz);
        const double dn = sqrt_up(d2);
        double angle = 3.14159265358979324; // the bound when nothing better is proven
        if (dn < 2.0) {
            const Interval q = Interval(1.0) - Interval(dn) * Interval(dn) / Interval(4.0);
            if (q.lo > 0.0) angle = std::min(angle, (Interval(dn) / num::sqrt(q)).hi);
        }
        double pos = 0.0;
        {
            const Interval p2 = Interval(eps[0]) * Interval(eps[0]) + Interval(eps[1]) * Interval(eps[1]) +
                                Interval(eps[2]) * Interval(eps[2]);
            pos = sqrt_up(p2);
        }
        out.gap = std::max(out.gap, pos);
        out.normal_angle = std::max(out.normal_angle, angle);
        out.db_tt = std::max(out.db_tt, bound(form2(fb, tB, tB)));
        out.db_tn = std::max(out.db_tn, bound(form2(fb, tB, nB)));
        out.db_nn = std::max(out.db_nn, bound(form2(fb, nB, nB)));
        out.dC_ttt = std::max(out.dC_ttt, bound(form3(fb, tB, tB, tB)));
        out.dC_ttn = std::max(out.dC_ttn, bound(form3(fb, tB, tB, nB)));
        out.dC_tnn = std::max(out.dC_tnn, bound(form3(fb, tB, nB, nB)));
        out.dC_nnn = std::max(out.dC_nnn, bound(form3(fb, nB, nB, nB)));
    }
    out.pieces = pieces;
    out.ok = true;
    return out;
}

inline SeamCertificate certify_seam(const PolyPatchD& A, const PolyPatchD& B, int side, int pieces = 8) {
    return certify_seam(to_interval(A), to_interval(B), side, pieces);
}

/// The polynomial patch P (every weight exactly 1) restricted to
/// [u0,u1] x [v0,v1], in interval arithmetic. The weights of a restriction of
/// a polynomial are exactly 1, and they are set so, rather than left as the
/// intervals that de Casteljau steps on them produce. A rational P is
/// restricted as it is (certify_seam then refuses it).
inline PolyPatchI restricted_polynomial(const PolyPatchD& P, double u0, double u1, double v0, double v1) {
    PolyPatchI r = to_interval(P).restricted(u0, u1, v0, v1);
    bool polynomial = true;
    for (const auto& cp : P.c) polynomial = polynomial && cp.w == 1.0;
    if (polynomial)
        for (auto& cp : r.c) cp.w = num::Interval(1.0);
    return r;
}

} // namespace ak::geom
