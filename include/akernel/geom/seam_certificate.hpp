// akernel — certified bounds on the geometric-continuity jumps across a seam
// between two Bezier patches, polynomial or rational (docs/BLUEPRINT.md
// §1.7; RESULTS.md §4, "Certified seam bounds").
//
// geom/continuity.hpp measures the jumps of the normal, of the second
// fundamental form b and of its covariant derivative C at sample points.
// Here they are bounded over the whole seam, rigorously.
//
// Setting. A and B are patches in one common parameterisation: the seam is the same side of [0,1]^2 for both (for
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
//     polynomial in the seam parameter: of the surface for polynomial
//     patches (weights exactly 1), of its homogeneous form for rational
//     ones, whose projection then belongs to H below. Rational patches must
//     therefore agree in homogeneous form up to the residual (pieces of one
//     rational parameterisation do). Its Bernstein coefficients on the
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
//  4. Parts. The enclosure of G over a part shrinks with the part, so more
//     parts give tighter bounds; through the projection of a rational patch
//     it shrinks more slowly, and its bounds are looser. A part on which the
//     orientation, a weight or the metric is not proven is halved.
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

using Curve = std::vector<Vec4<Interval>>; ///< homogeneous Bernstein coefficients in the seam parameter

/// Homogeneous jet entry (k, l) of P along `side`, as a curve in the seam parameter.
inline Curve jet_curve(const PolyPatchI& P, int k, int l, int side) {
    PolyPatchI d = P;
    for (int i = 0; i < k; ++i) d = d.deriv_u();
    for (int i = 0; i < l; ++i) d = d.deriv_v();
    Curve c;
    if (side == kSeamV0 || side == kSeamV1) {
        const int j = side == kSeamV0 ? 0 : d.q;
        for (int i = 0; i <= d.p; ++i) c.push_back(d.at(i, j));
    } else {
        const int i = side == kSeamU0 ? 0 : d.p;
        for (int j = 0; j <= d.q; ++j) c.push_back(d.at(i, j));
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
    detail::restrict1d<Vec4<Interval>, Interval>(c.data(), static_cast<int>(c.size()) - 1, a, b);
    return c;
}

inline const Interval& coord(const Vec4<Interval>& p, int c) { return c == 0 ? p.x : c == 1 ? p.y : c == 2 ? p.z : p.w; }

inline Interval hull_of(const Curve& c, int k) {
    Interval h = coord(c[0], k);
    for (const auto& p : c) h = num::hull(h, coord(p, k));
    return h;
}

/// Derivatives of the projected surface from those of its homogeneous form
/// (Piegl & Tiller, eq. 4.20; geom::rational_jet), generic over the scalar.
template <class T>
std::array<Vec3<T>, kJetSize> project_jet(const std::array<Vec4<T>, kJetSize>& A) {
    static constexpr double binom[4][4] = {{1, 0, 0, 0}, {1, 1, 0, 0}, {1, 2, 1, 0}, {1, 3, 3, 1}};
    std::array<Vec3<T>, kJetSize> S{};
    const T w = A[0].w;
    for (int k = 0; k <= 3; ++k)
        for (int l = 0; l <= 3 - k; ++l) {
            Vec3<T> v = xyz(A[static_cast<std::size_t>(jet_index(k, l))]);
            for (int j = 1; j <= l; ++j)
                v -= S[static_cast<std::size_t>(jet_index(k, l - j))] * (T(binom[l][j]) * A[static_cast<std::size_t>(jet_index(0, j))].w);
            for (int i = 1; i <= k; ++i) {
                v -= S[static_cast<std::size_t>(jet_index(k - i, l))] * (T(binom[k][i]) * A[static_cast<std::size_t>(jet_index(i, 0))].w);
                Vec3<T> v2{};
                for (int j = 1; j <= l; ++j)
                    v2 += S[static_cast<std::size_t>(jet_index(k - i, l - j))] *
                          (T(binom[l][j]) * A[static_cast<std::size_t>(jet_index(i, j))].w);
                v -= v2 * T(binom[k][i]);
            }
            S[static_cast<std::size_t>(jet_index(k, l))] = v / w;
        }
    return S;
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
    int pieces = 0;            ///< parts of the seam certified (after any bisection)
};

namespace seam_detail {

/// The nine bounds of one part, from boxes of A's (homogeneous) jet, of the
/// hull of both, and the residual bounds eps.
struct PartBounds {
    bool oriented = false; ///< the normals are proven to agree
    bool finite = false;
    std::array<double, 9> b{}; ///< gap, angle, db_tt, db_tn, db_nn, dC_ttt, dC_ttn, dC_tnn, dC_nnn
};

template <bool kRational, int N>
PartBounds part_bounds(const std::array<Vec4<Interval>, kJetSize>& HA, const std::array<Vec4<Interval>, kJetSize>& HH,
                       const std::array<double, N>& eps, std::size_t tangent_entry) {
    constexpr int C = kRational ? 4 : 3;
    PartBounds r;
    // A's jet and frame over the part, and the seam frame (tau, nu).
    std::array<Vec3<Interval>, kJetSize> JA;
    if constexpr (kRational) JA = project_jet(HA);
    else
        for (int e = 0; e < kJetSize; ++e) JA[static_cast<std::size_t>(e)] = xyz(HA[static_cast<std::size_t>(e)]);
    const Frame<Interval> fa = frame(JA);
    const Vec3<Interval> T = JA[tangent_entry];
    Vec3<Interval> tau = T - fa.n * dot(fa.n, T);
    tau = tau / sqrt(dot(tau, tau));
    const Vec3<Interval> nu = cross(fa.n, tau);
    // B's (homogeneous) jet as the variables, over the hull of both.
    std::array<Vec4<Dual<N>>, kJetSize> HB;
    for (int e = 0; e < kJetSize; ++e) {
        const auto ue = static_cast<std::size_t>(e);
        Vec4<Dual<N>>& v = HB[ue];
        v = {Dual<N>(HH[ue].x), Dual<N>(HH[ue].y), Dual<N>(HH[ue].z), Dual<N>(HH[ue].w)};
        v.x.d[static_cast<std::size_t>(C * e)] = Interval(1.0);
        v.y.d[static_cast<std::size_t>(C * e + 1)] = Interval(1.0);
        v.z.d[static_cast<std::size_t>(C * e + 2)] = Interval(1.0);
        if constexpr (kRational) v.w.d[static_cast<std::size_t>(C * e + 3)] = Interval(1.0);
    }
    std::array<Vec3<Dual<N>>, kJetSize> JB;
    if constexpr (kRational) JB = project_jet(HB);
    else
        for (int e = 0; e < kJetSize; ++e) JB[static_cast<std::size_t>(e)] = xyz(HB[static_cast<std::size_t>(e)]);
    const Frame<Dual<N>> fb = frame(JB);
    // Orientation: the normals must agree (no flipped parameterisation).
    const Interval c_nn = dot(fa.n, Vec3<Interval>{fb.n.x.v, fb.n.y.v, fb.n.z.v});
    r.oriented = c_nn.lo > 0.0;
    if (!r.oriented) return r;
    const auto tB = to_param<Dual<N>>(fb, tau), nB = to_param<Dual<N>>(fb, nu);
    // Each jump is (B's value as a Dual) - (A's value, a constant here), and
    // vanishes when B's jet is A's: it is bounded by |gradient| . eps.
    auto bound = [&](const Dual<N>& hb) {
        Interval s(0.0);
        for (int i = 0; i < N; ++i)
            s = s + Interval(num::mag(hb.d[static_cast<std::size_t>(i)])) * Interval(eps[static_cast<std::size_t>(i)]);
        return std::isfinite(s.hi) ? s.hi : std::numeric_limits<double>::infinity();
    };
    auto length = [](double x, double y, double z) {
        return sqrt_up(Interval(x) * Interval(x) + Interval(y) * Interval(y) + Interval(z) * Interval(z));
    };
    // Normal: |n_B - n_A| per coordinate, then the angle.
    const double dn = length(bound(fb.n.x), bound(fb.n.y), bound(fb.n.z));
    double angle = 3.14159265358979324; // the bound when nothing better is proven
    if (dn < 2.0) {
        const Interval q = Interval(1.0) - Interval(dn) * Interval(dn) / Interval(4.0);
        if (q.lo > 0.0) angle = std::min(angle, (Interval(dn) / num::sqrt(q)).hi);
    }
    r.b = {length(bound(JB[0].x), bound(JB[0].y), bound(JB[0].z)),
           angle,
           bound(form2(fb, tB, tB)),
           bound(form2(fb, tB, nB)),
           bound(form2(fb, nB, nB)),
           bound(form3(fb, tB, tB, tB)),
           bound(form3(fb, tB, tB, nB)),
           bound(form3(fb, tB, nB, nB)),
           bound(form3(fb, nB, nB, nB))};
    r.finite = std::isfinite(dn);
    for (double x : r.b) r.finite = r.finite && std::isfinite(x);
    return r;
}

/// The certificate, with the residuals taken in homogeneous coordinates
/// (kRational: 4 coordinates per jet entry, and the projection to the
/// surface is part of the differentiated function) or, for polynomial
/// patches, in the surface's own coordinates (3 per entry, weights exactly 1).
template <bool kRational>
SeamCertificate certify_core(const PolyPatchI& A, const PolyPatchI& B, int side, int pieces) {
    constexpr int C = kRational ? 4 : 3;
    constexpr int N = C * kJetSize;
    SeamCertificate out;
    // Jet entries along the seam, raised to a common degree per entry.
    std::array<Curve, kJetSize> ca, cb;
    for (int k = 0; k <= 3; ++k)
        for (int l = 0; l <= 3 - k; ++l) {
            const auto e = static_cast<std::size_t>(jet_index(k, l));
            ca[e] = jet_curve(A, k, l, side);
            cb[e] = jet_curve(B, k, l, side);
            const int n = static_cast<int>(std::max(ca[e].size(), cb[e].size())) - 1;
            ca[e] = elevate(ca[e], n);
            cb[e] = elevate(cb[e], n);
        }
    // The residual curves themselves, restricted to each part below. (Taking
    // the difference before restricting keeps rows the two patches share at
    // exactly zero, instead of the sum of two restriction roundings.)
    // For patches of equal degrees the difference is taken on the control
    // points, before any derivative, so that it is exactly zero where they
    // agree.
    std::array<Curve, kJetSize> cd;
    if (A.p == B.p && A.q == B.q) {
        PolyPatchI D(A.p, A.q);
        for (std::size_t i = 0; i < D.c.size(); ++i) D.c[i] = B.c[i] - A.c[i];
        for (int k = 0; k <= 3; ++k)
            for (int l = 0; l <= 3 - k; ++l) cd[static_cast<std::size_t>(jet_index(k, l))] = jet_curve(D, k, l, side);
    } else {
        for (int e = 0; e < kJetSize; ++e) {
            const auto ue = static_cast<std::size_t>(e);
            cd[ue].resize(ca[ue].size());
            for (std::size_t i = 0; i < ca[ue].size(); ++i) cd[ue][i] = cb[ue][i] - ca[ue][i];
        }
    }
    const bool along_u = side == kSeamV0 || side == kSeamV1;
    const auto tangent_entry = static_cast<std::size_t>(along_u ? jet_index(1, 0) : jet_index(0, 1));
    // Parts of the seam still to do, with their bisection depth. A part is
    // halved where something is not proven on it (the orientation, a weight
    // or the metric away from zero: its enclosures are too wide).
    constexpr int kMaxDepth = 10, kMaxParts = 4096;
    std::vector<std::array<double, 3>> work;
    for (int m = pieces - 1; m >= 0; --m)
        work.push_back({static_cast<double>(m) / pieces, static_cast<double>(m + 1) / pieces, 0.0});
    int done = 0, started = pieces;
    while (!work.empty()) {
        const double t0 = work.back()[0], t1 = work.back()[1];
        const int depth = static_cast<int>(work.back()[2]);
        work.pop_back();
        const bool may_split = depth < kMaxDepth && started < kMaxParts;
        auto split = [&]() {
            const double tm = 0.5 * (t0 + t1);
            work.push_back({tm, t1, depth + 1.0});
            work.push_back({t0, tm, depth + 1.0});
            ++started;
        };
        std::array<Vec4<Interval>, kJetSize> HA, HH; // A's box over the part, and the hull of both
        std::array<double, N> eps{};
        for (int e = 0; e < kJetSize; ++e) {
            const auto ue = static_cast<std::size_t>(e);
            const Curve a = restrict_curve(ca[ue], t0, t1);
            const Curve b = restrict_curve(cb[ue], t0, t1);
            const Curve d = restrict_curve(cd[ue], t0, t1);
            Interval box_a[4], box_h[4];
            for (int c = 0; c < 4; ++c) {
                // Polynomial: the weight is 1 and its derivatives 0, exactly.
                if (!kRational && c == 3) {
                    box_a[c] = box_h[c] = Interval(e == 0 ? 1.0 : 0.0);
                    continue;
                }
                box_a[c] = hull_of(a, c);
                box_h[c] = num::hull(box_a[c], hull_of(b, c));
                if (c >= C) continue;
                double r = 0.0;
                for (std::size_t i = 0; i < d.size(); ++i) r = std::max(r, num::mag(coord(d[i], c)));
                eps[static_cast<std::size_t>(C * e + c)] = r;
            }
            HA[ue] = {box_a[0], box_a[1], box_a[2], box_a[3]};
            HH[ue] = {box_h[0], box_h[1], box_h[2], box_h[3]};
        }
        const PartBounds r = part_bounds<kRational, N>(HA, HH, eps, tangent_entry);
        if (!r.oriented || !r.finite) {
            if (!may_split) {
                out.why = std::string(!r.oriented ? "the normals are not proven to agree in orientation"
                                                  : "a bound is not finite (a weight or the metric is not proven away from zero)") +
                          " on the part [" + std::to_string(t0) + ", " + std::to_string(t1) + "] of the seam";
                return out;
            }
            split();
            continue;
        }
        out.gap = std::max(out.gap, r.b[0]);
        out.normal_angle = std::max(out.normal_angle, r.b[1]);
        out.db_tt = std::max(out.db_tt, r.b[2]);
        out.db_tn = std::max(out.db_tn, r.b[3]);
        out.db_nn = std::max(out.db_nn, r.b[4]);
        out.dC_ttt = std::max(out.dC_ttt, r.b[5]);
        out.dC_ttn = std::max(out.dC_ttn, r.b[6]);
        out.dC_tnn = std::max(out.dC_tnn, r.b[7]);
        out.dC_nnn = std::max(out.dC_nnn, r.b[8]);
        ++done;
    }
    out.pieces = done;
    out.ok = true;
    return out;
}

} // namespace seam_detail

/// Bounds for the seam `side` of A and B, which share the parameterisation.
/// Polynomial patches (every weight exactly 1) are compared in their own
/// coordinates; rational ones in homogeneous coordinates, which must then
/// agree along the seam up to the residual (as for pieces of one rational
/// parameterisation).
/// `pieces` equal parts: more parts give tighter bounds (the enclosure of
/// the gradient over a part shrinks with it), at proportional cost.
inline SeamCertificate certify_seam(const PolyPatchI& A, const PolyPatchI& B, int side, int pieces = 32) {
    bool polynomial = true;
    for (const PolyPatchI* P : {&A, &B})
        for (const auto& cp : P->c) polynomial = polynomial && cp.w.lo == 1.0 && cp.w.hi == 1.0;
    return polynomial ? seam_detail::certify_core<false>(A, B, side, pieces)
                      : seam_detail::certify_core<true>(A, B, side, pieces);
}

inline SeamCertificate certify_seam(const PolyPatchD& A, const PolyPatchD& B, int side, int pieces = 32) {
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
