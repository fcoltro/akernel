// akernel — certified local offset validity and curvature bounds for a
// rational Bezier patch (docs/BLUEPRINT.md §1.3; docs/FEATURES.md, step F2).
//
// The offset S_d = S + d n of a regular surface is regular exactly where
// 1 - d kappa_i > 0 for both principal curvatures (with b = n . S_ab, so a
// sphere with outward normal has kappa = -1/R). Offsets, shells and
// rolling-ball blends share this predicate: the spine of a ball of radius r
// rolling on a face lies on the face's offset at distance r, which must be
// regular there (a ball that does not fit in a concave region gives a
// self-intersecting blend, the classic "blend too large" failure).
//
// Invariant form. kappa_i are the eigenvalues of g^-1 b (g positive
// definite, b symmetric, so they are real), and 1 - d kappa_i are those of
// g^-1 (g - d b). Both are positive iff
//     det(g - d b) > 0   and   tr(g^-1 (g - d b)) > 0.
// With N = S_u x S_v, L = |N| and the unnormalised b~ = N . S_ab = L b,
// multiply through by positive factors:
//     A = L g - d b~,    det A > 0,    g22 A11 + g11 A22 - 2 g12 A12 > 0.
// Everything is polynomial in the jet except L = sqrt(N . N).
//
// Certification. On a parameter box every partial up to order 3 of the
// homogeneous form is enclosed by the hull of its Bezier coefficients on the
// box (exact interval restriction, convex-hull property). The quantities
// are enclosed twice: directly from those hulls, and by the mean-value form
// q(c) + grad q(box) . (x - c), with q(c) evaluated at the box centre and
// grad q enclosed by forward differentiation over the box. The direct form
// alone suffers badly from dependency in the rational quotient rule; the
// mean-value form shrinks like O(h) |grad q|. The intersection is used. If
// both expressions are proven positive the offset is valid on the box; if
// at the box centre (a point interval) one is proven <= 0, the offset is
// invalid there and that point is the witness. Otherwise the box is split
// in four, down to a minimum size, beyond which the answer is Undecided:
// never a guess.
#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <vector>

#include <akernel/geom/bernstein.hpp>
#include <akernel/num/interval.hpp>

namespace ak::geom {

enum class OffsetVerdict { Valid, Invalid, Undecided };

struct OffsetCheck {
    OffsetVerdict verdict = OffsetVerdict::Undecided;
    int boxes = 0;            ///< parameter boxes examined
    double u = 0.0, v = 0.0;  ///< witness (Invalid) or unresolved box centre (Undecided)
};

namespace offset_detail {

using num::Interval;

/// Interval value with interval partial derivatives in (u, v): forward-mode
/// differentiation used for the mean-value form below.
struct D2 {
    Interval v{}, du{}, dv{};
    D2() = default;
    D2(double x) : v(x) {}
    D2(const Interval& x) : v(x) {}
    D2(const Interval& x, const Interval& a, const Interval& b) : v(x), du(a), dv(b) {}
};
inline D2 operator+(const D2& a, const D2& b) { return {a.v + b.v, a.du + b.du, a.dv + b.dv}; }
inline D2 operator-(const D2& a, const D2& b) { return {a.v - b.v, a.du - b.du, a.dv - b.dv}; }
inline D2 operator-(const D2& a) { return {-a.v, -a.du, -a.dv}; }
inline D2 operator*(const D2& a, const D2& b) {
    return {a.v * b.v, a.du * b.v + a.v * b.du, a.dv * b.v + a.v * b.dv};
}
inline D2 operator/(const D2& a, const D2& b) {
    const Interval q = a.v / b.v;
    return {q, (a.du - q * b.du) / b.v, (a.dv - q * b.dv) / b.v};
}
inline D2 sqrt_(const D2& a) {
    const Interval r = num::sqrt(a.v), t = Interval(2.0) * r;
    return {r, a.du / t, a.dv / t};
}
inline Interval sqrt_(const Interval& a) { return num::sqrt(a); }
inline D2 sqr_(const D2& a) { return a * a; }
inline Interval sqr_(const Interval& a) { return num::sqr(a); }
inline Interval value(const Interval& a) { return a; }
inline Interval value(const D2& a) { return a.v; }

/// Quantities decided on a box, in this order.
enum Q : int { kNN = 0, kDet, kTr, kH, kDisc, kQCount };

/// Homogeneous partials of order <= 3 as interval patches, indexed by jet_index.
struct Partials {
    std::array<PolyPatchI, kJetSize> d;
    Partials() = default;
    explicit Partials(const PolyPatchI& p) {
        d[0] = p;
        for (int n = 1; n <= 3; ++n)
            for (int l = 0; l <= n; ++l) {
                const int k = n - l;
                d[jet_index(k, l)] = k > 0 ? d[jet_index(k - 1, l)].deriv_u()
                                           : d[jet_index(k, l - 1)].deriv_v();
            }
    }
};

/// From homogeneous partials X[jet_index(k,l)], k + l <= 2, compute
/// |N|^2, det A, the trace expression, H and H^2 - K. `ok` is false when
/// the weight is not proven positive.
template <class T>
std::array<T, kQCount> quantities(const std::array<Vec4<T>, 6>& X, double d, bool& ok) {
    using V = Vec3<T>;
    std::array<T, kQCount> q{};
    ok = value(X[0].w).lo > 0.0;
    if (!ok) return q;
    const T w = X[0].w, two(2.0);
    const V s = xyz(X[0]) / w;
    const V su = (xyz(X[1]) - s * X[1].w) / w;
    const V sv = (xyz(X[2]) - s * X[2].w) / w;
    const V suu = (xyz(X[3]) - su * (two * X[1].w) - s * X[3].w) / w;
    const V suv = (xyz(X[4]) - su * X[2].w - sv * X[1].w - s * X[4].w) / w;
    const V svv = (xyz(X[5]) - sv * (two * X[2].w) - s * X[5].w) / w;
    const V N = cross(su, sv);
    const T G = dot(N, N), L = sqrt_(G);
    const T g11 = dot(su, su), g12 = dot(su, sv), g22 = dot(sv, sv);
    const T b11 = dot(N, suu), b12 = dot(N, suv), b22 = dot(N, svv);
    const T D(d);
    // A = L g - d b~;  det A and g22 A11 + g11 A22 - 2 g12 A12 = 2 L G - d T.
    const T a11 = L * g11 - D * b11, a12 = L * g12 - D * b12, a22 = L * g22 - D * b22;
    const T Tr = g22 * b11 + g11 * b22 - two * g12 * b12;
    q[kNN] = G;
    q[kDet] = a11 * a22 - sqr_(a12);
    q[kTr] = two * L * G - D * Tr;
    // H = T / (2 G L), K = det b~ / G^2.
    q[kH] = Tr / (two * G * L);
    q[kDisc] = sqr_(q[kH]) - (b11 * b22 - sqr_(b12)) / sqr_(G);
    return q;
}

/// The four quarters of a box's restricted partials, by exact halving
/// (de Casteljau at 1/2 in interval arithmetic). Order: (lo u, lo v),
/// (hi u, lo v), (lo u, hi v), (hi u, hi v).
inline std::array<Partials, 4> quarters(const Partials& B) {
    std::array<Partials, 4> r;
    const Interval h(0.5);
    for (int i = 0; i < kJetSize; ++i) {
        PolyPatchI lo, hi;
        B.d[i].split_u(h, lo, hi);
        lo.split_v(h, r[0].d[i], r[2].d[i]);
        hi.split_v(h, r[1].d[i], r[3].d[i]);
    }
    return r;
}

inline Interval meet(const Interval& a, const Interval& b) {
    return {std::max(a.lo, b.lo), std::min(a.hi, b.hi)};
}

/// Enclosures of the quantities over the parameter box, by the mean-value
/// form q(c) + q_u(box) (u - c_u) + q_v(box) (v - c_v), intersected with
/// the direct enclosure. The direct enclosure uses the hulls of the
/// partials' Bezier coefficients on the box (convex-hull property) and
/// suffers from dependency; the mean-value form has width O(h) |grad q|,
/// which is what makes boxes near a focal point decidable.
/// `B` holds the partials restricted to the box (values unscaled: the
/// restriction only reparameterises), so their coefficient hulls enclose
/// the global partials over the box.
inline std::array<Interval, kQCount> enclose_box(const Partials& P, const Partials& B, double u0,
                                                 double u1, double v0, double v1, double d,
                                                 bool& ok) {
    std::array<Vec4<Interval>, kJetSize> H;
    for (int i = 0; i < kJetSize; ++i) H[i] = hull(B.d[i]);
    std::array<Vec4<Interval>, 6> Xi;
    std::array<Vec4<D2>, 6> Xd;
    for (int k = 0; k <= 2; ++k)
        for (int l = 0; l <= 2 - k; ++l) {
            const int i = jet_index(k, l), iu = jet_index(k + 1, l), iv = jet_index(k, l + 1);
            Xi[i] = H[i];
            Xd[i] = {D2(H[i].x, H[iu].x, H[iv].x), D2(H[i].y, H[iu].y, H[iv].y),
                     D2(H[i].z, H[iu].z, H[iv].z), D2(H[i].w, H[iu].w, H[iv].w)};
        }
    bool ok1 = false, ok2 = false, ok3 = false;
    const auto direct = quantities(Xi, d, ok1);
    const auto grad = quantities(Xd, d, ok2);
    const double uc = 0.5 * u0 + 0.5 * u1, vc = 0.5 * v0 + 0.5 * v1;
    std::array<Vec4<Interval>, 6> Xc;
    for (int i = 0; i < 6; ++i) Xc[i] = P.d[i].eval(Interval(uc), Interval(vc));
    const auto centre = quantities(Xc, d, ok3);
    ok = ok1 && ok2 && ok3;
    std::array<Interval, kQCount> out{};
    if (!ok) return out;
    const Interval du = Interval(u0) - Interval(uc), dv = Interval(v0) - Interval(vc);
    const Interval eu = num::hull(du, -du), ev = num::hull(dv, -dv);
    for (int i = 0; i < kQCount; ++i) {
        Interval mv = centre[i] + grad[i].du * eu + grad[i].dv * ev;
        if (!num::is_finite(mv)) mv = Interval::entire();
        out[i] = meet(direct[i], mv);
        if (out[i].lo > out[i].hi) out[i] = num::hull(direct[i], mv); // cannot happen; stay safe
    }
    return out;
}

/// The quantities at one point (point intervals), for witnesses.
inline std::array<Interval, kQCount> at_point(const Partials& P, double u, double v, double d,
                                              bool& ok) {
    std::array<Vec4<Interval>, 6> X;
    for (int i = 0; i < 6; ++i) X[i] = P.d[i].eval(Interval(u), Interval(v));
    return quantities(X, d, ok);
}

/// Upper bound on max |kappa_i| from enclosures of H and H^2 - K.
inline double curvature_magnitude_bound(const std::array<Interval, kQCount>& q) {
    if (!(q[kNN].lo > 0.0) || !num::is_finite(q[kH]) || !num::is_finite(q[kDisc]))
        return std::numeric_limits<double>::infinity();
    const Interval disc = q[kDisc]; // = ((k1 - k2) / 2)^2 >= 0 in exact terms
    const double root = disc.hi <= 0.0 ? 0.0 : num::upper(num::sqrt(Interval(0.0, disc.hi)));
    const double b = num::detail::step_up(num::mag(q[kH]) + root);
    return std::isfinite(b) ? b : std::numeric_limits<double>::infinity();
}

} // namespace offset_detail

/// Is the offset S + d n of `patch` over [0,1]^2 regular everywhere
/// (1 - d kappa_i > 0 for both principal curvatures)? Boxes smaller than
/// 2^-max_depth in a side, or more than max_boxes boxes in total, give
/// Undecided. A patch that is not proven regular (|S_u x S_v| > 0) on a
/// small box is Undecided there too: a singular point has no offset normal.
inline OffsetCheck check_offset(const PolyPatchD& patch, double d, int max_depth = 14,
                                int max_boxes = 200000) {
    using namespace offset_detail;
    struct Box { double u0, u1, v0, v1; int depth; Partials B; };
    const Partials P(to_interval(patch));
    OffsetCheck out;
    std::vector<Box> stack;
    stack.push_back({0.0, 1.0, 0.0, 1.0, 0, P});
    while (!stack.empty()) {
        const Box b = std::move(stack.back());
        stack.pop_back();
        const double um = 0.5 * (b.u0 + b.u1), vm = 0.5 * (b.v0 + b.v1);
        if (++out.boxes > max_boxes) {
            out.verdict = OffsetVerdict::Undecided;
            out.u = um;
            out.v = vm;
            return out;
        }
        bool ok = false;
        const auto q = enclose_box(P, b.B, b.u0, b.u1, b.v0, b.v1, d, ok);
        if (ok && q[kNN].lo > 0.0 && q[kDet].lo > 0.0 && q[kTr].lo > 0.0) continue; // valid here
        const auto c = at_point(P, um, vm, d, ok);
        if (ok && c[kNN].lo > 0.0 && (c[kDet].hi <= 0.0 || c[kTr].hi <= 0.0)) {
            out.verdict = OffsetVerdict::Invalid;
            out.u = um;
            out.v = vm;
            return out;
        }
        if (b.depth >= max_depth) {
            out.verdict = OffsetVerdict::Undecided;
            out.u = um;
            out.v = vm;
            return out;
        }
        auto q4 = quarters(b.B);
        stack.push_back({b.u0, um, b.v0, vm, b.depth + 1, std::move(q4[0])});
        stack.push_back({um, b.u1, b.v0, vm, b.depth + 1, std::move(q4[1])});
        stack.push_back({b.u0, um, vm, b.v1, b.depth + 1, std::move(q4[2])});
        stack.push_back({um, b.u1, vm, b.v1, b.depth + 1, std::move(q4[3])});
    }
    out.verdict = OffsetVerdict::Valid;
    return out;
}

/// Certified upper bound on max |kappa_i| over the patch, from a uniform
/// grid of 2^depth x 2^depth boxes. Every offset with |d| below 1 / bound is
/// locally valid, on either side. Infinite if regularity is not proven on
/// some box. The bound converges to the true maximum as depth grows.
inline double curvature_bound(const PolyPatchD& patch, int depth = 5) {
    using namespace offset_detail;
    const Partials P(to_interval(patch));
    const int n = 1 << depth;
    double best = 0.0;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) {
            const double u0 = static_cast<double>(i) / n, u1 = static_cast<double>(i + 1) / n;
            const double v0 = static_cast<double>(j) / n, v1 = static_cast<double>(j + 1) / n;
            Partials B;
            for (int k = 0; k < kJetSize; ++k) B.d[k] = P.d[k].restricted(u0, u1, v0, v1);
            bool ok = false;
            const auto q = enclose_box(P, B, u0, u1, v0, v1, 0.0, ok);
            if (!ok) return std::numeric_limits<double>::infinity();
            best = std::max(best, curvature_magnitude_bound(q));
            if (!std::isfinite(best)) return best;
        }
    return best;
}

} // namespace ak::geom
