// akernel — certified enclosure of the volume of a Brep.
//
// op::volume (op/brep.hpp) measures the volume by Gauss quadrature in
// doubles: a measurement, not a bound (issue #5). Here the same divergence-
// theorem integral is evaluated so that the result is an interval proven to
// contain the exact volume of the solid the Brep describes.
//
// Volume = (1/3) * sum over faces of the flux of S . (Su x Sv) over the
// trimmed parameter domain. On a polynomial patch (all weights equal) the
// density f = S . (Su x Sv) is a polynomial of degree (3p - 1, 3q - 1); with
// G(u, v) = integral of f from 0 to u, Green's theorem turns the flux into the
// loop integral of G dv. Along a straight segment of the parameter plane
// G dv is a polynomial in the segment parameter, integrated exactly. Every
// coefficient is an interval and every vertex enters by its certified
// enclosure, so the sum encloses the exact volume. No quadrature, no error
// term: the only widening is rounding.
//
// Covered: faces on polynomial carriers whose loops are made of patch sides
// and straight edges in flat faces, including intersection edges of two flat
// carriers (boxes, boxes with polynomial tops, and every Boolean of
// flat-faced solids). Not covered yet, so reported as uncertified: rational
// carriers, intersection curves on a curved carrier, and procedural faces. Those need a remainder bound for the quadrature
// (BLUEPRINT §5.4).
#pragma once

#include <optional>
#include <string>
#include <vector>

#include <akernel/num/interval.hpp>
#include <akernel/op/brep.hpp>

namespace ak::op {

namespace cvol_detail {

/// A bivariate polynomial sum a[i][j] u^i v^j with interval coefficients.
struct Poly2 {
    int du = 0, dv = 0; ///< degrees
    std::vector<Interval> a;
    Poly2() = default;
    Poly2(int du_, int dv_) : du(du_), dv(dv_), a(static_cast<std::size_t>((du_ + 1) * (dv_ + 1)), Interval(0.0)) {}
    Interval& at(int i, int j) { return a[static_cast<std::size_t>(i * (dv + 1) + j)]; }
    const Interval& at(int i, int j) const { return a[static_cast<std::size_t>(i * (dv + 1) + j)]; }
};

inline Poly2 operator*(const Poly2& x, const Poly2& y) {
    Poly2 r(x.du + y.du, x.dv + y.dv);
    for (int i = 0; i <= x.du; ++i)
        for (int j = 0; j <= x.dv; ++j)
            for (int k = 0; k <= y.du; ++k)
                for (int l = 0; l <= y.dv; ++l) r.at(i + k, j + l) = r.at(i + k, j + l) + x.at(i, j) * y.at(k, l);
    return r;
}

inline Poly2 operator+(const Poly2& x, const Poly2& y) {
    Poly2 r(std::max(x.du, y.du), std::max(x.dv, y.dv));
    for (int i = 0; i <= x.du; ++i)
        for (int j = 0; j <= x.dv; ++j) r.at(i, j) = r.at(i, j) + x.at(i, j);
    for (int i = 0; i <= y.du; ++i)
        for (int j = 0; j <= y.dv; ++j) r.at(i, j) = r.at(i, j) + y.at(i, j);
    return r;
}

inline Poly2 operator-(const Poly2& x, const Poly2& y) {
    Poly2 n = y;
    for (Interval& c : n.a) c = -c;
    return x + n;
}

inline Poly2 d_du(const Poly2& x) {
    if (x.du == 0) return Poly2(0, x.dv);
    Poly2 r(x.du - 1, x.dv);
    for (int i = 1; i <= x.du; ++i)
        for (int j = 0; j <= x.dv; ++j) r.at(i - 1, j) = x.at(i, j) * Interval(static_cast<double>(i));
    return r;
}

inline Poly2 d_dv(const Poly2& x) {
    if (x.dv == 0) return Poly2(x.du, 0);
    Poly2 r(x.du, x.dv - 1);
    for (int i = 0; i <= x.du; ++i)
        for (int j = 1; j <= x.dv; ++j) r.at(i, j - 1) = x.at(i, j) * Interval(static_cast<double>(j));
    return r;
}

/// Integral in u from 0: a[i][j] u^(i+1) / (i+1) v^j.
inline Poly2 int_du(const Poly2& x) {
    Poly2 r(x.du + 1, x.dv);
    for (int i = 0; i <= x.du; ++i)
        for (int j = 0; j <= x.dv; ++j) r.at(i + 1, j) = x.at(i, j) / Interval(static_cast<double>(i + 1));
    return r;
}

/// Binomial coefficient, exact in doubles for the degrees used (n <= 30).
inline double binom(int n, int k) {
    double c = 1.0;
    for (int i = 1; i <= k; ++i) c = c * (n - k + i) / i;
    return c;
}

/// The monomial form of one coordinate of a polynomial Bernstein patch,
/// divided by the common weight w.
inline Poly2 monomial(const PolyPatchD& P, int coord, const Interval& w) {
    Poly2 r(P.p, P.q);
    // B_i^p(u) = sum over k >= i of (-1)^(k-i) C(p,k) C(k,i) u^k.
    for (int i = 0; i <= P.p; ++i)
        for (int j = 0; j <= P.q; ++j) {
            const Vec4d& c = P.at(i, j);
            const double x = coord == 0 ? c.x : coord == 1 ? c.y : c.z;
            for (int k = i; k <= P.p; ++k)
                for (int l = j; l <= P.q; ++l) {
                    const double s = ((k - i + l - j) % 2 ? -1.0 : 1.0) * binom(P.p, k) * binom(k, i) * binom(P.q, l) * binom(l, j);
                    r.at(k, l) = r.at(k, l) + Interval(x) * Interval(s);
                }
        }
    for (Interval& a : r.a) a = a / w;
    return r;
}

/// Integral over s in [0, 1] of G(ua + s (ub - ua), va + s (vb - va)) (vb - va).
inline Interval segment_integral(const Poly2& G, const Interval& ua, const Interval& va, const Interval& ub, const Interval& vb) {
    const Interval du = ub - ua, dv = vb - va;
    // Powers of (ua + s du) and (va + s dv) as polynomials in s.
    auto powers = [](const Interval& a, const Interval& d, int n) {
        std::vector<std::vector<Interval>> P(static_cast<std::size_t>(n + 1));
        P[0] = {Interval(1.0)};
        for (int k = 1; k <= n; ++k) {
            const auto& prev = P[static_cast<std::size_t>(k - 1)];
            std::vector<Interval> cur(prev.size() + 1, Interval(0.0));
            for (std::size_t m = 0; m < prev.size(); ++m) {
                cur[m] = cur[m] + prev[m] * a;
                cur[m + 1] = cur[m + 1] + prev[m] * d;
            }
            P[static_cast<std::size_t>(k)] = std::move(cur);
        }
        return P;
    };
    const auto PU = powers(ua, du, G.du), PV = powers(va, dv, G.dv);
    Interval total(0.0);
    for (int i = 0; i <= G.du; ++i)
        for (int j = 0; j <= G.dv; ++j) {
            const Interval& c = G.at(i, j);
            if (c.lo == 0.0 && c.hi == 0.0) continue;
            const auto& pu = PU[static_cast<std::size_t>(i)];
            const auto& pv = PV[static_cast<std::size_t>(j)];
            // Integral of pu(s) pv(s) over [0, 1].
            Interval acc(0.0);
            for (std::size_t m = 0; m < pu.size(); ++m)
                for (std::size_t n = 0; n < pv.size(); ++n)
                    acc = acc + pu[m] * pv[n] / Interval(static_cast<double>(m + n + 1));
            total = total + c * acc;
        }
    return total * dv;
}

} // namespace cvol_detail

/// A certified enclosure of the volume, or nothing (with the reason) when
/// some face is outside the covered class (see the file comment).
inline std::optional<Interval> certified_volume(const Brep& b, std::string* why = nullptr) {
    using namespace cvol_detail;
    auto fail = [&](const char* w) -> std::optional<Interval> {
        if (why) *why = w;
        return std::nullopt;
    };
    Interval total(0.0);
    for (const BFace& face : b.faces) {
        if (face.proc >= 0) return fail("a procedural face");
        const PolyPatchD& P = b.carriers[static_cast<std::size_t>(face.carrier)];
        for (const Vec4d& c : P.c)
            if (c.w != P.c.front().w) return fail("a rational carrier");
        const Interval w(P.c.front().w);
        const Poly2 X = monomial(P, 0, w), Y = monomial(P, 1, w), Z = monomial(P, 2, w);
        const Poly2 Xu = d_du(X), Yu = d_du(Y), Zu = d_du(Z), Xv = d_dv(X), Yv = d_dv(Y), Zv = d_dv(Z);
        // f = S . (Su x Sv)
        const Poly2 f = X * (Yu * Zv - Zu * Yv) + Y * (Zu * Xv - Xu * Zv) + Z * (Xu * Yv - Yu * Xv);
        const Poly2 G = int_du(f);
        Interval flux(0.0);
        for (const auto& loop : face.loops)
            for (const BCoedge& ce : loop) {
                const BEdge& e = b.edges[static_cast<std::size_t>(ce.edge)];
                if (e.proc) return fail("an edge on a procedural face");
                // An intersection curve of two flat carriers is a straight
                // line, and so is its image in an affine flat face: it is
                // handled like any straight edge below. On a curved carrier
                // the traced curve would need a remainder bound.
                if (e.arc && !(planar_affine(b.carriers[static_cast<std::size_t>(e.carrier[0])]) &&
                               planar_affine(b.carriers[static_cast<std::size_t>(e.carrier[1])])))
                    return fail("an edge along a traced intersection curve");
                Interval ua, va, ub, vb;
                const BEdge::SideRef* ref = e.arc ? nullptr : side_on(e, face.carrier);
                if (ref) {
                    const Interval ta = ce.reversed ? e.t1 : e.t0, tb = ce.reversed ? e.t0 : e.t1;
                    side_uv(ref->side, ref->reversed ? Interval(1.0) - ta : ta, ua, va);
                    side_uv(ref->side, ref->reversed ? Interval(1.0) - tb : tb, ub, vb);
                } else {
                    // A straight edge inside a flat face: a segment between the
                    // vertices' certified parameters.
                    if (!planar_affine(P)) return fail("an edge inside a curved face");
                    const Chart* a = b.vertices[static_cast<std::size_t>(ce.reversed ? e.v1 : e.v0)].chart(face.carrier);
                    const Chart* c = b.vertices[static_cast<std::size_t>(ce.reversed ? e.v0 : e.v1)].chart(face.carrier);
                    if (!a || !c) return fail("a vertex without parameters on its face");
                    ua = a->u; va = a->v;
                    ub = c->u; vb = c->v;
                }
                flux = flux + segment_integral(G, ua, va, ub, vb);
            }
        total = total + (face.flipped ? -flux : flux);
    }
    return total / Interval(3.0);
}

} // namespace ak::op
