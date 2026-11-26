// akernel — certified tangent continuity (G1) across the edges of a Brep.
//
// geom/seam_certificate.hpp bounds the continuity jumps across a seam between
// two patches in one common parameterisation (the variational fill). The
// edges of a Brep are different: the two faces along a side edge are
// parameterised independently (a fillet against the plate it is tangent to,
// a cylinder against its cap), and only the curve is shared. G1 depends on
// the point alone, so it can be certified there without matching
// parameterisations.
//
// The angle between the outward normals of the two faces along an edge is
// enclosed over the whole edge:
//  1. Along a side the homogeneous point H, and its partials Hu and Hv, are
//     polynomials in the side parameter. The products
//     Pu = xyz(Hu) w - xyz(H) w_u (and Pv likewise) are w^2 times the true
//     partials, so N = Pu x Pv is a positive multiple of Su x Sv. All of them
//     are formed exactly as polynomials, in Bernstein form with interval
//     coefficients: no division, and no dependency between factors.
//  2. With the edge parameter t mapped to each carrier's side parameter
//     (reversed where the side runs the other way), the cross product
//     X = NA x NB and the dot product D = NA . NB of the two faces' outward
//     normals are polynomials in t. On each part of the edge the hull of
//     their Bernstein coefficients encloses them (convex hull property).
//  3. The angle is atan2(|X|, D), which does not depend on the scale of NA
//     and NB. Over a box of (|X|, D) that does not contain the origin, its
//     extremes are at the corners. A part whose box does contain the origin
//     (a normal not proven nonzero) is halved, up to a depth limit.
// Because X is formed as a polynomial before it is bounded, an edge where
// the faces are tangent by construction gets a bound at the level of
// rounding (about 1e-15 rad), however curved the faces are. A sharp edge gets
// a positive lower bound: proven not G1.
//
// The angle uses num::det::atan2, which is within 2 ulp of the true value
// (fdlibm's bound is below 1 ulp; tests/num/detmath_test.cpp); each end is
// moved outward by 2^-48 relative, plus one smallest normal number.
//
// Covered: side edges between two faces on polynomial or rational carriers.
// Refused, with a reason: edges along traced intersection curves, edges on
// procedural faces, and edges that do not have exactly two faces.
#pragma once

#include <algorithm>
#include <limits>
#include <numbers>
#include <string>
#include <vector>

#include <akernel/geom/seam_certificate.hpp>
#include <akernel/num/detmath.hpp>
#include <akernel/num/interval.hpp>
#include <akernel/op/brep.hpp>

namespace ak::op {

/// Enclosure of the angle between the outward normals of the two faces along
/// an edge, over the whole edge: 0 where the faces are tangent (G1), pi/2
/// for a box edge, pi where the surface folds back on itself.
struct EdgeAngle {
    bool ok = false;
    std::string why;
    double lo = 0.0, hi = 0.0; ///< radians, lo <= angle <= hi at every point of the edge
    int parts = 0;              ///< parts the edge was split into
};

namespace g1_detail {

using Poly = std::vector<Interval>; ///< Bernstein coefficients in [0, 1]

inline double binom(int n, int k) {
    double c = 1.0;
    for (int i = 1; i <= k; ++i) c = c * (n - k + i) / i;
    return c;
}

/// Product of two Bernstein polynomials, exactly up to rounding.
inline Poly mul(const Poly& a, const Poly& b) {
    const int m = static_cast<int>(a.size()) - 1, n = static_cast<int>(b.size()) - 1;
    Poly r(static_cast<std::size_t>(m + n + 1), Interval(0.0));
    for (int i = 0; i <= m; ++i)
        for (int j = 0; j <= n; ++j)
            r[static_cast<std::size_t>(i + j)] +=
                a[static_cast<std::size_t>(i)] * b[static_cast<std::size_t>(j)] *
                (Interval(binom(m, i)) * Interval(binom(n, j)) / Interval(binom(m + n, i + j)));
    return r;
}

/// Degree elevation, so that two polynomials can be added.
inline Poly elevate(Poly c, int n) {
    while (static_cast<int>(c.size()) - 1 < n) {
        const int m = static_cast<int>(c.size());
        Poly e(static_cast<std::size_t>(m + 1));
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

inline Poly add(Poly a, Poly b, double sb = 1.0) {
    const int n = static_cast<int>(std::max(a.size(), b.size())) - 1;
    a = elevate(std::move(a), n);
    b = elevate(std::move(b), n);
    for (std::size_t i = 0; i < a.size(); ++i) a[i] = a[i] + b[i] * Interval(sb);
    return a;
}

inline Interval hull(const Poly& p) {
    Interval h = p[0];
    for (const Interval& c : p) h = num::hull(h, c);
    return h;
}

struct PVec { Poly x, y, z; };

inline PVec cross(const PVec& a, const PVec& b) {
    return {add(mul(a.y, b.z), mul(a.z, b.y), -1.0), add(mul(a.z, b.x), mul(a.x, b.z), -1.0),
            add(mul(a.x, b.y), mul(a.y, b.x), -1.0)};
}
inline Poly dot(const PVec& a, const PVec& b) { return add(add(mul(a.x, b.x), mul(a.y, b.y)), mul(a.z, b.z)); }

/// brep side numbering (patch_solid.hpp: 0 v = 0, 1 u = 1, 2 v = 1, 3 u = 0)
/// to the seam numbering of seam_certificate.hpp.
inline int seam_side(int side) {
    static constexpr int map[4] = {geom::kSeamV0, geom::kSeamU1, geom::kSeamV1, geom::kSeamU0};
    return map[side];
}

/// One face along an edge: its homogeneous jet curves along the side, and
/// how the edge parameter maps to the side parameter.
struct FaceSide {
    geom::seam_detail::Curve H, Hu, Hv;
    bool reversed = false;
    double sign = 1.0; ///< -1 for a flipped face: outward normal -(Su x Sv)
};

/// The outward normal, up to a positive factor, on edge parameters [a, b].
inline PVec normal_on(const FaceSide& f, double a, double b) {
    using geom::seam_detail::restrict_curve;
    // Coordinates and weight of a jet curve on the part, in edge order.
    auto piece = [&](const geom::seam_detail::Curve& c, PVec& v, Poly& wt) {
        auto r = f.reversed ? restrict_curve(c, 1.0 - b, 1.0 - a) : restrict_curve(c, a, b);
        if (f.reversed) std::reverse(r.begin(), r.end());
        for (const auto& q : r) {
            v.x.push_back(q.x);
            v.y.push_back(q.y);
            v.z.push_back(q.z);
            wt.push_back(q.w);
        }
    };
    PVec X, Xu, Xv;
    Poly w, wu, wv;
    piece(f.H, X, w);
    piece(f.Hu, Xu, wu);
    piece(f.Hv, Xv, wv);
    auto partial = [&](const PVec& D, const Poly& wd) {
        return PVec{add(mul(D.x, w), mul(X.x, wd), -1.0), add(mul(D.y, w), mul(X.y, wd), -1.0),
                    add(mul(D.z, w), mul(X.z, wd), -1.0)};
    };
    PVec N = cross(partial(Xu, wu), partial(Xv, wv));
    if (f.sign < 0)
        for (Poly* p : {&N.x, &N.y, &N.z})
            for (Interval& c : *p) c = -c;
    return N;
}

/// det::atan2 moved outward to a rigorous bound.
inline double atan2_down(double y, double x) {
    const double a = num::det::atan2(y, x);
    return std::max(0.0, a - std::fabs(a) * 0x1p-48 - std::numeric_limits<double>::min());
}
inline double atan2_up(double y, double x) {
    const double a = num::det::atan2(y, x);
    return std::min(std::numbers::pi * (1.0 + 0x1p-48), a + std::fabs(a) * 0x1p-48 + std::numeric_limits<double>::min());
}

} // namespace g1_detail

/// The certified angle between the outward normals of the two faces along
/// side edge `edge` (see the file comment). `pieces` equal parts, each
/// halved up to `max_depth` times where a normal is not proven nonzero.
inline EdgeAngle certify_edge_angle(const Brep& b, int edge, int pieces = 8, int max_depth = 12) {
    using namespace g1_detail;
    EdgeAngle res;
    auto fail = [&](const char* w) {
        res.ok = false;
        res.why = w;
        return res;
    };
    const BEdge& e = b.edges[static_cast<std::size_t>(edge)];
    if (e.proc) return fail("an edge on a procedural face");
    if (e.arc) return fail("an edge along a traced intersection curve");
    // The faces that use the edge.
    std::vector<const BFace*> faces;
    for (const BFace& f : b.faces)
        for (const auto& loop : f.loops)
            for (const BCoedge& ce : loop)
                if (ce.edge == edge) faces.push_back(&f);
    if (faces.size() != 2) return fail("an edge without exactly two faces");
    FaceSide fs[2];
    for (int k = 0; k < 2; ++k) {
        const BFace& f = *faces[static_cast<std::size_t>(k)];
        if (f.proc >= 0) return fail("an edge on a procedural face");
        const BEdge::SideRef* ref = side_on(e, f.carrier);
        if (!ref) return fail("a face whose carrier does not hold the edge as a side");
        const geom::PolyPatchI P = geom::to_interval(b.carriers[static_cast<std::size_t>(f.carrier)]);
        for (const auto& c : P.c)
            if (!(c.w.lo > 0.0)) return fail("a carrier with a weight not proven positive");
        const int s = seam_side(ref->side);
        fs[k].H = geom::seam_detail::jet_curve(P, 0, 0, s);
        fs[k].Hu = geom::seam_detail::jet_curve(P, 1, 0, s);
        fs[k].Hv = geom::seam_detail::jet_curve(P, 0, 1, s);
        fs[k].reversed = ref->reversed;
        fs[k].sign = f.flipped ? -1.0 : 1.0;
    }
    // The edge's parameter range, widened to cover its enclosures.
    const double t0 = std::max(0.0, std::min(e.t0.lo, e.t1.lo)), t1 = std::min(1.0, std::max(e.t0.hi, e.t1.hi));
    res.lo = std::numeric_limits<double>::infinity();
    res.hi = 0.0;
    struct Part { double a, b; int depth; };
    std::vector<Part> todo;
    for (int i = pieces - 1; i >= 0; --i)
        todo.push_back({t0 + (t1 - t0) * i / pieces, i + 1 == pieces ? t1 : t0 + (t1 - t0) * (i + 1) / pieces, 0});
    while (!todo.empty()) {
        const Part p = todo.back();
        todo.pop_back();
        const PVec NA = normal_on(fs[0], p.a, p.b), NB = normal_on(fs[1], p.a, p.b);
        const PVec X = cross(NA, NB);
        const Interval D = hull(dot(NA, NB));
        // |X| enclosed from its components' hulls.
        Interval x2(0.0);
        for (const Poly* c : {&X.x, &X.y, &X.z}) x2 = x2 + num::sqr(hull(*c));
        const Interval xn = num::sqrt(Interval(std::max(0.0, x2.lo), std::max(0.0, x2.hi)));
        if (xn.lo <= 0.0 && D.lo <= 0.0 && D.hi >= 0.0) {
            if (p.depth >= max_depth) return fail("a normal not proven nonzero");
            const double m = 0.5 * (p.a + p.b);
            todo.push_back({m, p.b, p.depth + 1});
            todo.push_back({p.a, m, p.depth + 1});
            continue;
        }
        ++res.parts;
        // Extremes of atan2(y, x) over the box are at its corners.
        for (double yy : {xn.lo, xn.hi})
            for (double dd : {D.lo, D.hi}) {
                res.lo = std::min(res.lo, atan2_down(yy, dd));
                res.hi = std::max(res.hi, atan2_up(yy, dd));
            }
    }
    res.ok = true;
    return res;
}

} // namespace ak::op
