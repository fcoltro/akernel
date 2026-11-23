// akernel — rational B-spline curves in homogeneous form: evaluation, knot
// insertion, Bezier extraction, degree elevation and compatibility.
//
// Curves are the profiles, rails and sections of the generative operations in
// geom/sweep.hpp (extrude, revolve, translational sweep, ruled and skinned
// lofts). Like NurbsSurface, every operation here is done in doubles: the
// result differs from the exact refinement by rounding only, and from here on
// the result is the definition (thesis T1: no tolerance is attached). The
// operations that only copy or scale control points by powers of two are
// exact (the exact circle, revolution by quarter turns).
#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

#include <akernel/geom/bspline.hpp>
#include <akernel/geom/nurbs_surface.hpp>
#include <akernel/geom/vec.hpp>

namespace ak::geom {

class NurbsCurve {
public:
    int p = 0;                ///< degree
    std::vector<double> U;    ///< clamped knot vector, size cp.size() + p + 1
    std::vector<Vec4d> cp;    ///< homogeneous (w*x, w*y, w*z, w)

    NurbsCurve() = default;
    NurbsCurve(int p_, std::vector<double> U_, std::vector<Vec4d> cp_)
        : p(p_), U(std::move(U_)), cp(std::move(cp_)) {}

    int n() const { return static_cast<int>(cp.size()) - 1; }
    double t0() const { return U[static_cast<std::size_t>(p)]; }
    double t1() const { return U[static_cast<std::size_t>(n() + 1)]; }

    Vec4d homogeneous(double t) const {
        const int s = find_span(n(), p, t, U);
        const auto N = ders_basis_funs(s, t, p, 0, U);
        Vec4d h{};
        for (int j = 0; j <= p; ++j) h = h + cp[static_cast<std::size_t>(s - p + j)] * N[0][j];
        return h;
    }
    Vec3d point(double t) const {
        const Vec4d h = homogeneous(t);
        return {h.x / h.w, h.y / h.w, h.z / h.w};
    }
    /// Smallest and largest control weight (all must be positive for the
    /// convex-hull bounds used by the kernel's certified routines).
    std::array<double, 2> weight_range() const {
        double lo = cp[0].w, hi = cp[0].w;
        for (const auto& c : cp) { lo = std::min(lo, c.w); hi = std::max(hi, c.w); }
        return {lo, hi};
    }
};

/// Polynomial Bezier segment through homogeneous control points, as a curve on [0,1].
inline NurbsCurve bezier_curve(std::vector<Vec4d> cp) {
    const int p = static_cast<int>(cp.size()) - 1;
    std::vector<double> U(static_cast<std::size_t>(2 * p + 2), 0.0);
    for (int i = p + 1; i < 2 * p + 2; ++i) U[static_cast<std::size_t>(i)] = 1.0;
    return NurbsCurve(p, std::move(U), std::move(cp));
}

/// Polyline-free convenience: a non-rational curve from 3D points.
inline NurbsCurve polynomial_curve(int p, std::vector<double> U, const std::vector<Vec3d>& pts) {
    std::vector<Vec4d> cp;
    cp.reserve(pts.size());
    for (const auto& x : pts) cp.push_back({x.x, x.y, x.z, 1.0});
    return NurbsCurve(p, std::move(U), std::move(cp));
}

/// The exact rational circle of radius r about the z axis in the plane z = h,
/// from angle 0 through `quarters` quarter turns (1..4).
///
/// Each quarter is the rational quadratic (1 - t^2, 2t, 1 + t^2) / (1 + t^2)
/// on t in [0,1], whose control data are exact in doubles when r and h are.
/// Quarter k is scaled by 2^k so that consecutive quarters share their joint
/// homogeneous point bitwise, which a single NURBS needs; the scaling is
/// exact. Weights run 1 .. 16 for a full circle; every point lies on the
/// circle exactly in real arithmetic.
inline NurbsCurve exact_circle(double r, double h, int quarters) {
    quarters = std::clamp(quarters, 1, 4);
    // Unit-radius homogeneous controls (cx, cy, w) of the four scaled quarters.
    static constexpr double arc[9][3] = {{1, 0, 1},   {1, 1, 1},   {0, 2, 2},
                                         {-2, 2, 2},  {-4, 0, 4},  {-4, -4, 4},
                                         {0, -8, 8},  {8, -8, 8},  {16, 0, 16}};
    std::vector<Vec4d> cp;
    for (int i = 0; i <= 2 * quarters; ++i)
        cp.push_back({arc[i][0] * r, arc[i][1] * r, arc[i][2] * h, arc[i][2]});
    std::vector<double> U = {0, 0, 0};
    for (int k = 1; k < quarters; ++k) {
        const double t = static_cast<double>(k) / quarters;
        U.push_back(t);
        U.push_back(t);
    }
    U.insert(U.end(), {1, 1, 1});
    return NurbsCurve(2, std::move(U), std::move(cp));
}

namespace curve_detail {

inline void insert_knot(NurbsCurve& c, double t) { detail::insert_knot(c.U, c.cp, c.p, t); }

/// Number of times t occurs in U.
inline int multiplicity(const std::vector<double>& U, double t) {
    return static_cast<int>(std::count(U.begin(), U.end(), t));
}

/// Raise a Bezier segment (homogeneous controls) by one degree.
inline std::vector<Vec4d> elevate_bezier(const std::vector<Vec4d>& P) {
    const int n = static_cast<int>(P.size()) - 1;
    std::vector<Vec4d> Q(static_cast<std::size_t>(n + 2));
    Q[0] = P[0];
    Q[static_cast<std::size_t>(n + 1)] = P[static_cast<std::size_t>(n)];
    for (int i = 1; i <= n; ++i) {
        const double a = static_cast<double>(i) / (n + 1);
        Q[static_cast<std::size_t>(i)] =
            P[static_cast<std::size_t>(i - 1)] * a + P[static_cast<std::size_t>(i)] * (1.0 - a);
    }
    return Q;
}

} // namespace curve_detail

/// Affine change of parameter so that the curve runs over [0, 1].
inline NurbsCurve normalized_domain(NurbsCurve c) {
    const double a = c.t0(), b = c.t1();
    if (a == 0.0 && b == 1.0) return c;
    for (auto& u : c.U) u = (u - a) / (b - a);
    c.U.front() = 0.0;
    for (int i = 0; i <= c.p; ++i) c.U[static_cast<std::size_t>(i)] = 0.0;
    for (int i = 0; i <= c.p; ++i) c.U[c.U.size() - 1 - static_cast<std::size_t>(i)] = 1.0;
    return c;
}

/// Bezier segments (homogeneous controls) of a curve, one per knot span.
inline std::vector<std::vector<Vec4d>> bezier_segments(const NurbsCurve& c) {
    NurbsCurve r = c;
    for (double t : detail::bezier_refinement(c.U, c.p)) curve_detail::insert_knot(r, t);
    const auto spans = detail::distinct_knots(c.U);
    std::vector<std::vector<Vec4d>> segs;
    for (std::size_t s = 0; s + 1 < spans.size(); ++s) {
        std::vector<Vec4d> seg(r.cp.begin() + static_cast<std::ptrdiff_t>(s * r.p),
                               r.cp.begin() + static_cast<std::ptrdiff_t>(s * r.p + r.p + 1));
        segs.push_back(std::move(seg));
    }
    return segs;
}

/// The same curve at degree p + k. The curve is split into Bezier segments,
/// each is elevated, and they are reassembled with interior knots of full
/// multiplicity (C0 in the knot vector, the same point set and parameterisation).
inline NurbsCurve elevate_degree(const NurbsCurve& c, int k) {
    if (k <= 0) return c;
    auto segs = bezier_segments(c);
    for (auto& s : segs)
        for (int i = 0; i < k; ++i) s = curve_detail::elevate_bezier(s);
    const int q = c.p + k;
    const auto spans = detail::distinct_knots(c.U);
    std::vector<double> U(static_cast<std::size_t>(q + 1), spans.front());
    std::vector<Vec4d> cp;
    for (std::size_t s = 0; s < segs.size(); ++s) {
        for (std::size_t i = (s == 0 ? 0 : 1); i < segs[s].size(); ++i) cp.push_back(segs[s][i]);
        if (s + 1 < segs.size())
            for (int m = 0; m < q; ++m) U.push_back(spans[s + 1]);
    }
    for (int m = 0; m <= q; ++m) U.push_back(spans.back());
    return NurbsCurve(q, std::move(U), std::move(cp));
}

/// Bring curves to one degree and one knot vector on [0, 1] (Piegl & Tiller
/// §10.3): normalise the domains, elevate to the highest degree, then insert
/// into each curve every knot it lacks, at the highest multiplicity found.
/// Shapes and parameterisations are unchanged up to rounding.
inline std::vector<NurbsCurve> make_compatible(std::vector<NurbsCurve> cs) {
    int q = 0;
    for (auto& c : cs) {
        c = normalized_domain(std::move(c));
        q = std::max(q, c.p);
    }
    for (auto& c : cs) c = elevate_degree(c, q - c.p);
    std::vector<double> all;
    for (const auto& c : cs)
        for (double t : detail::distinct_knots(c.U))
            if (t > 0.0 && t < 1.0) all.push_back(t);
    std::sort(all.begin(), all.end());
    all.erase(std::unique(all.begin(), all.end()), all.end());
    for (double t : all) {
        int m = 0;
        for (const auto& c : cs) m = std::max(m, curve_detail::multiplicity(c.U, t));
        for (auto& c : cs)
            for (int have = curve_detail::multiplicity(c.U, t); have < m; ++have)
                curve_detail::insert_knot(c, t);
    }
    return cs;
}

} // namespace ak::geom
