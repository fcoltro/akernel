// akernel — certified curve/surface intersection.
//
// A rational Bezier curve C(t) and a rational Bezier patch S(u,v) meet where
//     G(t,u,v) = C(t) w_S(u,v) - S(u,v) w_C(t) = 0     (homogeneous form),
// three polynomial equations in three unknowns.
//
// Every returned root comes with
//   - a tight enclosure of the exact root, and
//   - a box in which it is proven to be the only root (interval Krawczyk test).
// A regular root means the curve crosses the surface transversally there.
//
// Completeness: the solver subdivides the whole domain and shows that every
// sub-box is either free of roots or lies inside a root's uniqueness box.
// What it cannot show (a tangential contact, for instance) is reported.
//
// In a Boolean operation these roots are the vertices where an edge of one
// solid pierces a face of the other; a regular root also proves that the
// intersection curve ending there crosses the face boundary.
#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <vector>

#include <akernel/solve/ssi.hpp>

namespace ak::solve {

using Pt3 = std::array<double, 3>;
using Box3 = std::array<Interval, 3>;

/// A rational Bezier curve, stored as a patch of degree (n, 0) so that the
/// patch machinery (subdivision, restriction, hulls) applies unchanged. The
/// second parameter is ignored.
inline PolyPatchD make_curve(const std::vector<Vec4d>& cp) {
    PolyPatchD c(static_cast<int>(cp.size()) - 1, 0);
    for (std::size_t i = 0; i < cp.size(); ++i) c.c[i] = cp[i];
    return c;
}
inline PolyPatchD make_line(const Vec3d& a, const Vec3d& b) {
    return make_curve({{a.x, a.y, a.z, 1.0}, {b.x, b.y, b.z, 1.0}});
}

inline bool invert3(const double A[3][3], double inv[3][3]) {
    const double d = det3(A);
    if (!(std::fabs(d) > 1e-300) || !std::isfinite(d)) return false;
    inv[0][0] = (A[1][1] * A[2][2] - A[1][2] * A[2][1]) / d;
    inv[0][1] = (A[0][2] * A[2][1] - A[0][1] * A[2][2]) / d;
    inv[0][2] = (A[0][1] * A[1][2] - A[0][2] * A[1][1]) / d;
    inv[1][0] = (A[1][2] * A[2][0] - A[1][0] * A[2][2]) / d;
    inv[1][1] = (A[0][0] * A[2][2] - A[0][2] * A[2][0]) / d;
    inv[1][2] = (A[0][2] * A[1][0] - A[0][0] * A[1][2]) / d;
    inv[2][0] = (A[1][0] * A[2][1] - A[1][1] * A[2][0]) / d;
    inv[2][1] = (A[0][1] * A[2][0] - A[0][0] * A[2][1]) / d;
    inv[2][2] = (A[0][0] * A[1][1] - A[0][1] * A[1][0]) / d;
    return true;
}

class CsiSystem {
    PolyPatchD c_, ct_, s_, su_, sv_;
    PolyPatchI ci_, cti_, si_, sui_, svi_;

public:
    CsiSystem(const PolyPatchD& curve, const PolyPatchD& surface)
        : c_(curve), ct_(curve.deriv_u()), s_(surface), su_(surface.deriv_u()), sv_(surface.deriv_v()),
          ci_(geom::to_interval(c_)), cti_(geom::to_interval(ct_)), si_(geom::to_interval(s_)),
          sui_(geom::to_interval(su_)), svi_(geom::to_interval(sv_)) {}

    const PolyPatchI& curve_i() const { return ci_; }
    const PolyPatchI& surface_i() const { return si_; }

    /// G and its Jacobian at x = (t, u, v), in doubles.
    void eval(const Pt3& x, double G[3], double J[3][3]) const {
        const Vec4d C = c_.eval(x[0], 0.0), Ct = ct_.eval(x[0], 0.0);
        const Vec4d S = s_.eval(x[1], x[2]), Su = su_.eval(x[1], x[2]), Sv = sv_.eval(x[1], x[2]);
        const double c[3] = {C.x, C.y, C.z}, ct[3] = {Ct.x, Ct.y, Ct.z};
        const double s[3] = {S.x, S.y, S.z}, su[3] = {Su.x, Su.y, Su.z}, sv[3] = {Sv.x, Sv.y, Sv.z};
        for (int k = 0; k < 3; ++k) {
            G[k] = c[k] * S.w - s[k] * C.w;
            J[k][0] = ct[k] * S.w - s[k] * Ct.w;
            J[k][1] = c[k] * Su.w - su[k] * C.w;
            J[k][2] = c[k] * Sv.w - sv[k] * C.w;
        }
    }

    /// Rigorous enclosure of G at a point.
    std::array<Interval, 3> G_at(const Pt3& x) const {
        const auto C = ci_.eval(Interval(x[0]), Interval(0.0));
        const auto S = si_.eval(Interval(x[1]), Interval(x[2]));
        return {C.x * S.w - S.x * C.w, C.y * S.w - S.y * C.w, C.z * S.w - S.z * C.w};
    }

    /// Rigorous enclosure of the Jacobian over a box.
    void J_over(const Box3& X, Interval J[3][3]) const {
        const double t0 = X[0].lo, t1 = X[0].hi, u0 = X[1].lo, u1 = X[1].hi, v0 = X[2].lo, v1 = X[2].hi;
        const auto C = geom::hull(ci_.restricted(t0, t1, 0.0, 1.0));
        const auto Ct = geom::hull(cti_.restricted(t0, t1, 0.0, 1.0));
        const auto S = geom::hull(si_.restricted(u0, u1, v0, v1));
        const auto Su = geom::hull(sui_.restricted(u0, u1, v0, v1));
        const auto Sv = geom::hull(svi_.restricted(u0, u1, v0, v1));
        const Interval c[3] = {C.x, C.y, C.z}, ct[3] = {Ct.x, Ct.y, Ct.z};
        const Interval s[3] = {S.x, S.y, S.z}, su[3] = {Su.x, Su.y, Su.z}, sv[3] = {Sv.x, Sv.y, Sv.z};
        for (int k = 0; k < 3; ++k) {
            J[k][0] = ct[k] * S.w - s[k] * Ct.w;
            J[k][1] = c[k] * Su.w - su[k] * C.w;
            J[k][2] = c[k] * Sv.w - sv[k] * C.w;
        }
    }

    /// Rigorous enclosure of the curve point for a parameter interval.
    Vec3<Interval> curve_point(const Interval& t) const {
        const auto h = ci_.eval(t, Interval(0.0));
        return {h.x / h.w, h.y / h.w, h.z / h.w};
    }
};

/// Newton iteration in doubles. Returns false if it does not settle.
inline bool csi_newton(const CsiSystem& sys, Pt3& x) {
    for (int it = 0; it < 30; ++it) {
        double G[3], J[3][3], inv[3][3];
        sys.eval(x, G, J);
        if (!invert3(J, inv)) return false;
        double step = 0.0;
        for (int i = 0; i < 3; ++i) {
            double dx = 0.0;
            for (int k = 0; k < 3; ++k) dx -= inv[i][k] * G[k];
            x[i] += dx;
            step = std::max(step, std::fabs(dx));
        }
        if (!std::isfinite(step) || step > 4.0) return false;
        if (step < 1e-13) return true;
    }
    return false;
}

/// Krawczyk test. On success K is strictly inside X and encloses the unique
/// root in X.
inline bool csi_krawczyk(const CsiSystem& sys, const Pt3& xhat, const Box3& X, Box3& K) {
    double G[3], J[3][3], C[3][3];
    sys.eval(xhat, G, J);
    if (!invert3(J, C)) return false;
    const auto g = sys.G_at(xhat);
    Interval JX[3][3];
    sys.J_over(X, JX);
    for (int i = 0; i < 3; ++i) {
        Interval ki(xhat[i]);
        for (int k = 0; k < 3; ++k) ki -= Interval(C[i][k]) * g[k];
        for (int j = 0; j < 3; ++j) {
            Interval mij(i == j ? 1.0 : 0.0);
            for (int k = 0; k < 3; ++k) mij -= Interval(C[i][k]) * JX[k][j];
            ki += mij * (X[j] - Interval(xhat[j]));
        }
        K[i] = ki;
        if (!num::strictly_inside(ki, X[i])) return false;
    }
    return true;
}

struct CsiRoot {
    Pt3 x{};        ///< approximate root (t, u, v)
    Box3 enc{};     ///< certified enclosure of the exact root
    Box3 unique{};  ///< the root is the only one in this box
    bool interior = false;    ///< the enclosure lies strictly inside the searched domain
    bool uv_interior = false; ///< its surface parameters lie strictly inside (0,1)^2
};

struct CsiOptions {
    int max_depth = 48;
    std::size_t node_limit = 200000;
};

struct CsiResult {
    std::vector<CsiRoot> roots;       ///< roots not proven to lie outside the domain
    bool complete = false;            ///< proven: no other root in [0,1]^3
    std::size_t unverified_boxes = 0;
    Box3 unverified_hull{};
    std::size_t nodes = 0;
};

namespace detail {

inline bool subset3(const Box3& a, const Box3& b) {
    return num::subset(a[0], b[0]) && num::subset(a[1], b[1]) && num::subset(a[2], b[2]);
}

struct CsiSolver {
    const CsiSystem& sys;
    const CsiOptions& opt;
    std::vector<CsiRoot> found;
    std::vector<Box3> covers; ///< boxes proven to hold exactly one root, already in `found`
    std::vector<Box3> unverified;
    std::size_t nodes = 0;

    /// Enlarge the region around a proven root in which it is the only root
    /// (boxes of doubling size, each proven by the Krawczyk test), so that
    /// boxes beside a root where the curve leaves the surface at a small
    /// angle need not be refined to the depth limit. Sound: a cover holds no
    /// root other than the one already found.
    void grow_cover(const Pt3& x, const Box3& X) {
        double rad[3];
        for (int i = 0; i < 3; ++i) rad[i] = std::max(0.5 * (X[i].hi - X[i].lo), 1e-12);
        for (int step = 0; step < 24; ++step) {
            Box3 G, K;
            bool small = true;
            for (int i = 0; i < 3; ++i) {
                rad[i] *= 2.0;
                small = small && rad[i] <= 0.25;
                G[i] = Interval(x[i] - rad[i], x[i] + rad[i]);
            }
            if (!small || !csi_krawczyk(sys, x, G, K)) return;
            covers.push_back(G);
        }
    }

    void run(const PolyPatchI& C, const PolyPatchI& S, const Box3& box, int depth) {
        ++nodes;
        if (boxes_disjoint(geom::projected_bounds(C), geom::projected_bounds(S))) return;
        for (const Box3& cover : covers)
            if (subset3(box, cover)) return;

        // The curve piece may run close to the surface piece without meeting
        // it; along the surface normal the surface piece is thin.
        {
            const Vec3d n = corner_normal(S);
            const double len = geom::norm(n);
            if (len > 0.0 && std::isfinite(len) &&
                num::disjoint(extent_along(C, n / len), extent_along(S, n / len)))
                return;
        }

        // Look for a root from the box centre and try to prove it unique in
        // a box that covers this one.
        Pt3 x = {num::mid(box[0]), num::mid(box[1]), num::mid(box[2])};
        if (csi_newton(sys, x)) {
            bool near = true, known = false;
            for (int i = 0; i < 3; ++i) {
                const double w = box[i].hi - box[i].lo;
                if (x[i] < box[i].lo - w || x[i] > box[i].hi + w) near = false;
            }
            for (const CsiRoot& r : found)
                if (num::contains(r.unique[0], x[0]) && num::contains(r.unique[1], x[1]) &&
                    num::contains(r.unique[2], x[2]))
                    known = true;
            if (near) {
                Box3 X, K;
                for (int i = 0; i < 3; ++i) {
                    const double pad = 0.05 * (box[i].hi - box[i].lo) + 1e-9;
                    X[i] = Interval(std::min(box[i].lo, x[i]) - pad, std::max(box[i].hi, x[i]) + pad);
                }
                if (csi_krawczyk(sys, x, X, K)) {
                    // This box is settled either way: X holds exactly one
                    // root, which is new unless it lies where an earlier one
                    // was proven unique.
                    covers.push_back(X);
                    if (!known) {
                        CsiRoot r;
                        r.x = x;
                        r.unique = X;
                        r.enc = K;
                        for (double rad = 1e-14; rad <= 1e-7; rad *= 10.0) {
                            Box3 Xs, Ks;
                            for (int i = 0; i < 3; ++i) Xs[i] = Interval(x[i] - rad, x[i] + rad);
                            if (subset3(Xs, X) && csi_krawczyk(sys, x, Xs, Ks)) { r.enc = Ks; break; }
                        }
                        found.push_back(r);
                        grow_cover(x, X);
                    }
                    return;
                }
            }
        }

        if (depth >= opt.max_depth || nodes >= opt.node_limit) { unverified.push_back(box); return; }

        int dim = 0;
        double wmax = 0.0;
        for (int i = 0; i < 3; ++i) {
            const double w = box[i].hi - box[i].lo;
            if (w > wmax) { wmax = w; dim = i; }
        }
        const double midv = 0.5 * (box[dim].lo + box[dim].hi);
        Box3 lo = box, hi = box;
        lo[dim].hi = midv;
        hi[dim].lo = midv;
        PolyPatchI Pl, Ph;
        if (dim == 0) {
            C.split_u(Interval(0.5), Pl, Ph);
            run(Pl, S, lo, depth + 1);
            run(Ph, S, hi, depth + 1);
        } else {
            if (dim == 1) S.split_u(Interval(0.5), Pl, Ph); else S.split_v(Interval(0.5), Pl, Ph);
            run(C, Pl, lo, depth + 1);
            run(C, Ph, hi, depth + 1);
        }
    }
};

} // namespace detail

/// All intersections of the curve (t in [t_lo, t_hi], by default the whole
/// curve) with the patch ((u,v) in [0,1]^2).
inline CsiResult intersect_curve_surface(const CsiSystem& sys, const CsiOptions& opt = {}, double t_lo = 0.0,
                                         double t_hi = 1.0) {
    detail::CsiSolver solver{sys, opt, {}, {}, {}, 0};
    const Box3 domain = {Interval(t_lo, t_hi), Interval(0.0, 1.0), Interval(0.0, 1.0)};
    const bool whole = t_lo == 0.0 && t_hi == 1.0;
    solver.run(whole ? sys.curve_i() : sys.curve_i().restricted(t_lo, t_hi, 0.0, 1.0), sys.surface_i(), domain, 0);

    CsiResult res;
    res.nodes = solver.nodes;
    for (const CsiRoot& f : solver.found) {
        bool outside = f.enc[0].hi < t_lo || f.enc[0].lo > t_hi;
        bool interior = f.enc[0].lo > t_lo && f.enc[0].hi < t_hi;
        bool uv_interior = true;
        for (int i = 1; i < 3; ++i) {
            if (f.enc[i].hi < 0.0 || f.enc[i].lo > 1.0) outside = true;
            if (!(f.enc[i].lo > 0.0 && f.enc[i].hi < 1.0)) uv_interior = false;
        }
        if (outside) continue;
        CsiRoot r = f;
        r.uv_interior = uv_interior;
        r.interior = interior && uv_interior;
        res.roots.push_back(r);
    }
    res.unverified_boxes = solver.unverified.size();
    for (std::size_t i = 0; i < solver.unverified.size(); ++i) {
        if (i == 0) res.unverified_hull = solver.unverified[i];
        else
            for (int k = 0; k < 3; ++k)
                res.unverified_hull[k] = num::hull(res.unverified_hull[k], solver.unverified[i][k]);
    }
    res.complete = solver.unverified.empty();
    return res;
}

} // namespace ak::solve
