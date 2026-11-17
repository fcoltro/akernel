// akernel — certified three-surface intersection points.
//
// After a Boolean, some edges of a solid are intersection curves of two
// carriers S1 and S2. Where a third surface S3 cuts such an edge, the new
// vertex lies on all three surfaces:
//     F12(u1,v1,u2,v2) = 0      (S1 = S2, three equations)
//     F13(u1,v1,s,t)   = 0      (S1 = S3, three equations)
// six polynomial equations in the six unknowns y = (u1,v1,u2,v2,s,t).
//
// The search runs along the edge's certified tube: a node is a slice range
// of one tube segment together with a rectangle of S3. A node is discarded
// if the enclosure of that piece of the edge is separated from that piece
// of S3, or if it lies in a box where a root was proven unique (Krawczyk
// test in six dimensions). What cannot be decided is reported.
#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <map>
#include <tuple>
#include <vector>

#include <akernel/op/trim.hpp>
#include <akernel/solve/csi.hpp>

namespace ak::op {

using geom::PolyPatchD;
using geom::PolyPatchI;
using geom::Vec3;
using geom::Vec3d;

using Pt6 = std::array<double, 6>;
using Box6 = std::array<Interval, 6>;

inline bool invert6(const double A[6][6], double inv[6][6]) {
    double a[6][12];
    for (int i = 0; i < 6; ++i)
        for (int j = 0; j < 6; ++j) {
            a[i][j] = A[i][j];
            a[i][j + 6] = i == j ? 1.0 : 0.0;
        }
    for (int k = 0; k < 6; ++k) {
        int piv = k;
        for (int i = k + 1; i < 6; ++i)
            if (std::fabs(a[i][k]) > std::fabs(a[piv][k])) piv = i;
        if (!(std::fabs(a[piv][k]) > 1e-300)) return false;
        if (piv != k)
            for (int j = 0; j < 12; ++j) std::swap(a[k][j], a[piv][j]);
        const double d = a[k][k];
        for (int j = 0; j < 12; ++j) a[k][j] /= d;
        for (int i = 0; i < 6; ++i) {
            if (i == k) continue;
            const double f = a[i][k];
            if (f == 0.0) continue;
            for (int j = 0; j < 12; ++j) a[i][j] -= f * a[k][j];
        }
    }
    for (int i = 0; i < 6; ++i)
        for (int j = 0; j < 6; ++j) {
            inv[i][j] = a[i][j + 6];
            if (!std::isfinite(inv[i][j])) return false;
        }
    return true;
}

class TsiSystem {
    const SsiSystem& edge_; ///< S1 against S2: unknowns y0..y3
    SsiSystem third_;       ///< S1 against S3: unknowns (y0, y1, y4, y5)

public:
    TsiSystem(const SsiSystem& edge, const PolyPatchD& third) : edge_(edge), third_(edge.patch_a(), third) {}
    const SsiSystem& edge() const { return edge_; }

    void eval(const Pt6& y, double G[6], double J[6][6]) const {
        double F[3], J4[3][4];
        edge_.eval({y[0], y[1], y[2], y[3]}, F, J4);
        for (int r = 0; r < 3; ++r) {
            G[r] = F[r];
            for (int c = 0; c < 4; ++c) J[r][c] = J4[r][c];
            J[r][4] = J[r][5] = 0.0;
        }
        third_.eval({y[0], y[1], y[4], y[5]}, F, J4);
        for (int r = 0; r < 3; ++r) {
            G[3 + r] = F[r];
            J[3 + r][0] = J4[r][0];
            J[3 + r][1] = J4[r][1];
            J[3 + r][2] = J[3 + r][3] = 0.0;
            J[3 + r][4] = J4[r][2];
            J[3 + r][5] = J4[r][3];
        }
    }

    std::array<Interval, 6> G_at(const Pt6& y) const {
        const auto a = edge_.F_at({y[0], y[1], y[2], y[3]});
        const auto b = third_.F_at({y[0], y[1], y[4], y[5]});
        return {a[0], a[1], a[2], b[0], b[1], b[2]};
    }

    void J_over(const Box6& Y, Interval J[6][6]) const {
        Interval J4[3][4];
        edge_.J_over({Y[0], Y[1], Y[2], Y[3]}, J4);
        for (int r = 0; r < 3; ++r) {
            for (int c = 0; c < 4; ++c) J[r][c] = J4[r][c];
            J[r][4] = J[r][5] = Interval(0.0);
        }
        third_.J_over({Y[0], Y[1], Y[4], Y[5]}, J4);
        for (int r = 0; r < 3; ++r) {
            J[3 + r][0] = J4[r][0];
            J[3 + r][1] = J4[r][1];
            J[3 + r][2] = J[3 + r][3] = Interval(0.0);
            J[3 + r][4] = J4[r][2];
            J[3 + r][5] = J4[r][3];
        }
    }
};

inline bool tsi_newton(const TsiSystem& sys, Pt6& y) {
    for (int it = 0; it < 30; ++it) {
        double G[6], J[6][6], inv[6][6];
        sys.eval(y, G, J);
        if (!invert6(J, inv)) return false;
        double step = 0.0;
        for (int i = 0; i < 6; ++i) {
            double dy = 0.0;
            for (int k = 0; k < 6; ++k) dy -= inv[i][k] * G[k];
            y[i] += dy;
            step = std::max(step, std::fabs(dy));
        }
        if (!std::isfinite(step) || step > 4.0) return false;
        if (step < 1e-13) return true;
    }
    return false;
}

/// Krawczyk test. On success K is strictly inside Y and encloses the unique
/// root in Y.
inline bool tsi_krawczyk(const TsiSystem& sys, const Pt6& yhat, const Box6& Y, Box6& K) {
    double G[6], J[6][6], C[6][6];
    sys.eval(yhat, G, J);
    if (!invert6(J, C)) return false;
    const auto g = sys.G_at(yhat);
    Interval JY[6][6];
    sys.J_over(Y, JY);
    for (int i = 0; i < 6; ++i) {
        Interval ki(yhat[i]);
        for (int k = 0; k < 6; ++k) ki -= Interval(C[i][k]) * g[static_cast<std::size_t>(k)];
        for (int j = 0; j < 6; ++j) {
            Interval mij(i == j ? 1.0 : 0.0);
            for (int k = 0; k < 6; ++k) mij -= Interval(C[i][k]) * JY[k][j];
            ki += mij * (Y[j] - Interval(yhat[j]));
        }
        K[i] = ki;
        if (!num::strictly_inside(ki, Y[i])) return false;
    }
    return true;
}

struct TsiRoot {
    Pt6 y{};       ///< approximate root (u1, v1, u2, v2, s, t)
    Box6 enc{};    ///< certified enclosure of the exact root
    Box6 unique{}; ///< the root is the only one in this box
};

struct TsiResult {
    std::vector<TsiRoot> roots; ///< roots whose (s, t) is not proven outside [0,1]^2
    bool complete = false;      ///< proven: no other root on the searched part of the edge
    std::size_t unverified_boxes = 0;
    Box6 unverified_hull{};     ///< hull of the boxes left (curve 4, then surface 2), when there are any
    std::size_t nodes = 0;
};

/// Enclosures of pieces of an arc's tube, by segment and slice range. They
/// depend only on the arc, not on the surface it is searched against, and
/// the search bisects at deterministic midpoints, so one edge searched
/// against several surfaces asks for the same enclosures again. Kept per
/// edge; the values are exactly those computed without it.
struct ArcEnclosureCache {
    std::map<std::tuple<const TubeSegment*, double, double>, Box4> boxes;
    std::size_t hits = 0, misses = 0;
    Box4 get(const SsiSystem& sys, const TubeSegment& seg, double c0, double c1) {
        const auto key = std::make_tuple(&seg, c0, c1);
        const auto it = boxes.find(key);
        if (it != boxes.end()) { ++hits; return it->second; }
        ++misses;
        return boxes[key] = arc_enclosure(sys, seg, c0, c1);
    }
};

namespace detail {

/// Is a box of space separated from patch S along the direction of S's
/// corner normal? (The extents of both along it are disjoint.)
inline bool separated_along_normal(const Vec3<Interval>& pa, const PolyPatchI& S) {
    const Vec3d n = solve::detail::corner_normal(S);
    const double len = geom::norm(n);
    if (!(len > 0.0) || !std::isfinite(len)) return false;
    const Vec3d d = n / len;
    const Interval e = Interval(d.x) * pa.x + Interval(d.y) * pa.y + Interval(d.z) * pa.z;
    return num::disjoint(e, solve::detail::extent_along(S, d));
}

/// A box around the arc of `seg` over the slice range [c0, c1], given a box X
/// that contains it: X cut down by the slab c0 <= m . (x - a) <= c1, one
/// coordinate at a time. Cheap and valid, but looser than arc_enclosure.
inline Box4 slab_box(const TubeSegment& seg, const Box4& X, double c0, double c1) {
    Box4 B = X;
    const Interval c(c0, c1);
    for (int i = 0; i < 4; ++i) {
        if (seg.m[i] == 0.0) continue;
        Interval rest(0.0);
        for (int j = 0; j < 4; ++j)
            if (j != i) rest += Interval(seg.m[j]) * (B[j] - Interval(seg.a[j]));
        const Interval xi = (c - rest) / Interval(seg.m[i]) + Interval(seg.a[i]);
        const double lo = std::max(B[i].lo, xi.lo), hi = std::min(B[i].hi, xi.hi);
        if (lo <= hi) B[i] = Interval(lo, hi);
    }
    return B;
}

/// Can a piece of the edge inside box X not meet patch S? The two cheap
/// tests of the search: bounding boxes, and extents along S's normal.
inline bool cheap_exclusion(const SsiSystem& edge, const Box4& X, const PolyPatchI& S) {
    const Vec3<Interval> pa = edge.point_a(X);
    return solve::detail::boxes_disjoint(pa, geom::projected_bounds(S)) || separated_along_normal(pa, S);
}

struct TsiSolver {
    const TsiSystem& sys;
    const solve::CsiOptions& opt;
    ArcEnclosureCache& cache;
    std::vector<TsiRoot> found;
    std::vector<Box6> covers; ///< boxes proven to hold exactly one root, already in `found`
    std::size_t unverified = 0, nodes = 0;
    Box6 hull{};

    /// Enlarge the region around a proven root in which it is the only root:
    /// boxes of doubling size around y, each proven by the Krawczyk test, are
    /// added as covers. Near a root where the edge leaves the surface at a
    /// small angle, the boxes next to the root's own box can be neither
    /// excluded nor certified until they are very small; a larger proven
    /// region settles them at once. Search nodes inside a cover hold no root
    /// other than one already found, so skipping them loses nothing.
    void grow_cover(const Pt6& y, const Box6& Y) {
        double rad[6];
        for (int i = 0; i < 6; ++i) rad[i] = std::max(0.5 * (Y[i].hi - Y[i].lo), 1e-12);
        for (int step = 0; step < 24; ++step) {
            Box6 G, K;
            bool small = true;
            for (int i = 0; i < 6; ++i) {
                rad[i] *= 2.0;
                small = small && rad[i] <= 0.25;
                G[i] = Interval(y[i] - rad[i], y[i] + rad[i]);
            }
            if (!small || !tsi_krawczyk(sys, y, G, K)) return;
            covers.push_back(G);
        }
    }

    /// X must enclose the arc of `seg` over the slice range [c0, c1].
    void run(const TubeSegment& seg, double c0, double c1, const Box4& X, const PolyPatchI& S, const Interval& s,
             const Interval& t, int depth) {
        ++nodes;
        const Vec3<Interval> pa = sys.edge().point_a(X);
        const Vec3<Interval> sb = geom::projected_bounds(S);
        if (solve::detail::boxes_disjoint(pa, sb)) return;
        const Box6 cur = {X[0], X[1], X[2], X[3], s, t};
        for (const Box6& cover : covers) {
            bool inside = true;
            for (int i = 0; i < 6 && inside; ++i) inside = num::subset(cur[i], cover[i]);
            if (inside) return;
        }
        if (separated_along_normal(pa, S)) return;

        const double cm = 0.5 * (c0 + c1);
        Pt4 x;
        for (int i = 0; i < 4; ++i) x[i] = seg.a[i] + (cm / seg.len) * (seg.b[i] - seg.a[i]);
        {
            Pt4 z = x;
            bool ok = solve::newton_slice(sys.edge(), seg.m, seg.a, cm, z);
            for (int i = 0; i < 4 && ok; ++i) ok = std::isfinite(z[i]);
            if (ok) x = z;
        }
        Pt6 y = {x[0], x[1], x[2], x[3], num::mid(s), num::mid(t)};
        if (tsi_newton(sys, y)) {
            bool near = true, known = false;
            for (int i = 0; i < 6; ++i) {
                const double w = cur[i].hi - cur[i].lo;
                if (y[i] < cur[i].lo - w || y[i] > cur[i].hi + w) near = false;
            }
            // A root is the same as an earlier one if it lies in a box where
            // that one was proven unique.
            for (const TsiRoot& r : found) {
                bool in = true;
                for (int i = 0; i < 6 && in; ++i) in = num::contains(r.unique[i], y[i]);
                if (in) known = true;
            }
            if (near) {
                Box6 Y, K;
                for (int i = 0; i < 6; ++i) {
                    const double pad = 0.05 * (cur[i].hi - cur[i].lo) + 1e-9;
                    Y[i] = Interval(std::min(cur[i].lo, y[i]) - pad, std::max(cur[i].hi, y[i]) + pad);
                }
                if (tsi_krawczyk(sys, y, Y, K)) {
                    // This node is settled either way: Y holds exactly one root.
                    covers.push_back(Y);
                    if (!known) {
                        TsiRoot r;
                        r.y = y;
                        r.unique = Y;
                        r.enc = K;
                        for (double rad = 1e-14; rad <= 1e-7; rad *= 10.0) {
                            Box6 Ys, Ks;
                            bool sub = true;
                            for (int i = 0; i < 6; ++i) {
                                Ys[i] = Interval(y[i] - rad, y[i] + rad);
                                sub = sub && num::subset(Ys[i], Y[i]);
                            }
                            if (sub && tsi_krawczyk(sys, y, Ys, Ks)) { r.enc = Ks; break; }
                        }
                        found.push_back(r);
                        grow_cover(y, Y);
                    }
                    return;
                }
            }
        }

        if (depth >= opt.max_depth || nodes >= opt.node_limit) {
            for (int i = 0; i < 6; ++i) hull[static_cast<std::size_t>(i)] = unverified == 0 ? cur[i] : num::hull(hull[static_cast<std::size_t>(i)], cur[i]);
            ++unverified;
            return;
        }

        // Split whichever is larger in space: the piece of the edge or the
        // piece of the surface. Splitting the surface cannot separate it
        // from a piece of edge that is bigger than it.
        const double wc = c1 - c0, ws = s.hi - s.lo, wt = t.hi - t.lo;
        const double size_edge = std::max(num::width(pa.x), std::max(num::width(pa.y), num::width(pa.z)));
        const double size_surf = std::max(num::width(sb.x), std::max(num::width(sb.y), num::width(sb.z)));
        if (size_edge >= 0.5 * size_surf && wc > 1e-13) {
            for (int half = 0; half < 2; ++half) {
                const double lo = half == 0 ? c0 : cm, hi = half == 0 ? cm : c1;
                // Excluded with a cheap box around the half: its tight
                // enclosure lies inside, so the half's search would stop at once.
                if (cheap_exclusion(sys.edge(), slab_box(seg, X, lo, hi), S)) continue;
                Box4 Xh = cache.get(sys.edge(), seg, lo, hi);
                for (int i = 0; i < 4; ++i) {
                    const double l = std::max(Xh[i].lo, X[i].lo), h = std::min(Xh[i].hi, X[i].hi);
                    if (l <= h) Xh[i] = Interval(l, h);
                }
                run(seg, lo, hi, Xh, S, s, t, depth + 1);
            }
        } else {
            PolyPatchI Pl, Ph;
            if (ws >= wt) {
                const double m = 0.5 * (s.lo + s.hi);
                S.split_u(Interval(0.5), Pl, Ph);
                run(seg, c0, c1, X, Pl, Interval(s.lo, m), t, depth + 1);
                run(seg, c0, c1, X, Ph, Interval(m, s.hi), t, depth + 1);
            } else {
                const double m = 0.5 * (t.lo + t.hi);
                S.split_v(Interval(0.5), Pl, Ph);
                run(seg, c0, c1, X, Pl, s, Interval(t.lo, m), depth + 1);
                run(seg, c0, c1, X, Ph, s, Interval(m, t.hi), depth + 1);
            }
        }
    }
};

} // namespace detail

/// All points where the arc (an edge lying on two carriers) meets a third
/// surface.
inline TsiResult intersect_arc_surface(const ArcGeom& arc, const PolyPatchD& third, const solve::CsiOptions& opt,
                                       ArcEnclosureCache& cache) {
    const SsiSystem& esys = *arc.base->sys;
    const TsiSystem sys(esys, third);
    detail::TsiSolver solver{sys, opt, cache, {}, {}, 0, 0};
    const PolyPatchI S = geom::to_interval(third);
    const Vec3<Interval> sb = geom::projected_bounds(S);
    for (const Span& sp : arc.spans) {
        const TubeSegment& seg = arc.base->curve->segs[static_cast<std::size_t>(sp.seg)];
        const double lo = std::max(std::min(sp.ca, sp.cb) - kSlicePad, seg.c_lo);
        const double hi = std::min(std::max(sp.ca, sp.cb) + kSlicePad, seg.c_hi);
        if (!(lo <= hi)) continue;
        // Cheap tests with the segment's own box first; tighten only if
        // needed. (Whatever they exclude, the search below would exclude at
        // its first node: the tight enclosure lies inside this box.)
        const Vec3<Interval> loose = esys.point_a(seg.X);
        if (solve::detail::boxes_disjoint(loose, sb)) continue;
        if (detail::separated_along_normal(loose, S)) continue;
        if (lo > seg.c_lo || hi < seg.c_hi)
            if (detail::cheap_exclusion(esys, detail::slab_box(seg, seg.X, lo, hi), S)) continue;
        solver.run(seg, lo, hi, cache.get(esys, seg, lo, hi), S, Interval(0.0, 1.0), Interval(0.0, 1.0), 0);
    }
    TsiResult res;
    res.nodes = solver.nodes;
    for (const TsiRoot& r : solver.found) {
        if (r.enc[4].hi < 0.0 || r.enc[4].lo > 1.0 || r.enc[5].hi < 0.0 || r.enc[5].lo > 1.0) continue;
        res.roots.push_back(r);
    }
    res.unverified_boxes = solver.unverified;
    res.unverified_hull = solver.hull;
    res.complete = solver.unverified == 0;
    return res;
}
inline TsiResult intersect_arc_surface(const ArcGeom& arc, const PolyPatchD& third,
                                       const solve::CsiOptions& opt = {}) {
    ArcEnclosureCache cache;
    return intersect_arc_surface(arc, third, opt, cache);
}

} // namespace ak::op
