// akernel — certified surface/surface intersection (docs/BLUEPRINT.md §1.4).
//
// Two rational Bezier patches S1(u,v) and S2(s,t) intersect in the variety
//     E = { x = (u,v,s,t) : F(x) = P1(u,v) w2(s,t) - P2(s,t) w1(u,v) = 0 },
// where (P, w) is the homogeneous form. F is polynomial, so rational patches
// need no special case.
//
// The result is a set of curves, each a chain of TubeSegments. A segment is
//   - a box X in R^4,
//   - a family of parallel hyperplanes  m.(x - a) = c,  c in [c_lo, c_hi],
// together with a proof (interval Krawczyk test) that every hyperplane of the
// family meets E in exactly one point inside X. So E restricted to the
// segment's region is exactly one regular arc. Consecutive segments overlap,
// and the overlap is checked with a certified point, so the arcs join into
// one curve.
//
// Completeness: a second pass subdivides the whole parameter domain and
// shows that every sub-box is either free of intersection (bounding boxes of
// the two sub-patches are disjoint) or lies inside a segment's region. What
// cannot be shown is reported; the kernel never claims completeness it did
// not prove.
#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <cstddef>
#include <vector>

#include <akernel/geom/bernstein.hpp>
#include <akernel/num/interval.hpp>

namespace ak::solve {

using geom::PolyPatchD;
using geom::PolyPatchI;
using geom::Vec3;
using geom::Vec3d;
using geom::Vec4;
using geom::Vec4d;
using num::Interval;

using Pt4 = std::array<double, 4>;
using Box4 = std::array<Interval, 4>;

// --- small fixed-size linear algebra ---------------------------------------

/// Inverse of a 4x4 matrix by Gauss-Jordan with partial pivoting.
inline bool invert4(const double A[4][4], double inv[4][4]) {
    double a[4][8];
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) {
            a[i][j] = A[i][j];
            a[i][j + 4] = i == j ? 1.0 : 0.0;
        }
    for (int k = 0; k < 4; ++k) {
        int piv = k;
        for (int i = k + 1; i < 4; ++i)
            if (std::fabs(a[i][k]) > std::fabs(a[piv][k])) piv = i;
        if (!(std::fabs(a[piv][k]) > 1e-300)) return false;
        if (piv != k)
            for (int j = 0; j < 8; ++j) std::swap(a[k][j], a[piv][j]);
        const double d = a[k][k];
        for (int j = 0; j < 8; ++j) a[k][j] /= d;
        for (int i = 0; i < 4; ++i) {
            if (i == k) continue;
            const double f = a[i][k];
            if (f == 0.0) continue;
            for (int j = 0; j < 8; ++j) a[i][j] -= f * a[k][j];
        }
    }
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) inv[i][j] = a[i][j + 4];
    return true;
}

inline double det3(const double m[3][3]) {
    return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) -
           m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
           m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
}

// --- the intersection system -------------------------------------------------

class SsiSystem {
    PolyPatchD a_, au_, av_, b_, bs_, bt_;
    PolyPatchI ai_, aui_, avi_, bi_, bsi_, bti_;

    static Vec3<Interval> xyzI(const Vec4<Interval>& h) { return {h.x, h.y, h.z}; }

public:
    SsiSystem(const PolyPatchD& A, const PolyPatchD& B)
        : a_(A), au_(A.deriv_u()), av_(A.deriv_v()), b_(B), bs_(B.deriv_u()), bt_(B.deriv_v()),
          ai_(geom::to_interval(a_)), aui_(geom::to_interval(au_)), avi_(geom::to_interval(av_)),
          bi_(geom::to_interval(b_)), bsi_(geom::to_interval(bs_)), bti_(geom::to_interval(bt_)) {}
    /// A system whose second surface is known only as an enclosure Bi (every
    /// control point an interval holding the exact one); Bmid, inside it, is
    /// used where doubles serve as guesses (Newton steps, seeds). Every
    /// certificate uses the interval patches, so it holds for the exact
    /// surface.
    SsiSystem(const PolyPatchD& A, const PolyPatchD& Bmid, const PolyPatchI& Bi)
        : a_(A), au_(A.deriv_u()), av_(A.deriv_v()), b_(Bmid), bs_(Bmid.deriv_u()), bt_(Bmid.deriv_v()),
          ai_(geom::to_interval(a_)), aui_(geom::to_interval(au_)), avi_(geom::to_interval(av_)), bi_(Bi), bsi_(Bi.deriv_u()),
          bti_(Bi.deriv_v()) {}

    const PolyPatchD& patch_a() const { return a_; }
    const PolyPatchD& patch_b() const { return b_; }
    const PolyPatchI& patch_ai() const { return ai_; }
    const PolyPatchI& patch_bi() const { return bi_; }

    /// F and its 3x4 Jacobian at a point, in doubles.
    void eval(const Pt4& x, double F[3], double J[3][4]) const {
        const Vec4d A = a_.eval(x[0], x[1]), Au = au_.eval(x[0], x[1]), Av = av_.eval(x[0], x[1]);
        const Vec4d B = b_.eval(x[2], x[3]), Bs = bs_.eval(x[2], x[3]), Bt = bt_.eval(x[2], x[3]);
        const double a[3] = {A.x, A.y, A.z}, au[3] = {Au.x, Au.y, Au.z}, av[3] = {Av.x, Av.y, Av.z};
        const double b[3] = {B.x, B.y, B.z}, bs[3] = {Bs.x, Bs.y, Bs.z}, bt[3] = {Bt.x, Bt.y, Bt.z};
        for (int k = 0; k < 3; ++k) {
            F[k] = a[k] * B.w - b[k] * A.w;
            J[k][0] = au[k] * B.w - b[k] * Au.w;
            J[k][1] = av[k] * B.w - b[k] * Av.w;
            J[k][2] = a[k] * Bs.w - bs[k] * A.w;
            J[k][3] = a[k] * Bt.w - bt[k] * A.w;
        }
    }

    /// Rigorous enclosure of F at a point.
    std::array<Interval, 3> F_at(const Pt4& x) const {
        const auto A = ai_.eval(Interval(x[0]), Interval(x[1]));
        const auto B = bi_.eval(Interval(x[2]), Interval(x[3]));
        return {A.x * B.w - B.x * A.w, A.y * B.w - B.y * A.w, A.z * B.w - B.z * A.w};
    }

    /// Rigorous enclosure of the Jacobian over a box.
    void J_over(const Box4& X, Interval J[3][4]) const {
        const double u0 = X[0].lo, u1 = X[0].hi, v0 = X[1].lo, v1 = X[1].hi;
        const double s0 = X[2].lo, s1 = X[2].hi, t0 = X[3].lo, t1 = X[3].hi;
        const auto A = geom::hull(ai_.restricted(u0, u1, v0, v1));
        const auto Au = geom::hull(aui_.restricted(u0, u1, v0, v1));
        const auto Av = geom::hull(avi_.restricted(u0, u1, v0, v1));
        const auto B = geom::hull(bi_.restricted(s0, s1, t0, t1));
        const auto Bs = geom::hull(bsi_.restricted(s0, s1, t0, t1));
        const auto Bt = geom::hull(bti_.restricted(s0, s1, t0, t1));
        const Interval a[3] = {A.x, A.y, A.z}, au[3] = {Au.x, Au.y, Au.z}, av[3] = {Av.x, Av.y, Av.z};
        const Interval b[3] = {B.x, B.y, B.z}, bs[3] = {Bs.x, Bs.y, Bs.z}, bt[3] = {Bt.x, Bt.y, Bt.z};
        for (int k = 0; k < 3; ++k) {
            J[k][0] = au[k] * B.w - b[k] * Au.w;
            J[k][1] = av[k] * B.w - b[k] * Av.w;
            J[k][2] = a[k] * Bs.w - bs[k] * A.w;
            J[k][3] = a[k] * Bt.w - bt[k] * A.w;
        }
    }

    /// Rigorous enclosures of the two surface points for a parameter box.
    Vec3<Interval> point_a(const Box4& X) const {
        const auto h = ai_.eval(X[0], X[1]);
        return {h.x / h.w, h.y / h.w, h.z / h.w};
    }
    Vec3<Interval> point_b(const Box4& X) const {
        const auto h = bi_.eval(X[2], X[3]);
        return {h.x / h.w, h.y / h.w, h.z / h.w};
    }

    /// Rigorous enclosures over X of the first derivatives of both projected
    /// surfaces. Unbounded where a weight enclosure contains zero.
    void frames(const Box4& X, Vec3<Interval>& Su, Vec3<Interval>& Sv, Vec3<Interval>& Bs,
                Vec3<Interval>& Bt) const {
        auto rat = [](const Vec4<Interval>& P, const Vec4<Interval>& D) {
            const Vec3<Interval> p{P.x / P.w, P.y / P.w, P.z / P.w};
            return Vec3<Interval>{(D.x - p.x * D.w) / P.w, (D.y - p.y * D.w) / P.w,
                                  (D.z - p.z * D.w) / P.w};
        };
        const auto A = ai_.eval(X[0], X[1]);
        const auto B = bi_.eval(X[2], X[3]);
        Su = rat(A, aui_.eval(X[0], X[1]));
        Sv = rat(A, avi_.eval(X[0], X[1]));
        Bs = rat(B, bsi_.eval(X[2], X[3]));
        Bt = rat(B, bti_.eval(X[2], X[3]));
    }

    /// Sine of the angle between the two surface normals: the conditioning
    /// indicator of blueprint §2.2. Zero means tangential intersection.
    double sin_angle(const Pt4& x) const {
        const Vec4d A = a_.eval(x[0], x[1]), Au = au_.eval(x[0], x[1]), Av = av_.eval(x[0], x[1]);
        const Vec4d B = b_.eval(x[2], x[3]), Bs = bs_.eval(x[2], x[3]), Bt = bt_.eval(x[2], x[3]);
        auto rat = [](const Vec4d& P, const Vec4d& D) {
            return Vec3d{(D.x - P.x / P.w * D.w) / P.w, (D.y - P.y / P.w * D.w) / P.w,
                         (D.z - P.z / P.w * D.w) / P.w};
        };
        const Vec3d n1 = cross(rat(A, Au), rat(A, Av)), n2 = cross(rat(B, Bs), rat(B, Bt));
        const double d = geom::norm(n1) * geom::norm(n2);
        return d > 0.0 ? geom::norm(cross(n1, n2)) / d : 0.0;
    }

    /// Unit tangent of E at x (null vector of J), with a continuous orientation.
    bool tangent(const Pt4& x, Pt4& t) const {
        double F[3], J[3][4];
        eval(x, F, J);
        double n2 = 0.0;
        for (int i = 0; i < 4; ++i) {
            double m[3][3];
            for (int r = 0; r < 3; ++r)
                for (int c = 0, cc = 0; c < 4; ++c)
                    if (c != i) m[r][cc++] = J[r][c];
            t[i] = (i % 2 == 0 ? 1.0 : -1.0) * det3(m);
            n2 += t[i] * t[i];
        }
        if (!(n2 > 0.0) || !std::isfinite(n2)) return false;
        const double inv = 1.0 / std::sqrt(n2);
        for (double& ti : t) ti *= inv;
        return true;
    }
};

// --- Newton and Krawczyk on the sliced system ---------------------------------
//
// G(x) = ( F(x), m.(x - a) - c ).  With c an interval the Krawczyk test
// proves existence and uniqueness for every slice value in c at once.

inline double dot4(const Pt4& a, const Pt4& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2] + a[3] * b[3]; }

/// Newton iteration in doubles for G(x) = 0. Returns false if it does not settle.
inline bool newton_slice(const SsiSystem& sys, const Pt4& m, const Pt4& a, double c, Pt4& x) {
    double last = std::numeric_limits<double>::infinity();
    for (int it = 0; it < 25; ++it) {
        double F[3], J[3][4], A[4][4], inv[4][4];
        sys.eval(x, F, J);
        for (int r = 0; r < 3; ++r)
            for (int k = 0; k < 4; ++k) A[r][k] = J[r][k];
        for (int k = 0; k < 4; ++k) A[3][k] = m[k];
        if (!invert4(A, inv)) return false;
        Pt4 d{};
        for (int k = 0; k < 4; ++k) d[k] = x[k] - a[k];
        const double g[4] = {F[0], F[1], F[2], dot4(m, d) - c};
        double step = 0.0;
        for (int i = 0; i < 4; ++i) {
            double dx = 0.0;
            for (int k = 0; k < 4; ++k) dx -= inv[i][k] * g[k];
            x[i] += dx;
            step = std::max(step, std::fabs(dx));
        }
        if (!std::isfinite(step) || step > 10.0) return false;
        if (step < 1e-14) return true;
        last = step;
    }
    // Converged to the level of rounding noise without a step below 1e-14
    // (small patches amplify the noise through the inverse Jacobian). The
    // point is only a guess for the Krawczyk test that follows, so accept it.
    return last < 1e-12;
}

/// Krawczyk operator of the sliced system over X, centred at xhat (which must
/// lie in X). K encloses every zero of G in X for every slice value in c,
/// whether or not K itself lies inside X, so it can be used as a contractor.
/// Returns false only if the preconditioner could not be formed.
inline bool krawczyk_operator(const SsiSystem& sys, const Pt4& m, const Pt4& a, const Interval& c,
                              const Pt4& xhat, const Box4& X, Box4& K) {
    double F[3], J[3][4], A[4][4], C[4][4];
    sys.eval(xhat, F, J);
    for (int r = 0; r < 3; ++r)
        for (int k = 0; k < 4; ++k) A[r][k] = J[r][k];
    for (int k = 0; k < 4; ++k) A[3][k] = m[k];
    if (!invert4(A, C)) return false;

    Interval g[4];
    const auto Fi = sys.F_at(xhat);
    g[0] = Fi[0]; g[1] = Fi[1]; g[2] = Fi[2];
    g[3] = -c;
    for (int k = 0; k < 4; ++k) g[3] += Interval(m[k]) * (Interval(xhat[k]) - Interval(a[k]));

    Interval JX[4][4];
    sys.J_over(X, JX);
    for (int k = 0; k < 4; ++k) JX[3][k] = Interval(m[k]);

    for (int i = 0; i < 4; ++i) {
        Interval ki(xhat[i]);
        for (int k = 0; k < 4; ++k) ki -= Interval(C[i][k]) * g[k];
        for (int j = 0; j < 4; ++j) {
            Interval mij(i == j ? 1.0 : 0.0);
            for (int k = 0; k < 4; ++k) mij -= Interval(C[i][k]) * JX[k][j];
            ki += mij * (X[j] - Interval(xhat[j]));
        }
        K[i] = ki;
    }
    return true;
}

/// Krawczyk test. On success K is strictly inside X and encloses the unique
/// zero of G in X for every slice value in c.
inline bool krawczyk(const SsiSystem& sys, const Pt4& m, const Pt4& a, const Interval& c,
                     const Pt4& xhat, const Box4& X, Box4& K) {
    if (!krawczyk_operator(sys, m, a, c, xhat, X, K)) return false;
    for (int i = 0; i < 4; ++i)
        if (!num::strictly_inside(K[i], X[i])) return false;
    return true;
}

/// Refine x onto the slice m.(x-a) = c and return a certified enclosure of
/// the exact intersection point.
inline bool certify_point(const SsiSystem& sys, const Pt4& m, const Pt4& a, double c, Pt4& x, Box4& enc) {
    if (!newton_slice(sys, m, a, c, x)) return false;
    for (double r = 1e-13; r <= 1e-7; r *= 10.0) {
        Box4 X;
        for (int i = 0; i < 4; ++i) X[i] = Interval(x[i] - r, x[i] + r);
        if (krawczyk(sys, m, a, Interval(c), x, X, enc)) return true;
    }
    return false;
}

// --- tubes ---------------------------------------------------------------------

struct TubeSegment {
    Pt4 a{}, b{};     ///< approximate end points; a, m and the c-range define the region exactly
    Pt4 m{};          ///< slicing direction
    double len = 0;   ///< nominal slice value at b
    double c_lo = 0, c_hi = 0;
    Box4 X{};         ///< every slice c in [c_lo, c_hi] meets E exactly once inside X
    bool reversed = false; ///< curve orientation runs from b to a
};

/// Is the box provably inside the segment's certified region?
inline bool in_region(const Box4& P, const TubeSegment& s) {
    Interval c(0.0);
    for (int i = 0; i < 4; ++i) {
        if (!num::subset(P[i], s.X[i])) return false;
        c += Interval(s.m[i]) * (P[i] - Interval(s.a[i]));
    }
    return s.c_lo <= c.lo && c.hi <= s.c_hi;
}

enum class EndKind {
    Boundary, ///< reached the edge of a patch domain
    Closed,   ///< the curve returned to its start
    Singular, ///< certification failed at the smallest step: tangency suspected
};

struct SsiCurve {
    std::vector<TubeSegment> segs;
    bool known = false;  ///< supplied by the caller (for example a side shared by construction), not traced
    bool closed = false;
    EndKind end_front = EndKind::Boundary, end_back = EndKind::Boundary;
    double min_sin_angle = 1.0; ///< worst conditioning seen along the curve
};

struct SsiOptions {
    double h0 = 0.05;      ///< initial step in parameter space
    double h_max = 0.1;
    double h_min = 1e-7;   ///< below this the curve is declared singular
    double overlap = 0.1;  ///< slab overlap between consecutive segments (fraction of length)
    double inflate = 0.3;  ///< box inflation (fraction of length)
    /// Lower bounds on the overlap and inflation, tried first. They keep the
    /// certified region of a very short segment (typically the last step onto
    /// a domain boundary) large enough for the completeness pass to use.
    double min_overlap = 0.004;
    double min_inflate = 0.012;
    int seed_depth = 3;    ///< seeding subdivision: 2^depth cells per parameter
    int verify_depth = 64; ///< maximum bisections in the completeness pass
    std::size_t verify_node_limit = 4000000;
    int reseed_rounds = 3;
    int trace_steps = 50000; ///< accepted steps allowed per traced branch before it is declared singular
    /// The parameter box searched, (u, v) of the first patch and (s, t) of the
    /// second; [0,1]^4 by default. A caller that only needs the part of the
    /// intersection over certain regions of the patches (the faces on them)
    /// can give a smaller box: completeness is then proven over that box, and
    /// curves stop where they leave it.
    std::array<double, 4> dom_lo{0.0, 0.0, 0.0, 0.0};
    std::array<double, 4> dom_hi{1.0, 1.0, 1.0, 1.0};
    bool full_domain() const {
        for (int i = 0; i < 4; ++i)
            if (dom_lo[static_cast<std::size_t>(i)] != 0.0 || dom_hi[static_cast<std::size_t>(i)] != 1.0) return false;
        return true;
    }
};

struct SsiStats {
    std::size_t steps_accepted = 0, steps_rejected = 0;
    std::size_t verify_nodes = 0, excluded = 0, covered = 0;
};

struct SsiResult {
    std::vector<SsiCurve> curves;
    bool complete = false;            ///< proven: every intersection point lies on a returned curve
    std::size_t unverified_boxes = 0; ///< sub-boxes neither excluded nor covered
    double unverified_volume = 0.0;   ///< their total 4D parameter volume
    Box4 unverified_hull{};           ///< where they are
    SsiStats stats;
};

namespace detail {

inline bool inside_domain(const Pt4& x, double tol, const SsiOptions& opt) {
    for (std::size_t i = 0; i < 4; ++i)
        if (x[i] < opt.dom_lo[i] - tol || x[i] > opt.dom_hi[i] + tol) return false;
    return true;
}

inline TubeSegment make_segment(const Pt4& a, const Pt4& b, const SsiOptions& opt, bool generous) {
    TubeSegment s;
    s.a = a;
    s.b = b;
    double n2 = 0.0;
    for (int i = 0; i < 4; ++i) n2 += (b[i] - a[i]) * (b[i] - a[i]);
    const double n = std::sqrt(n2);
    for (int i = 0; i < 4; ++i) s.m[i] = (b[i] - a[i]) / n;
    s.len = 0.0;
    for (int i = 0; i < 4; ++i) s.len += s.m[i] * (b[i] - a[i]);
    const double over = generous ? std::max(opt.overlap * s.len, opt.min_overlap) : opt.overlap * s.len;
    const double r = generous ? std::max(opt.inflate * s.len, opt.min_inflate) : opt.inflate * s.len;
    s.c_lo = -over;
    s.c_hi = s.len + over;
    for (int i = 0; i < 4; ++i)
        s.X[i] = Interval(std::min(a[i], b[i]) - r, std::max(a[i], b[i]) + r);
    return s;
}

inline bool certify_segment(const SsiSystem& sys, const TubeSegment& s) {
    Pt4 mid;
    for (int i = 0; i < 4; ++i) mid[i] = 0.5 * (s.a[i] + s.b[i]);
    Box4 K;
    return krawczyk(sys, s.m, s.a, Interval(s.c_lo, s.c_hi), mid, s.X, K);
}

struct Branch {
    std::vector<TubeSegment> segs;
    EndKind end = EndKind::Boundary;
    double min_sin = 1.0;
};

/// Trace from a certified start point in one direction (dir = +1 or -1).
/// `start_enc` is a certified enclosure of the start point. With
/// `detect_closure` the branch stops when it has left its first segment's
/// region and a certified point of it lies in that region again: by
/// uniqueness of the arc there, the curve has closed.
inline Branch trace(const SsiSystem& sys, const Pt4& start, const Box4& start_enc, double dir,
                    const SsiOptions& opt, SsiStats& stats, bool detect_closure) {
    Branch br;
    Pt4 x = start;
    Box4 junction = start_enc; // certified point shared with the previous segment
    double h = opt.h0;
    TubeSegment first;
    bool have_first = false, left_first = false;

    for (int guard = 0; guard < opt.trace_steps; ++guard) {
        Pt4 t;
        if (!sys.tangent(x, t)) { br.end = EndKind::Singular; return br; }
        for (double& ti : t) ti *= dir;
        br.min_sin = std::min(br.min_sin, sys.sin_angle(x));

        bool accepted = false;
        while (h >= opt.h_min) {
            // Predict along the tangent, correct onto the curve.
            Pt4 xp, y;
            for (int i = 0; i < 4; ++i) xp[i] = x[i] + h * t[i];
            y = xp;
            if (!newton_slice(sys, t, xp, 0.0, y)) { h *= 0.5; ++stats.steps_rejected; continue; }

            // Leaving the domain: stop exactly on the boundary instead.
            bool at_boundary = false;
            if (!inside_domain(y, 0.0, opt)) {
                bool found = false;
                double best = 0.0;
                Pt4 best_z{};
                for (int i = 0; i < 4; ++i) {
                    const double lo = opt.dom_lo[static_cast<std::size_t>(i)], hi = opt.dom_hi[static_cast<std::size_t>(i)];
                    if (y[i] >= lo && y[i] <= hi) continue;
                    const double bound = y[i] < lo ? lo : hi;
                    Pt4 e{};
                    e[i] = 1.0;
                    const double f = (bound - x[i]) / (y[i] - x[i]);
                    Pt4 z;
                    for (int k = 0; k < 4; ++k) z[k] = x[k] + f * (y[k] - x[k]);
                    if (!newton_slice(sys, e, Pt4{}, bound, z)) continue;
                    z[i] = bound;
                    if (!inside_domain(z, 1e-12, opt)) continue;
                    double along = 0.0;
                    for (int k = 0; k < 4; ++k) along += t[k] * (z[k] - x[k]);
                    if (!found || along < best) { found = true; best = along; best_z = z; }
                }
                if (!found) { h *= 0.5; ++stats.steps_rejected; continue; }
                if (best < 1e-10) { br.end = EndKind::Boundary; return br; } // already on the boundary
                y = best_z;
                at_boundary = true;
            }

            // Certify the segment, preferring the larger region.
            TubeSegment seg;
            Box4 end_enc;
            bool ok = false;
            for (bool generous : {true, false}) {
                seg = make_segment(x, y, opt, generous);
                Pt4 y_ref = y;
                ok = seg.len > 0.0 && in_region(junction, seg) && certify_segment(sys, seg) &&
                     certify_point(sys, seg.m, seg.a, seg.len, y_ref, end_enc) &&
                     in_region(end_enc, seg);
                if (ok) break;
            }
            if (!ok) { h *= 0.5; ++stats.steps_rejected; continue; }

            ++stats.steps_accepted;
            seg.reversed = dir < 0.0;
            br.segs.push_back(seg);
            junction = end_enc;
            x = y;
            accepted = true;
            if (at_boundary) { br.end = EndKind::Boundary; return br; }

            double h_next = std::min(h * 1.5, opt.h_max);
            if (detect_closure) {
                if (!have_first) { first = seg; have_first = true; }
                else if (!left_first) { left_first = !in_region(end_enc, first); }
                else if (in_region(end_enc, first)) { br.end = EndKind::Closed; return br; }
                // Coming back around: keep steps shorter than the first
                // region so an end point must land inside it.
                if (left_first) {
                    double d2 = 0.0;
                    for (int i = 0; i < 4; ++i) d2 += (x[i] - start[i]) * (x[i] - start[i]);
                    if (std::sqrt(d2) < 3.0 * opt.h_max)
                        h_next = std::min(h_next, std::max(0.5 * first.len, opt.h_min));
                }
            }
            h = h_next;
            break;
        }
        if (!accepted) { br.end = EndKind::Singular; return br; }
    }
    br.end = EndKind::Singular;
    return br;
}

/// Gauss-Newton projection of a guess onto E (minimum-norm steps), in doubles.
inline bool project_to_curve(const SsiSystem& sys, Pt4& x) {
    for (int it = 0; it < 40; ++it) {
        double F[3], J[3][4];
        sys.eval(x, F, J);
        double M[3][3];
        for (int r = 0; r < 3; ++r)
            for (int c = 0; c < 3; ++c) {
                M[r][c] = 0.0;
                for (int k = 0; k < 4; ++k) M[r][c] += J[r][k] * J[c][k];
            }
        const double d = det3(M);
        if (!(std::fabs(d) > 1e-300)) return false;
        // Solve M y = F by Cramer's rule.
        double y[3];
        for (int c = 0; c < 3; ++c) {
            double Mc[3][3];
            for (int r = 0; r < 3; ++r)
                for (int k = 0; k < 3; ++k) Mc[r][k] = k == c ? F[r] : M[r][k];
            y[c] = det3(Mc) / d;
        }
        double step = 0.0;
        for (int k = 0; k < 4; ++k) {
            double dx = 0.0;
            for (int r = 0; r < 3; ++r) dx -= J[r][k] * y[r];
            x[k] += dx;
            step = std::max(step, std::fabs(dx));
        }
        if (!std::isfinite(step) || step > 2.0) return false;
        if (step < 1e-14) return true;
    }
    return false;
}

/// Collect seed guesses: centres of sub-box pairs whose bounding boxes overlap.
inline void collect_seeds(const PolyPatchD& A, double u0, double u1, double v0, double v1,
                          const PolyPatchD& B, double s0, double s1, double t0, double t1, int depth,
                          std::vector<Pt4>& seeds) {
    Vec3d alo, ahi, blo, bhi;
    geom::projected_bounds(A, alo, ahi);
    geom::projected_bounds(B, blo, bhi);
    if (alo.x > bhi.x || blo.x > ahi.x || alo.y > bhi.y || blo.y > ahi.y || alo.z > bhi.z || blo.z > ahi.z)
        return;
    if (depth == 0) {
        seeds.push_back({0.5 * (u0 + u1), 0.5 * (v0 + v1), 0.5 * (s0 + s1), 0.5 * (t0 + t1)});
        return;
    }
    PolyPatchD Al, Ah, All, Alh, Ahl, Ahh, Bl, Bh, Bll, Blh, Bhl, Bhh;
    A.split_u(0.5, Al, Ah);
    Al.split_v(0.5, All, Alh);
    Ah.split_v(0.5, Ahl, Ahh);
    B.split_u(0.5, Bl, Bh);
    Bl.split_v(0.5, Bll, Blh);
    Bh.split_v(0.5, Bhl, Bhh);
    const double um = 0.5 * (u0 + u1), vm = 0.5 * (v0 + v1), sm = 0.5 * (s0 + s1), tm = 0.5 * (t0 + t1);
    const PolyPatchD* As[4] = {&All, &Alh, &Ahl, &Ahh};
    const PolyPatchD* Bs[4] = {&Bll, &Blh, &Bhl, &Bhh};
    const double ab[4][4] = {{u0, um, v0, vm}, {u0, um, vm, v1}, {um, u1, v0, vm}, {um, u1, vm, v1}};
    const double bb[4][4] = {{s0, sm, t0, tm}, {s0, sm, tm, t1}, {sm, s1, t0, tm}, {sm, s1, tm, t1}};
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            collect_seeds(*As[i], ab[i][0], ab[i][1], ab[i][2], ab[i][3], *Bs[j], bb[j][0], bb[j][1],
                          bb[j][2], bb[j][3], depth - 1, seeds);
}

inline bool boxes_disjoint(const Vec3<Interval>& a, const Vec3<Interval>& b) {
    return num::disjoint(a.x, b.x) || num::disjoint(a.y, b.y) || num::disjoint(a.z, b.z);
}

/// Approximate normal of a patch from its corner control points (doubles;
/// only used to choose a direction, so it need not be rigorous).
inline Vec3d corner_normal(const PolyPatchI& P) {
    auto pt = [&](int i, int j) {
        const auto& c = P.at(i, j);
        const double w = num::mid(c.w);
        return Vec3d{num::mid(c.x) / w, num::mid(c.y) / w, num::mid(c.z) / w};
    };
    const Vec3d d1 = (pt(P.p, 0) - pt(0, 0)) + (pt(P.p, P.q) - pt(0, P.q));
    const Vec3d d2 = (pt(0, P.q) - pt(0, 0)) + (pt(P.p, P.q) - pt(P.p, 0));
    return cross(d1, d2);
}

/// Rigorous range of n . S over the patch, for any fixed direction n: the
/// surface point is a convex combination of the projected control points.
inline Interval extent_along(const PolyPatchI& P, const Vec3d& n) {
    Interval r;
    bool first = true;
    for (const auto& c : P.c) {
        if (!(c.w.lo > 0.0)) return Interval::entire();
        const Interval d = (Interval(n.x) * c.x + Interval(n.y) * c.y + Interval(n.z) * c.z) / c.w;
        r = first ? d : num::hull(r, d);
        first = false;
    }
    return r;
}

/// Separating-axis test along the two patch normals. Two nearly parallel
/// pieces that are close together have overlapping axis-aligned boxes, but
/// along their common normal each piece is thin, so this separates them at a
/// much coarser subdivision level.
inline bool separated_along_normals(const PolyPatchI& A, const PolyPatchI& B) {
    for (const PolyPatchI* P : {&A, &B}) {
        const Vec3d n = corner_normal(*P);
        const double len = geom::norm(n);
        if (!(len > 0.0) || !std::isfinite(len)) continue;
        const Vec3d dir = n / len;
        if (num::disjoint(extent_along(A, dir), extent_along(B, dir))) return true;
    }
    return false;
}

struct Verifier {
    const std::vector<const TubeSegment*>& all;
    const SsiOptions& opt;
    SsiResult& res;
    std::vector<Box4> unverified;
    bool node_limit_hit = false;

    void run(const PolyPatchI& A, const Vec3<Interval>& boundsA, const PolyPatchI& B,
             const Vec3<Interval>& boundsB, const Box4& box, const std::vector<int>& cand, int depth) {
        ++res.stats.verify_nodes;
        if (boxes_disjoint(boundsA, boundsB)) { ++res.stats.excluded; return; }

        std::vector<int> near;
        near.reserve(cand.size());
        for (int k : cand) {
            const TubeSegment& s = *all[static_cast<std::size_t>(k)];
            bool overlap = true;
            for (int i = 0; i < 4 && overlap; ++i) overlap = !num::disjoint(box[i], s.X[i]);
            if (!overlap) continue;
            if (in_region(box, s)) { ++res.stats.covered; return; }
            near.push_back(k);
        }

        if (separated_along_normals(A, B)) { ++res.stats.excluded; return; }

        if (depth >= opt.verify_depth || res.stats.verify_nodes >= opt.verify_node_limit) {
            if (res.stats.verify_nodes >= opt.verify_node_limit) node_limit_hit = true;
            unverified.push_back(box);
            return;
        }

        // Bisect the widest parameter.
        int dim = 0;
        double wmax = 0.0;
        for (int i = 0; i < 4; ++i) {
            const double w = box[i].hi - box[i].lo;
            if (w > wmax) { wmax = w; dim = i; }
        }
        const double midv = 0.5 * (box[dim].lo + box[dim].hi);
        Box4 lo = box, hi = box;
        lo[dim].hi = midv;
        hi[dim].lo = midv;
        PolyPatchI Pl, Ph;
        if (dim < 2) {
            if (dim == 0) A.split_u(Interval(0.5), Pl, Ph); else A.split_v(Interval(0.5), Pl, Ph);
            run(Pl, geom::projected_bounds(Pl), B, boundsB, lo, near, depth + 1);
            run(Ph, geom::projected_bounds(Ph), B, boundsB, hi, near, depth + 1);
        } else {
            if (dim == 2) B.split_u(Interval(0.5), Pl, Ph); else B.split_v(Interval(0.5), Pl, Ph);
            run(A, boundsA, Pl, geom::projected_bounds(Pl), lo, near, depth + 1);
            run(A, boundsA, Ph, geom::projected_bounds(Ph), hi, near, depth + 1);
        }
    }
};

} // namespace detail

/// Intersect two rational Bezier patches over [0,1]^2 x [0,1]^2.
///
/// `known` are curves of the intersection that the caller already has as
/// certified segment chains, typically a side the two patches share by
/// construction. Such a curve lies on the boundary of both domains, where
/// no subdivision can separate the patches; its certified tube covers that
/// neighbourhood in the completeness proof instead. Known curves are
/// returned with `known` set and are not traced again.
inline SsiResult intersect(const SsiSystem& sys, const SsiOptions& opt = {},
                           const std::vector<SsiCurve>& known = {}) {
    SsiResult res;
    res.curves = known;
    for (SsiCurve& c : res.curves) c.known = true;

    auto already_traced = [&](const Box4& enc) {
        for (const auto& c : res.curves)
            for (const auto& s : c.segs)
                if (in_region(enc, s)) return true;
        return false;
    };

    auto trace_from = [&](Pt4 x) {
        if (!detail::project_to_curve(sys, x) || !detail::inside_domain(x, 0.0, opt)) return;
        Pt4 t;
        if (!sys.tangent(x, t)) return;
        Box4 enc;
        const Pt4 anchor = x;
        if (!certify_point(sys, t, anchor, 0.0, x, enc)) return;
        if (!detail::inside_domain(x, 0.0, opt)) return;
        if (already_traced(enc)) return;

        detail::Branch fwd = detail::trace(sys, x, enc, +1.0, opt, res.stats, true);
        SsiCurve curve;
        curve.min_sin_angle = fwd.min_sin;
        if (fwd.end == EndKind::Closed) {
            curve.segs = std::move(fwd.segs);
            curve.closed = true;
            curve.end_front = curve.end_back = EndKind::Closed;
        } else {
            detail::Branch back = detail::trace(sys, x, enc, -1.0, opt, res.stats, false);
            curve.min_sin_angle = std::min(curve.min_sin_angle, back.min_sin);
            std::reverse(back.segs.begin(), back.segs.end());
            curve.segs = std::move(back.segs);
            curve.segs.insert(curve.segs.end(), fwd.segs.begin(), fwd.segs.end());
            curve.end_front = back.end;
            curve.end_back = fwd.end;
        }
        if (!curve.segs.empty()) res.curves.push_back(std::move(curve));
    };

    const auto& lo = opt.dom_lo;
    const auto& hi = opt.dom_hi;
    const bool full = opt.full_domain();
    std::vector<Pt4> seeds;
    if (full) {
        detail::collect_seeds(sys.patch_a(), 0, 1, 0, 1, sys.patch_b(), 0, 1, 0, 1, opt.seed_depth, seeds);
    } else {
        detail::collect_seeds(sys.patch_a().restricted(lo[0], hi[0], lo[1], hi[1]), lo[0], hi[0], lo[1], hi[1],
                              sys.patch_b().restricted(lo[2], hi[2], lo[3], hi[3]), lo[2], hi[2], lo[3], hi[3],
                              opt.seed_depth, seeds);
    }
    for (const Pt4& s : seeds) trace_from(s);

    // Completeness pass, with re-seeding from whatever it could not resolve.
    // Over a smaller domain the verifier starts from the patches restricted to
    // it (interval restriction: enclosures of the exact pieces).
    const Box4 domain = {Interval(lo[0], hi[0]), Interval(lo[1], hi[1]), Interval(lo[2], hi[2]), Interval(lo[3], hi[3])};
    const PolyPatchI Ai = full ? sys.patch_ai() : sys.patch_ai().restricted(lo[0], hi[0], lo[1], hi[1]);
    const PolyPatchI Bi = full ? sys.patch_bi() : sys.patch_bi().restricted(lo[2], hi[2], lo[3], hi[3]);
    for (int round = 0;; ++round) {
        std::vector<const TubeSegment*> all;
        for (const auto& c : res.curves)
            for (const auto& s : c.segs) all.push_back(&s);
        std::vector<int> cand(all.size());
        for (std::size_t i = 0; i < all.size(); ++i) cand[i] = static_cast<int>(i);

        res.stats.verify_nodes = res.stats.excluded = res.stats.covered = 0;
        detail::Verifier ver{all, opt, res, {}, false};
        ver.run(Ai, geom::projected_bounds(Ai), Bi, geom::projected_bounds(Bi), domain, cand, 0);

        res.unverified_boxes = ver.unverified.size();
        res.unverified_volume = 0.0;
        for (std::size_t i = 0; i < ver.unverified.size(); ++i) {
            const Box4& b = ver.unverified[i];
            double vol = 1.0;
            for (int k = 0; k < 4; ++k) vol *= b[k].hi - b[k].lo;
            res.unverified_volume += vol;
            if (i == 0) res.unverified_hull = b;
            else
                for (int k = 0; k < 4; ++k) res.unverified_hull[k] = num::hull(res.unverified_hull[k], b[k]);
        }
        res.complete = ver.unverified.empty();
        if (res.complete || round >= opt.reseed_rounds) break;

        const std::size_t before = res.curves.size();
        const std::size_t stride = std::max<std::size_t>(1, ver.unverified.size() / 64);
        for (std::size_t i = 0; i < ver.unverified.size(); i += stride) {
            const Box4& b = ver.unverified[i];
            trace_from({num::mid(b[0]), num::mid(b[1]), num::mid(b[2]), num::mid(b[3])});
        }
        if (res.curves.size() == before) break; // nothing new to learn
    }
    return res;
}

/// Certified point on a segment at local parameter s in [0,1] along the curve.
inline bool curve_point(const SsiSystem& sys, const TubeSegment& seg, double s, Box4& enc) {
    const double f = seg.reversed ? 1.0 - s : s;
    Pt4 x;
    for (int i = 0; i < 4; ++i) x[i] = seg.a[i] + f * (seg.b[i] - seg.a[i]);
    return certify_point(sys, seg.m, seg.a, f * seg.len, x, enc) && in_region(enc, seg);
}

} // namespace ak::solve
