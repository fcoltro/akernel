// akernel — certified intersection of procedural carriers, and the spine of
// a rolling-ball blend (docs/FEATURES.md, steps F3b and §5).
//
// CarrierPairSystem is the surface/surface system of solve/ssi.hpp for two
// carriers (geom/carrier.hpp): F(u, v, s, t) = A(u, v) - B(s, t), with the
// Jacobian [A_u A_v -B_s -B_t]. It supplies what the templated Newton,
// Krawczyk and tracing routines of ssi.hpp use, so the curves it traces
// carry the same certificates: every slice of a tube segment meets the
// curve exactly once inside the segment's box.
//
// The spine. A ball of radius r rolling between faces A and B (normals
// pointing to the ball's side) has its centre on both offsets at distance
// r: the spine is the intersection of OffsetCarrier(A, r) and
// OffsetCarrier(B, r), a curve in (u, v, s, t). At each spine point the
// contacts are A(u, v) and B(s, t), and the blend's cross-section is the
// circular arc between them (geom/blend.hpp, fillet_arc). The rolling-ball
// blend is therefore defined by the certified spine and that arc, with no
// fitted surface: the procedural blend carrier.
//
// intersect_carriers() finds the curves without a guess and proves that no
// intersection point in the domain is missed (or reports the boxes it could
// not resolve), like ssi.hpp's intersect() for patches; find_spines()
// applies it to two offsets. trace_spine() follows one branch from a guess.
#pragma once

#include <algorithm>
#include <optional>
#include <vector>

#include <akernel/geom/blend.hpp>
#include <akernel/geom/carrier.hpp>
#include <akernel/solve/ssi.hpp>

namespace ak::solve {

class CarrierPairSystem {
    const geom::Carrier& a_;
    const geom::Carrier& b_;

public:
    CarrierPairSystem(const geom::Carrier& a, const geom::Carrier& b) : a_(a), b_(b) {}

    void eval(const Pt4& x, double F[3], double J[3][4]) const {
        Vec3d A, Au, Av, B, Bs, Bt;
        a_.eval(x[0], x[1], A, Au, Av);
        b_.eval(x[2], x[3], B, Bs, Bt);
        const double a[3] = {A.x, A.y, A.z}, au[3] = {Au.x, Au.y, Au.z}, av[3] = {Av.x, Av.y, Av.z};
        const double b[3] = {B.x, B.y, B.z}, bs[3] = {Bs.x, Bs.y, Bs.z}, bt[3] = {Bt.x, Bt.y, Bt.z};
        for (int k = 0; k < 3; ++k) {
            F[k] = a[k] - b[k];
            J[k][0] = au[k];
            J[k][1] = av[k];
            J[k][2] = -bs[k];
            J[k][3] = -bt[k];
        }
    }

    std::array<Interval, 3> F_at(const Pt4& x) const {
        const auto A = a_.point_at(x[0], x[1]);
        const auto B = b_.point_at(x[2], x[3]);
        return {A.x - B.x, A.y - B.y, A.z - B.z};
    }

    void J_over(const Box4& X, Interval J[3][4]) const {
        geom::IVec3 A, Au, Av, B, Bs, Bt;
        const bool ok = a_.jet1_over(X[0].lo, X[0].hi, X[1].lo, X[1].hi, A, Au, Av) &&
                        b_.jet1_over(X[2].lo, X[2].hi, X[3].lo, X[3].hi, B, Bs, Bt);
        for (int k = 0; k < 3; ++k)
            for (int j = 0; j < 4; ++j) J[k][j] = Interval::entire();
        if (!ok) return; // unbounded Jacobian: every Krawczyk test on X fails
        const Interval au[3] = {Au.x, Au.y, Au.z}, av[3] = {Av.x, Av.y, Av.z};
        const Interval bs[3] = {Bs.x, Bs.y, Bs.z}, bt[3] = {Bt.x, Bt.y, Bt.z};
        for (int k = 0; k < 3; ++k) {
            J[k][0] = au[k];
            J[k][1] = av[k];
            J[k][2] = -bs[k];
            J[k][3] = -bt[k];
        }
    }

    double sin_angle(const Pt4& x) const {
        Vec3d A, Au, Av, B, Bs, Bt;
        a_.eval(x[0], x[1], A, Au, Av);
        b_.eval(x[2], x[3], B, Bs, Bt);
        const Vec3d n1 = cross(Au, Av), n2 = cross(Bs, Bt);
        const double d = geom::norm(n1) * geom::norm(n2);
        return d > 0.0 ? geom::norm(cross(n1, n2)) / d : 0.0;
    }

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

namespace carrier_detail {

inline bool disjoint3(const geom::IVec3& a, const geom::IVec3& b) {
    return num::disjoint(a.x, b.x) || num::disjoint(a.y, b.y) || num::disjoint(a.z, b.z);
}

/// Rigorous enclosure of a carrier's points over a parameter box.
inline bool point_over(const geom::Carrier& c, const Interval& u, const Interval& v, geom::IVec3& P,
                       geom::IVec3& Pu, geom::IVec3& Pv) {
    return c.jet1_over(u.lo, u.hi, v.lo, v.hi, P, Pu, Pv);
}

/// Range of dir . P over the box, from the enclosures of P (rigorous).
inline Interval along(const geom::IVec3& P, const Vec3d& dir) {
    return Interval(dir.x) * P.x + Interval(dir.y) * P.y + Interval(dir.z) * P.z;
}

/// Separation along an approximate normal of either piece: the mean-value
/// form gives dir . P over the box as dir . P(c) + (dir . P_u) du + ..., and
/// near a normal direction dir . P_u and dir . P_v are small, so thin pieces
/// separate at a coarse level (the analogue of ssi.hpp's
/// separated_along_normals).
inline bool separated_along_normals(const geom::Carrier& a, const Box4& X, const geom::Carrier& b) {
    const double uc = num::mid(X[0]), vc = num::mid(X[1]), sc = num::mid(X[2]), tc = num::mid(X[3]);
    Vec3d pa, pau, pav, pb, pbs, pbt;
    a.eval(uc, vc, pa, pau, pav);
    b.eval(sc, tc, pb, pbs, pbt);
    geom::IVec3 A, Au, Av, B, Bs, Bt;
    if (!a.jet1_over(X[0].lo, X[0].hi, X[1].lo, X[1].hi, A, Au, Av)) return false;
    if (!b.jet1_over(X[2].lo, X[2].hi, X[3].lo, X[3].hi, B, Bs, Bt)) return false;
    for (const Vec3d& n : {cross(pau, pav), cross(pbs, pbt)}) {
        const double len = geom::norm(n);
        if (!(len > 0.0) || !std::isfinite(len)) continue;
        const Vec3d dir = n / len;
        // Mean-value ranges of dir . A and dir . B around the box centres.
        auto mv = [&](const geom::Carrier& c, double uu, double vv, const geom::IVec3& Du,
                      const geom::IVec3& Dv, const Interval& U, const Interval& V) {
            const Interval c0 = along(c.point_at(uu, vv), dir);
            const Interval eu = U - Interval(uu), ev = V - Interval(vv);
            return c0 + along(Du, dir) * eu + along(Dv, dir) * ev;
        };
        const Interval ra = mv(a, uc, vc, Au, Av, X[0], X[1]);
        const Interval rb = mv(b, sc, tc, Bs, Bt, X[2], X[3]);
        if (num::is_finite(ra) && num::is_finite(rb) && num::disjoint(ra, rb)) return true;
    }
    return false;
}

inline void collect_seeds(const geom::Carrier& a, const geom::Carrier& b, const Box4& X, int depth,
                          std::vector<Pt4>& seeds) {
    geom::IVec3 A, Au, Av, B, Bs, Bt;
    if (a.jet1_over(X[0].lo, X[0].hi, X[1].lo, X[1].hi, A, Au, Av) &&
        b.jet1_over(X[2].lo, X[2].hi, X[3].lo, X[3].hi, B, Bs, Bt) && disjoint3(A, B))
        return;
    if (depth == 0) {
        seeds.push_back({num::mid(X[0]), num::mid(X[1]), num::mid(X[2]), num::mid(X[3])});
        return;
    }
    for (int m = 0; m < 16; ++m) {
        Box4 Y = X;
        for (int i = 0; i < 4; ++i) {
            const double mid = num::mid(X[i]);
            Y[i] = (m >> i) & 1 ? Interval(mid, X[i].hi) : Interval(X[i].lo, mid);
        }
        collect_seeds(a, b, Y, depth - 1, seeds);
    }
}

struct Verifier {
    const geom::Carrier& a;
    const geom::Carrier& b;
    const std::vector<const TubeSegment*>& all;
    const SsiOptions& opt;
    SsiResult& res;
    std::vector<Box4> unverified;

    void run(const Box4& box, const std::vector<int>& cand, int depth) {
        ++res.stats.verify_nodes;
        geom::IVec3 A, Au, Av, B, Bs, Bt;
        const bool ok = a.jet1_over(box[0].lo, box[0].hi, box[1].lo, box[1].hi, A, Au, Av) &&
                        b.jet1_over(box[2].lo, box[2].hi, box[3].lo, box[3].hi, B, Bs, Bt);
        if (ok && disjoint3(A, B)) { ++res.stats.excluded; return; }

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
        if (ok && separated_along_normals(a, box, b)) { ++res.stats.excluded; return; }

        if (depth >= opt.verify_depth || res.stats.verify_nodes >= opt.verify_node_limit) {
            unverified.push_back(box);
            return;
        }
        int dim = 0;
        double wmax = 0.0;
        for (int i = 0; i < 4; ++i) {
            const double w = box[i].hi - box[i].lo;
            if (w > wmax) { wmax = w; dim = i; }
        }
        const double midv = num::mid(box[dim]);
        Box4 lo = box, hi = box;
        lo[dim].hi = midv;
        hi[dim].lo = midv;
        run(lo, near, depth + 1);
        run(hi, near, depth + 1);
    }
};

} // namespace carrier_detail

/// Intersect two carriers over the domain of `opt` (default [0,1]^4): seeds
/// from a subdivision of the parameter box, certified tracing (as in
/// ssi.hpp), then a completeness pass that proves every intersection point
/// in the domain lies on a returned curve, or reports the boxes it could not
/// resolve. Same result type and guarantees as intersect() for patches.
inline SsiResult intersect_carriers(const geom::Carrier& a, const geom::Carrier& b,
                                    const SsiOptions& opt = {}) {
    const CarrierPairSystem sys(a, b);
    SsiResult res;
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
        if (!detail::inside_domain(x, 0.0, opt) || already_traced(enc)) return;
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

    const Box4 domain = {Interval(opt.dom_lo[0], opt.dom_hi[0]), Interval(opt.dom_lo[1], opt.dom_hi[1]),
                         Interval(opt.dom_lo[2], opt.dom_hi[2]), Interval(opt.dom_lo[3], opt.dom_hi[3])};
    std::vector<Pt4> seeds;
    carrier_detail::collect_seeds(a, b, domain, opt.seed_depth, seeds);
    for (const Pt4& s : seeds) trace_from(s);

    for (int round = 0;; ++round) {
        std::vector<const TubeSegment*> all;
        for (const auto& c : res.curves)
            for (const auto& s : c.segs) all.push_back(&s);
        std::vector<int> cand(all.size());
        for (std::size_t i = 0; i < all.size(); ++i) cand[i] = static_cast<int>(i);
        res.stats.verify_nodes = res.stats.excluded = res.stats.covered = 0;
        carrier_detail::Verifier ver{a, b, all, opt, res, {}};
        ver.run(domain, cand, 0);
        res.unverified_boxes = ver.unverified.size();
        res.unverified_volume = 0.0;
        for (std::size_t i = 0; i < ver.unverified.size(); ++i) {
            const Box4& bx = ver.unverified[i];
            double vol = 1.0;
            for (int k = 0; k < 4; ++k) vol *= bx[k].hi - bx[k].lo;
            res.unverified_volume += vol;
            if (i == 0) res.unverified_hull = bx;
            else
                for (int k = 0; k < 4; ++k) res.unverified_hull[k] = num::hull(res.unverified_hull[k], bx[k]);
        }
        res.complete = ver.unverified.empty();
        if (res.complete || round >= opt.reseed_rounds) break;
        const std::size_t before = res.curves.size();
        const std::size_t stride = std::max<std::size_t>(1, ver.unverified.size() / 64);
        for (std::size_t i = 0; i < ver.unverified.size(); i += stride) {
            const Box4& bx = ver.unverified[i];
            trace_from({num::mid(bx[0]), num::mid(bx[1]), num::mid(bx[2]), num::mid(bx[3])});
        }
        if (res.curves.size() == before) break;
    }
    return res;
}

/// Every spine of a ball of radius r between patches A and B (normals
/// pointing to the ball's side), found without a guess, with the
/// completeness proof of intersect_carriers.
inline SsiResult find_spines(const PolyPatchD& A, const PolyPatchD& B, double r,
                             const SsiOptions& opt = {}) {
    const geom::OffsetCarrier oa(A, r), ob(B, r);
    return intersect_carriers(oa, ob, opt);
}

/// A certified piece of a rolling-ball spine.
struct Spine {
    std::vector<TubeSegment> segs;  ///< in (u, v, s, t); each certified as in ssi.hpp
    EndKind end_front = EndKind::Singular, end_back = EndKind::Singular;
    double min_sin_angle = 1.0;     ///< worst angle between the two offsets along the spine
};

/// Trace the spine of a ball of radius r between patches A and B whose
/// normals (S_u x S_v) point to the ball's side, from a guess x0 in
/// (u, v, s, t), in both directions until it leaves `opt`'s domain. Returns
/// nothing if the guess does not converge to a certified spine point.
inline std::optional<Spine> trace_spine(const PolyPatchD& A, const PolyPatchD& B, double r,
                                        Pt4 x0, const SsiOptions& opt = {}) {
    const geom::OffsetCarrier oa(A, r), ob(B, r);
    const CarrierPairSystem sys(oa, ob);
    if (!detail::project_to_curve(sys, x0)) return std::nullopt;
    Pt4 t;
    if (!sys.tangent(x0, t)) return std::nullopt;
    Pt4 x = x0;
    Box4 enc;
    if (!certify_point(sys, t, x0, 0.0, x, enc)) return std::nullopt;
    SsiStats stats;
    const auto fwd = detail::trace(sys, x, enc, +1.0, opt, stats, false);
    const auto back = detail::trace(sys, x, enc, -1.0, opt, stats, false);
    Spine sp;
    sp.segs.assign(back.segs.rbegin(), back.segs.rend());
    sp.segs.insert(sp.segs.end(), fwd.segs.begin(), fwd.segs.end());
    sp.end_back = back.end;
    sp.end_front = fwd.end;
    sp.min_sin_angle = std::min(fwd.min_sin, back.min_sin);
    return sp;
}

/// The blend's cross-section at a spine point: the arc from the contact on
/// A to the contact on B about the ball's centre (in doubles; the spine
/// point's certificate is in the tube segment that holds it).
inline geom::FilletArc spine_section(const PolyPatchD& A, const PolyPatchD& B, double r,
                                     const Pt4& x) {
    const geom::PatchCarrier pa(A), pb(B);
    Vec3d a, au, av, b, bs, bt;
    pa.eval(x[0], x[1], a, au, av);
    pb.eval(x[2], x[3], b, bs, bt);
    const Vec3d na = geom::normalized(cross(au, av)), nb = geom::normalized(cross(bs, bt));
    const Vec3d c = a + na * r;
    return geom::fillet_arc(c, r, -na, -nb);
}

} // namespace ak::solve
