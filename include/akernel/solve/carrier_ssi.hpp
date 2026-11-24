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
#include <atomic>
#include <cstdint>
#include <map>
#include <memory>
#include <mutex>
#include <optional>
#include <tuple>
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
/// The enclosures over X may be passed in when the caller has them. The
/// direction between the two centre points is tried as well: where the
/// pieces are apart but neither is flat across the gap, it is the one that
/// separates them.
inline bool separated_along_normals(const geom::Carrier& a, const Box4& X, const geom::Carrier& b,
                                    const geom::IVec3* jets = nullptr) {
    const double uc = num::mid(X[0]), vc = num::mid(X[1]), sc = num::mid(X[2]), tc = num::mid(X[3]);
    Vec3d pa, pau, pav, pb, pbs, pbt;
    a.eval(uc, vc, pa, pau, pav);
    b.eval(sc, tc, pb, pbs, pbt);
    geom::IVec3 Au, Av, Bs, Bt;
    if (jets) {
        Au = jets[0]; Av = jets[1]; Bs = jets[2]; Bt = jets[3];
    } else {
        geom::IVec3 A, B;
        if (!a.jet1_over(X[0].lo, X[0].hi, X[1].lo, X[1].hi, A, Au, Av)) return false;
        if (!b.jet1_over(X[2].lo, X[2].hi, X[3].lo, X[3].hi, B, Bs, Bt)) return false;
    }
    for (const Vec3d& n : {cross(pau, pav), cross(pbs, pbt), pa - pb}) {
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
        const geom::IVec3 jets[4] = {Au, Av, Bs, Bt};
        if (ok && separated_along_normals(a, box, b, jets)) { ++res.stats.excluded; return; }

        if (depth >= opt.verify_depth || res.stats.verify_nodes >= opt.verify_node_limit) {
            unverified.push_back(box);
            return;
        }
        // Split where the box is widest in space: parameter width times the
        // largest speed in that direction (plain width if not bounded).
        double speed[4] = {1.0, 1.0, 1.0, 1.0};
        if (ok) {
            auto mag = [](const geom::IVec3& d) {
                const double x = num::mag(d.x), y = num::mag(d.y), z = num::mag(d.z);
                return std::sqrt(x * x + y * y + z * z);
            };
            const double sp[4] = {mag(Au), mag(Av), mag(Bs), mag(Bt)};
            bool finite = true;
            for (double v : sp) finite = finite && std::isfinite(v) && v > 0.0;
            if (finite)
                for (int i = 0; i < 4; ++i) speed[i] = sp[i];
        }
        int dim = 0;
        double wmax = 0.0;
        for (int i = 0; i < 4; ++i) {
            const double w = (box[i].hi - box[i].lo) * speed[i];
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

// --- the rolling-ball blend as a carrier ----------------------------------------
//
// BlendCarrier is the blend surface of one certified spine: sigma in [0, 1]
// runs along the spine (segment k of n covers [k/n, (k+1)/n], its slice
// value linear in sigma), theta in [0, 1] across the section, from the
// contact on A (theta = 0) to the contact on B (theta = 1). A point is the
// section arc of spine_section() at theta. Nothing is fitted: the surface is
// defined by the two supports, the radius and the spine's certificate.
//
// Enclosures over a (sigma, theta) box, per spine segment it overlaps:
//  1. Spine parameters. The Krawczyk operator of the segment over its box
//     for the slice values of the sigma range gives K, holding the spine
//     points of those slices (ssi.hpp, krawczyk_operator).
//  2. Spine derivative. x' = dx/dc solves M x' = e4, M = [J(x); m]. With C
//     approximately inv(mid M(K)) and E = I - C M(K), if q = |E|_inf < 1 then
//     |x'|_inf <= b = |C e4|_inf / (1 - q) and x' lies in C e4 + E [-b, b]^4.
//  3. Contacts, centre and their sigma-derivatives from the patch and offset
//     carriers' enclosures over K by the chain rule; the arc in forward-mode
//     interval arithmetic (value, d/dsigma, d/dtheta).
//  4. Tightening. Each quantity q of step 3 also lies in q(sigma_mid) +
//     q' [-h, h] (mean-value form in sigma); the box result likewise around
//     the box centre in (sigma, theta). This keeps the correlation between
//     contacts and centre that the direct form loses.
// The theta-independent part (steps 1 to 4) is cached per sigma range, since
// subdivision revisits each range with many theta boxes.
namespace carrier_detail {

using geom::carrier_detail::meet;
using D2 = geom::offset_detail::D2;
using V2 = Vec3<D2>;

inline V2 lift(const geom::IVec3& v, const geom::IVec3& dv_dsigma) {
    return {D2(v.x, dv_dsigma.x, Interval(0.0)), D2(v.y, dv_dsigma.y, Interval(0.0)),
            D2(v.z, dv_dsigma.z, Interval(0.0))};
}

/// Rigorous enclosure of dx/dc over K for the sliced system (see above).
inline bool spine_derivative(const CarrierPairSystem& sys, const Pt4& m, const Box4& K,
                             std::array<Interval, 4>& xp) {
    Interval M[4][4];
    sys.J_over(K, M);
    for (int k = 0; k < 4; ++k) M[3][k] = Interval(m[k]);
    double Mm[4][4], C[4][4];
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) {
            if (!num::is_finite(M[i][j])) return false;
            Mm[i][j] = num::mid(M[i][j]);
        }
    if (!invert4(Mm, C)) return false;
    Interval E[4][4];
    Interval q(0.0);
    double qmax = 0.0;
    for (int i = 0; i < 4; ++i) {
        Interval row(0.0);
        for (int j = 0; j < 4; ++j) {
            Interval e(i == j ? 1.0 : 0.0);
            for (int k = 0; k < 4; ++k) e -= Interval(C[i][k]) * M[k][j];
            E[i][j] = e;
            row += Interval(num::mag(e));
        }
        qmax = std::max(qmax, row.hi);
    }
    (void)q;
    if (!(qmax < 1.0)) return false;
    double ce = 0.0;
    for (int i = 0; i < 4; ++i) ce = std::max(ce, std::fabs(C[i][3]));
    const double b = (Interval(ce) / (Interval(1.0) - Interval(qmax))).hi;
    const Interval Y(-b, b);
    for (int i = 0; i < 4; ++i) {
        Interval v(C[i][3]);
        for (int j = 0; j < 4; ++j) v += E[i][j] * Y;
        xp[static_cast<std::size_t>(i)] = v;
    }
    return true;
}

/// The arc point and its derivatives from lifted contacts, centre and the
/// unit directions uu = (a - c) / r, vv = (b - c) / r (passed separately so
/// that their enclosures can be tighter than the difference's).
inline V2 arc_point(const V2& a, const V2& c, const V2& b, const V2& uu, const V2& vv, double r,
                    const D2& theta, bool& ok) {
    const D2 R(r), one(1.0), two(2.0);
    const D2 k = dot(uu, vv);
    const D2 onek = one + k;
    ok = onek.v.lo > 0.0;
    if (!ok) return {};
    const V2 corner = c + (uu + vv) * (R / onek);
    const D2 w2 = two / onek;
    const D2 s = one - theta;
    const D2 B0 = s * s, B1 = two * s * theta, B2 = theta * theta * w2;
    return (a * B0 + corner * B1 + b * B2) / (B0 + B1 + B2);
}

} // namespace carrier_detail

class BlendCarrier final : public geom::Carrier {
    geom::PatchCarrier pa_, pb_;
    geom::OffsetCarrier oa_, ob_;
    double r_;
    std::vector<TubeSegment> segs_; ///< in curve order
    /// sigma at the segments' ends: brk_[k] to brk_[k + 1] is segment k, a
    /// share of [0, 1] proportional to the length of its centre path, so
    /// the speed in sigma is nearly even along the blend (a traced curve can
    /// end in a very short segment).
    std::vector<double> brk_;

    CarrierPairSystem sys() const { return CarrierPairSystem(oa_, ob_); }
    int count() const { return static_cast<int>(segs_.size()); }
    /// Slice value of segment k at local fraction f in [0, 1] along the curve.
    double slice(int k, double f) const {
        const TubeSegment& s = segs_[static_cast<std::size_t>(k)];
        return (s.reversed ? 1.0 - f : f) * s.len;
    }
    int segment_of(double sigma) const {
        const auto it = std::upper_bound(brk_.begin() + 1, brk_.end() - 1, sigma);
        return static_cast<int>(it - brk_.begin()) - 1;
    }
    /// Local fraction of sigma along segment k (doubles: the map sigma -> slice
    /// is defined by these operations, and enclosures are taken over the
    /// slices it gives).
    double local(int k, double sigma) const {
        const double a = brk_[static_cast<std::size_t>(k)], b = brk_[static_cast<std::size_t>(k) + 1];
        return std::clamp((sigma - a) / (b - a), 0.0, 1.0);
    }
    void make_breaks() {
        const int n = count();
        std::vector<double> w(static_cast<std::size_t>(n), 1.0);
        double total = 0.0;
        for (int k = 0; k < n; ++k) {
            const TubeSegment& seg = segs_[static_cast<std::size_t>(k)];
            Vec3d prev{};
            double len = 0.0;
            bool ok = true;
            for (int i = 0; i <= 4 && ok; ++i) {
                const double c = slice(k, i / 4.0);
                Pt4 x;
                for (int d = 0; d < 4; ++d) x[d] = seg.a[d] + c * seg.m[d];
                ok = newton_slice(sys(), seg.m, seg.a, c, x);
                Vec3d C, Cu, Cv;
                oa_.eval(x[0], x[1], C, Cu, Cv);
                if (i > 0) len += norm(C - prev);
                prev = C;
            }
            if (ok && std::isfinite(len)) w[static_cast<std::size_t>(k)] = len;
            total += w[static_cast<std::size_t>(k)];
        }
        brk_.assign(static_cast<std::size_t>(n) + 1, 0.0);
        const bool even = !(total > 0.0);
        // A floor keeps every segment a positive share of sigma.
        double sum = 0.0;
        for (int k = 0; k < n; ++k) sum += even ? 1.0 : std::max(w[static_cast<std::size_t>(k)], 1e-6 * total);
        double acc = 0.0;
        for (int k = 0; k < n; ++k) {
            acc += even ? 1.0 : std::max(w[static_cast<std::size_t>(k)], 1e-6 * total);
            brk_[static_cast<std::size_t>(k) + 1] = acc / sum;
        }
        brk_.back() = 1.0;
    }

    /// Contacts and centre over a sigma range, lifted with their sigma-derivatives.
    struct SpineJet {
        carrier_detail::V2 a, c, b; ///< contacts and centre
        carrier_detail::V2 u, v;    ///< (a - c) / r and (b - c) / r
    };
    struct Cache {
        std::mutex mutex;
        std::map<std::tuple<int, double, double>, std::optional<SpineJet>> map;
    };
    /// Shared by copies (same geometry), so the carrier stays copyable.
    std::shared_ptr<Cache> cache_ = std::make_shared<Cache>();
    /// Names the geometry for the per-thread section memo (never reused).
    std::uint64_t serial_ = next_serial();
    static std::uint64_t next_serial() {
        static std::atomic<std::uint64_t> n{0};
        return ++n;
    }

    /// The theta-independent part of the enclosures for sigma in [s0, s1]
    /// within segment k. Subdivision revisits the same sigma ranges with many
    /// theta (and other-carrier) boxes, so results are cached.
    std::optional<SpineJet> spine_jet(int k, double s0, double s1) const {
        const auto key = std::make_tuple(k, s0, s1);
        {
            std::lock_guard<std::mutex> lock(cache_->mutex);
            if (auto it = cache_->map.find(key); it != cache_->map.end()) return it->second;
        }
        auto res = compute_spine_jet(k, s0, s1);
        std::lock_guard<std::mutex> lock(cache_->mutex);
        if (cache_->map.size() > 200000) cache_->map.clear();
        cache_->map.emplace(key, res);
        return res;
    }

    std::optional<SpineJet> compute_spine_jet(int k, double s0, double s1) const {
        using namespace carrier_detail;
        const TubeSegment& seg = segs_[static_cast<std::size_t>(k)];
        const double c0 = slice(k, local(k, s0)), c1 = slice(k, local(k, s1));
        const Interval Cs(std::min(c0, c1), std::max(c0, c1));
        const auto sy = sys();
        Pt4 xhat;
        const double cm = num::mid(Cs);
        for (int i = 0; i < 4; ++i) xhat[i] = seg.a[i] + cm * seg.m[i];
        if (!newton_slice(sy, seg.m, seg.a, cm, xhat)) return std::nullopt;
        // K(X) for the segment's certified box X holds the spine points of the
        // slices Cs; iterating K <- K(K) ∩ K keeps that and contracts towards
        // them (tight for a single slice, as point_at needs).
        Box4 K = seg.X;
        for (int it = 0; it < 8; ++it) {
            // The mean-value step behind K needs xhat in the box.
            bool in = true;
            for (int i = 0; i < 4; ++i) in = in && K[i].lo <= xhat[i] && xhat[i] <= K[i].hi;
            if (!in) break;
            Box4 Kn;
            if (!krawczyk_operator(sy, seg.m, seg.a, Cs, xhat, K, Kn)) {
                if (it == 0) return std::nullopt;
                break;
            }
            double before = 0.0, after = 0.0;
            for (int i = 0; i < 4; ++i) {
                const Interval m = meet(Kn[i], K[i]);
                before = std::max(before, K[i].hi - K[i].lo);
                after = std::max(after, m.hi - m.lo);
                K[i] = m;
                if (!num::is_finite(K[i])) return std::nullopt;
            }
            if (after > 0.5 * before || after < 1e-15) break; // stalled or tight
        }
        std::array<Interval, 4> xp;
        if (!spine_derivative(sy, seg.m, K, xp)) return std::nullopt;
        const Interval wk = Interval(brk_[static_cast<std::size_t>(k) + 1]) - Interval(brk_[static_cast<std::size_t>(k)]);
        const Interval dcds = Interval(seg.reversed ? -seg.len : seg.len) / wk;
        for (auto& d : xp) d = d * dcds;

        geom::IVec3 SA, SAu, SAv, OA, OAu, OAv, SB, SBs, SBt;
        geom::IVec3 OB, OBs, OBt;
        // Each support and its offset from one restriction of the patch.
        if (!oa_.jet1_over_with_patch(K[0].lo, K[0].hi, K[1].lo, K[1].hi, SA, SAu, SAv, OA, OAu, OAv)) return std::nullopt;
        if (!ob_.jet1_over_with_patch(K[2].lo, K[2].hi, K[3].lo, K[3].hi, SB, SBs, SBt, OB, OBs, OBt)) return std::nullopt;
        SpineJet j;
        j.a = lift(SA, SAu * xp[0] + SAv * xp[1]);
        // The centre is on both offsets: A's enclosure meets B's.
        j.c = lift(geom::carrier_detail::meet(OA, OB), OAu * xp[0] + OAv * xp[1]);
        j.b = lift(SB, SBs * xp[2] + SBt * xp[3]);
        const D2 R(r_);
        j.u = (j.a - j.c) / R;
        j.v = (j.b - j.c) / R;
        // Mean-value form in sigma around the range's middle: a quantity q
        // lies in q(mid) + q' [-h, h]. The direct values lose the correlation
        // between contacts and centre (dependency); this keeps it, so
        // u . v stays tight and the section arc stays defined on wide ranges.
        if (s1 > s0) {
            const double sm = 0.5 * s0 + 0.5 * s1;
            if (const auto m = spine_jet(k, sm, sm)) {
                const Interval h = Interval(s1) - Interval(sm);
                const Interval e = num::hull(h, -h);
                auto tighten = [&](V2& q, const V2& qm) {
                    auto one = [&](D2& x, const D2& xm) {
                        const Interval mv = xm.v + x.du * e;
                        if (num::is_finite(mv)) x.v = meet(x.v, mv);
                    };
                    one(q.x, qm.x); one(q.y, qm.y); one(q.z, qm.z);
                };
                tighten(j.a, m->a); tighten(j.c, m->c); tighten(j.b, m->b);
                tighten(j.u, m->u); tighten(j.v, m->v);
            }
        }
        return j;
    }

    /// Enclosures for sigma in [s0, s1] within segment k and theta in T.
    bool segment_jet(int k, double s0, double s1, const Interval& T, geom::IVec3& P,
                     geom::IVec3& Ps, geom::IVec3& Pt) const {
        using namespace carrier_detail;
        const auto j = spine_jet(k, s0, s1);
        if (!j) return false;
        // The arc's rational form loses the correlation between its basis
        // functions over a theta interval; pieces of theta, hulled, keep
        // the derivative enclosures within a few times their true spread.
        // (Measured on the rounded block: 40 times the box width before, 16
        // times with 8 pieces; pieces of sigma would cost a spine jet each.)
        const int pieces = T.hi > T.lo ? 8 : 1;
        bool first = true;
        for (int q = 0; q < pieces; ++q) {
            const double t0 = T.lo + (T.hi - T.lo) * q / pieces, t1 = q + 1 == pieces ? T.hi : T.lo + (T.hi - T.lo) * (q + 1) / pieces;
            bool ok = false;
            const V2 p = arc_point(j->a, j->c, j->b, j->u, j->v, r_, D2(Interval(t0, t1), Interval(0.0), Interval(1.0)), ok);
            if (!ok) return false;
            const geom::IVec3 p0{p.x.v, p.y.v, p.z.v}, p1{p.x.du, p.y.du, p.z.du}, p2{p.x.dv, p.y.dv, p.z.dv};
            auto h = [](geom::IVec3& x, const geom::IVec3& y) { x = {num::hull(x.x, y.x), num::hull(x.y, y.y), num::hull(x.z, y.z)}; };
            if (first) { P = p0; Ps = p1; Pt = p2; first = false; }
            else { h(P, p0); h(Ps, p1); h(Pt, p2); }
        }
        return num::is_finite(P.x) && num::is_finite(Ps.x) && num::is_finite(Pt.x) &&
               num::is_finite(P.y) && num::is_finite(Ps.y) && num::is_finite(Pt.y) &&
               num::is_finite(P.z) && num::is_finite(Ps.z) && num::is_finite(Pt.z);
    }

public:
    /// `spine` is one curve of find_spines(A, B, r) (or trace_spine's segments).
    /// `side` = -1 puts the ball on the side opposite the supports' normals
    /// (outward normals of a solid and a convex edge): the spine is then that
    /// of OffsetCarrier(A, -r) and OffsetCarrier(B, -r).
    BlendCarrier(const PolyPatchD& A, const PolyPatchD& B, double r, std::vector<TubeSegment> spine, int side = 1)
        : pa_(A), pb_(B), oa_(A, side * r), ob_(B, side * r), r_(r), segs_(std::move(spine)) {
        make_breaks();
    }

    std::size_t segments() const { return segs_.size(); }
    /// sigma is smooth on each spine segment, kinked between them.
    std::vector<double> breaks_u() const override { return brk_; }

    /// Spine parameters (u, v, s, t) at sigma, in doubles.
    bool spine_at(double sigma, Pt4& x) const {
        const int k = segment_of(sigma);
        const TubeSegment& seg = segs_[static_cast<std::size_t>(k)];
        const double c = slice(k, local(k, sigma));
        for (int i = 0; i < 4; ++i) x[i] = seg.a[i] + c * seg.m[i];
        return newton_slice(sys(), seg.m, seg.a, c, x);
    }

    /// The section arc at sigma, as the homogeneous control points of its
    /// rational quadratic (or the centre if no ball fits). Quadrature and
    /// sampling evaluate many theta at one sigma, so the last few sections
    /// are kept per thread.
    struct Section {
        bool ok = false;
        geom::Vec4d cp[3]{};
        Vec3d centre{};
    };
    Section section(double sg) const {
        struct Entry {
            std::uint64_t id = 0;
            double sigma = -1.0;
            Section sec;
        };
        thread_local std::array<Entry, 16> memo{};
        thread_local std::size_t next = 0;
        const std::uint64_t id = serial_; // copies share geometry, and so sections
        for (const Entry& e : memo)
            if (e.id == id && e.sigma == sg) return e.sec;
        Section S;
        Pt4 x;
        spine_at(std::clamp(sg, 0.0, 1.0), x);
        Vec3d a, au, av, b, bs, bt, c, cu, cv;
        pa_.eval(x[0], x[1], a, au, av);
        pb_.eval(x[2], x[3], b, bs, bt);
        oa_.eval(x[0], x[1], c, cu, cv);
        S.centre = c;
        const auto f = geom::fillet_arc(c, r_, (a - c) / r_, (b - c) / r_);
        if (f.status == geom::BlendStatus::Ok) {
            S.ok = true;
            for (int i = 0; i < 3; ++i) S.cp[i] = f.arc.cp[static_cast<std::size_t>(i)];
        }
        memo[next] = {id, sg, S};
        next = (next + 1) % memo.size();
        return S;
    }

    void eval(double sigma, double theta, Vec3d& P, Vec3d& Ps, Vec3d& Pt) const override {
        // Point from the section; derivatives by central differences (doubles
        // are only guesses for Newton steps and quadrature nodes; certified
        // values come from the enclosures).
        auto at = [&](double sg, double th) {
            const Section S = section(sg);
            if (!S.ok) return S.centre;
            const double u = 1.0 - th;
            const geom::Vec4d h = S.cp[0] * (u * u) + S.cp[1] * (2.0 * u * th) + S.cp[2] * (th * th);
            return Vec3d{h.x / h.w, h.y / h.w, h.z / h.w};
        };
        P = at(sigma, theta);
        // Differences in sigma stay inside sigma's spine segment (the speed
        // may jump at its ends): one-sided, second order, near them.
        const int k = segment_of(sigma);
        const double lo = brk_[static_cast<std::size_t>(k)], hi = brk_[static_cast<std::size_t>(k) + 1];
        const double h = std::min(1e-6, 0.25 * (hi - lo));
        if (sigma - h < lo)
            Ps = (at(sigma + h, theta) * 4.0 - at(sigma + 2 * h, theta) - P * 3.0) / (2 * h);
        else if (sigma + h > hi)
            Ps = (P * 3.0 - at(sigma - h, theta) * 4.0 + at(sigma - 2 * h, theta)) / (2 * h);
        else
            Ps = (at(sigma + h, theta) - at(sigma - h, theta)) / (2 * h);
        Pt = (at(sigma, theta + 1e-6) - at(sigma, theta - 1e-6)) / 2e-6;
    }

    geom::IVec3 point_at(double sigma, double theta) const override {
        geom::IVec3 P, Ps, Pt;
        const int k = segment_of(sigma);
        if (!segment_jet(k, sigma, sigma, Interval(theta), P, Ps, Pt))
            return {Interval::entire(), Interval::entire(), Interval::entire()};
        return P;
    }

    bool jet1_over(double s0, double s1, double t0, double t1, geom::IVec3& P, geom::IVec3& Ps,
                   geom::IVec3& Pt) const override {
        const int k0 = segment_of(s0), k1 = segment_of(s1);
        bool first = true;
        for (int k = k0; k <= k1; ++k) {
            const double a = std::max(s0, brk_[static_cast<std::size_t>(k)]);
            const double b = std::min(s1, brk_[static_cast<std::size_t>(k) + 1]);
            geom::IVec3 p, ps, pt;
            if (!segment_jet(k, a, b, Interval(t0, t1), p, ps, pt)) return false;
            if (first) { P = p; Ps = ps; Pt = pt; first = false; }
            else {
                auto h = [](geom::IVec3& x, const geom::IVec3& y) {
                    x = {num::hull(x.x, y.x), num::hull(x.y, y.y), num::hull(x.z, y.z)};
                };
                h(P, p); h(Ps, ps); h(Pt, pt);
            }
        }
        if (first) return false;
        // Mean-value form around the box centre, intersected with the direct
        // enclosure: P(s, t) in P(sc, tc) + Ps [s - sc] + Pt [t - tc] since the
        // segment from the centre stays in the box. It converges where the
        // direct form does not tighten.
        const double sc = 0.5 * s0 + 0.5 * s1, tc = 0.5 * t0 + 0.5 * t1;
        geom::IVec3 C, Cs_, Ct_;
        if (segment_jet(segment_of(sc), sc, sc, Interval(tc), C, Cs_, Ct_)) {
            const Interval ds = Interval(s0) - Interval(sc), dt = Interval(t0) - Interval(tc);
            const Interval es = num::hull(ds, -ds), et = num::hull(dt, -dt);
            auto mv = [&](const Interval& c, const Interval& gs, const Interval& gt) {
                return c + gs * es + gt * et;
            };
            const geom::IVec3 M{mv(C.x, Ps.x, Pt.x), mv(C.y, Ps.y, Pt.y), mv(C.z, Ps.z, Pt.z)};
            P = geom::carrier_detail::meet(P, M);
        }
        return true;
    }
};

} // namespace ak::solve
