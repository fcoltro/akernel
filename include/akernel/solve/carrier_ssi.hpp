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
// What is not here yet: the completeness proof of ssi.hpp's intersect()
// (which subdivides Bezier patches) for carriers, so the spine is traced
// from given seeds and is complete only along the traced branch.
#pragma once

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
