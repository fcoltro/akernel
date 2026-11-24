// akernel — procedural carriers: surfaces kept as their defining equation
// and evaluated with certified enclosures (docs/FEATURES.md §5, decided
// 2026-11-23: procedural carriers are the definition of offsets, blends and
// sweeps that are not rational).
//
// A carrier gives, for parameters (u, v):
//   - eval:      the point and its first partials in doubles (Newton, seeds),
//   - point_at:  a rigorous enclosure of the point at exact (u, v),
//   - jet1_over: rigorous enclosures of the point and its first partials over
//                a parameter box (Krawczyk Jacobians, broad phase).
// That is what the certified intersection machinery needs
// (solve/carrier_ssi.hpp), so any carrier can be intersected with any
// other, with the same certificates as two rational patches.
//
// Two carriers here:
//   PatchCarrier   a rational Bezier patch (the kernel's existing carrier),
//   OffsetCarrier  S + d n for a rational Bezier patch S, n = S_u x S_v / |..|.
// The offset of a general patch is not rational (BLUEPRINT §1.6), so it has
// no Bezier form; its enclosures come from the patch's homogeneous partials
// up to order 3 through the normal's formula, n_u = (N_u - n (n . N_u)) / |N|
// with N_u = S_uu x S_v + S_u x S_uv. As in geom/offset_validity.hpp each
// box enclosure is the intersection of the direct interval evaluation and
// the mean-value form around the box centre, which converges as the box
// shrinks even where the direct form is useless.
#pragma once

#include <array>
#include <memory>
#include <vector>

#include <akernel/geom/bernstein.hpp>
#include <akernel/geom/offset_validity.hpp>
#include <akernel/num/interval.hpp>

namespace ak::geom {

using IVec3 = Vec3<num::Interval>;

class Carrier {
public:
    virtual ~Carrier() = default;
    virtual void eval(double u, double v, Vec3d& P, Vec3d& Pu, Vec3d& Pv) const = 0;
    virtual IVec3 point_at(double u, double v) const = 0;
    /// Enclosures over [u0,u1] x [v0,v1]; `false` if not proven finite
    /// (a weight or |S_u x S_v| not proven positive on the box).
    virtual bool jet1_over(double u0, double u1, double v0, double v1, IVec3& P, IVec3& Pu,
                           IVec3& Pv) const = 0;
    /// Break points 0 = b0 < ... < bn = 1 in u: the carrier is smooth on each
    /// piece, its u-derivative may jump between them (e.g. at the segments
    /// of a traced spine). Quadrature and sampling align with them.
    virtual std::vector<double> breaks_u() const { return {0.0, 1.0}; }
};

namespace carrier_detail {

using num::Interval;
using offset_detail::D2;

template <class T>
struct Jet1 {
    Vec3<T> p, pu, pv;
    bool ok = false;
};

/// Point and first partials of the projected patch, or of its offset by d,
/// from the homogeneous partials X[jet_index(k, l)], k + l <= 2.
template <class T>
Jet1<T> first_jet(const std::array<Vec4<T>, 6>& X, double d, bool offset) {
    using offset_detail::sqrt_;
    using offset_detail::value;
    using V = Vec3<T>;
    Jet1<T> r;
    if (!(value(X[0].w).lo > 0.0)) return r;
    const T w = X[0].w, two(2.0);
    const V s = xyz(X[0]) / w;
    const V su = (xyz(X[1]) - s * X[1].w) / w;
    const V sv = (xyz(X[2]) - s * X[2].w) / w;
    if (!offset) {
        r.p = s;
        r.pu = su;
        r.pv = sv;
        r.ok = true;
        return r;
    }
    const V suu = (xyz(X[3]) - su * (two * X[1].w) - s * X[3].w) / w;
    const V suv = (xyz(X[4]) - su * X[2].w - sv * X[1].w - s * X[4].w) / w;
    const V svv = (xyz(X[5]) - sv * (two * X[2].w) - s * X[5].w) / w;
    const V N = cross(su, sv);
    const T G = dot(N, N);
    if (!(value(G).lo > 0.0)) return r;
    const T L = sqrt_(G);
    const V n = N / L;
    const V Nu = cross(suu, sv) + cross(su, suv);
    const V Nv = cross(suv, sv) + cross(su, svv);
    const V nu = (Nu - n * dot(n, Nu)) / L;
    const V nv = (Nv - n * dot(n, Nv)) / L;
    const T D(d);
    r.p = s + n * D;
    r.pu = su + nu * D;
    r.pv = sv + nv * D;
    r.ok = true;
    return r;
}

inline Interval meet(const Interval& a, const Interval& b) {
    const Interval m{std::max(a.lo, b.lo), std::min(a.hi, b.hi)};
    return m.lo <= m.hi ? m : num::hull(a, b); // empty only through a bug; stay safe
}
inline IVec3 meet(const IVec3& a, const IVec3& b) { return {meet(a.x, b.x), meet(a.y, b.y), meet(a.z, b.z)}; }

/// c + g_u e_u + g_v e_v componentwise, for the mean-value form.
inline IVec3 mean_value(const IVec3& c, const Vec3<D2>& g, const Interval& eu, const Interval& ev) {
    auto one = [&](const Interval& cc, const D2& gg) {
        const Interval m = cc + gg.du * eu + gg.dv * ev;
        return num::is_finite(m) ? m : Interval::entire();
    };
    return {one(c.x, g.x), one(c.y, g.y), one(c.z, g.z)};
}
inline IVec3 values(const Vec3<D2>& a) { return {a.x.v, a.y.v, a.z.v}; }

/// Shared implementation for patch and offset carriers.
class PatchBased : public Carrier {
protected:
    offset_detail::Partials P_;  ///< interval homogeneous partials to order 3
    std::array<PolyPatchD, 6> D_; ///< double homogeneous partials to order 2
    double d_ = 0.0;
    bool offset_ = false;

    PatchBased(const PolyPatchD& patch, double d, bool offset)
        : P_(to_interval(patch)), d_(d), offset_(offset) {
        D_[0] = patch;
        D_[1] = patch.deriv_u();
        D_[2] = patch.deriv_v();
        D_[3] = D_[1].deriv_u();
        D_[4] = D_[1].deriv_v();
        D_[5] = D_[2].deriv_v();
    }

public:
    void eval(double u, double v, Vec3d& P, Vec3d& Pu, Vec3d& Pv) const override {
        std::array<Vec4d, 6> X;
        for (int i = 0; i < 6; ++i) X[i] = D_[i].eval(u, v);
        // Same formula in doubles (Interval-free instantiation through point intervals
        // would be slower; doubles are only guesses here).
        const Vec3d s = xyz(X[0]) / X[0].w;
        const Vec3d su = (xyz(X[1]) - s * X[1].w) / X[0].w;
        const Vec3d sv = (xyz(X[2]) - s * X[2].w) / X[0].w;
        if (!offset_) { P = s; Pu = su; Pv = sv; return; }
        const Vec3d suu = (xyz(X[3]) - su * (2.0 * X[1].w) - s * X[3].w) / X[0].w;
        const Vec3d suv = (xyz(X[4]) - su * X[2].w - sv * X[1].w - s * X[4].w) / X[0].w;
        const Vec3d svv = (xyz(X[5]) - sv * (2.0 * X[2].w) - s * X[5].w) / X[0].w;
        const Vec3d N = cross(su, sv);
        const double L = norm(N);
        const Vec3d n = N / L;
        const Vec3d Nu = cross(suu, sv) + cross(su, suv), Nv = cross(suv, sv) + cross(su, svv);
        P = s + n * d_;
        Pu = su + (Nu - n * dot(n, Nu)) * (d_ / L);
        Pv = sv + (Nv - n * dot(n, Nv)) * (d_ / L);
    }

    IVec3 point_at(double u, double v) const override {
        std::array<Vec4<Interval>, 6> X;
        for (int i = 0; i < (offset_ ? 6 : 3); ++i) X[i] = P_.d[i].eval(Interval(u), Interval(v));
        const auto j = first_jet(X, d_, offset_);
        if (!j.ok) return {Interval::entire(), Interval::entire(), Interval::entire()};
        return j.p;
    }

    bool jet1_over(double u0, double u1, double v0, double v1, IVec3& P, IVec3& Pu,
                   IVec3& Pv) const override {
        return jets(u0, u1, v0, v1, P, Pu, Pv, nullptr);
    }
    /// For an offset carrier: the enclosures of the offset and of the patch
    /// itself over one box, from one restriction of the patch's partials
    /// (a blend needs both for each support).
    bool jet1_over_with_patch(double u0, double u1, double v0, double v1, IVec3& S, IVec3& Su, IVec3& Sv,
                              IVec3& P, IVec3& Pu, IVec3& Pv) const {
        IVec3 patch[3];
        if (!jets(u0, u1, v0, v1, P, Pu, Pv, offset_ ? patch : nullptr)) return false;
        if (!offset_) { S = P; Su = Pu; Sv = Pv; return true; }
        S = patch[0];
        Su = patch[1];
        Sv = patch[2];
        return true;
    }

private:
    bool jets(double u0, double u1, double v0, double v1, IVec3& P, IVec3& Pu, IVec3& Pv, IVec3* patch) const {
        // The offset needs the patch's partials to order 2 with their own
        // derivatives (order 3); the patch itself only to order 1 (order 2).
        const int order = offset_ ? 2 : 1;
        const int count = (order + 2) * (order + 3) / 2;
        std::array<Vec4<Interval>, kJetSize> H;
        for (int i = 0; i < count; ++i) H[i] = hull(P_.d[i].restricted(u0, u1, v0, v1));
        std::array<Vec4<D2>, 6> Xd;
        for (int k = 0; k <= order; ++k)
            for (int l = 0; l <= order - k; ++l) {
                const int i = jet_index(k, l), iu = jet_index(k + 1, l), iv = jet_index(k, l + 1);
                Xd[i] = {D2(H[i].x, H[iu].x, H[iv].x), D2(H[i].y, H[iu].y, H[iv].y),
                         D2(H[i].z, H[iu].z, H[iv].z), D2(H[i].w, H[iu].w, H[iv].w)};
            }
        const double uc = 0.5 * u0 + 0.5 * u1, vc = 0.5 * v0 + 0.5 * v1;
        std::array<Vec4<Interval>, 6> Xc;
        for (int i = 0; i < (offset_ ? 6 : 3); ++i) Xc[i] = P_.d[i].eval(Interval(uc), Interval(vc));
        const Interval du = Interval(u0) - Interval(uc), dv = Interval(v0) - Interval(vc);
        const Interval eu = num::hull(du, -du), ev = num::hull(dv, -dv);
        auto finish = [&](double d, bool off, IVec3& Q, IVec3& Qu, IVec3& Qv) {
            const auto g = first_jet(Xd, d, off);
            const auto c = first_jet(Xc, d, off);
            if (!g.ok || !c.ok) return false;
            Q = meet(values(g.p), mean_value(c.p, g.p, eu, ev));
            Qu = meet(values(g.pu), mean_value(c.pu, g.pu, eu, ev));
            Qv = meet(values(g.pv), mean_value(c.pv, g.pv, eu, ev));
            return true;
        };
        if (!finish(d_, offset_, P, Pu, Pv)) return false;
        if (patch && !finish(0.0, false, patch[0], patch[1], patch[2])) return false;
        return true;
    }
};

} // namespace carrier_detail

/// A rational Bezier patch seen through the carrier interface.
class PatchCarrier final : public carrier_detail::PatchBased {
public:
    explicit PatchCarrier(const PolyPatchD& patch) : PatchBased(patch, 0.0, false) {}
};

/// The offset S + d n of a rational Bezier patch, n = (S_u x S_v) / |S_u x S_v|.
/// Where the offset is not regular (1 - d kappa_i <= 0, see
/// geom/offset_validity.hpp) the map is still defined but folds; callers
/// check validity first. Where S is singular the enclosures report failure.
class OffsetCarrier final : public carrier_detail::PatchBased {
public:
    OffsetCarrier(const PolyPatchD& patch, double d) : PatchBased(patch, d, true) {}
    double distance() const { return d_; }
};

} // namespace ak::geom
