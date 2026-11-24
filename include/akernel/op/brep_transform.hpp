// akernel — placing solids: exact rigid maps of a Brep, including its
// procedural faces and edges, and re-indexing of procedural edges when two
// Breps are merged (docs/FEATURES.md §5.1).
//
// A SignedPermutation map x'_i = s_i x_{p(i)} + t_i (quarter turns, mirror
// images and translations) is exact on doubles for dyadic t of moderate
// size: rational carriers keep exact control points, and procedural ones
// are wrapped so that every evaluation, enclosure and edge point is mapped
// the same way (components permuted and negated, which is exact in
// intervals; the translation is added with outward rounding).
#pragma once

#include <array>
#include <memory>
#include <vector>

#include <akernel/geom/carrier.hpp>
#include <akernel/op/brep.hpp>

namespace ak::op {

/// x'_i = sign[i] * x_{perm[i]} + shift_i.
struct SignedPermutation {
    std::array<int, 3> perm{0, 1, 2};
    std::array<double, 3> sign{1.0, 1.0, 1.0};
    Vec3d shift{0, 0, 0};

    /// Identity.
    static SignedPermutation translation(const Vec3d& t) {
        SignedPermutation m;
        m.shift = t;
        return m;
    }
    /// A quarter turn about the z axis, (x, y, z) -> (-y, x, z), then t.
    static SignedPermutation quarter_turn_z(const Vec3d& t = {0, 0, 0}) {
        SignedPermutation m;
        m.perm = {1, 0, 2};
        m.sign = {-1.0, 1.0, 1.0};
        m.shift = t;
        return m;
    }
    double det() const {
        // Sign of the permutation times the product of the signs.
        int inv = 0;
        for (int i = 0; i < 3; ++i)
            for (int j = i + 1; j < 3; ++j) inv += perm[static_cast<std::size_t>(i)] > perm[static_cast<std::size_t>(j)];
        return (inv % 2 ? -1.0 : 1.0) * sign[0] * sign[1] * sign[2];
    }
    template <class V>
    V linear(const V& x) const {
        // Negation is exact for doubles and intervals alike.
        const auto c = [&](int i) {
            const int k = perm[static_cast<std::size_t>(i)];
            const auto& a = k == 0 ? x.x : k == 1 ? x.y : x.z;
            return sign[static_cast<std::size_t>(i)] < 0.0 ? -a : a;
        };
        return {c(0), c(1), c(2)};
    }
    Vec3d apply(const Vec3d& x) const { return linear(x) + shift; }
    geom::IVec3 apply(const geom::IVec3& x) const {
        const geom::IVec3 l = linear(x);
        return {l.x + num::Interval(shift.x), l.y + num::Interval(shift.y), l.z + num::Interval(shift.z)};
    }
    geom::Vec4d apply_h(const geom::Vec4d& h) const {
        const Vec3d l = linear(Vec3d{h.x, h.y, h.z});
        return {l.x + shift.x * h.w, l.y + shift.y * h.w, l.z + shift.z * h.w, h.w};
    }
};

namespace transform_detail {

/// A carrier seen through a map.
class MappedCarrier final : public geom::Carrier {
    std::shared_ptr<const geom::Carrier> inner_;
    SignedPermutation m_;

public:
    MappedCarrier(std::shared_ptr<const geom::Carrier> inner, const SignedPermutation& m) : inner_(std::move(inner)), m_(m) {}
    void eval(double u, double v, Vec3d& P, Vec3d& Pu, Vec3d& Pv) const override {
        inner_->eval(u, v, P, Pu, Pv);
        P = m_.apply(P);
        Pu = m_.linear(Pu);
        Pv = m_.linear(Pv);
    }
    geom::IVec3 point_at(double u, double v) const override { return m_.apply(inner_->point_at(u, v)); }
    bool jet1_over(double u0, double u1, double v0, double v1, geom::IVec3& P, geom::IVec3& Pu,
                   geom::IVec3& Pv) const override {
        if (!inner_->jet1_over(u0, u1, v0, v1, P, Pu, Pv)) return false;
        P = m_.apply(P);
        Pu = m_.linear(Pu);
        Pv = m_.linear(Pv);
        return true;
    }
    std::vector<double> breaks_u() const override { return inner_->breaks_u(); }
};

/// A procedural edge through a map, with its carrier and proc indices
/// shifted (the Brep it lives in was appended to another).
class MappedEdge final : public ProcEdge {
    std::shared_ptr<const ProcEdge> inner_;
    SignedPermutation m_;
    int dc_, dp_;

public:
    MappedEdge(std::shared_ptr<const ProcEdge> inner, const SignedPermutation& m, int carrier_shift, int proc_shift)
        : inner_(std::move(inner)), m_(m), dc_(carrier_shift), dp_(proc_shift) {}
    Vec3d point(double tau) const override { return m_.apply(inner_->point(tau)); }
    bool on_carrier(int c, double tau, double& u, double& v, double& du, double& dv) const override {
        return c - dc_ >= 0 && inner_->on_carrier(c - dc_, tau, u, v, du, dv);
    }
    bool on_proc(int p, double tau, double& s, double& t) const override {
        return p - dp_ >= 0 && inner_->on_proc(p - dp_, tau, s, t);
    }
    std::vector<double> breaks() const override { return inner_->breaks(); }
    std::shared_ptr<const ProcEdge> part(double t0, double t1) const override {
        auto p = inner_->part(t0, t1);
        if (!p) return nullptr;
        return std::make_shared<MappedEdge>(std::move(p), m_, dc_, dp_);
    }
};

} // namespace transform_detail

/// The procedural edge e read in a Brep whose carriers and procs were
/// appended after carrier_shift carriers and proc_shift procs.
inline std::shared_ptr<const ProcEdge> reindexed(std::shared_ptr<const ProcEdge> e, int carrier_shift, int proc_shift) {
    if (!e || (carrier_shift == 0 && proc_shift == 0)) return e;
    return std::make_shared<transform_detail::MappedEdge>(std::move(e), SignedPermutation{}, carrier_shift, proc_shift);
}

/// B moved by m. A mirror map (det < 0) turns every face over, so that
/// outward stays outward.
inline Brep transformed(const Brep& B, const SignedPermutation& m) {
    Brep R = B;
    for (PolyPatchD& C : R.carriers)
        for (geom::Vec4d& c : C.c) c = m.apply_h(c);
    for (auto& p : R.procs) p = std::make_shared<transform_detail::MappedCarrier>(p, m);
    for (BVertex& v : R.vertices) {
        v.point = m.apply(v.point);
        if (v.exact) v.xp = m.apply_h(v.xp);
    }
    for (BEdge& e : R.edges)
        if (e.proc) e.proc = std::make_shared<transform_detail::MappedEdge>(e.proc, m, 0, 0);
    if (m.det() < 0.0)
        for (BFace& f : R.faces) f.flipped = !f.flipped;
    R.outward_known = false;
    return R;
}

} // namespace ak::op
