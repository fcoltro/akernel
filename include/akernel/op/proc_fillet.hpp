// akernel — a freeform fillet inside a solid: the blend is a procedural face
// (solve/carrier_ssi.hpp, BlendCarrier) and the supports are trimmed along
// the certified contact curves (docs/FEATURES.md §5.1, first target).
//
// round_block_edge() takes the block x in [0, a], y in [0, L], 0 <= z <=
// top(x, y), whose top is a polynomial graph patch over the rectangle, and
// rounds the convex edge between the top and the side x = a with radius r.
// The top's dihedral angle with the side varies along the edge, so the
// rolling-ball blend is not a rational surface: it is kept as its defining
// equation (the certified spine plus the section arc).
//
// Faces (each outward with Su x Sv, so none is flipped except possibly the
// blend):
//   bottom   (L u, a v) at z = 0          side x = 0   (0, L v, u top(0, v))
//   side B   (a, L u, v top(1, u))        top A        (a u, L v, top(u, v))
//   end y=0  (a u, 0, v top(u, 0))        end y=L      (a v, L, u top(v, 1))
//   blend    sigma along the spine, theta from the contact on the top to the
//            contact on side B.
//
// The ends. Where the top has zero slope across an end plane (its first two
// control rows in y equal, or its last two), the top's normal there has no
// y component, the ball centre and both contacts lie in the end plane, and
// the blend's end section is a planar arc in the end face: no trim is
// needed. At any other end (a sloped top, or a cut) the supports are
// extended a quarter of their length past the end (the polynomial
// continued), the spine is traced over the extended supports, and the blend
// is trimmed by its certified intersection curve with the end face
// (op/proc_edge.hpp); the procedural face is then less than its square.
//
// Edges are side edges (parts of the patches' sides, as in op/brep.hpp),
// or procedural: the two contact curves and the two end sections, each a
// side of the blend's square, read on the supports through the spine's own
// coordinates and on the end faces by inverting the planar patch.
#pragma once

#include <algorithm>
#include <cmath>
#include <functional>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include <akernel/num/detmath.hpp>
#include <akernel/op/brep.hpp>
#include <akernel/op/proc_edge.hpp>
#include <akernel/solve/blend_corner.hpp>
#include <akernel/solve/carrier_ssi.hpp>

namespace ak::op {

namespace proc_fillet_detail {

using geom::Vec4d;

/// Parameters of a point on a patch by Newton from a guess (doubles: these
/// feed quadrature and sampling, not certificates).
inline bool invert(const PolyPatchD& P, const Vec3d& q, double& u, double& v) {
    const geom::BezierSurface S(P);
    for (int it = 0; it < 50; ++it) {
        const geom::Jet3 J = S.jet(u, v);
        const Vec3d d = q - J[0];
        const double a11 = dot(J[1], J[1]), a12 = dot(J[1], J[2]), a22 = dot(J[2], J[2]);
        const double b1 = dot(J[1], d), b2 = dot(J[2], d), det = a11 * a22 - a12 * a12;
        if (!(std::fabs(det) > 0.0)) return false;
        const double du = (a22 * b1 - a12 * b2) / det, dv = (a11 * b2 - a12 * b1) / det;
        u += du;
        v += dv;
        if (std::fabs(du) + std::fabs(dv) < 1e-15) return true;
    }
    return true;
}

/// A side of the blend's square, or part of one: a contact curve (theta = 0
/// or 1, sigma from `from` to `to`) or an end section (sigma fixed, theta
/// from `from` to `to`); tau runs from the edge's v0 to v1. On a support the
/// contact curve is read through the spine's coordinates, divided by
/// `scale` (a support restricted to [0, scale] in that coordinate).
class BlendSide final : public ProcEdge {
public:
    std::shared_ptr<const solve::BlendCarrier> blend;
    int proc = -1;
    bool contact = true;   ///< contact curve (theta fixed) or end section (sigma fixed)
    double fixed = 0.0;    ///< the fixed theta or sigma
    double from = 0.0, to = 1.0;
    int support = -1;      ///< contact: the carrier the spine coordinates refer to
    int coord = 0;         ///< 0: (u, v) on A, 2: (s, t) on B
    double scale[2] = {1.0, 1.0};
    double shift[2] = {0.0, 0.0}; ///< face coordinate = (spine coordinate - shift) / scale
    int plane = -1;        ///< end section: the end face's carrier
    PolyPatchD plane_patch;

    double param(double tau) const { return from + tau * (to - from); }
    void st(double tau, double& s, double& t) const {
        if (contact) { s = param(tau); t = fixed; }
        else { s = fixed; t = param(tau); }
    }
    Vec3d point(double tau) const override {
        double s, t;
        st(tau, s, t);
        Vec3d P, Ps, Pt;
        blend->eval(s, t, P, Ps, Pt);
        return P;
    }
    bool on_proc(int p, double tau, double& s, double& t) const override {
        if (p != proc) return false;
        st(tau, s, t);
        return true;
    }
    bool uv(int carrier, double tau, double& u, double& v) const {
        if (contact) {
            if (carrier != support) return false;
            solve::Pt4 x;
            if (!blend->spine_at(param(tau), x)) return false;
            u = (x[static_cast<std::size_t>(coord)] - shift[0]) / scale[0];
            v = (x[static_cast<std::size_t>(coord + 1)] - shift[1]) / scale[1];
            return true;
        }
        if (carrier != plane) return false;
        u = v = 0.5;
        return invert(plane_patch, point(tau), u, v);
    }
    bool on_carrier(int carrier, double tau, double& u, double& v, double& du, double& dv) const override {
        if (!uv(carrier, tau, u, v)) return false;
        return piece_derivative(breaks(), tau, [&](double t, double& p, double& q) { return uv(carrier, t, p, q); }, du, dv);
    }
    std::shared_ptr<const ProcEdge> part(double t0, double t1) const override {
        auto e = std::make_shared<BlendSide>(*this);
        e->from = param(t0);
        e->to = param(t1);
        return e;
    }
    std::vector<double> breaks() const override {
        if (!contact) return ProcEdge::breaks();
        const double lo = std::min(from, to), hi = std::max(from, to);
        std::vector<double> t{0.0};
        for (double s : blend->breaks_u())
            if (s > lo && s < hi) t.push_back((s - from) / (to - from));
        t.push_back(1.0);
        std::sort(t.begin(), t.end());
        return t;
    }
};

/// A point's (sigma, theta) on a blend, from a guess of sigma: its footpoint
/// on the spine (Newton on (P - c) . c', c' by differences), then theta on
/// that section's arc (from the angle, refined by Newton). Doubles.
/// The footpoint of P on a blend's spine (centre curve), by Newton on
/// (P - c) . c' with c' from the spine's own derivative enclosure (its
/// middle): differences would leave an error of their step in sigma.
inline double spine_foot(const solve::BlendCarrier& b, const Vec3d& P, double sg) {
    for (int it = 0; it < 80; ++it) {
        solve::BlendCarrier::SpineBox j;
        if (!b.spine_box(sg, sg, j)) break;
        const Vec3d c{num::mid(j.c.x), num::mid(j.c.y), num::mid(j.c.z)};
        const Vec3d cp{num::mid(j.c1.x), num::mid(j.c1.y), num::mid(j.c1.z)};
        const Vec3d cpp{num::mid(j.c2.x), num::mid(j.c2.y), num::mid(j.c2.z)};
        const double f = dot(P - c, cp), g = dot(P - c, cpp) - dot(cp, cp);
        if (!(std::fabs(g) > 0.0)) break;
        const double step = -f / g;
        sg = std::clamp(sg + step, 0.0, 1.0);
        if (std::fabs(step) < 1e-17) break;
    }
    return sg;
}

inline bool blend_coords(const solve::BlendCarrier& b, const Vec3d& P, double sg, double& s, double& t) {
    sg = spine_foot(b, P, sg);
    const auto S = b.section(sg);
    if (!S.ok) return false;
    const Vec3d a{S.cp[0].x / S.cp[0].w, S.cp[0].y / S.cp[0].w, S.cp[0].z / S.cp[0].w};
    const Vec3d e{S.cp[2].x / S.cp[2].w, S.cp[2].y / S.cp[2].w, S.cp[2].z / S.cp[2].w};
    const Vec3d ua = a - S.centre, ue = e - S.centre, up = P - S.centre;
    auto ang = [](const Vec3d& x, const Vec3d& y) { return num::det::atan2(norm(cross(x, y)), dot(x, y)); };
    const double phi = ang(ua, ue), psi = ang(ua, up);
    // theta of the angle psi on a rational quadratic arc of angle phi.
    const double tp = num::det::tan(0.5 * psi), tq = num::det::tan(0.5 * (phi - psi));
    double th = tp / (tp + tq);
    for (int it = 0; it < 30; ++it) {
        Vec3d Q, Qs, Qt;
        b.eval(sg, th, Q, Qs, Qt);
        const double f = dot(Q - P, Qt), g = dot(Qt, Qt);
        if (!(g > 0.0)) break;
        const double step = f / g;
        th -= step;
        if (std::fabs(step) < 1e-16) break;
    }
    s = sg;
    t = th;
    return std::isfinite(s) && std::isfinite(t);
}

/// Whether P lies in the wedge of the blend's section at its spine foot:
/// between the rays from the centre through the arc's two ends.
inline bool in_wedge(const solve::BlendCarrier& b, const Vec3d& P, double sg) {
    sg = spine_foot(b, P, sg);
    const auto S = b.section(sg);
    if (!S.ok) return false;
    const Vec3d a{S.cp[0].x / S.cp[0].w, S.cp[0].y / S.cp[0].w, S.cp[0].z / S.cp[0].w};
    const Vec3d e{S.cp[2].x / S.cp[2].w, S.cp[2].y / S.cp[2].w, S.cp[2].z / S.cp[2].w};
    const Vec3d ua = a - S.centre, ue = e - S.centre, up = P - S.centre;
    auto ang = [](const Vec3d& x, const Vec3d& y) { return num::det::atan2(norm(cross(x, y)), dot(x, y)); };
    return ang(ua, up) + ang(up, ue) < ang(ua, ue) + 1e-9;
}

/// The piece of a mitre curve next to the touching point A of two blends
/// (solve/blend_corner.hpp): the proven branch on the first blend's arc,
/// in the second blend's scaled coordinates (a, b) for b from 0 to b_end,
/// a found by Newton on h (distance to the first spine minus r1).
class CornerEdge final : public ProcEdge {
public:
    std::shared_ptr<const solve::BlendCarrier> b1, b2;
    solve::ContactOnFace f1, f2;
    solve::CornerLocal L;
    int p1 = -1, p2 = -1;
    double a_end = 0.0, b_end = 0.0;

    void ab(double tau, double& a, double& bb) const {
        bb = tau * b_end;
        a = tau * a_end;
        if (tau == 0.0) return;
        const double s2 = num::mid(L.A.s2);
        auto h = [&](double x) {
            Vec3d P, Ps, Pt;
            b2->eval(s2 + L.k_s * x, f2.theta_c + L.k_t * bb, P, Ps, Pt);
            const double sg = spine_foot(*b1, P, num::mid(L.A.s1));
            return norm(P - b1->section(sg).centre) - b1->radius();
        };
        for (int it = 0; it < 60; ++it) {
            const double e = 1e-3 * std::max(std::fabs(bb), 1e-12), f = h(a), g = (h(a + e) - h(a - e)) / (2 * e);
            if (!(std::fabs(g) > 0.0)) break;
            const double step = f / g;
            a -= step;
            if (std::fabs(step) < 1e-17) break;
        }
    }
    void st2(double tau, double& s, double& t) const {
        double a, bb;
        ab(tau, a, bb);
        s = num::mid(L.A.s2) + L.k_s * a;
        t = f2.theta_c + L.k_t * bb;
    }
    Vec3d point(double tau) const override {
        double s, t;
        st2(tau, s, t);
        Vec3d P, Ps, Pt;
        b2->eval(s, t, P, Ps, Pt);
        return P;
    }
    bool on_carrier(int, double, double&, double&, double&, double&) const override { return false; }
    bool on_proc(int p, double tau, double& s, double& t) const override {
        if (p == p2) {
            st2(tau, s, t);
            return true;
        }
        if (p == p1) return blend_coords(*b1, point(tau), num::mid(L.A.s1), s, t);
        return false;
    }
    std::vector<double> breaks() const override { return {0.0, 0.25, 0.5, 0.75, 1.0}; }
    /// Derivatives from a five-point Lagrange stencil of step 1e-3 inside the
    /// piece (error about 1e-12): the coordinates come from Newton solves,
    /// whose last-digit noise a two-point difference of step 1e-6 turns into
    /// 1e-10, and the volume weighs ds by the face's flux function.
    bool on_proc_jet(int p, double tau, double& s, double& t, double& ds, double& dt) const override {
        if (!on_proc(p, tau, s, t)) return false;
        const std::vector<double> br = breaks();
        const auto it = std::upper_bound(br.begin() + 1, br.end() - 1, tau);
        const double lo = *(it - 1), hi = *it, h = 1e-3;
        const double c = std::min(std::max(tau, lo + 2 * h), hi - 2 * h);
        double x[5], fs[5], ft[5];
        for (int k = 0; k < 5; ++k) {
            x[k] = c + (k - 2) * h;
            if (!on_proc(p, x[k], fs[k], ft[k])) return false;
        }
        ds = dt = 0.0;
        for (int k = 0; k < 5; ++k) {
            double den = 1.0, num_ = 0.0;
            for (int n = 0; n < 5; ++n)
                if (n != k) den *= x[k] - x[n];
            for (int m = 0; m < 5; ++m) {
                if (m == k) continue;
                double pr = 1.0;
                for (int n = 0; n < 5; ++n)
                    if (n != k && n != m) pr *= tau - x[n];
                num_ += pr;
            }
            ds += fs[k] * num_ / den;
            dt += ft[k] * num_ / den;
        }
        if (p == p2) dt = L.k_t * b_end; // theta is linear in tau
        return true;
    }
};

/// Remove the edges no loop uses and the vertices no edge uses, renumbering
/// the rest (a corner whose fillet swallows an earlier end curve leaves
/// both behind).
inline void drop_unused(Brep& b) {
    std::vector<int> emap(b.edges.size(), -1), vmap(b.vertices.size(), -1);
    for (const BFace& F : b.faces)
        for (const auto& L : F.loops)
            for (const BCoedge& c : L) emap[static_cast<std::size_t>(c.edge)] = 0;
    std::vector<BEdge> edges;
    for (std::size_t i = 0; i < b.edges.size(); ++i)
        if (emap[i] == 0) {
            emap[i] = static_cast<int>(edges.size());
            edges.push_back(b.edges[i]);
            vmap[static_cast<std::size_t>(b.edges[i].v0)] = 0;
            vmap[static_cast<std::size_t>(b.edges[i].v1)] = 0;
        }
    std::vector<BVertex> vertices;
    for (std::size_t i = 0; i < b.vertices.size(); ++i)
        if (vmap[i] == 0) {
            vmap[i] = static_cast<int>(vertices.size());
            vertices.push_back(b.vertices[i]);
        }
    for (BEdge& e : edges) {
        e.v0 = vmap[static_cast<std::size_t>(e.v0)];
        e.v1 = vmap[static_cast<std::size_t>(e.v1)];
    }
    for (BFace& F : b.faces)
        for (auto& L : F.loops)
            for (BCoedge& c : L) c.edge = emap[static_cast<std::size_t>(c.edge)];
    b.edges = std::move(edges);
    b.vertices = std::move(vertices);
}

} // namespace proc_fillet_detail

struct ProcFilletResult {
    std::optional<Brep> solid;
    std::shared_ptr<const solve::BlendCarrier> blend;
    std::string error;
    /// round_edge: the supports' carriers (theta = 0 on the first, 1 on the
    /// second) and their maps, face coordinate k = alpha[i][k] x[2 i + k] +
    /// beta[i][k] of the spine coordinates x.
    int support[2] = {-1, -1};
    double alpha[2][2] = {{1.0, 1.0}, {1.0, 1.0}}, beta[2][2] = {{0.0, 0.0}, {0.0, 0.0}};
};

/// The block under `top` (degree (p, q), controls (a i/p, L j/q, z_ij, 1))
/// before rounding, as a patch solid (built: its edges and vertices are
/// found from the shared sides).
inline PatchSolid freeform_block(const PolyPatchD& top) {
    using geom::Vec4d;
    const int p = top.p, q = top.q;
    auto at = [&](int i, int j) { return top.at(i, j); };
    auto flat = [](Vec4d h) { h.z = 0.0; return h; };
    PatchSolid S;
    // Each face is built with outward Su x Sv (see the header).
    PolyPatchD bottom(q, p), x0(1, q), xb(q, 1), endA(p, 1), endL(1, p);
    for (int j = 0; j <= q; ++j)
        for (int i = 0; i <= p; ++i) bottom.at(j, i) = flat(at(i, j));
    for (int j = 0; j <= q; ++j) {
        x0.at(0, j) = flat(at(0, j));
        x0.at(1, j) = at(0, j);
        xb.at(j, 0) = flat(at(p, j));
        xb.at(j, 1) = at(p, j);
    }
    for (int i = 0; i <= p; ++i) {
        endA.at(i, 0) = flat(at(i, 0));
        endA.at(i, 1) = at(i, 0);
        endL.at(0, i) = flat(at(i, q));
        endL.at(1, i) = at(i, q);
    }
    S.faces = {bottom, x0, xb, top, endA, endL};
    std::string err;
    S.build(err);
    return S;
}

/// Round the edge between the top and the side x = a by radius r.
///
/// `cut` in (0, 1]: the block ends at y = cut L instead of y = L. That end,
/// like any end where the top is not flat across the end plane, is trimmed
/// by a traced curve (see the file comment). `cut` should be
/// dyadic so that the restricted faces are exact.
inline ProcFilletResult round_block_edge(const PolyPatchD& top, double r, const solve::SsiOptions& opt = {},
                                         double cut = 1.0) {
    using namespace proc_fillet_detail;
    ProcFilletResult res;
    const int p = top.p, q = top.q;
    if (!(cut > 0.0 && cut <= 1.0)) { res.error = "the cut must be in (0, 1]"; return res; }
    if (q < 1) { res.error = "the top needs degree >= 1 across the ends"; return res; }
    // An end where the top has zero slope across the end plane has the
    // ball's end section in that plane: a planar arc ends the blend there.
    // At any other end (or at a cut) the supports run on past the end and
    // the blend is trimmed by its certified curve with the end face.
    bool flat0 = q >= 3, flat1 = q >= 3;
    for (int i = 0; i <= p; ++i) {
        flat0 = flat0 && top.at(i, 0).z == top.at(i, 1).z;
        flat1 = flat1 && top.at(i, q).z == top.at(i, q - 1).z;
    }
    const bool trim0 = !flat0, trim1 = cut < 1.0 || !flat1;
    // Supports: the top over v in [-e0, 1 + e1] (the polynomial continued),
    // with the solid's top at w in [c0, c1] of the support's coordinate w.
    const double e0 = trim0 ? 0.25 : 0.0, e1 = !flat1 ? 0.25 : 0.0;
    const PolyPatchD A = (e0 > 0.0 || e1 > 0.0) ? top.restricted(0.0, 1.0, -e0, 1.0 + e1) : top;
    const double span = 1.0 + e0 + e1, c0 = e0 / span, c1 = (e0 + cut) / span, cw = c1 - c0;
    const PolyPatchD B = freeform_block(A).faces[2];
    const PatchSolid block = freeform_block(cut < 1.0 ? top.restricted(0.0, 1.0, 0.0, cut) : top);
    enum { kBottom, kX0, kB, kTop, kEnd0, kEndL };

    // The spine: the ball inside the material, opposite the outward normals.
    const geom::OffsetCarrier oa(A, -r), ob(B, -r);
    const solve::SsiResult sp = solve::intersect_carriers(oa, ob, opt);
    if (!sp.complete || sp.curves.size() != 1) {
        res.error = "expected exactly one spine with a completeness proof";
        return res;
    }
    auto blend = std::make_shared<solve::BlendCarrier>(A, B, r, sp.curves[0].segs, -1);
    res.blend = blend;
    solve::Pt4 x0, x1;
    if (!blend->spine_at(0.0, x0) || !blend->spine_at(1.0, x1)) { res.error = "spine ends"; return res; }
    // Orient sigma from y = 0 to y = L (support coordinate w = x[1]).
    const bool sigma_up = x0[1] < x1[1];
    const solve::Pt4 sp0 = sigma_up ? x0 : x1, spL = sigma_up ? x1 : x0;
    if (!trim0 && !(std::fabs(sp0[1]) < 1e-9 && std::fabs(sp0[2]) < 1e-9)) {
        res.error = "the spine does not end on the end plane y = 0";
        return res;
    }
    if (!trim1 && !(std::fabs(spL[1] - 1.0) < 1e-9 && std::fabs(spL[2] - 1.0) < 1e-9)) {
        res.error = "the spine does not end on the end plane y = L";
        return res;
    }

    Brep b;
    b.carriers = block.faces;
    b.procs.push_back(blend);
    // Each end: sigma of the contacts there (A on the top, B on the side)
    // and, where trimmed, the curve with the end face, run from theta = 0 to 1.
    auto end_cut = [&](int plane_face, double& sA_, double& sB_, std::shared_ptr<TracedEdge>& curve) {
        auto plane = std::make_shared<geom::PatchCarrier>(b.carriers[static_cast<std::size_t>(plane_face)]);
        const solve::SsiResult cr = solve::intersect_carriers(*blend, *plane, opt);
        if (!cr.complete || cr.curves.size() != 1) {
            res.error = "expected exactly one curve where the blend meets an end face, with a completeness proof";
            return false;
        }
        curve = std::make_shared<TracedEdge>();
        curve->a = blend;
        curve->b = plane;
        curve->proc[0] = 0;
        curve->carrier[1] = plane_face;
        curve->path.segs = cr.curves[0].segs;
        curve->path.make_breaks(curve->sys(), *blend);
        curve->path.find_kinks(curve->sys(), *curve->a, *curve->b);
        solve::Pt4 q0, q1;
        curve->from = 0.0;
        curve->to = 1.0;
        if (!curve->at(0.0, q0) || !curve->at(1.0, q1)) { res.error = "trimming curve ends"; return false; }
        if (q0[1] > q1[1]) {
            std::swap(q0, q1);
            curve->from = 1.0;
            curve->to = 0.0;
        }
        if (!(std::fabs(q0[1]) < 1e-9 && std::fabs(q1[1] - 1.0) < 1e-9)) {
            res.error = "the trimming curve does not run from contact to contact";
            return false;
        }
        sA_ = q0[0];
        sB_ = q1[0];
        return true;
    };
    const double s0 = sigma_up ? 0.0 : 1.0;
    double sA0 = s0, sB0 = s0, sA = 1.0 - s0, sB = 1.0 - s0;
    std::shared_ptr<TracedEdge> start_curve, end_curve;
    if (trim0 && !end_cut(kEnd0, sA0, sB0, start_curve)) return res;
    if (trim1 && !end_cut(kEndL, sA, sB, end_curve)) return res;

    auto vertex = [&](const Vec3d& pt) {
        b.vertices.push_back({});
        b.vertices.back().point = pt;
        return static_cast<int>(b.vertices.size()) - 1;
    };
    auto corner = [&](int face, double u, double v) {
        return geom::BezierSurface(b.carriers[static_cast<std::size_t>(face)]).point(u, v);
    };
    auto chart = [&](int vtx, int face, double u, double v) {
        b.vertices[static_cast<std::size_t>(vtx)].charts.push_back({face, Interval(u), Interval(v)});
    };
    // Vertices.
    const int P000 = vertex(corner(kBottom, 0, 0)), P0L0 = vertex(corner(kBottom, 1, 0)), PaL0 = vertex(corner(kBottom, 1, 1)),
              Pa00 = vertex(corner(kBottom, 0, 1)), T00 = vertex(corner(kTop, 0, 0)), T0L = vertex(corner(kTop, 0, 1));
    Vec3d P, Ps, Pt;
    blend->eval(sA0, 0, P, Ps, Pt);
    const int A0 = vertex(P);
    blend->eval(sA, 0, P, Ps, Pt);
    const int AL = vertex(P);
    blend->eval(sB0, 1, P, Ps, Pt);
    const int B0 = vertex(P);
    blend->eval(sB, 1, P, Ps, Pt);
    const int BL = vertex(P);
    solve::Pt4 xA0, xB0, xA, xB;
    blend->spine_at(sA0, xA0);
    blend->spine_at(sB0, xB0);
    blend->spine_at(sA, xA);
    blend->spine_at(sB, xB);
    const double uA0 = xA0[0], uAL = xA[0], vB0 = xB0[3], vBL = xB[3];
    chart(P000, kBottom, 0, 0); chart(P000, kX0, 0, 0); chart(P000, kEnd0, 0, 0);
    chart(P0L0, kBottom, 1, 0); chart(P0L0, kX0, 0, 1); chart(P0L0, kEndL, 0, 0);
    chart(PaL0, kBottom, 1, 1); chart(PaL0, kB, 1, 0);  chart(PaL0, kEndL, 0, 1);
    chart(Pa00, kBottom, 0, 1); chart(Pa00, kB, 0, 0);  chart(Pa00, kEnd0, 1, 0);
    chart(T00, kX0, 1, 0);      chart(T00, kTop, 0, 0); chart(T00, kEnd0, 0, 1);
    chart(T0L, kX0, 1, 1);      chart(T0L, kTop, 0, 1); chart(T0L, kEndL, 1, 0);
    chart(A0, kTop, uA0, 0); chart(A0, kEnd0, uA0, 1);
    chart(B0, kB, 0, vB0);   chart(B0, kEnd0, 1, vB0);
    // At the far end, on the end face (u: height fraction, v: x / a): A is
    // on its side u = 1 and B on its side v = 1; the other coordinate is
    // that of the contact on the top or the side.
    chart(AL, kTop, uAL, 1); chart(AL, kEndL, 1, uAL);
    chart(BL, kB, 1, vBL);   chart(BL, kEndL, vBL, 1);

    // Side edges: (carrier, side) pairs, edge parameter = side parameter on both.
    auto side_edge = [&](int v0, int v1, int f0, int sd0, int f1, int sd1, double t0, double t1) {
        BEdge e;
        e.v0 = v0;
        e.v1 = v1;
        e.on[0] = {f0, sd0, false};
        e.on[1] = {f1, sd1, false};
        e.t0 = Interval(t0);
        e.t1 = Interval(t1);
        b.edges.push_back(e);
        return static_cast<int>(b.edges.size()) - 1;
    };
    const int E1 = side_edge(P000, P0L0, kBottom, 0, kX0, 3, 0, 1);
    const int E2 = side_edge(P0L0, PaL0, kBottom, 1, kEndL, 3, 0, 1);
    const int E3 = side_edge(Pa00, PaL0, kBottom, 2, kB, 0, 0, 1);
    const int E4 = side_edge(P000, Pa00, kBottom, 3, kEnd0, 0, 0, 1);
    const int E5 = side_edge(P000, T00, kX0, 0, kEnd0, 3, 0, 1);
    const int E6 = side_edge(T00, T0L, kX0, 1, kTop, 3, 0, 1);
    const int E7 = side_edge(P0L0, T0L, kX0, 2, kEndL, 0, 0, 1);
    const int E8 = side_edge(PaL0, BL, kB, 1, kEndL, 2, 0, vBL);
    const int E9 = side_edge(Pa00, B0, kB, 3, kEnd0, 1, 0, vB0);
    const int E10 = side_edge(T00, A0, kTop, 0, kEnd0, 2, 0, uA0);
    const int E11 = side_edge(T0L, AL, kTop, 2, kEndL, 1, 0, uAL);
    // Procedural edges: parts of the blend's square sides, and the trimming curve.
    auto add_edge = [&](int v0, int v1, std::shared_ptr<const ProcEdge> pe) {
        BEdge e;
        e.v0 = v0;
        e.v1 = v1;
        e.proc = std::move(pe);
        b.edges.push_back(e);
        return static_cast<int>(b.edges.size()) - 1;
    };
    auto contact = [&](int v0, int v1, double theta, double from, double to, int support, int coord, double sc0, double sc1,
                       double sh0, double sh1) {
        auto pe = std::make_shared<BlendSide>();
        pe->blend = blend;
        pe->proc = 0;
        pe->contact = true;
        pe->fixed = theta;
        pe->from = from;
        pe->to = to;
        pe->support = support;
        pe->coord = coord;
        pe->scale[0] = sc0;
        pe->scale[1] = sc1;
        pe->shift[0] = sh0;
        pe->shift[1] = sh1;
        return add_edge(v0, v1, pe);
    };
    auto section = [&](int v0, int v1, double sigma, int plane) {
        auto pe = std::make_shared<BlendSide>();
        pe->blend = blend;
        pe->proc = 0;
        pe->contact = false;
        pe->fixed = sigma;
        pe->plane = plane;
        pe->plane_patch = b.carriers[static_cast<std::size_t>(plane)];
        return add_edge(v0, v1, pe);
    };
    // The supports' coordinates on the solid's faces: (w - c0) / (c1 - c0)
    // for v on the top and s on the side.
    const int EcA = contact(A0, AL, 0.0, sA0, sA, kTop, 0, 1.0, cw, 0.0, c0);
    const int EcB = contact(B0, BL, 1.0, sB0, sB, kB, 2, cw, 1.0, c0, 0.0);
    const int Es0 = trim0 ? add_edge(A0, B0, start_curve) : section(A0, B0, s0, kEnd0);
    const int EsL = trim1 ? add_edge(AL, BL, end_curve) : section(AL, BL, sA, kEndL);

    // Loops: counter-clockwise in each face's parameters; a coedge is
    // reversed when the loop meets the edge's v1 first.
    auto loop = [&](std::initializer_list<std::pair<int, int>> path) { // (edge, from vertex)
        std::vector<BCoedge> L;
        for (const auto& [e, from] : path) L.push_back({e, b.edges[static_cast<std::size_t>(e)].v0 != from});
        return L;
    };
    auto face = [&](int carrier, std::vector<BCoedge> L) {
        BFace f;
        f.carrier = carrier;
        f.loops.push_back(std::move(L));
        b.faces.push_back(std::move(f));
    };
    face(kBottom, loop({{E1, P000}, {E2, P0L0}, {E3, PaL0}, {E4, Pa00}}));
    face(kX0, loop({{E5, P000}, {E6, T00}, {E7, T0L}, {E1, P0L0}}));
    face(kB, loop({{E3, Pa00}, {E8, PaL0}, {EcB, BL}, {E9, B0}}));
    face(kTop, loop({{E10, T00}, {EcA, A0}, {E11, AL}, {E6, T0L}}));
    face(kEnd0, loop({{E4, P000}, {E9, Pa00}, {Es0, B0}, {E10, A0}, {E5, T00}}));
    face(kEndL, loop({{E7, P0L0}, {E11, T0L}, {EsL, AL}, {E8, BL}, {E2, PaL0}}));
    // The blend: counter-clockwise in (sigma, theta). With sigma rising from
    // the end y = 0 the region lies at larger sigma than that end.
    {
        BFace f;
        f.proc = 0;
        if (sigma_up) f.loops.push_back(loop({{EcA, A0}, {EsL, AL}, {EcB, BL}, {Es0, B0}}));
        else f.loops.push_back(loop({{EcA, AL}, {Es0, A0}, {EcB, B0}, {EsL, BL}}));
        // Outward: away from the ball centre.
        const double sm = 0.25 * (sA0 + sB0 + sA + sB);
        Vec3d Q, Qs, Qt;
        blend->eval(sm, 0.5, Q, Qs, Qt);
        solve::Pt4 x;
        blend->spine_at(sm, x);
        Vec3d c, cu, cv;
        geom::OffsetCarrier(A, -r).eval(x[0], x[1], c, cu, cv);
        f.flipped = dot(cross(Qs, Qt), Q - c) < 0.0;
        b.faces.push_back(std::move(f));
    }
    res.solid = std::move(b);
    return res;
}

/// Round a convex side edge of any Brep by radius r (docs/FEATURES.md §5.1,
/// fillets on general solids).
///
/// The edge must be a side edge between two faces with one loop each, and
/// each end vertex must have exactly three faces: the two supports and one
/// end face, met through two side edges (or, at a mitre, see below). The
/// supports are the two faces' carriers, turned so that Su x Sv is outward
/// and continued a quarter of their length past both ends along the edge;
/// the spine is their certified -r offset intersection. At
/// each end the blend is trimmed by its certified curve with the end face's
/// carrier, so the end face may be curved and meet the edge at any angle.
///
/// The edge becomes the contact curve on the first face, the end vertex
/// becomes the contact point there, and the edges along the end faces are
/// shortened to the contacts. Checked: convexity at the edge's midpoint, one
/// spine and one trimming curve per end with completeness proofs, the
/// contacts inside the edges they split and (on samples) the contact curves
/// inside their faces' parameter squares. Not checked yet: the ball meeting
/// faces other than the three at each end (a radius too large for the
/// solid).
///
/// An end face may also be an earlier blend whose contact runs along one of
/// the two faces (a corner of two fillets). The new blend then meets it in
/// a mitre through the point A where both touch that face. A and the curve's
/// two branches at A are proven by solve/blend_corner.hpp (a polar blow-up:
/// the blends touch at A, so the ordinary tracer cannot prove the curve
/// there). The rest of the curve is traced on the new blend with a small
/// rectangle around A removed. Where it ends depends on the radii. A
/// smaller new radius: on the new blend's far contact, inside the earlier
/// blend's end curve, which is split there. Equal radii: at the end curve's
/// far vertex, and the end curve goes. A larger new radius: at P on the
/// earlier blend's far contact (split there); the new blend then meets that
/// contact's face along a traced curve from P to Q on the face's edge with
/// the end curve's face (split there), and the end curve and its far vertex
/// go. The earlier contact on the shared face is split at A. Three blends:
/// when the earlier end curve was itself cut short by a second earlier
/// blend b2, and the new blend's far contact crosses b2's contact at a
/// second touching point B, the new blend runs from A along its mitre with
/// the first earlier blend to the point T where all three meet, then along
/// its mitre with b2 to B (both touching points proven, both mitres traced
/// on the region between their rectangles).
inline ProcFilletResult round_edge(const Brep& S, int edge, double r, const solve::SsiOptions& opt = {}) {
    using namespace proc_fillet_detail;
    ProcFilletResult res;
    auto fail = [&](const char* m) { res.error = m; return res; };
    if (edge < 0 || edge >= static_cast<int>(S.edges.size())) return fail("no such edge");
    if (!(r > 0.0)) return fail("the radius must be positive");
    const BEdge E = S.edges[static_cast<std::size_t>(edge)];
    if (E.proc || E.arc) return fail("the edge must be a side edge");
    auto face_of = [&](int carrier) {
        int f = -1, n = 0;
        for (std::size_t i = 0; i < S.faces.size(); ++i)
            if (S.faces[i].proc < 0 && S.faces[i].carrier == carrier) { f = static_cast<int>(i); ++n; }
        return n == 1 ? f : -1;
    };
    auto one_loop = [&](int f) { return S.faces[static_cast<std::size_t>(f)].loops.size() == 1; };
    auto faces_with = [&](int e) {
        std::vector<int> out;
        for (std::size_t i = 0; i < S.faces.size(); ++i)
            for (const auto& L : S.faces[i].loops)
                for (const BCoedge& c : L)
                    if (c.edge == e) out.push_back(static_cast<int>(i));
        return out;
    };
    const int CA = E.on[0].carrier, CB = E.on[1].carrier, FA = face_of(CA), FB = face_of(CB);
    if (FA < 0 || FB < 0 || !one_loop(FA) || !one_loop(FB)) return fail("both faces at the edge must be one carrier with one loop");

    // A support: the face's carrier turned outward and continued along the
    // edge; face coordinate k = alpha[k] w + beta[k] of support coordinate w.
    struct Support {
        PolyPatchD P;
        double alpha[2] = {1.0, 1.0}, beta[2] = {0.0, 0.0};
    };
    const double ext = 0.25, span = 1.0 + 2.0 * ext;
    auto support = [&](int f, int side) {
        Support s;
        const BFace& F = S.faces[static_cast<std::size_t>(f)];
        const PolyPatchD& C = S.carriers[static_cast<std::size_t>(F.carrier)];
        s.P = C;
        if (F.flipped) {
            for (int i = 0; i <= C.p; ++i)
                for (int j = 0; j <= C.q; ++j) s.P.at(i, j) = C.at(C.p - i, j);
            s.alpha[0] = -1.0;
            s.beta[0] = 1.0;
        }
        const int along = (side == 0 || side == 2) ? 0 : 1;
        s.P = along == 0 ? s.P.restricted(-ext, 1.0 + ext, 0.0, 1.0) : s.P.restricted(0.0, 1.0, -ext, 1.0 + ext);
        s.beta[along] -= s.alpha[along] * ext;
        s.alpha[along] *= span;
        return s;
    };
    const Support SA = support(FA, E.on[0].side), SB = support(FB, E.on[1].side);

    // Convex: at the edge's midpoint, the direction into the first face
    // points against the second face's outward normal.
    {
        const double tm = 0.5 * (E.t0.lo + E.t1.hi);
        auto frame = [&](int f, const BEdge::SideRef& ref, Vec3d& n, Vec3d& in) {
            double u, v;
            side_uv(ref.side, ref.reversed ? 1.0 - tm : tm, u, v);
            const BFace& F = S.faces[static_cast<std::size_t>(f)];
            const geom::Jet3 J = geom::BezierSurface(S.carriers[static_cast<std::size_t>(F.carrier)]).jet(u, v);
            n = cross(J[1], J[2]) * (F.flipped ? -1.0 : 1.0);
            const int sd = ref.side;
            in = sd == 0 ? J[2] : sd == 1 ? J[1] * -1.0 : sd == 2 ? J[2] * -1.0 : J[1];
        };
        Vec3d nA, inA, nB, inB;
        frame(FA, E.on[0], nA, inA);
        frame(FB, E.on[1], nB, inB);
        if (!(dot(inA, nB) < 0.0 && dot(inB, nA) < 0.0)) return fail("the edge is not convex");
    }

    // The spine: the ball inside the material.
    const geom::OffsetCarrier oa(SA.P, -r), ob(SB.P, -r);
    const solve::SsiResult sp = solve::intersect_carriers(oa, ob, opt);
    if (!sp.complete || sp.curves.size() != 1) return fail("expected exactly one spine with a completeness proof");
    auto blend = std::make_shared<solve::BlendCarrier>(SA.P, SB.P, r, sp.curves[0].segs, -1);
    res.blend = blend;

    Brep b = S;
    b.outward_known = false;
    const int pi = static_cast<int>(b.procs.size());
    b.procs.push_back(blend);

    // The ends.
    struct End {
        int V = -1, eA = -1, eB = -1, CE = -1, FE = -1, B = -1;
        double sA = 0.0, sB = 0.0;
        std::shared_ptr<TracedEdge> curve;
        bool corner = false; ///< the end face is an earlier blend: a mitre
        std::vector<BCoedge> chain; ///< edges from vertex V (theta 0) to B (theta 1)
        std::vector<BCoedge> fe_chain; ///< what goes into the end face's loop, if not `chain`
        int drop = -1;              ///< an edge the end removes from every loop
    };
    End end[2];
    auto other_edge = [&](int f, int V) {
        int found = -1, n = 0;
        for (const BCoedge& c : S.faces[static_cast<std::size_t>(f)].loops[0]) {
            const BEdge& e = S.edges[static_cast<std::size_t>(c.edge)];
            if (c.edge != edge && (e.v0 == V || e.v1 == V)) { found = c.edge; ++n; }
        }
        return n == 1 ? found : -1;
    };
    auto other_carrier = [&](int e, int c) {
        const BEdge& x = S.edges[static_cast<std::size_t>(e)];
        return x.on[0].carrier == c ? x.on[1].carrier : x.on[0].carrier;
    };
    auto reversed_chain = [](const std::vector<BCoedge>& c) {
        std::vector<BCoedge> out(c.rbegin(), c.rend());
        for (BCoedge& x : out) x.reversed = !x.reversed;
        return out;
    };
    // Split side edge e (one end at V) at the contact, which is at (u, v) on
    // carrier c; the end at V becomes vertex W.
    auto split = [&](int e, int V, int c, double u, double v, int W, int CE) {
        BEdge& x = b.edges[static_cast<std::size_t>(e)];
        const BEdge::SideRef* ref = side_on(x, c);
        const double sp_ = (ref->side == 0 || ref->side == 2) ? u : v;
        const double t = ref->reversed ? 1.0 - sp_ : sp_;
        if (!(t > std::min(x.t0.hi, x.t1.hi) && t < std::max(x.t0.lo, x.t1.lo))) return false; // strictly inside
        if (x.v0 == V) { x.v0 = W; x.t0 = Interval(t); }
        else { x.v1 = W; x.t1 = Interval(t); }
        const BEdge::SideRef* re = side_on(x, CE);
        double ue, ve;
        side_uv(re->side, re->reversed ? 1.0 - t : t, ue, ve);
        BVertex& w = b.vertices[static_cast<std::size_t>(W)];
        w.charts.push_back({c, Interval(u), Interval(v)});
        w.charts.push_back({CE, Interval(ue), Interval(ve)});
        return true;
    };
    auto support_of = [&](int j) -> const Support& { return j == 0 ? SA : SB; };
    // An end at an earlier blend's face: the new blend meets it in a mitre
    // through the point A where both touch the shared face (solve/
    // blend_corner.hpp). Returns an error or null.
    
    auto corner_end = [&](End& en, int k) -> const char* {
        (void)k;
        const std::vector<int> fa = faces_with(en.eA), fb = faces_with(en.eB);
        int FE = -1;
        for (int f : fa)
            if (f != FA && std::find(fb.begin(), fb.end(), f) != fb.end()) FE = f;
        if (FE < 0) return "an end vertex must have exactly three faces";
        const BFace& Fe = S.faces[static_cast<std::size_t>(FE)];
        if (Fe.proc < 0) return "an end face with procedural edges must be a blend";
        auto b1 = std::dynamic_pointer_cast<const solve::BlendCarrier>(S.procs[static_cast<std::size_t>(Fe.proc)]);
        if (!b1) return "the end face is not a blend";
        en.FE = FE;
        // The shared face: the one whose other edge at V is the earlier
        // blend's contact; the other face's edge is its end curve.
        const int eSh[2] = {en.eA, en.eB}, Csh[2] = {CA, CB};
        int j = -1;
        std::shared_ptr<const BlendSide> bs;
        for (int i = 0; i < 2; ++i) {
            auto c = std::dynamic_pointer_cast<const BlendSide>(S.edges[static_cast<std::size_t>(eSh[i])].proc);
            if (c && c->contact && c->support == Csh[i]) {
                if (j >= 0) return "both faces at the end hold contacts of the earlier blend";
                j = i;
                bs = c;
            }
        }
        if (j < 0) return "no face at the end holds a contact of the earlier blend";
        const int eS = eSh[j], eO = eSh[1 - j], CO = Csh[1 - j];
        auto other = std::dynamic_pointer_cast<const TracedEdge>(S.edges[static_cast<std::size_t>(eO)].proc);
        if (!other || other->proc[0] != Fe.proc || other->carrier[1] != CO)
            return "the other face's edge at the end is not the earlier blend's end curve";
        solve::ContactOnFace f1, f2;
        f1.blend = b1.get();
        f1.coord = bs->coord;
        f1.theta_c = bs->fixed;
        for (int i = 0; i < 2; ++i) {
            f1.alpha[i] = 1.0 / bs->scale[i];
            f1.beta[i] = -bs->shift[i] / bs->scale[i];
        }
        f2.blend = blend.get();
        f2.coord = 2 * j;
        f2.theta_c = j;
        for (int i = 0; i < 2; ++i) {
            f2.alpha[i] = support_of(j).alpha[i];
            f2.beta[i] = support_of(j).beta[i];
        }
        // A: guesses from V on the earlier contact and the nearest sample
        // of the new contact in the shared face.
        const BEdge& xs = S.edges[static_cast<std::size_t>(eS)];
        const double g1 = xs.v0 == en.V ? bs->from : bs->to;
        double w1[2], g2 = 0.5, best = 1e300;
        if (!solve::corner_detail::face_w(f1, g1, w1)) return "the earlier contact at the end";
        for (int i = 0; i <= 400; ++i) {
            double w2[2];
            if (!solve::corner_detail::face_w(f2, i / 400.0, w2)) continue;
            const double d = std::sqrt((w2[0] - w1[0]) * (w2[0] - w1[0]) + (w2[1] - w1[1]) * (w2[1] - w1[1]));
            if (d < best) { best = d; g2 = i / 400.0; }
        }
        const auto A = solve::corner_point(f1, f2, g1, g2);
        if (!A) return "the contacts on the shared face do not cross transversally";
        const double s1A = num::mid(A->s1), s2A = num::mid(A->s2);
        if (!(s1A > std::min(bs->from, bs->to) && s1A < std::max(bs->from, bs->to)))
            return "the new contact meets the earlier one outside its edge";
        // A rectangle around a touching point P of an earlier blend (bp) and
        // the new one: the local branch on the earlier blend, the
        // rectangle's size, and the strips beside it proven free of bp.
        struct Rect {
            solve::CornerLocal L;
            double rin = 0.0, sc0 = 0.0, sc1 = 0.0, tb = 0.0;
        };
        auto rect_at = [&](const solve::ContactOnFace& ca, const solve::ContactOnFace& cb, const solve::CornerPoint& P,
                           const geom::Carrier& bp, double rho) -> std::optional<Rect> {
            const auto Lp = solve::corner_local(ca, cb, P, rho);
            if (!Lp) return std::nullopt;
            const solve::CornerBranch* br = nullptr;
            int n_on = 0;
            for (const auto& x : Lp->branches)
                if (x.on_first) { br = &x; ++n_on; }
            if (n_on != 1) return std::nullopt;
            const double tau = num::mid(br->m);
            const double da = br->chart == 0 ? 1.0 - tau * tau : tau * tau - 1.0, db = 2.0 * tau;
            if (!(db > 0.0)) return std::nullopt;
            const double slope = std::max(std::fabs(da / db), 1e-3);
            const double w = 0.95 * rho / std::sqrt(1.0 + 0.81 / (slope * slope));
            Rect q;
            q.L = *Lp;
            q.rin = 0.9 * w / slope;
            if (!(q.rin > 1e-5 * rho)) return std::nullopt;
            const double sP = num::mid(P.s2);
            q.sc0 = sP - Lp->k_s * w;
            q.sc1 = sP + Lp->k_s * w;
            if (q.sc0 > q.sc1) std::swap(q.sc0, q.sc1);
            if (!(q.sc0 > 0.0 && q.sc1 < 1.0)) return std::nullopt;
            q.tb = cb.theta_c + Lp->k_t * q.rin;
            const double sl = std::min(cb.theta_c, q.tb), sh = std::max(cb.theta_c, q.tb);
            const solve::RectCarrier left(blend, 0.0, q.sc0, sl, sh), right(blend, q.sc1, 1.0, sl, sh);
            const solve::SsiResult cl = solve::intersect_carriers(left, bp, opt);
            const solve::SsiResult cr = solve::intersect_carriers(right, bp, opt);
            if (!cl.complete || !cl.curves.empty() || !cr.complete || !cr.curves.empty()) return std::nullopt;
            return q;
        };
        // Three blends at the corner: the earlier end curve (eO) was cut
        // short by a second earlier blend b2, whose mitre with b1 (eM) and
        // contact on CO (eC2) meet at eO's far vertex W. When the new
        // blend's far contact crosses b2's contact on CO at B, the new blend
        // runs A -> mitre with b1 -> T -> mitre with b2 -> B, T where all
        // three blends meet; eO, W and the pieces of eM and eC2 past T and B
        // go.
        {
            const BEdge& xo = S.edges[static_cast<std::size_t>(eO)];
            const int W = xo.v0 == en.V ? xo.v1 : xo.v0;
            auto neighbour = [&](int f, int e, int Vx) {
                for (const BCoedge& c : S.faces[static_cast<std::size_t>(f)].loops[0]) {
                    const BEdge& x = S.edges[static_cast<std::size_t>(c.edge)];
                    if (c.edge != e && (x.v0 == Vx || x.v1 == Vx)) return c.edge;
                }
                return -1;
            };
            int FO = -1;
            for (int f : faces_with(eO))
                if (f != FE) FO = f;
            const int eM = FO >= 0 ? neighbour(FE, eO, W) : -1, eC2 = FO >= 0 ? neighbour(FO, eO, W) : -1;
            int FE2 = -1;
            if (eM >= 0)
                for (int f : faces_with(eM))
                    if (f != FE) FE2 = f;
            std::shared_ptr<const solve::BlendCarrier> b2;
            std::shared_ptr<const BlendSide> bc;
            if (FE2 >= 0 && eC2 >= 0 && S.faces[static_cast<std::size_t>(FE2)].proc >= 0) {
                const int p2 = S.faces[static_cast<std::size_t>(FE2)].proc;
                b2 = std::dynamic_pointer_cast<const solve::BlendCarrier>(S.procs[static_cast<std::size_t>(p2)]);
                bc = std::dynamic_pointer_cast<const BlendSide>(S.edges[static_cast<std::size_t>(eC2)].proc);
                if (bc && !(bc->contact && bc->proc == p2 && bc->support == CO)) bc = nullptr;
            }
            std::optional<solve::CornerPoint> Bp;
            solve::ContactOnFace k1, k2;
            if (b2 && bc) {
                k1.blend = b2.get();
                k1.coord = bc->coord;
                k1.theta_c = bc->fixed;
                for (int i = 0; i < 2; ++i) {
                    k1.alpha[i] = 1.0 / bc->scale[i];
                    k1.beta[i] = -bc->shift[i] / bc->scale[i];
                }
                k2.blend = blend.get();
                k2.coord = 2 * (1 - j);
                k2.theta_c = 1 - j;
                for (int i = 0; i < 2; ++i) {
                    k2.alpha[i] = support_of(1 - j).alpha[i];
                    k2.beta[i] = support_of(1 - j).beta[i];
                }
                const BEdge& xc = S.edges[static_cast<std::size_t>(eC2)];
                const double h1 = xc.v0 == W ? bc->from : bc->to;
                double wc[2], h2 = 0.5, bst = 1e300;
                if (solve::corner_detail::face_w(k1, h1, wc)) {
                    for (int i = 0; i <= 400; ++i) {
                        double w2[2];
                        if (!solve::corner_detail::face_w(k2, i / 400.0, w2)) continue;
                        const double d = std::sqrt((w2[0] - wc[0]) * (w2[0] - wc[0]) + (w2[1] - wc[1]) * (w2[1] - wc[1]));
                        if (d < bst) { bst = d; h2 = i / 400.0; }
                    }
                    Bp = solve::corner_point(k1, k2, h1, h2);
                }
                if (Bp && !(num::mid(Bp->s1) > std::min(bc->from, bc->to) && num::mid(Bp->s1) < std::max(bc->from, bc->to)))
                    Bp.reset();
            }
            if (Bp) {
                const int p2 = S.faces[static_cast<std::size_t>(FE2)].proc;
                const double s1B = num::mid(Bp->s1), s3B = num::mid(Bp->s2);
                // h_b(P) = distance to b's spine - its radius.
                auto hdist = [&](const solve::BlendCarrier& bb, const Vec3d& P, double& sg) {
                    sg = spine_foot(bb, P, sg);
                    return norm(P - bb.section(sg).centre) - bb.radius();
                };
                // The parameter in [0, 1] of a curve where it crosses blend bb
                // inside bb's theta range: the sign changes of h by samples,
                // each refined by bisection; exactly one must lie on bb.
                auto crossing = [&](const std::function<Vec3d(double)>& pt, const solve::BlendCarrier& bb, double sg0,
                                    double& tx) {
                    const int n = 200;
                    double sg = sg0, prev = hdist(bb, pt(0.0), sg);
                    int found = 0;
                    for (int i = 1; i <= n; ++i) {
                        const double h = hdist(bb, pt(static_cast<double>(i) / n), sg);
                        if ((h > 0.0) != (prev > 0.0)) {
                            double lo = static_cast<double>(i - 1) / n, hi = static_cast<double>(i) / n, sl = sg;
                            const double flo = hdist(bb, pt(lo), sl);
                            for (int it = 0; it < 60; ++it) {
                                const double m = 0.5 * (lo + hi);
                                double sm = sl;
                                if ((hdist(bb, pt(m), sm) > 0.0) == (flo > 0.0)) lo = m;
                                else hi = m;
                            }
                            if (in_wedge(bb, pt(0.5 * (lo + hi)), sl)) {
                                tx = 0.5 * (lo + hi);
                                ++found;
                            }
                        }
                        prev = h;
                    }
                    return found == 1;
                };
                // The rectangles: the largest radius that works at each point.
                std::optional<Rect> RA, RB;
                for (double rho : {1e-3, 5e-4, 2.5e-4, 1.25e-4})
                    if (!RA) RA = rect_at(f1, f2, *A, *b1, rho);
                for (double rho : {1e-3, 5e-4, 2.5e-4, 1.25e-4})
                    if (RA && !RB) RB = rect_at(k1, k2, *Bp, *b2, rho);
                for (int once = 0; once < 1 && RA && RB; ++once) {
                    const double rl = std::min(RA->tb, RB->tb), rh = std::max(RA->tb, RB->tb);
                    // The rectangles at A and B are free of the other earlier blend.
                    {
                        const solve::RectCarrier sa(blend, RA->sc0, RA->sc1, std::min(f2.theta_c, RA->tb), std::max(f2.theta_c, RA->tb));
                        const solve::RectCarrier sb(blend, RB->sc0, RB->sc1, std::min(k2.theta_c, RB->tb), std::max(k2.theta_c, RB->tb));
                        const solve::SsiResult ca = solve::intersect_carriers(sa, *b2, opt);
                        const solve::SsiResult cb = solve::intersect_carriers(sb, *b1, opt);
                        if (!ca.complete || !ca.curves.empty() || !cb.complete || !cb.curves.empty()) continue;
                    }
                    auto R = std::make_shared<solve::RectCarrier>(blend, 0.0, 1.0, rl, rh);
                    // The traced mitre with earlier blend bp, starting on the
                    // strip side v = vs of R.
                    auto traced = [&](const std::shared_ptr<const solve::BlendCarrier>& bp, int pp, double vs,
                                      const Rect& q) -> std::shared_ptr<TracedEdge> {
                        const solve::SsiResult cm = solve::intersect_carriers(*R, *bp, opt);
                        if (!cm.complete || cm.curves.size() != 1) return nullptr;
                        auto t = std::make_shared<TracedEdge>();
                        t->a = R;
                        t->b = bp;
                        t->proc[0] = pi;
                        t->proc[1] = pp;
                        t->off[1] = rl;
                        t->scl[1] = rh - rl;
                        t->path.segs = cm.curves[0].segs;
                        t->path.make_breaks(t->sys(), *t->a);
                        t->path.find_kinks(t->sys(), *t->a, *t->b);
                        solve::Pt4 a0, a1;
                        if (!t->at(0.0, a0) || !t->at(1.0, a1)) return nullptr;
                        if (std::fabs(a0[1] - vs) > std::fabs(a1[1] - vs)) {
                            std::swap(a0, a1);
                            t->from = 1.0;
                            t->to = 0.0;
                        }
                        if (!(std::fabs(a0[1] - vs) < 1e-9 && a0[0] > q.sc0 && a0[0] < q.sc1)) return nullptr;
                        return t;
                    };
                    const double vA = RA->tb == rl ? 0.0 : 1.0, vB = 1.0 - vA;
                    const auto t1 = traced(b1, Fe.proc, vA, *RA);
                    const auto t2 = t1 ? traced(b2, p2, vB, *RB) : nullptr;
                    if (!t1 || !t2) continue;
                    // T on both mitres and on the earlier mitre eM.
                    const auto& xm = S.edges[static_cast<std::size_t>(eM)];
                    double tT1 = 0.0, tT2 = 0.0, tTM = 0.0;
                    if (!crossing([&](double x) { return t1->point(x); }, *b2, s1B, tT1)) continue;
                    if (!crossing([&](double x) { return t2->point(x); }, *b1, s1A, tT2)) continue;
                    if (!crossing([&](double x) { return xm.proc->point(x); }, *blend, s3B, tTM)) continue;
                    const Vec3d PT = t1->point(tT1);
                    if (!(norm(t2->point(tT2) - PT) < 1e-9 && norm(xm.proc->point(tTM) - PT) < 1e-9)) continue;
                    // The corner pieces.
                    auto corner_piece = [&](const std::shared_ptr<const solve::BlendCarrier>& bp, int pp,
                                            const solve::ContactOnFace& ka, const solve::ContactOnFace& kb, const Rect& q,
                                            const TracedEdge& t, double sP) {
                        auto c = std::make_shared<CornerEdge>();
                        c->b1 = bp;
                        c->b2 = blend;
                        c->f1 = ka;
                        c->f2 = kb;
                        c->L = q.L;
                        c->p1 = pp;
                        c->p2 = pi;
                        solve::Pt4 a0;
                        t.at(0.0, a0);
                        c->a_end = (a0[0] - sP) / q.L.k_s;
                        c->b_end = q.rin;
                        return c;
                    };
                    const auto ceA = corner_piece(b1, Fe.proc, f1, f2, *RA, *t1, s2A);
                    const auto ceB = corner_piece(b2, p2, k1, k2, *RB, *t2, s3B);
                    if (!(norm(ceA->point(1.0) - t1->point(0.0)) < 1e-9 && norm(ceB->point(1.0) - t2->point(0.0)) < 1e-9))
                        continue;
                    // Commit. Vertices: side 0 of the new blend keeps V.
                    auto push_vertex = [&]() {
                        b.vertices.push_back({});
                        return static_cast<int>(b.vertices.size()) - 1;
                    };
                    const int vA_ = j == 0 ? en.V : push_vertex();
                    const int vB_ = j == 0 ? push_vertex() : en.V;
                    en.B = j == 0 ? vB_ : vA_;
                    auto contact_vertex = [&](int side, double s3, int carrier) {
                        BVertex v;
                        solve::Pt4 x;
                        Vec3d P, Ps, Pt;
                        blend->eval(s3, static_cast<double>(side), P, Ps, Pt);
                        v.point = P;
                        if (blend->spine_at(s3, x)) {
                            const Support& sp_ = support_of(side);
                            v.charts.push_back({carrier, Interval(sp_.alpha[0] * x[static_cast<std::size_t>(2 * side)] + sp_.beta[0]),
                                                Interval(sp_.alpha[1] * x[static_cast<std::size_t>(2 * side + 1)] + sp_.beta[1])});
                        }
                        return v;
                    };
                    b.vertices[static_cast<std::size_t>(vA_)] = contact_vertex(j, s2A, Csh[j]);
                    b.vertices[static_cast<std::size_t>(vB_)] = contact_vertex(1 - j, s3B, CO);
                    const int vXA = push_vertex(), vXB = push_vertex(), vT = push_vertex();
                    b.vertices[static_cast<std::size_t>(vXA)].point = t1->point(0.0);
                    b.vertices[static_cast<std::size_t>(vXB)].point = t2->point(0.0);
                    b.vertices[static_cast<std::size_t>(vT)].point = PT;
                    // The earlier edges: b1's contact at A, b2's contact at B,
                    // the earlier mitre at T.
                    const double tS = (s1A - bs->from) / (bs->to - bs->from);
                    BEdge& ns = b.edges[static_cast<std::size_t>(eS)];
                    if (ns.v0 == en.V) { ns.proc = bs->part(tS, 1.0); ns.v0 = vA_; }
                    else { ns.proc = bs->part(0.0, tS); ns.v1 = vA_; }
                    const double tB = (s1B - bc->from) / (bc->to - bc->from);
                    BEdge& nc = b.edges[static_cast<std::size_t>(eC2)];
                    if (nc.v0 == W) { nc.proc = bc->part(tB, 1.0); nc.v0 = vB_; }
                    else { nc.proc = bc->part(0.0, tB); nc.v1 = vB_; }
                    BEdge& nm = b.edges[static_cast<std::size_t>(eM)];
                    if (nm.v0 == W) { nm.proc = xm.proc->part(tTM, 1.0); nm.v0 = vT; }
                    else { nm.proc = xm.proc->part(0.0, tTM); nm.v1 = vT; }
                    if (!ns.proc || !nc.proc || !nm.proc) return "an edge at the end cannot be split";
                    en.drop = eO;
                    // The new edges: A -> XA -> T and B -> XB -> T.
                    auto add_edge = [&](int v0, int v1, std::shared_ptr<const ProcEdge> pe) {
                        BEdge e;
                        e.v0 = v0;
                        e.v1 = v1;
                        e.proc = std::move(pe);
                        b.edges.push_back(e);
                        return static_cast<int>(b.edges.size()) - 1;
                    };
                    const int iA = add_edge(vA_, vXA, ceA), i1 = add_edge(vXA, vT, t1->part(0.0, tT1));
                    const int iB = add_edge(vB_, vXB, ceB), i2 = add_edge(vXB, vT, t2->part(0.0, tT2));
                    const std::vector<BCoedge> a2t = {{iA, false}, {i1, false}}, b2t = {{iB, false}, {i2, false}};
                    std::vector<BCoedge> a2b = a2t;
                    for (const BCoedge& c : reversed_chain(b2t)) a2b.push_back(c);
                    en.chain = j == 0 ? a2b : reversed_chain(a2b);
                    en.fe_chain = j == 0 ? a2t : reversed_chain(a2t);
                    // b2's loop: T -> B between the shortened eM and eC2.
                    std::vector<BCoedge>& L2 = b.faces[static_cast<std::size_t>(FE2)].loops[0];
                    const int n2 = static_cast<int>(L2.size());
                    int iM = -1, iC = -1;
                    for (int i = 0; i < n2; ++i) {
                        if (L2[static_cast<std::size_t>(i)].edge == eM) iM = i;
                        if (L2[static_cast<std::size_t>(i)].edge == eC2) iC = i;
                    }
                    const std::vector<BCoedge> t2b = reversed_chain(b2t);
                    if (iM >= 0 && iC >= 0 && (iM + 1) % n2 == iC) L2.insert(L2.begin() + iM + 1, t2b.begin(), t2b.end());
                    else if (iM >= 0 && iC >= 0 && (iC + 1) % n2 == iM) L2.insert(L2.begin() + iC + 1, b2t.begin(), b2t.end());
                    else return "the second earlier blend's edges at the corner are not consecutive";
                    en.sA = j == 0 ? s2A : s3B;
                    en.sB = j == 0 ? s3B : s2A;
                    return nullptr;
                }
                return "the corner of three blends could not be proven";
            }
        }
        // The local branches, then the regular rest around them.
        for (double rho : {1e-3, 5e-4, 2.5e-4}) {
            const auto L = solve::corner_local(f1, f2, *A, rho);
            if (!L) continue;
            const solve::CornerBranch* br = nullptr;
            int n_on = 0;
            for (const auto& x : L->branches)
                if (x.on_first) { br = &x; ++n_on; }
            if (n_on != 1) continue;
            const double tau = num::mid(br->m);
            const double da = br->chart == 0 ? 1.0 - tau * tau : tau * tau - 1.0, db = 2.0 * tau;
            // The central rectangle [-w, w] x [0, rin] lies in the half-disc
            // of radius rho, and the branch leaves it through its far side.
            if (!(db > 0.0)) continue;
            const double slope = std::max(std::fabs(da / db), 1e-3);
            const double w = 0.95 * rho / std::sqrt(1.0 + 0.81 / (slope * slope)), rin = 0.9 * w / slope;
            if (!(rin > 1e-5 * rho)) continue;
            const double ks = L->k_s, kt = L->k_t;
            const double sc0 = s2A - ks * w, sc1 = s2A + ks * w;
            if (!(sc0 > 0.0 && sc1 < 1.0)) continue;
            const double tb = f2.theta_c + kt * rin; // the strip's far side
            const double sl = std::min(f2.theta_c, tb), sh = std::max(f2.theta_c, tb);
            const double rl = f2.theta_c == 0.0 ? tb : 0.0, rh = f2.theta_c == 0.0 ? 1.0 : tb;
            auto R = std::make_shared<solve::RectCarrier>(blend, 0.0, 1.0, rl, rh);
            const solve::RectCarrier left(blend, 0.0, sc0, sl, sh), right(blend, sc1, 1.0, sl, sh);
            const solve::SsiResult cl = solve::intersect_carriers(left, *b1, opt);
            const solve::SsiResult cr = solve::intersect_carriers(right, *b1, opt);
            if (!cl.complete || !cl.curves.empty() || !cr.complete || !cr.curves.empty()) continue;
            const solve::SsiResult cm = solve::intersect_carriers(*R, *b1, opt);
            if (!cm.complete || cm.curves.size() != 1) continue;
            auto te = std::make_shared<TracedEdge>();
            te->a = R;
            te->b = b1;
            te->proc[0] = pi;
            te->proc[1] = Fe.proc;
            te->off[1] = rl;
            te->scl[1] = rh - rl;
            te->path.segs = cm.curves[0].segs;
            te->path.make_breaks(te->sys(), *te->a);
            te->path.find_kinks(te->sys(), *te->a, *te->b);
            solve::Pt4 q0, q1;
            if (!te->at(0.0, q0) || !te->at(1.0, q1)) continue;
            const double vj = f2.theta_c == 0.0 ? 0.0 : 1.0; // the strip side in R's v
            if (std::fabs(q0[1] - vj) > std::fabs(q1[1] - vj)) {
                std::swap(q0, q1);
                te->from = 1.0;
                te->to = 0.0;
            }
            // Where the mitre ends. On the new blend's far contact, inside the
            // earlier blend's end curve: the smaller radius second. On the
            // earlier blend's far contact: the larger radius second, and the
            // new blend goes on to meet that contact's face. On both: equal
            // radii, at the earlier blend's end vertex.
            if (!(std::fabs(q0[1] - vj) < 1e-9)) continue;
            const bool far2 = std::fabs(q1[1] - (1.0 - vj)) < 1e-9;
            const bool far1 = std::fabs(q1[3] - (1.0 - f1.theta_c)) < 1e-9;
            if (!far1 && !far2) continue;
            if (!(q0[0] > sc0 && q0[0] < sc1)) continue;
            // The pieces, and the vertices: X joins them; T (touching) and F
            // (far end) are the chain's ends.
            auto ce = std::make_shared<CornerEdge>();
            ce->b1 = b1;
            ce->b2 = blend;
            ce->f1 = f1;
            ce->f2 = f2;
            ce->L = *L;
            ce->p1 = Fe.proc;
            ce->p2 = pi;
            ce->a_end = (q0[0] - s2A) / ks;
            ce->b_end = rin;
            const BEdge& xo = S.edges[static_cast<std::size_t>(eO)];
            const int A0 = xo.v0 == en.V ? xo.v1 : xo.v0; // the end curve's far vertex
            const double tS = (s1A - bs->from) / (bs->to - bs->from);
            BVertex vt, vf;
            vt.point = ce->point(0.0);
            {
                solve::Pt4 x;
                if (!blend->spine_at(s2A, x)) continue;
                const Support& sp_ = support_of(j);
                const double u = sp_.alpha[0] * x[static_cast<std::size_t>(2 * j)] + sp_.beta[0];
                const double v = sp_.alpha[1] * x[static_cast<std::size_t>(2 * j + 1)] + sp_.beta[1];
                vt.charts.push_back({Csh[j], Interval(u), Interval(v)});
            }
            double s2F = q1[0];
            // The smaller radius second: the end curve split where the mitre
            // meets it.
            double tO = -1.0;
            if (far2 && !far1) {
                const double t1B = q1[3];
                double lo = 0.0, hi = 1.0, slo = 0.0, shi = 0.0, tmp;
                other->on_proc(Fe.proc, 0.0, tmp, slo);
                other->on_proc(Fe.proc, 1.0, tmp, shi);
                if (!((slo - t1B) * (shi - t1B) < 0.0)) continue;
                for (int it = 0; it < 80; ++it) {
                    const double m = 0.5 * (lo + hi);
                    double sm, tm;
                    other->on_proc(Fe.proc, m, sm, tm);
                    if ((tm - t1B) * (slo - t1B) > 0.0) lo = m;
                    else hi = m;
                }
                tO = 0.5 * (lo + hi);
                vf.point = te->point(1.0);
                double uo, vo, duo, dvo;
                if (other->on_carrier(CO, tO, uo, vo, duo, dvo)) vf.charts.push_back({CO, Interval(uo), Interval(vo)});
            }
            // Equal radii: the mitre ends at the end curve's far vertex.
            if (far1 && far2 && !(norm(te->point(1.0) - S.vertices[static_cast<std::size_t>(A0)].point) < 1e-9)) continue;
            // The larger radius second: at P on the earlier blend's far
            // contact (edge eT, on carrier CT), the new blend goes on along
            // CT to Q on the edge eFT that CT shares with the end curve's
            // face. The end curve and its far vertex go.
            int eT = -1, eFT = -1, CT = -1, FT = -1;
            std::shared_ptr<const BlendSide> bsT;
            std::shared_ptr<const ProcEdge> te2;
            double tPT = 0.0, uQ = 0.0, vQ = 0.0;
            if (far1 && !far2) {
                auto neighbour = [&](int f, int e, int V) {
                    for (const BCoedge& c : S.faces[static_cast<std::size_t>(f)].loops[0]) {
                        const BEdge& x = S.edges[static_cast<std::size_t>(c.edge)];
                        if (c.edge != e && (x.v0 == V || x.v1 == V)) return c.edge;
                    }
                    return -1;
                };
                int FO = -1;
                for (int f : faces_with(eO))
                    if (f != FE) FO = f;
                if (FO < 0) return "the earlier blend's end curve has no other face";
                eT = neighbour(FE, eO, A0);
                eFT = neighbour(FO, eO, A0);
                if (eT < 0 || eFT < 0) return "the earlier blend's end vertex";
                bsT = std::dynamic_pointer_cast<const BlendSide>(S.edges[static_cast<std::size_t>(eT)].proc);
                if (!bsT || !bsT->contact || bsT->proc != Fe.proc) return "the earlier blend's far contact";
                CT = bsT->support;
                FT = face_of(CT);
                const BEdge& xft = S.edges[static_cast<std::size_t>(eFT)];
                if (FT < 0 || xft.proc || xft.arc || !side_on(xft, CT) || !side_on(xft, CO))
                    return "the larger radius second needs a plain edge between the faces past the earlier blend";
                tPT = (q1[2] - bsT->from) / (bsT->to - bsT->from);
                if (!(tPT > 0.0 && tPT < 1.0)) continue;
                auto plane = std::make_shared<geom::PatchCarrier>(S.carriers[static_cast<std::size_t>(CT)]);
                const solve::SsiResult cp = solve::intersect_carriers(*R, *plane, opt);
                if (!cp.complete || cp.curves.size() != 1) continue;
                auto tc = std::make_shared<TracedEdge>();
                tc->a = R;
                tc->b = plane;
                tc->proc[0] = pi;
                tc->carrier[1] = CT;
                tc->off[1] = rl;
                tc->scl[1] = rh - rl;
                tc->path.segs = cp.curves[0].segs;
                tc->path.make_breaks(tc->sys(), *tc->a);
                tc->path.find_kinks(tc->sys(), *tc->a, *tc->b);
                solve::Pt4 c0, c1;
                if (!tc->at(0.0, c0) || !tc->at(1.0, c1)) continue;
                // Q: the end on the new blend's far contact; P where R's v
                // is the mitre's.
                const double tQ = std::fabs(c1[1] - (1.0 - vj)) < 1e-9 ? 1.0 : 0.0;
                const solve::Pt4& cq = tQ == 1.0 ? c1 : c0;
                if (!(std::fabs(cq[1] - (1.0 - vj)) < 1e-9)) continue;
                const double vP = q1[1];
                double lo = 0.0, hi = 1.0;
                auto vof = [&](double t) {
                    solve::Pt4 x;
                    tc->at(t, x);
                    return x[1] - vP;
                };
                const double flo = vof(lo);
                if (!(flo * vof(hi) < 0.0)) continue;
                for (int it = 0; it < 80; ++it) {
                    const double m = 0.5 * (lo + hi);
                    if (vof(m) * flo > 0.0) lo = m;
                    else hi = m;
                }
                const double tP = 0.5 * (lo + hi);
                if (!(norm(tc->point(tP) - te->point(1.0)) < 1e-9)) continue;
                te2 = tc->part(tP, tQ);
                double du, dv, sq, tq;
                if (!te2->on_carrier(CT, 1.0, uQ, vQ, du, dv) || !te2->on_proc(pi, 1.0, sq, tq)) continue;
                s2F = sq;
                vf.point = te2->point(1.0);
            }
            // Commit: vertices. Side 0 of the new blend keeps V.
            int vT, vF;
            const bool equal = far1 && far2;
            auto push_vertex = [&]() {
                b.vertices.push_back({});
                return static_cast<int>(b.vertices.size()) - 1;
            };
            if (j == 0) {
                vT = en.V;
                vF = equal ? A0 : push_vertex();
                en.B = vF;
            } else {
                vF = en.V;
                vT = push_vertex();
                en.B = vT;
            }
            b.vertices[static_cast<std::size_t>(vT)] = vt;
            if (!equal) b.vertices[static_cast<std::size_t>(vF)] = vf;
            BEdge& ns = b.edges[static_cast<std::size_t>(eS)];
            if (ns.v0 == en.V) { ns.proc = bs->part(tS, 1.0); ns.v0 = vT; }
            else { ns.proc = bs->part(0.0, tS); ns.v1 = vT; }
            if (!ns.proc) return "an edge at the end cannot be split";
            if (tO >= 0.0) {
                BEdge& no = b.edges[static_cast<std::size_t>(eO)];
                if (no.v0 == en.V) { no.proc = other->part(tO, 1.0); no.v0 = vF; }
                else { no.proc = other->part(0.0, tO); no.v1 = vF; }
                if (!no.proc) return "an edge at the end cannot be split";
            } else en.drop = eO;
            BVertex vx;
            vx.point = te->point(0.0);
            b.vertices.push_back(vx);
            const int vX = static_cast<int>(b.vertices.size()) - 1;
            const int vP = te2 ? push_vertex() : vF;
            BEdge e1, e2;
            e1.v0 = vT; e1.v1 = vX; e1.proc = ce;
            e2.v0 = vX; e2.v1 = vP; e2.proc = te;
            b.edges.push_back(e1);
            const int i1 = static_cast<int>(b.edges.size()) - 1;
            b.edges.push_back(e2);
            const int i2 = static_cast<int>(b.edges.size()) - 1;
            std::vector<BCoedge> t2f = {{i1, false}, {i2, false}};
            en.fe_chain = t2f;
            if (te2) {
                // P on the earlier far contact; the contact, and the edge
                // past it on CT, shortened; the new piece in CT's face.
                BVertex vp;
                vp.point = te->point(1.0);
                double up, vp_, dup, dvp;
                if (bsT->on_carrier(CT, tPT, up, vp_, dup, dvp)) vp.charts.push_back({CT, Interval(up), Interval(vp_)});
                b.vertices[static_cast<std::size_t>(vP)] = vp;
                BEdge& nt = b.edges[static_cast<std::size_t>(eT)];
                if (nt.v0 == A0) { nt.proc = bsT->part(tPT, 1.0); nt.v0 = vP; }
                else { nt.proc = bsT->part(0.0, tPT); nt.v1 = vP; }
                if (!split(eFT, A0, CT, uQ, vQ, vF, CO)) return "the new blend meets the edge past the earlier blend outside it";
                BEdge e3;
                e3.v0 = vP; e3.v1 = vF; e3.proc = te2;
                b.edges.push_back(e3);
                const int i3 = static_cast<int>(b.edges.size()) - 1;
                t2f.push_back({i3, false});
                std::vector<BCoedge>& LT = b.faces[static_cast<std::size_t>(FT)].loops[0];
                const int n = static_cast<int>(LT.size());
                int iT = -1, iF = -1;
                for (int i = 0; i < n; ++i) {
                    if (LT[static_cast<std::size_t>(i)].edge == eT) iT = i;
                    if (LT[static_cast<std::size_t>(i)].edge == eFT) iF = i;
                }
                if (iT < 0 || iF < 0) return "the face past the earlier blend does not hold both edges";
                if ((iT + 1) % n == iF) LT.insert(LT.begin() + iT + 1, BCoedge{i3, false});
                else if ((iF + 1) % n == iT) LT.insert(LT.begin() + iF + 1, BCoedge{i3, true});
                else return "the face past the earlier blend: its edges at the vertex are not consecutive";
            }
            if (equal && j == 1) {
                // V takes the far vertex's place.
                b.vertices[static_cast<std::size_t>(en.V)] = S.vertices[static_cast<std::size_t>(A0)];
                for (BEdge& x : b.edges) {
                    if (x.v0 == A0) x.v0 = en.V;
                    if (x.v1 == A0) x.v1 = en.V;
                }
            }
            if (!te2) en.fe_chain = t2f;
            en.chain = j == 0 ? t2f : reversed_chain(t2f);
            if (j == 1) en.fe_chain = reversed_chain(en.fe_chain);
            en.sA = j == 0 ? s2A : s2F;
            en.sB = j == 0 ? s2F : s2A;
            return nullptr;
        }
        return "the mitre at the end could not be proven";
    };
    for (int k = 0; k < 2; ++k) {
        End& en = end[k];
        en.V = k == 0 ? E.v0 : E.v1;
        en.eA = other_edge(FA, en.V);
        en.eB = other_edge(FB, en.V);
        if (en.eA < 0 || en.eB < 0) return fail("an end vertex needs one more edge on each face");
        const BEdge& xa = S.edges[static_cast<std::size_t>(en.eA)];
        const BEdge& xb = S.edges[static_cast<std::size_t>(en.eB)];
        if (xa.proc || xb.proc) {
            en.corner = true;
            if (const char* m = corner_end(en, k)) return fail(m);
            continue;
        }
        if (xa.arc || xb.arc) return fail("arc edges at an end vertex are not handled");
        en.CE = other_carrier(en.eA, CA);
        if (other_carrier(en.eB, CB) != en.CE) return fail("an end vertex must have exactly three faces");
        en.FE = face_of(en.CE);
        if (en.FE < 0) return fail("the end face must be one carrier");
        auto plane = std::make_shared<geom::PatchCarrier>(b.carriers[static_cast<std::size_t>(en.CE)]);
        const solve::SsiResult cr = solve::intersect_carriers(*blend, *plane, opt);
        if (!cr.complete || cr.curves.size() != 1)
            return fail("expected exactly one curve where the blend meets an end face, with a completeness proof");
        auto curve = std::make_shared<TracedEdge>();
        curve->a = blend;
        curve->b = plane;
        curve->proc[0] = pi;
        curve->carrier[1] = en.CE;
        curve->path.segs = cr.curves[0].segs;
        curve->path.make_breaks(curve->sys(), *blend);
        curve->path.find_kinks(curve->sys(), *curve->a, *curve->b);
        solve::Pt4 q0, q1;
        if (!curve->at(0.0, q0) || !curve->at(1.0, q1)) return fail("trimming curve ends");
        if (q0[1] > q1[1]) {
            std::swap(q0, q1);
            curve->from = 1.0;
            curve->to = 0.0;
        }
        if (!(std::fabs(q0[1]) < 1e-9 && std::fabs(q1[1] - 1.0) < 1e-9))
            return fail("the trimming curve does not run from contact to contact");
        en.sA = q0[0];
        en.sB = q1[0];
        en.curve = curve;
    }
    const bool sigma_up = end[0].sA < end[1].sA;
    if (sigma_up != (end[0].sB < end[1].sB)) return fail("the trimming curves cross");

    // Contacts: face coordinates, split parameters, vertices.
    auto face_uv = [&](const Support& s, int coord, double sigma, double& u, double& v) {
        solve::Pt4 x;
        if (!blend->spine_at(sigma, x)) return false;
        u = s.alpha[0] * x[static_cast<std::size_t>(coord)] + s.beta[0];
        v = s.alpha[1] * x[static_cast<std::size_t>(coord + 1)] + s.beta[1];
        return true;
    };
    for (End& en : end) {
        if (en.corner) continue;
        Vec3d P, Ps, Pt;
        double u, v;
        // A: the end vertex itself, moved to the contact on the first face.
        BVertex& a = b.vertices[static_cast<std::size_t>(en.V)];
        blend->eval(en.sA, 0.0, P, Ps, Pt);
        a = BVertex{};
        a.point = P;
        if (!face_uv(SA, 0, en.sA, u, v) || !split(en.eA, en.V, CA, u, v, en.V, en.CE))
            return fail("a contact is not inside the edge it splits");
        blend->eval(en.sB, 1.0, P, Ps, Pt);
        b.vertices.push_back({});
        b.vertices.back().point = P;
        en.B = static_cast<int>(b.vertices.size()) - 1;
        if (!face_uv(SB, 2, en.sB, u, v) || !split(en.eB, en.V, CB, u, v, en.B, en.CE))
            return fail("a contact is not inside the edge it splits");
    }
    // The contact curves inside their faces (samples).
    for (int side = 0; side < 2; ++side) {
        const Support& s = side == 0 ? SA : SB;
        const double s0 = side == 0 ? end[0].sA : end[0].sB, s1 = side == 0 ? end[1].sA : end[1].sB;
        for (int i = 0; i <= 64; ++i) {
            double u, v;
            if (!face_uv(s, 2 * side, s0 + (s1 - s0) * i / 64.0, u, v) || !(u > -1e-12 && u < 1 + 1e-12 && v > -1e-12 && v < 1 + 1e-12))
                return fail("a contact curve leaves its face: the radius is too large");
        }
    }

    // Edges: the old edge becomes the contact on the first face; a new one
    // on the second face; the trimming curves.
    auto contact = [&](double theta, double from, double to, int carrier, int coord, const Support& s) {
        auto pe = std::make_shared<BlendSide>();
        pe->blend = blend;
        pe->proc = pi;
        pe->contact = true;
        pe->fixed = theta;
        pe->from = from;
        pe->to = to;
        pe->support = carrier;
        pe->coord = coord;
        for (int k = 0; k < 2; ++k) {
            pe->scale[k] = 1.0 / s.alpha[k];
            pe->shift[k] = -s.beta[k] / s.alpha[k];
        }
        return pe;
    };
    {
        BEdge e;
        e.v0 = E.v0;
        e.v1 = E.v1;
        e.proc = contact(0.0, end[0].sA, end[1].sA, CA, 0, SA);
        b.edges[static_cast<std::size_t>(edge)] = e;
    }
    BEdge eb;
    eb.v0 = end[0].B;
    eb.v1 = end[1].B;
    eb.proc = contact(1.0, end[0].sB, end[1].sB, CB, 2, SB);
    b.edges.push_back(eb);
    const int EcB = static_cast<int>(b.edges.size()) - 1;
    for (BCoedge& c : b.faces[static_cast<std::size_t>(FB)].loops[0])
        if (c.edge == edge) c.edge = EcB;
    for (int k = 0; k < 2; ++k) {
        End& en = end[k];
        if (!en.corner) {
            BEdge e;
            e.v0 = en.V;
            e.v1 = en.B;
            e.proc = en.curve;
            b.edges.push_back(e);
            en.chain = {BCoedge{static_cast<int>(b.edges.size()) - 1, false}};
        }
        // In the end face's loop, between the two edges it shortened.
        std::vector<BCoedge>& L = b.faces[static_cast<std::size_t>(en.FE)].loops[0];
        const int n = static_cast<int>(L.size());
        int iA = -1, iB = -1;
        for (int i = 0; i < n; ++i) {
            if (L[static_cast<std::size_t>(i)].edge == en.eA) iA = i;
            if (L[static_cast<std::size_t>(i)].edge == en.eB) iB = i;
        }
        if (iA < 0 || iB < 0) return fail("the end face does not hold both edges at the vertex");
        const std::vector<BCoedge>& ch = en.fe_chain.empty() ? en.chain : en.fe_chain;
        if ((iA + 1) % n == iB) L.insert(L.begin() + iA + 1, ch.begin(), ch.end());
        else if ((iB + 1) % n == iA) {
            const std::vector<BCoedge> rc = reversed_chain(ch);
            L.insert(L.begin() + iB + 1, rc.begin(), rc.end());
        } else return fail("the end face's edges at the vertex are not consecutive");
        if (en.drop >= 0)
            for (BFace& F : b.faces)
                for (auto& lp : F.loops)
                    lp.erase(std::remove_if(lp.begin(), lp.end(), [&](const BCoedge& c) { return c.edge == en.drop; }),
                             lp.end());
    }
    // The blend: counter-clockwise in (sigma, theta).
    {
        BFace f;
        f.proc = pi;
        std::vector<BCoedge> L;
        auto add = [&](const std::vector<BCoedge>& c) { L.insert(L.end(), c.begin(), c.end()); };
        if (sigma_up) {
            L.push_back({edge, false});
            add(end[1].chain);
            L.push_back({EcB, true});
            add(reversed_chain(end[0].chain));
        } else {
            L.push_back({edge, true});
            add(end[0].chain);
            L.push_back({EcB, false});
            add(reversed_chain(end[1].chain));
        }
        f.loops.push_back(std::move(L));
        const double sm = 0.25 * (end[0].sA + end[0].sB + end[1].sA + end[1].sB);
        Vec3d Q, Qs, Qt;
        blend->eval(sm, 0.5, Q, Qs, Qt);
        solve::Pt4 x;
        blend->spine_at(sm, x);
        Vec3d c, cu, cv;
        oa.eval(x[0], x[1], c, cu, cv);
        f.flipped = dot(cross(Qs, Qt), Q - c) < 0.0;
        b.faces.push_back(std::move(f));
    }
    // Clearance: the blend, over its sigma span, meets no face of the solid
    // other than those at its own vertices (the supports, the end faces, and
    // at a corner the faces at the earlier blend's end vertex). Faces whose
    // enclosure misses the blend's are skipped; the rest must have no
    // intersection curve with it, with a completeness proof. Conservative:
    // a face is refused when its carrier, not only its trimmed region, meets
    // the blend.
    {
        const int nb = static_cast<int>(b.faces.size()) - 1;
        std::vector<char> at_blend(b.vertices.size(), 0);
        for (const BCoedge& c : b.faces[static_cast<std::size_t>(nb)].loops[0]) {
            at_blend[static_cast<std::size_t>(b.edges[static_cast<std::size_t>(c.edge)].v0)] = 1;
            at_blend[static_cast<std::size_t>(b.edges[static_cast<std::size_t>(c.edge)].v1)] = 1;
        }
        double smin = 1e300, smax = -1e300;
        for (const End& en : end) {
            smin = std::min({smin, en.sA, en.sB});
            smax = std::max({smax, en.sA, en.sB});
        }
        auto boxes = [](const geom::Carrier& C, double u0, double u1, int nu, int nv, std::vector<geom::IVec3>& out) {
            for (int i = 0; i < nu; ++i)
                for (int k = 0; k < nv; ++k) {
                    geom::IVec3 P, Pu, Pv;
                    if (!C.jet1_over(u0 + (u1 - u0) * i / nu, u0 + (u1 - u0) * (i + 1) / nu, static_cast<double>(k) / nv,
                                     static_cast<double>(k + 1) / nv, P, Pu, Pv))
                        return false;
                    out.push_back(P);
                }
            return true;
        };
        std::vector<geom::IVec3> bb;
        const bool bb_ok = boxes(*blend, smin, smax, 16, 2, bb);
        auto rb = std::make_shared<solve::RectCarrier>(blend, smin, smax, 0.0, 1.0);
        for (int f = 0; f < nb; ++f) {
            const BFace& F = b.faces[static_cast<std::size_t>(f)];
            bool touches = false;
            for (const auto& L : F.loops)
                for (const BCoedge& c : L) {
                    const BEdge& x = b.edges[static_cast<std::size_t>(c.edge)];
                    touches = touches || at_blend[static_cast<std::size_t>(x.v0)] || at_blend[static_cast<std::size_t>(x.v1)];
                }
            if (touches) continue;
            std::shared_ptr<const geom::Carrier> surf;
            if (F.proc >= 0) surf = b.procs[static_cast<std::size_t>(F.proc)];
            else surf = std::make_shared<geom::PatchCarrier>(b.carriers[static_cast<std::size_t>(F.carrier)]);
            std::vector<geom::IVec3> fb;
            if (bb_ok && boxes(*surf, 0.0, 1.0, 4, 4, fb)) {
                bool near = false;
                auto meet = [](const num::Interval& x, const num::Interval& y) { return x.lo <= y.hi && y.lo <= x.hi; };
                for (const geom::IVec3& p : bb)
                    for (const geom::IVec3& q : fb) near = near || (meet(p.x, q.x) && meet(p.y, q.y) && meet(p.z, q.z));
                if (!near) continue;
            }
            const solve::SsiResult cr = solve::intersect_carriers(*rb, *surf, opt);
            if (!cr.complete || !cr.curves.empty()) return fail("the ball meets another face of the solid: the radius is too large here");
        }
    }
    // Every loop closes: each coedge ends where the next begins.
    for (const BFace& F : b.faces)
        for (const auto& L : F.loops)
            for (std::size_t i = 0; i < L.size(); ++i) {
                const BEdge& x = b.edges[static_cast<std::size_t>(L[i].edge)];
                const BEdge& y = b.edges[static_cast<std::size_t>(L[(i + 1) % L.size()].edge)];
                const int xe = L[i].reversed ? x.v0 : x.v1, ys = L[(i + 1) % L.size()].reversed ? y.v1 : y.v0;
                if (xe != ys) return fail("internal: a loop does not close");
            }
    drop_unused(b);
    res.support[0] = CA;
    res.support[1] = CB;
    for (int k = 0; k < 2; ++k) {
        res.alpha[0][k] = SA.alpha[k];
        res.beta[0][k] = SA.beta[k];
        res.alpha[1][k] = SB.alpha[k];
        res.beta[1][k] = SB.beta[k];
    }
    res.solid = std::move(b);
    return res;
}

} // namespace ak::op
