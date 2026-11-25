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
#include <memory>
#include <optional>
#include <string>
#include <vector>

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
    auto ang = [](const Vec3d& x, const Vec3d& y) { return std::atan2(norm(cross(x, y)), dot(x, y)); };
    const double phi = ang(ua, ue), psi = ang(ua, up);
    // theta of the angle psi on a rational quadratic arc of angle phi.
    const double tp = std::tan(0.5 * psi), tq = std::tan(0.5 * (phi - psi));
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
};

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
/// rectangle around A removed. The earlier blend's contact and end curve
/// are split where the mitre meets them. Handled: the new radius smaller
/// than the earlier one, so that its far contact meets the earlier end
/// curve.
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
        std::vector<BCoedge> r(c.rbegin(), c.rend());
        for (BCoedge& x : r) x.reversed = !x.reversed;
        return r;
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
            const double d = std::hypot(w2[0] - w1[0], w2[1] - w1[1]);
            if (d < best) { best = d; g2 = i / 400.0; }
        }
        const auto A = solve::corner_point(f1, f2, g1, g2);
        if (!A) return "the contacts on the shared face do not cross transversally";
        const double s1A = num::mid(A->s1), s2A = num::mid(A->s2);
        if (!(s1A > std::min(bs->from, bs->to) && s1A < std::max(bs->from, bs->to)))
            return "the new contact meets the earlier one outside its edge";
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
            if (!(std::fabs(q0[1] - vj) < 1e-9 && std::fabs(q1[1] - (1.0 - vj)) < 1e-9)) continue;
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
            const double s1B = q1[2], t1B = q1[3], s2B = q1[0];
            // Split the earlier contact at A and its end curve at the far end.
            const double tS = (s1A - bs->from) / (bs->to - bs->from);
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
            const double tO = 0.5 * (lo + hi);
            (void)s1B;
            // Vertices: side 0 is V, side 1 new.
            b.vertices.push_back({});
            en.B = static_cast<int>(b.vertices.size()) - 1;
            const int vT = j == 0 ? en.V : en.B, vF = j == 0 ? en.B : en.V;
            BVertex vt, vf;
            vt.point = ce->point(0.0);
            vf.point = te->point(1.0);
            {
                solve::Pt4 x;
                if (!blend->spine_at(s2A, x)) continue;
                const Support& sp_ = support_of(j);
                const double u = sp_.alpha[0] * x[static_cast<std::size_t>(2 * j)] + sp_.beta[0];
                const double v = sp_.alpha[1] * x[static_cast<std::size_t>(2 * j + 1)] + sp_.beta[1];
                vt.charts.push_back({Csh[j], Interval(u), Interval(v)});
                double uo, vo, duo, dvo;
                if (other->on_carrier(CO, tO, uo, vo, duo, dvo)) vf.charts.push_back({CO, Interval(uo), Interval(vo)});
            }
            b.vertices[static_cast<std::size_t>(vT)] = vt;
            b.vertices[static_cast<std::size_t>(vF)] = vf;
            BEdge& ns = b.edges[static_cast<std::size_t>(eS)];
            if (ns.v0 == en.V) { ns.proc = bs->part(tS, 1.0); ns.v0 = vT; }
            else { ns.proc = bs->part(0.0, tS); ns.v1 = vT; }
            BEdge& no = b.edges[static_cast<std::size_t>(eO)];
            if (no.v0 == en.V) { no.proc = other->part(tO, 1.0); no.v0 = vF; }
            else { no.proc = other->part(0.0, tO); no.v1 = vF; }
            if (!ns.proc || !no.proc) return "an edge at the end cannot be split";
            BVertex vx;
            vx.point = te->point(0.0);
            b.vertices.push_back(vx);
            const int vX = static_cast<int>(b.vertices.size()) - 1;
            BEdge e1, e2;
            e1.v0 = vT; e1.v1 = vX; e1.proc = ce;
            e2.v0 = vX; e2.v1 = vF; e2.proc = te;
            b.edges.push_back(e1);
            const int i1 = static_cast<int>(b.edges.size()) - 1;
            b.edges.push_back(e2);
            const int i2 = static_cast<int>(b.edges.size()) - 1;
            const std::vector<BCoedge> t2f = {{i1, false}, {i2, false}};
            en.chain = j == 0 ? t2f : reversed_chain(t2f);
            en.sA = j == 0 ? s2A : s2B;
            en.sB = j == 0 ? s2B : s2A;
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
        if ((iA + 1) % n == iB) L.insert(L.begin() + iA + 1, en.chain.begin(), en.chain.end());
        else if ((iB + 1) % n == iA) {
            const std::vector<BCoedge> r = reversed_chain(en.chain);
            L.insert(L.begin() + iB + 1, r.begin(), r.end());
        } else return fail("the end face's edges at the vertex are not consecutive");
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
    // Every loop closes: each coedge ends where the next begins.
    for (const BFace& F : b.faces)
        for (const auto& L : F.loops)
            for (std::size_t i = 0; i < L.size(); ++i) {
                const BEdge& x = b.edges[static_cast<std::size_t>(L[i].edge)];
                const BEdge& y = b.edges[static_cast<std::size_t>(L[(i + 1) % L.size()].edge)];
                const int xe = L[i].reversed ? x.v0 : x.v1, ys = L[(i + 1) % L.size()].reversed ? y.v1 : y.v0;
                if (xe != ys) return fail("internal: a loop does not close");
            }
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
