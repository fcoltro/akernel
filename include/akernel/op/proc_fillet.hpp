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
// The ends. The top must have zero slope across the end planes (its first
// two control rows in y equal, and its last two). Then at each end the
// top's normal has no y component, the ball centre and both contacts lie in
// the end plane, and the blend's end section is a planar arc in the end
// face: the end faces are trimmed by it, with no end cap needed.
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
            u = x[static_cast<std::size_t>(coord)] / scale[0];
            v = x[static_cast<std::size_t>(coord + 1)] / scale[1];
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

} // namespace proc_fillet_detail

struct ProcFilletResult {
    std::optional<Brep> solid;
    std::shared_ptr<const solve::BlendCarrier> blend;
    std::string error;
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
/// `cut` in (0, 1]: the block ends at y = cut L instead of y = L. There the
/// top need not be flat across the end, so the ball's end section leaves the
/// end plane: the blend (built on the whole supports) is trimmed by its
/// intersection with the end face, a traced curve (op/proc_edge.hpp), and
/// the procedural face is no longer its whole square. `cut` should be
/// dyadic so that the restricted faces are exact.
inline ProcFilletResult round_block_edge(const PolyPatchD& top, double r, const solve::SsiOptions& opt = {},
                                         double cut = 1.0) {
    using namespace proc_fillet_detail;
    ProcFilletResult res;
    const int p = top.p, q = top.q;
    const bool trimmed = cut < 1.0;
    if (!(cut > 0.0 && cut <= 1.0)) { res.error = "the cut must be in (0, 1]"; return res; }
    if (q < 3) { res.error = "the top needs degree >= 3 across the ends"; return res; }
    for (int i = 0; i <= p; ++i)
        if (top.at(i, 0).z != top.at(i, 1).z || (!trimmed && top.at(i, q).z != top.at(i, q - 1).z)) {
            res.error = "the top must have zero slope across the end planes";
            return res;
        }
    // The supports of the blend are the whole faces; the solid's faces are
    // the block under the top restricted to v in [0, cut].
    const PolyPatchD B = freeform_block(top).faces[2];
    const PatchSolid block = freeform_block(trimmed ? top.restricted(0.0, 1.0, 0.0, cut) : top);
    enum { kBottom, kX0, kB, kTop, kEnd0, kEndL };

    // The spine: the ball inside the material, opposite the outward normals.
    const geom::OffsetCarrier oa(top, -r), ob(B, -r);
    const solve::SsiResult sp = solve::intersect_carriers(oa, ob, opt);
    if (!sp.complete || sp.curves.size() != 1) {
        res.error = "expected exactly one spine with a completeness proof";
        return res;
    }
    auto blend = std::make_shared<solve::BlendCarrier>(top, B, r, sp.curves[0].segs, -1);
    res.blend = blend;
    solve::Pt4 x0, x1;
    if (!blend->spine_at(0.0, x0) || !blend->spine_at(1.0, x1)) { res.error = "spine ends"; return res; }
    // Orient sigma from y = 0 to y = L (top coordinate v = x[1]).
    const bool sigma_up = x0[1] < x1[1];
    const solve::Pt4 e0 = sigma_up ? x0 : x1, eL = sigma_up ? x1 : x0;
    if (!(std::fabs(e0[1]) < 1e-9 && std::fabs(e0[2]) < 1e-9)) {
        res.error = "the spine does not end on the end plane y = 0";
        return res;
    }
    if (!trimmed && !(std::fabs(eL[1] - 1.0) < 1e-9 && std::fabs(eL[2] - 1.0) < 1e-9)) {
        res.error = "the spine does not end on the end plane y = L";
        return res;
    }

    Brep b;
    b.carriers = block.faces;
    b.procs.push_back(blend);
    const double s0 = sigma_up ? 0.0 : 1.0;
    // The far end: sigma of the contacts there (A on the top, B on the side),
    // and, for a cut, the trimming curve.
    double sA = 1.0 - s0, sB = 1.0 - s0;
    std::shared_ptr<TracedEdge> end_curve;
    if (trimmed) {
        auto plane = std::make_shared<geom::PatchCarrier>(b.carriers[kEndL]);
        const solve::SsiResult cr = solve::intersect_carriers(*blend, *plane, opt);
        if (!cr.complete || cr.curves.size() != 1) {
            res.error = "expected exactly one curve where the blend meets the end face, with a completeness proof";
            return res;
        }
        end_curve = std::make_shared<TracedEdge>();
        end_curve->a = blend;
        end_curve->b = plane;
        end_curve->proc[0] = 0;
        end_curve->carrier[1] = kEndL;
        end_curve->path.segs = cr.curves[0].segs;
        end_curve->path.make_breaks(end_curve->sys(), *blend);
        // Run it from theta = 0 (on the top) to theta = 1 (on the side).
        solve::Pt4 c0, c1;
        end_curve->from = 0.0;
        end_curve->to = 1.0;
        if (!end_curve->at(0.0, c0) || !end_curve->at(1.0, c1)) { res.error = "trimming curve ends"; return res; }
        if (c0[1] > c1[1]) {
            std::swap(c0, c1);
            end_curve->from = 1.0;
            end_curve->to = 0.0;
        }
        if (!(std::fabs(c0[1]) < 1e-9 && std::fabs(c1[1] - 1.0) < 1e-9)) {
            res.error = "the trimming curve does not run from contact to contact";
            return res;
        }
        sA = c0[0];
        sB = c1[0];
    }

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
    blend->eval(s0, 0, P, Ps, Pt);
    const int A0 = vertex(P);
    blend->eval(sA, 0, P, Ps, Pt);
    const int AL = vertex(P);
    blend->eval(s0, 1, P, Ps, Pt);
    const int B0 = vertex(P);
    blend->eval(sB, 1, P, Ps, Pt);
    const int BL = vertex(P);
    solve::Pt4 xA, xB;
    blend->spine_at(sA, xA);
    blend->spine_at(sB, xB);
    const double uA0 = e0[0], uAL = xA[0], vB0 = e0[3], vBL = xB[3];
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
    auto contact = [&](int v0, int v1, double theta, double from, double to, int support, int coord, double sc0, double sc1) {
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
    // The supports' coordinates on the restricted faces: v / cut on the top,
    // s / cut on the side.
    const int EcA = contact(A0, AL, 0.0, s0, sA, kTop, 0, 1.0, cut);
    const int EcB = contact(B0, BL, 1.0, s0, sB, kB, 2, cut, 1.0);
    const int Es0 = section(A0, B0, s0, kEnd0);
    const int EsL = trimmed ? add_edge(AL, BL, end_curve) : section(AL, BL, sA, kEndL);

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
        const double sm = 0.5 * (s0 + 0.5 * (sA + sB));
        Vec3d Q, Qs, Qt;
        blend->eval(sm, 0.5, Q, Qs, Qt);
        solve::Pt4 x;
        blend->spine_at(sm, x);
        Vec3d c, cu, cv;
        geom::OffsetCarrier(top, -r).eval(x[0], x[1], c, cu, cv);
        f.flipped = dot(cross(Qs, Qt), Q - c) < 0.0;
        b.faces.push_back(std::move(f));
    }
    res.solid = std::move(b);
    return res;
}

} // namespace ak::op
