// akernel — boundary representation with trimmed faces.
//
// One type for the operands and the results of Boolean operations, so that
// results can be operands.
//
//   carrier  a rational Bezier patch
//   face     a carrier, a flag saying whether the outward normal is Su x Sv
//            or its opposite, and loops in the carrier's parameter square
//            (the face is on the left of every loop)
//   edge     one definition, used by two faces. Either
//              - a side edge: part of a boundary curve shared by two
//                carriers, given by a certified parameter range, or
//              - an arc edge: part of the certified intersection curve of
//                two carriers
//   vertex   a point with certified parameters ("charts") on every carrier
//            it lies on
//
// No entity stores a tolerance. What is stored instead are the certificates:
// parameter enclosures of vertices and tube segments of arc edges.
#pragma once

#include <algorithm>
#include <cmath>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <akernel/num/detmath.hpp>
#include <akernel/geom/carrier.hpp>
#include <akernel/op/contact.hpp>
#include <akernel/op/patch_solid.hpp>
#include <akernel/op/planar.hpp>
#include <akernel/op/surfaces.hpp>
#include <akernel/op/trim.hpp>
#include <akernel/topo/complex.hpp>

namespace ak::op {

using geom::Vec3;

struct Chart {
    int carrier = -1;
    Interval u, v;
};

struct BVertex {
    Vec3d point{}; ///< approximate position
    std::vector<Chart> charts;
    /// The vertex exactly, where an operation proved it: homogeneous
    /// coordinates that are doubles (a point with round coordinates, such as
    /// (0.75, 1, 0.5)). Its parameters on a surface may still be enclosures:
    /// they need not be dyadic where the point is.
    bool exact = false;
    geom::Vec4d xp{};
    const Chart* chart(int carrier) const {
        for (const Chart& c : charts)
            if (c.carrier == carrier) return &c;
        return nullptr;
    }
};

/// An edge defined procedurally: a side of a procedural face's parameter
/// square, read by every face that uses it in that face's own parameters
/// (docs/FEATURES.md §5.1). tau in [0, 1] runs from the edge's v0 to v1.
struct ProcEdge {
    virtual ~ProcEdge() = default;
    virtual Vec3d point(double tau) const = 0;
    /// Parameters (u, v) on Brep::carriers[carrier] at tau, and their
    /// derivatives in tau; false if the edge is not on that carrier.
    virtual bool on_carrier(int carrier, double tau, double& u, double& v, double& du, double& dv) const = 0;
    /// Parameters on Brep::procs[proc]; false if the edge is not on it.
    virtual bool on_proc(int proc, double tau, double& s, double& t) const = 0;
    /// Break points 0 = b0 < ... < bn = 1 of tau: the edge is smooth on each
    /// piece (e.g. a spine segment). Quadrature and sampling align with them.
    virtual std::vector<double> breaks() const { return {0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0}; }
    /// The part tau in [t0, t1] of this edge as an edge of its own (tau
    /// rescaled to [0, 1]); null if the edge cannot be split.
    virtual std::shared_ptr<const ProcEdge> part(double t0, double t1) const {
        (void)t0;
        (void)t1;
        return nullptr;
    }
    /// Parameters on Brep::procs[proc] and their derivatives in tau. By
    /// default from central differences inside tau's piece of breaks();
    /// edges that know their tangent override it (differences lose digits
    /// on short pieces).
    virtual bool on_proc_jet(int proc, double tau, double& s, double& t, double& ds, double& dt) const;
};

/// Central differences of a curve's parameters in tau, inside tau's piece of
/// `breaks` (one-sided at its ends): the curve may have a kink between pieces.
template <class F>
bool piece_derivative(const std::vector<double>& breaks, double tau, const F& uv, double& du, double& dv) {
    const auto it = std::upper_bound(breaks.begin() + 1, breaks.end() - 1, tau);
    const double lo = *(it - 1), hi = *it;
    const double h = std::min(1e-6, 0.25 * (hi - lo)), ta = std::max(lo, tau - h), tb = std::min(hi, tau + h);
    double ua, va, ub, vb;
    if (!uv(ta, ua, va) || !uv(tb, ub, vb)) return false;
    du = (ub - ua) / (tb - ta);
    dv = (vb - va) / (tb - ta);
    return true;
}

inline bool ProcEdge::on_proc_jet(int proc, double tau, double& s, double& t, double& ds, double& dt) const {
    if (!on_proc(proc, tau, s, t)) return false;
    auto st = [&](double x, double& a, double& b) { return on_proc(proc, x, a, b); };
    return piece_derivative(breaks(), tau, st, ds, dt);
}

struct BEdge {
    int v0 = -1, v1 = -1;
    bool arc = false;
    std::shared_ptr<const ProcEdge> proc; ///< set for a procedural edge (then arc is false)

    // Side edge: the curve is side `on[0].side` of carrier `on[0].carrier`,
    // for the edge parameter t in [t0, t1]; the same curve is a side of
    // `on[1].carrier`, possibly with reversed parameter.
    struct SideRef {
        int carrier = -1;
        int side = 0;
        bool reversed = false; ///< side parameter = 1 - edge parameter
    };
    SideRef on[2];
    Interval t0{0.0}, t1{1.0};

    // Arc edge: part of the intersection curve of carrier[0] (curve
    // coordinates 0, 1) and carrier[1] (coordinates 2, 3).
    int carrier[2] = {-1, -1};
    std::shared_ptr<const ArcGeom> geom;
};

struct BCoedge {
    int edge = -1;
    bool reversed = false;
};

struct BFace {
    int carrier = -1;
    int proc = -1;        ///< >= 0: the face is Brep::procs[proc] over its whole parameter square
    bool flipped = false; ///< outward normal is -(Su x Sv)
    std::vector<std::vector<BCoedge>> loops;
};

struct Brep {
    std::vector<PolyPatchD> carriers;
    /// Procedural carriers (offsets, blends): faces on them have proc >= 0.
    /// Booleans do not read them yet and refuse such operands.
    std::vector<std::shared_ptr<const geom::Carrier>> procs;
    std::vector<BVertex> vertices;
    std::vector<BEdge> edges;
    std::vector<BFace> faces;
    /// Set by a Boolean operation on its result. Its faces are regions of
    /// operands that were checked, kept with their orientation or flipped
    /// as the operation requires, so the outward normal is known and the
    /// inside-out check (a volume quadrature) can be skipped. Clear it after
    /// editing a Brep by hand.
    bool outward_known = false;
};

/// A patch solid as a Brep: every face is a whole carrier.
inline Brep to_brep(const PatchSolid& s) {
    Brep b;
    b.carriers = s.faces;
    b.vertices.resize(s.vertex_points.size());
    for (std::size_t i = 0; i < s.vertex_points.size(); ++i) b.vertices[i].point = s.vertex_points[i];
    static constexpr double corner_uv[4][2] = {{0.0, 0.0}, {1.0, 0.0}, {1.0, 1.0}, {0.0, 1.0}};
    for (std::size_t f = 0; f < s.faces.size(); ++f)
        for (int c = 0; c < 4; ++c)
            b.vertices[static_cast<std::size_t>(s.corners[f][static_cast<std::size_t>(c)])].charts.push_back(
                {static_cast<int>(f), Interval(corner_uv[c][0]), Interval(corner_uv[c][1])});
    for (const PatchSolid::ModelEdge& me : s.edges) {
        BEdge e;
        e.v0 = me.v0;
        e.v1 = me.v1;
        for (int k = 0; k < 2; ++k) e.on[k] = {me.face[k], me.side[k], me.rev[k]};
        b.edges.push_back(e);
    }
    for (std::size_t f = 0; f < s.faces.size(); ++f) {
        BFace face;
        face.carrier = static_cast<int>(f);
        std::vector<BCoedge> loop;
        for (int side = 0; side < 4; ++side) {
            const PatchSolid::SideRef ref = s.sides[f][static_cast<std::size_t>(side)];
            loop.push_back({ref.edge, side_ccw_forward(side) == ref.reversed});
        }
        face.loops.push_back(loop);
        b.faces.push_back(face);
    }
    return b;
}

/// How a side edge runs on one of its two carriers.
inline const BEdge::SideRef* side_on(const BEdge& e, int carrier) {
    if (e.on[0].carrier == carrier) return &e.on[0];
    if (e.on[1].carrier == carrier) return &e.on[1];
    return nullptr;
}
/// Curve coordinate offset (0 or 2) of an arc edge on one of its carriers, or -1.
inline int coord_on(const BEdge& e, int carrier) {
    if (e.carrier[0] == carrier) return 0;
    if (e.carrier[1] == carrier) return 2;
    return -1;
}

/// How to read an arc edge in the parameters of a face's carrier: the curve
/// coordinate offset, and a map if the edge lies on a different carrier of
/// the same surface (the same plane, or an affine reparameterisation; exact
/// tests). False if the edge does not lie on it.
inline bool arc_on(const Brep& b, const BEdge& e, int carrier, int& co, bool& mapped, Affine& map) {
    co = coord_on(e, carrier);
    mapped = false;
    if (co >= 0) return true;
    const PolyPatchD& F = b.carriers[static_cast<std::size_t>(carrier)];
    for (int k = 0; k < 2; ++k) {
        const PolyPatchD& C = b.carriers[static_cast<std::size_t>(e.carrier[k])];
        if (same_surface_map(C, F, map)) {
            co = 2 * k;
            mapped = true;
            return true;
        }
    }
    return false;
}

/// Chord-and-box geometry of the loops of a face, for point classification.
struct FaceGeom {
    std::vector<CycleGeom> loops;
};

inline bool face_geometry(const Brep& b, int f, FaceGeom& out, std::string& error) {
    const BFace& face = b.faces[static_cast<std::size_t>(f)];
    out.loops.clear();
    for (const auto& loop : face.loops) {
        CycleGeom cyc;
        for (const BCoedge& ce : loop) {
            const BEdge& e = b.edges[static_cast<std::size_t>(ce.edge)];
            const Chart* c0 = b.vertices[static_cast<std::size_t>(e.v0)].chart(face.carrier);
            const Chart* c1 = b.vertices[static_cast<std::size_t>(e.v1)].chart(face.carrier);
            if (!c0 || !c1) { error = "a vertex has no parameters on a face it bounds"; return false; }
            const double p0[2] = {num::mid(c0->u), num::mid(c0->v)}, p1[2] = {num::mid(c1->u), num::mid(c1->v)};
            const Interval b0[2] = {c0->u, c0->v}, b1[2] = {c1->u, c1->v};
            if (!e.arc) {
                if (ce.reversed) append_straight_link(cyc.links, p1, b1, p0, b0);
                else append_straight_link(cyc.links, p0, b0, p1, b1);
            } else {
                int co = -1;
                bool mapped = false;
                Affine map;
                if (!arc_on(b, e, face.carrier, co, mapped, map)) {
                    error = "an arc edge does not lie on the carrier of a face it bounds";
                    return false;
                }
                append_arc_links(cyc.links, *e.geom, co, !ce.reversed, p0, p1, mapped ? &map : nullptr);
            }
        }
        out.loops.push_back(std::move(cyc));
    }
    return true;
}

enum class Where { Inside, Outside, Boundary, Undecided };

/// Certified classification of a parameter box against a trimmed face.
/// Refines the loop geometry as a side effect.
inline Where locate(std::vector<CycleGeom>& loops, const Interval& u, const Interval& v, int max_depth = 40) {
    int total = 0;
    for (CycleGeom& cyc : loops) {
        const WindingResult w = winding(cyc, u, v, max_depth);
        if (!w.decided) return w.on_boundary ? Where::Boundary : Where::Undecided;
        total += w.wn;
    }
    return total == 1 ? Where::Inside : (total == 0 ? Where::Outside : Where::Undecided);
}
inline Where locate(FaceGeom& g, const Interval& u, const Interval& v, int max_depth = 40) {
    return locate(g.loops, u, v, max_depth);
}

/// The incidence structure of a Brep, for the integer audit.
inline topo::Complex to_complex(const Brep& b) {
    topo::Complex c;
    for (std::size_t i = 0; i < b.vertices.size(); ++i) c.add_vertex(static_cast<int>(i));
    for (std::size_t i = 0; i < b.edges.size(); ++i) c.add_edge(b.edges[i].v0, b.edges[i].v1, static_cast<int>(i));
    for (const BFace& f : b.faces) {
        std::vector<topo::Loop> loops;
        for (const auto& loop : f.loops) {
            topo::Loop L;
            for (const BCoedge& ce : loop) L.push_back({ce.edge, ce.reversed});
            if (f.flipped) {
                std::reverse(L.begin(), L.end());
                for (topo::Coedge& ce : L) ce.reversed = !ce.reversed;
            }
            loops.push_back(std::move(L));
        }
        c.add_face(f.proc >= 0 ? static_cast<int>(b.carriers.size()) + f.proc : f.carrier, std::move(loops));
    }
    return c;
}

inline topo::AuditReport audit(const Brep& b) { return to_complex(b).audit(); }

/// Volume enclosed by the Brep: (1/3) * sum over faces of the integral of
/// S . n over the trimmed face, turned into loop integrals by Green's
/// theorem. Numerical quadrature in doubles: a measurement, not a bound.
/// `quick`: on polynomial faces, Gauss rules just large enough to be exact
/// for the degree (two and three points on a flat bilinear face instead of
/// 16 and 12); used where only the sign matters.
inline double volume(const Brep& b, bool quick = false) {
    const UnitRule inner16(16), outer12(12);
    std::vector<std::unique_ptr<UnitRule>> rules(17);
    auto rule = [&](int n) -> const UnitRule& {
        if (!rules[static_cast<std::size_t>(n)]) rules[static_cast<std::size_t>(n)] = std::make_unique<UnitRule>(n);
        return *rules[static_cast<std::size_t>(n)];
    };
    double total = 0.0;
    for (const BFace& face : b.faces) {
        if (face.proc >= 0) {
            // Green's theorem in the carrier's (s, t): the flux of the trimmed
            // face is the loop integral of F ds, F(s, t) = -integral from 0 to
            // t of S . (Ss x St) (Gauss 8 on 2 cells in t). Edges are taken
            // piece by piece of their breaks, 4 cells each.
            const geom::Carrier& C = *b.procs[static_cast<std::size_t>(face.proc)];
            const UnitRule& g = rule(8);
            auto F = [&](double s, double t) {
                double acc = 0.0;
                for (int c = 0; c < 2; ++c)
                    for (std::size_t j = 0; j < g.x.size(); ++j) {
                        const double tt = t * (c + g.x[j]) / 2;
                        Vec3d P, Ps, Pt;
                        C.eval(s, tt, P, Ps, Pt);
                        acc += g.w[j] * dot(P, cross(Ps, Pt)) * t / 2;
                    }
                return -acc;
            };
            double flux = 0.0;
            for (const auto& loop : face.loops)
                for (const BCoedge& ce : loop) {
                    const BEdge& e = b.edges[static_cast<std::size_t>(ce.edge)];
                    if (!e.proc) continue; // a procedural face is bounded by procedural edges
                    const std::vector<double> bt = e.proc->breaks();
                    for (std::size_t piece = 0; piece + 1 < bt.size(); ++piece)
                        for (int k = 0; k < 4; ++k) {
                            const double dt = (bt[piece + 1] - bt[piece]) / 4, t0 = bt[piece] + k * dt;
                            for (std::size_t i = 0; i < outer12.x.size(); ++i) {
                                const double tau = t0 + outer12.x[i] * dt;
                                double s, t, ds, dtt;
                                if (!e.proc->on_proc_jet(face.proc, tau, s, t, ds, dtt)) continue;
                                if (ds == 0.0) continue;
                                flux += outer12.w[i] * F(s, t) * (ce.reversed ? -ds : ds) * dt;
                            }
                        }
                }
            total += face.flipped ? -flux : flux;
            continue;
        }
        const PolyPatchD& P = b.carriers[static_cast<std::size_t>(face.carrier)];
        const geom::BezierSurface S(P);
        // The density S . (Su x Sv) of a polynomial patch has degree 3p - 1 in
        // u and 3q - 1 in v; G below then has degree 3p + 3q - 1 along a
        // straight segment of the parameter plane.
        bool polynomial = quick;
        for (const Vec4d& c : P.c) polynomial = polynomial && c.w == P.c.front().w;
        const int n_in = (3 * P.p + 1) / 2, n_out = (3 * P.p + 3 * P.q + 1) / 2;
        const bool small = polynomial && n_in <= 16 && n_out <= 16;
        const UnitRule& inner = small ? rule(n_in) : inner16;
        const UnitRule& straight = small ? rule(n_out) : outer12;
        const UnitRule& outer = outer12;
        // G(u, v) = integral of the density from 0 to u; the area integral of
        // the density equals the loop integral of G dv.
        auto G = [&](double u, double v) {
            double s = 0.0;
            for (std::size_t i = 0; i < inner.x.size(); ++i) s += inner.w[i] * volume_density(S, u * inner.x[i], v);
            return s * u;
        };
        double flux = 0.0;
        for (const auto& loop : face.loops)
            for (const BCoedge& ce : loop) {
                const BEdge& e = b.edges[static_cast<std::size_t>(ce.edge)];
                if (e.proc) {
                    // The loop integral of G dv along the edge in this face's parameters.
                    const std::vector<double> bt = e.proc->breaks();
                    for (std::size_t piece = 0; piece + 1 < bt.size(); ++piece)
                        for (int k = 0; k < 4; ++k) {
                            const double dt = (bt[piece + 1] - bt[piece]) / 4, t0 = bt[piece] + k * dt;
                            // Along the edge's own tau (nodes aligned with its
                            // breaks); a reversed coedge changes the sign.
                            for (std::size_t i = 0; i < outer.x.size(); ++i) {
                                const double tau = t0 + outer.x[i] * dt;
                                double u, v, du, dv;
                                if (!e.proc->on_carrier(face.carrier, tau, u, v, du, dv)) continue;
                                flux += outer.w[i] * G(u, v) * (ce.reversed ? -dv : dv) * dt;
                            }
                        }
                    continue;
                }
                if (!e.arc) {
                    double ua, va, ub, vb;
                    if (const BEdge::SideRef* ref = side_on(e, face.carrier)) {
                        const double ta = num::mid(ce.reversed ? e.t1 : e.t0), tb = num::mid(ce.reversed ? e.t0 : e.t1);
                        side_uv(ref->side, ref->reversed ? 1.0 - ta : ta, ua, va);
                        side_uv(ref->side, ref->reversed ? 1.0 - tb : tb, ub, vb);
                    } else {
                        // A straight edge that is not a side of this carrier (it lies in
                        // its plane): a segment between the vertices' parameters.
                        const Chart* a = b.vertices[static_cast<std::size_t>(ce.reversed ? e.v1 : e.v0)].chart(face.carrier);
                        const Chart* c = b.vertices[static_cast<std::size_t>(ce.reversed ? e.v0 : e.v1)].chart(face.carrier);
                        if (!a || !c) continue;
                        ua = num::mid(a->u); va = num::mid(a->v);
                        ub = num::mid(c->u); vb = num::mid(c->v);
                    }
                    if (va == vb) continue;
                    for (std::size_t i = 0; i < straight.x.size(); ++i)
                        flux += straight.w[i] * (vb - va) * G(ua + straight.x[i] * (ub - ua), va + straight.x[i] * (vb - va));
                    continue;
                }
                int coord = -1;
                bool mapped = false;
                Affine map;
                if (!arc_on(b, e, face.carrier, coord, mapped, map)) continue;
                const std::size_t co = static_cast<std::size_t>(coord);
                const SsiSystem& sys = *e.geom->base->sys;
                for (const Span& sp : e.geom->spans) {
                    const TubeSegment& seg = e.geom->base->curve->segs[static_cast<std::size_t>(sp.seg)];
                    const double c0 = ce.reversed ? sp.cb : sp.ca, c1 = ce.reversed ? sp.ca : sp.cb;
                    for (std::size_t i = 0; i < outer.x.size(); ++i) {
                        const double c = c0 + outer.x[i] * (c1 - c0);
                        Pt4 x, t;
                        for (int k = 0; k < 4; ++k) x[k] = seg.a[k] + (c / seg.len) * (seg.b[k] - seg.a[k]);
                        solve::newton_slice(sys, seg.m, seg.a, c, x);
                        if (!sys.tangent(x, t)) continue;
                        double u = x[co], v = x[co + 1], dv = t[co + 1];
                        if (mapped) {
                            Interval mu, mv, du, dw;
                            map.apply(Interval(x[co]), Interval(x[co + 1]), mu, mv);
                            map.apply_linear(Interval(t[co]), Interval(t[co + 1]), du, dw);
                            u = num::mid(mu);
                            v = num::mid(mv);
                            dv = num::mid(dw);
                        }
                        flux += outer.w[i] * (c1 - c0) * G(u, v) * (dv / solve::dot4(seg.m, t));
                    }
                }
            }
        total += face.flipped ? -flux : flux;
    }
    return total / 3.0;
}

/// Largest certified distance between the two carriers of any arc edge,
/// taken at one certified point per span.
inline double edge_gap_bound(const Brep& b) {
    double worst = 0.0;
    for (const BEdge& e : b.edges) {
        if (!e.arc) continue;
        const SsiSystem& sys = *e.geom->base->sys;
        for (const Span& sp : e.geom->spans) {
            const TubeSegment& seg = e.geom->base->curve->segs[static_cast<std::size_t>(sp.seg)];
            const double c = 0.5 * (sp.ca + sp.cb);
            Pt4 x;
            for (int i = 0; i < 4; ++i) x[i] = seg.a[i] + (c / seg.len) * (seg.b[i] - seg.a[i]);
            Box4 enc;
            if (!solve::certify_point(sys, seg.m, seg.a, c, x, enc)) continue;
            const auto pa = sys.point_a(enc), pb = sys.point_b(enc);
            const Interval d2 = num::sqr(pa.x - pb.x) + num::sqr(pa.y - pb.y) + num::sqr(pa.z - pb.z);
            worst = std::max(worst, num::sqrt(d2).hi);
        }
    }
    return worst;
}

/// Geometry of the faces around every edge used by more than two faces (a
/// solid whose parts touch along an edge). At a point of the edge, each face
/// that uses it has an interior direction perpendicular to the edge:
/// (Su x Sv) x T, with T the edge's tangent in the coedge's direction (the
/// face lies on the left of its loops). Sorted by angle around the edge,
/// consecutive faces must bound material and empty space alternately, so
/// the effective orientations of the uses (coedge direction times the face's
/// flip) must alternate. Faces at the same angle touch tangentially along the
/// edge: a curved face and a flat one are ordered by the side of the plane on
/// which the curved face lies (proven exactly, op/contact.hpp); any other
/// such pair is refused. Doubles with a margin of 1e-9 rad: a check that can
/// refuse a valid solid, not one that decides a topology.
inline bool check_edge_fans(const Brep& b, std::string& why) {
    std::vector<std::vector<std::pair<int, bool>>> uses(b.edges.size()); // (face, coedge reversed)
    for (std::size_t f = 0; f < b.faces.size(); ++f)
        for (const auto& loop : b.faces[f].loops)
            for (const BCoedge& ce : loop) uses[static_cast<std::size_t>(ce.edge)].push_back({static_cast<int>(f), ce.reversed});
    const std::size_t iu = static_cast<std::size_t>(geom::jet_index(1, 0)), iv = static_cast<std::size_t>(geom::jet_index(0, 1));
    for (std::size_t ei = 0; ei < b.edges.size(); ++ei) {
        const auto& U = uses[ei];
        if (U.size() <= 2) continue;
        const BEdge& e = b.edges[ei];
        // Tangent line of the edge at one point, and that point's parameters
        // on the carriers that contain the edge.
        Vec3d T{};
        std::map<int, std::pair<double, double>> uv_on;
        if (!e.arc) {
            const double t = 0.5 * (num::mid(e.t0) + num::mid(e.t1));
            for (int k = 0; k < 2; ++k) {
                double u, v;
                side_uv(e.on[k].side, e.on[k].reversed ? 1.0 - t : t, u, v);
                uv_on[e.on[k].carrier] = {u, v};
            }
            const auto& p = uv_on[e.on[0].carrier];
            const geom::Jet3 J = geom::BezierSurface(b.carriers[static_cast<std::size_t>(e.on[0].carrier)]).jet(p.first, p.second);
            T = (e.on[0].side == 0 || e.on[0].side == 2) ? J[iu] : J[iv];
        } else {
            Box4 mid;
            if (!arc_point(*e.geom, mid)) { why = "no point could be placed on an edge with more than two faces"; return false; }
            Pt4 x{num::mid(mid[0]), num::mid(mid[1]), num::mid(mid[2]), num::mid(mid[3])}, t4;
            if (!e.geom->base->sys->tangent(x, t4)) { why = "no tangent on an edge with more than two faces"; return false; }
            uv_on[e.carrier[0]] = {x[0], x[1]};
            uv_on[e.carrier[1]] = {x[2], x[3]};
            const geom::Jet3 J = geom::BezierSurface(b.carriers[static_cast<std::size_t>(e.carrier[0])]).jet(x[0], x[1]);
            T = J[iu] * t4[0] + J[iv] * t4[1];
        }
        const double tl = norm(T);
        if (!(tl > 0.0)) { why = "an edge with more than two faces has no tangent"; return false; }
        T = T / tl;
        Vec3d A = std::fabs(T.x) < 0.9 ? cross(T, Vec3d{1, 0, 0}) : cross(T, Vec3d{0, 1, 0});
        A = A / norm(A);
        const Vec3d B = cross(T, A);
        struct Use {
            double angle;
            int orient;
            int carrier = -1;
            int tie = 0; ///< second-order offset at an equal angle: -1, 0 or +1
        };
        std::vector<Use> around;
        for (const auto& [f, rev] : U) {
            const BFace& face = b.faces[static_cast<std::size_t>(f)];
            // A face on a carrier that does not contain the edge lies on another
            // carrier of the same surface: map the point's parameters (for a
            // plane any point will do; its normal is the same everywhere).
            auto it = uv_on.find(face.carrier);
            double u = 0.5, v = 0.5;
            if (it != uv_on.end()) {
                u = it->second.first;
                v = it->second.second;
            } else {
                for (const auto& [c, p] : uv_on) {
                    Affine m;
                    if (!same_surface_map(b.carriers[static_cast<std::size_t>(c)], b.carriers[static_cast<std::size_t>(face.carrier)], m))
                        continue;
                    Interval mu, mv;
                    m.apply(Interval(p.first), Interval(p.second), mu, mv);
                    u = num::mid(mu);
                    v = num::mid(mv);
                    break;
                }
            }
            const geom::Jet3 J = geom::BezierSurface(b.carriers[static_cast<std::size_t>(face.carrier)]).jet(u, v);
            const Vec3d n = cross(J[iu], J[iv]);
            const Vec3d d = cross(n, rev ? -T : T);
            around.push_back({num::det::atan2(dot(d, B), dot(d, A)), (rev ? -1 : 1) * (face.flipped ? -1 : 1), face.carrier, 0});
        }
        // Faces at one angle touch tangentially along the edge. A curved face
        // whose side along this edge lies in the plane of a flat face there,
        // and that meets that plane nowhere else (decided exactly,
        // op/contact.hpp), lies strictly on one side of it: it leaves the edge
        // turned that way to second order, which orders the two.
        if (!e.arc)
            for (Use& cu : around)
                for (const Use& pu : around) {
                    if (&cu == &pu || std::fabs(cu.angle - pu.angle) >= 1e-9 || cu.tie != 0) continue;
                    const PolyPatchD& P = b.carriers[static_cast<std::size_t>(pu.carrier)];
                    const PolyPatchD& C = b.carriers[static_cast<std::size_t>(cu.carrier)];
                    if (!planar_affine(P) || planar_affine(C)) continue;
                    int side = -1;
                    for (int k = 0; k < 2; ++k)
                        if (e.on[k].carrier == cu.carrier) side = e.on[k].side;
                    if (side < 0 || !side_in_plane(C, side, P) || !side_contact(C, side, P).only_side) continue;
                    // Which side of the plane: at the patch's centre, which is off it.
                    const Vec4d c00 = corner_point(P, 0, 0), c10 = corner_point(P, 1, 0), c01 = corner_point(P, 0, 1);
                    const Vec3d nP = cross(Vec3d{c10.x - c00.x, c10.y - c00.y, c10.z - c00.z}, Vec3d{c01.x - c00.x, c01.y - c00.y, c01.z - c00.z});
                    const Vec4d h = C.eval(0.5, 0.5);
                    const double off = dot(nP, Vec3d{h.x / h.w - c00.x, h.y / h.w - c00.y, h.z / h.w - c00.z});
                    // Direction of increasing angle at this angle.
                    const Vec3d up = A * -num::det::sin(cu.angle) + B * num::det::cos(cu.angle);
                    const double turn = dot(nP, up) * (off > 0.0 ? 1.0 : -1.0);
                    if (off != 0.0 && turn != 0.0) cu.tie = turn > 0.0 ? 1 : -1;
                }
        std::sort(around.begin(), around.end(), [](const Use& a, const Use& c) { return a.angle < c.angle; });
        // Within a run of equal angles, order by the second-order offset.
        for (std::size_t i = 0; i + 1 < around.size(); ++i)
            for (std::size_t j = i + 1; j < around.size() && around[j].angle - around[i].angle < 1e-9; ++j)
                if (around[j].tie < around[i].tie) std::swap(around[i], around[j]);
        const std::size_t m = around.size();
        for (std::size_t i = 0; i < m; ++i) {
            const Use& a = around[i];
            const Use& c = around[(i + 1) % m];
            double gap = c.angle - a.angle;
            if (i + 1 == m) gap += 2.0 * 3.141592653589793;
            if (gap < 1e-9 && a.tie == c.tie) { why = "two faces meet tangentially along an edge with more than two faces"; return false; }
            if (a.orient == c.orient) {
                why = "the faces around an edge with more than two faces do not alternate between material and empty space";
                return false;
            }
        }
    }
    return true;
}

/// Widest certified parameter enclosure of any vertex.
inline double vertex_width(const Brep& b) {
    double worst = 0.0;
    for (const BVertex& v : b.vertices)
        for (const Chart& c : v.charts) worst = std::max(worst, std::max(num::width(c.u), num::width(c.v)));
    return worst;
}

} // namespace ak::op
