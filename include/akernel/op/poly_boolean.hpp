// akernel — Booleans between a solid that may have procedural faces and a
// solid whose faces are all rational patches, flat or curved
// (docs/FEATURES.md §5.1; decided 2026-11-24: cuts by planar solids first,
// curved faces after).
//
// The rational Boolean engine (op/boolean.hpp) reads polynomial carriers
// only. This one works on any surfaces of the first operand A (polynomial
// or procedural, through geom::Carrier) and needs the second operand P to
// have polynomial carriers, so that a ray against P is a certified
// line-surface intersection.
//
// 1. Every pair (face of A, face of P) whose enclosures meet is intersected
//    by intersect_carriers (certified, complete), over the whole surfaces.
// 2. New vertices: where a curve leaves one of its two faces it crosses an
//    edge of that face (an edge of A crossing a face of P, or an edge of P
//    piercing a face of A). The curves are complete, so walking them finds
//    every such point; each is refined onto its edge by Gauss-Newton.
// 3. The pieces of each curve inside both faces are traced edges
//    (op/proc_edge.hpp). Each face is split by the pieces lying on it: an
//    arrangement in its parameters, traversed with the outgoing edges at
//    each vertex sorted by angle. Cycles running the face's way round are
//    regions; cycles running the other way are holes in the region around
//    them.
// 4. Regions next to a curve are inside or outside the other solid by the
//    side its outward normal points to; the others inherit it across the
//    uncut edges they share, and a part of A no curve touches takes a ray
//    test against P (certified line-surface roots).
// 5. The result keeps A's regions outside P (union, difference) or inside
//    it (intersection), and P's regions outside A (union), inside A
//    (intersection) or inside A reversed (difference).
//
// Certified: the intersection curves and their completeness, and the ray
// roots. Measured in doubles, with refusals near a decision: where a curve
// leaves a face (winding on sampled loops, then refined), the order of edges
// around vertices. Refused with a message for now: arc edges, faces with
// inner loops in the operands, faces of the two solids on one surface,
// curves through vertices, tangent contact, curves closed inside both
// faces, and parts of the second solid that never meet the first's
// boundary.
#pragma once

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include <akernel/op/brep.hpp>
#include <akernel/op/brep_transform.hpp>
#include <akernel/op/plane_cut.hpp>
#include <akernel/op/proc_edge.hpp>
#include <akernel/solve/carrier_ssi.hpp>
#include <akernel/solve/csi.hpp>

namespace ak::op {

enum class PolyOp { Union, Intersection, Difference };

struct PolyBooleanResult {
    std::optional<Brep> solid;
    std::string error;
    int curves = 0; ///< pieces of intersection curve (new edges)
    double seconds = 0.0;
};

namespace poly_boolean_detail {

using plane_cut_detail::signed_area;
using plane_cut_detail::winding;

} // namespace poly_boolean_detail

inline PolyBooleanResult poly_boolean(const Brep& A, const Brep& P, PolyOp op, const solve::SsiOptions& opt = {}) {
    using namespace poly_boolean_detail;
    const auto t_start = std::chrono::steady_clock::now();
    PolyBooleanResult res;
    auto fail = [&](std::string msg) {
        res.error = std::move(msg);
        res.solid.reset();
        return res;
    };
    const double tol = 1e-9;
    for (const Brep* b : {&A, &P})
        for (const BEdge& e : b->edges)
            if (e.arc) return fail("arc edges are not supported by the polyhedral Boolean yet");

    // The merged model: A's carriers then P's, A's procs; vertices, edges and
    // faces of A then of P, with P's carrier indices shifted.
    Brep M;
    const int ncA = static_cast<int>(A.carriers.size());
    M.carriers = A.carriers;
    M.carriers.insert(M.carriers.end(), P.carriers.begin(), P.carriers.end());
    M.procs = A.procs;
    M.procs.insert(M.procs.end(), P.procs.begin(), P.procs.end());
    const int npA = static_cast<int>(A.procs.size());
    const int nvA = static_cast<int>(A.vertices.size()), neA = static_cast<int>(A.edges.size());
    const int nfA = static_cast<int>(A.faces.size());
    M.vertices = A.vertices;
    for (BVertex v : P.vertices) {
        for (Chart& c : v.charts) c.carrier += ncA;
        M.vertices.push_back(v);
    }
    M.edges = A.edges;
    for (BEdge e : P.edges) {
        e.v0 += nvA;
        e.v1 += nvA;
        for (auto& r : e.on)
            if (r.carrier >= 0) r.carrier += ncA;
        e.proc = reindexed(e.proc, ncA, npA);
        M.edges.push_back(e);
    }
    M.faces = A.faces;
    for (BFace f : P.faces) {
        if (f.proc >= 0) f.proc += npA;
        else f.carrier += ncA;
        for (auto& loop : f.loops)
            for (BCoedge& ce : loop) ce.edge += neA;
        M.faces.push_back(f);
    }
    const int nF = static_cast<int>(M.faces.size());
    auto solid_of_face = [&](int f) { return f < nfA ? 0 : 1; };

    // Surfaces of the faces.
    std::vector<std::shared_ptr<const geom::Carrier>> fsurf(static_cast<std::size_t>(nF));
    for (int f = 0; f < nF; ++f) {
        const BFace& F = M.faces[static_cast<std::size_t>(f)];
        if (F.proc >= 0) fsurf[static_cast<std::size_t>(f)] = M.procs[static_cast<std::size_t>(F.proc)];
        else fsurf[static_cast<std::size_t>(f)] = std::make_shared<geom::PatchCarrier>(M.carriers[static_cast<std::size_t>(F.carrier)]);
    }
    // Edge evaluation, and an edge read in a face's parameters.
    auto edge_point = [&](const BEdge& e, double tau) {
        if (e.proc) return e.proc->point(tau);
        const BEdge::SideRef& ref = e.on[0];
        const double t = num::mid(e.t0) + tau * (num::mid(e.t1) - num::mid(e.t0));
        double u, v;
        side_uv(ref.side, ref.reversed ? 1.0 - t : t, u, v);
        return geom::BezierSurface(M.carriers[static_cast<std::size_t>(ref.carrier)]).point(u, v);
    };
    auto edge_uv = [&](const BFace& face, const BEdge& e, double tau, double& u, double& v, double& du, double& dv) {
        if (face.proc >= 0) {
            if (!e.proc) return false;
            auto st = [&](double t, double& s, double& w) { return e.proc->on_proc(face.proc, t, s, w); };
            return st(tau, u, v) && piece_derivative(e.proc->breaks(), tau, st, du, dv);
        }
        if (e.proc) return e.proc->on_carrier(face.carrier, tau, u, v, du, dv);
        const BEdge::SideRef* ref = side_on(e, face.carrier);
        if (!ref) return false;
        const double a = num::mid(e.t0), b = num::mid(e.t1), t = a + tau * (b - a);
        side_uv(ref->side, ref->reversed ? 1.0 - t : t, u, v);
        const double dt = (b - a) * (ref->reversed ? -1.0 : 1.0);
        du = (ref->side == 0 || ref->side == 2) ? dt : 0.0;
        dv = (ref->side == 1 || ref->side == 3) ? dt : 0.0;
        return true;
    };
    // Loops of a face as polylines in its parameters.
    std::vector<std::vector<std::array<double, 2>>> outline(static_cast<std::size_t>(nF));
    for (int f = 0; f < nF; ++f) {
        const BFace& F = M.faces[static_cast<std::size_t>(f)];
        auto& poly = outline[static_cast<std::size_t>(f)];
        for (const auto& loop : F.loops) {
            // Several loops: joined into one polyline through a slit would be
            // wrong; the winding test takes the sum over loops instead.
            for (const BCoedge& ce : loop) {
                const BEdge& e = M.edges[static_cast<std::size_t>(ce.edge)];
                for (int k = 0; k < 64; ++k) {
                    const double tau = ce.reversed ? 1.0 - k / 64.0 : k / 64.0;
                    double u, v, du, dv;
                    if (!edge_uv(F, e, tau, u, v, du, dv)) return fail("an edge cannot be read in its face's parameters");
                    poly.push_back({u, v});
                }
            }
            if (F.loops.size() > 1) return fail("faces with inner loops are not supported by the polyhedral Boolean yet");
        }
    }
    auto inside_face = [&](int f, double u, double v) { return winding(outline[static_cast<std::size_t>(f)], u, v) != 0; };

    // Faces that use each edge.
    std::vector<std::vector<int>> edge_faces(M.edges.size());
    for (int f = 0; f < nF; ++f)
        for (const auto& loop : M.faces[static_cast<std::size_t>(f)].loops)
            for (const BCoedge& ce : loop) edge_faces[static_cast<std::size_t>(ce.edge)].push_back(f);

    // New vertices with their charts.
    auto add_vertex = [&](const Vec3d& X) {
        BVertex v;
        v.point = X;
        M.vertices.push_back(v);
        return static_cast<int>(M.vertices.size()) - 1;
    };
    auto add_chart = [&](int vx, int carrier, double u, double v) {
        BVertex& V = M.vertices[static_cast<std::size_t>(vx)];
        if (!V.chart(carrier)) V.charts.push_back({carrier, Interval(u), Interval(v)});
    };

    // 1. Curves of (A face, P face) pairs whose enclosures meet.
    auto enclosure = [&](int f, geom::IVec3& box) {
        const BFace& F = M.faces[static_cast<std::size_t>(f)];
        if (F.proc >= 0) {
            geom::IVec3 Pu, Pv;
            return M.procs[static_cast<std::size_t>(F.proc)]->jet1_over(0, 1, 0, 1, box, Pu, Pv);
        }
        const PolyPatchD& C = M.carriers[static_cast<std::size_t>(F.carrier)];
        bool first = true;
        for (int i = 0; i <= C.p; ++i)
            for (int j = 0; j <= C.q; ++j) {
                const geom::Vec4d h = C.at(i, j);
                if (!(h.w > 0.0)) return false;
                const geom::IVec3 c{Interval(h.x) / Interval(h.w), Interval(h.y) / Interval(h.w), Interval(h.z) / Interval(h.w)};
                if (first) { box = c; first = false; }
                else box = {num::hull(box.x, c.x), num::hull(box.y, c.y), num::hull(box.z, c.z)};
            }
        return true;
    };
    std::vector<geom::IVec3> fbox(static_cast<std::size_t>(nF));
    std::vector<char> fbox_ok(static_cast<std::size_t>(nF));
    for (int f = 0; f < nF; ++f) fbox_ok[static_cast<std::size_t>(f)] = enclosure(f, fbox[static_cast<std::size_t>(f)]);

    struct CurveRec {
        int fa, fp;
        std::shared_ptr<TracedEdge> te; ///< sigma in [0, 1] over the whole traced curve
        bool closed;
    };
    std::vector<CurveRec> curves;
    for (int fa = 0; fa < nfA; ++fa)
        for (int fp = nfA; fp < nF; ++fp) {
            if (fbox_ok[static_cast<std::size_t>(fa)] && fbox_ok[static_cast<std::size_t>(fp)]) {
                const geom::IVec3& b = fbox[static_cast<std::size_t>(fa)];
                const geom::IVec3& c = fbox[static_cast<std::size_t>(fp)];
                if (num::disjoint(b.x, c.x) || num::disjoint(b.y, c.y) || num::disjoint(b.z, c.z)) continue;
            }
            const solve::SsiResult sr = solve::intersect_carriers(*fsurf[static_cast<std::size_t>(fa)], *fsurf[static_cast<std::size_t>(fp)], opt);
            if (!sr.complete) return fail("an intersection of two faces could not be proven complete (faces on one surface are not supported yet)");
            for (const solve::SsiCurve& c : sr.curves) {
                auto te = std::make_shared<TracedEdge>();
                te->a = fsurf[static_cast<std::size_t>(fa)];
                te->b = fsurf[static_cast<std::size_t>(fp)];
                const BFace& FA = M.faces[static_cast<std::size_t>(fa)];
                if (FA.proc >= 0) te->proc[0] = FA.proc;
                else te->carrier[0] = FA.carrier;
                const BFace& FP = M.faces[static_cast<std::size_t>(fp)];
                if (FP.proc >= 0) te->proc[1] = FP.proc;
                else te->carrier[1] = FP.carrier;
                te->path.segs = c.segs;
                te->path.make_breaks(te->sys(), *te->a);
                te->path.find_kinks(te->sys(), *te->a, *te->b);
                curves.push_back({fa, fp, te, c.closed});
            }
        }

    // 2. New vertices. A curve is traced over whole surfaces; where it leaves
    // one of its two faces it crosses an edge of that face: an edge of A
    // crossing a face of P, or an edge of P piercing a face of A. Every such
    // point lies on a curve (the curves are complete), so walking the curves
    // finds them all. Inside-ness is sampled along the curve (measured); each
    // change is refined to the edge by Gauss-Newton on curve(sigma) = edge(tau).
    struct Split {
        double tau;
        int vertex;
    };
    std::vector<std::vector<Split>> splits(M.edges.size());
    struct Crossing {
        int edge;
        Vec3d X;
        int vertex;
    };
    std::vector<Crossing> crossings;
    std::vector<std::vector<std::pair<double, int>>> curve_on(curves.size());
    auto curve_at = [&](const CurveRec& cr, double sg, solve::Pt4& x, Vec3d& X) {
        if (!cr.te->path.point(cr.te->sys(), std::clamp(sg, 0.0, 1.0), x)) return false;
        Vec3d Xu, Xv;
        cr.te->a->eval(x[0], x[1], X, Xu, Xv);
        return true;
    };
    // Closest approach of curve and edge near (sg, tau); returns the distance.
    auto meet_edge = [&](const CurveRec& cr, int ei, double s_lo, double s_hi, double& sg, double& tau, Vec3d& X) {
        const BEdge& e = M.edges[static_cast<std::size_t>(ei)];
        solve::Pt4 x;
        auto C = [&](double s) {
            Vec3d Y;
            if (!curve_at(cr, s, x, Y)) Y = Vec3d{1e300, 1e300, 1e300};
            return Y;
        };
        auto Ed = [&](double t) { return edge_point(e, std::clamp(t, 0.0, 1.0)); };
        auto diff = [&](auto&& f, double t) {
            const double h = 1e-7;
            const double a = std::max(0.0, t - h), b = std::min(1.0, t + h);
            return (f(b) - f(a)) / (b - a);
        };
        double d = 1e300;
        for (int it = 0; it < 40; ++it) {
            const Vec3d r = C(sg) - Ed(tau);
            d = norm(r);
            if (d < 1e-15) break;
            const Vec3d Cs = diff(C, sg), Et = diff(Ed, tau);
            const double a11 = dot(Cs, Cs), a12 = -dot(Cs, Et), a22 = dot(Et, Et);
            const double b1 = -dot(Cs, r), b2 = dot(Et, r), det = a11 * a22 - a12 * a12;
            if (!(std::fabs(det) > 0.0)) break;
            const double ds = (a22 * b1 - a12 * b2) / det, dt = (a11 * b2 - a12 * b1) / det;
            sg = std::clamp(sg + ds, s_lo, s_hi);
            tau = std::clamp(tau + dt, 0.0, 1.0);
            if (std::fabs(ds) + std::fabs(dt) < 1e-16) break;
        }
        X = C(sg);
        return norm(X - Ed(tau));
    };
    // The edge of face f that the curve meets for sigma in [s_lo, s_hi]
    // (starting at sg); -1 if none.
    auto find_edge = [&](const CurveRec& cr, int f, double s_lo, double s_hi, double& sg, double& tau, Vec3d& X) {
        Vec3d P0;
        solve::Pt4 x0;
        if (!curve_at(cr, sg, x0, P0)) return -1;
        struct Cand {
            double d, tau;
            int edge;
        };
        std::vector<Cand> cand;
        for (const auto& loop : M.faces[static_cast<std::size_t>(f)].loops)
            for (const BCoedge& ce : loop) {
                const BEdge& e = M.edges[static_cast<std::size_t>(ce.edge)];
                Cand c{1e300, 0.0, ce.edge};
                for (int k = 0; k <= 64; ++k) {
                    const double t = k / 64.0, dd = norm(edge_point(e, t) - P0);
                    if (dd < c.d) { c.d = dd; c.tau = t; }
                }
                cand.push_back(c);
            }
        std::sort(cand.begin(), cand.end(), [](const Cand& a, const Cand& b) { return a.d < b.d; });
        int best = -1;
        double bd = 1e-8;
        for (std::size_t k = 0; k < std::min<std::size_t>(3, cand.size()); ++k) {
            double s = sg, t = cand[k].tau;
            Vec3d Y;
            const double dd = meet_edge(cr, cand[k].edge, s_lo, s_hi, s, t, Y);
            if (dd < bd) { bd = dd; best = cand[k].edge; sg = s; tau = t; X = Y; }
        }
        return best;
    };
    for (std::size_t ci = 0; ci < curves.size(); ++ci) {
        const CurveRec& cr = curves[ci];
        const int ns = 256;
        std::vector<char> in_a(ns + 1), in_p(ns + 1);
        for (int i = 0; i <= ns; ++i) {
            solve::Pt4 x;
            Vec3d X;
            if (!curve_at(cr, static_cast<double>(i) / ns, x, X)) return fail("an intersection curve could not be evaluated");
            in_a[static_cast<std::size_t>(i)] = inside_face(cr.fa, x[0], x[1]);
            in_p[static_cast<std::size_t>(i)] = inside_face(cr.fp, x[2], x[3]);
        }
        // Events: (sigma guess, face whose boundary is crossed).
        std::vector<std::pair<double, int>> events;
        for (int i = 0; i < ns; ++i) {
            const double sm = (i + 0.5) / ns;
            if (in_a[static_cast<std::size_t>(i)] != in_a[static_cast<std::size_t>(i) + 1]) events.push_back({sm, cr.fa});
            if (in_p[static_cast<std::size_t>(i)] != in_p[static_cast<std::size_t>(i) + 1]) events.push_back({sm, cr.fp});
        }
        // A curve's ends lie on the boundary of a surface's patch, which is
        // often an edge of the face (a box face is its whole patch).
        if (!cr.closed)
            for (double se : {0.0, 1.0})
                for (int f : {cr.fa, cr.fp}) events.push_back({se, f});
        for (const auto& [s_guess, f] : events) {
            double sg = s_guess, tau = 0.0;
            Vec3d X;
            const bool at_end = s_guess == 0.0 || s_guess == 1.0;
            const double win = 1.5 / ns;
            const int ei = find_edge(cr, f, std::max(0.0, s_guess - win), std::min(1.0, s_guess + win), sg, tau, X);
            if (ei < 0) {
                if (at_end) continue; // the end is off the face
                return fail("an intersection curve leaves a face but meets none of its edges");
            }
            if (!(tau > tol && tau < 1.0 - tol)) return fail("an intersection passes through a vertex; not supported yet");
            int vx = -1;
            for (const Crossing& c : crossings)
                if (c.edge == ei && norm(c.X - X) < 1e-7) vx = c.vertex;
            if (vx < 0) {
                vx = add_vertex(X);
                crossings.push_back({ei, X, vx});
                splits[static_cast<std::size_t>(ei)].push_back({tau, vx});
                for (int g : edge_faces[static_cast<std::size_t>(ei)]) {
                    const BFace& F = M.faces[static_cast<std::size_t>(g)];
                    if (F.proc >= 0) continue;
                    double u, v, du, dv;
                    if (!edge_uv(F, M.edges[static_cast<std::size_t>(ei)], tau, u, v, du, dv))
                        return fail("an edge cannot be read in its face's parameters");
                    add_chart(vx, F.carrier, u, v);
                }
            }
            solve::Pt4 x;
            Vec3d Y;
            curve_at(cr, sg, x, Y);
            const BFace& FA = M.faces[static_cast<std::size_t>(cr.fa)];
            if (FA.proc < 0) add_chart(vx, FA.carrier, x[0], x[1]);
            const BFace& FP = M.faces[static_cast<std::size_t>(cr.fp)];
            if (FP.proc < 0) add_chart(vx, FP.carrier, x[2], x[3]);
            bool seen = false;
            for (const auto& o : curve_on[ci]) seen = seen || o.second == vx;
            if (!seen) curve_on[ci].push_back({sg, vx});
        }
    }

    // Edge pieces.
    struct EPiece {
        double t0, t1;
        int edge; ///< in the output edge list (M.edges is reused: pieces are appended)
    };
    std::vector<BEdge> E; // output edges
    std::vector<std::vector<EPiece>> pieces(M.edges.size());
    for (std::size_t ei = 0; ei < M.edges.size(); ++ei) {
        const BEdge& e = M.edges[ei];
        auto& sp = splits[ei];
        std::sort(sp.begin(), sp.end(), [](const Split& a, const Split& b) { return a.tau < b.tau; });
        std::vector<std::pair<double, int>> ends{{0.0, e.v0}};
        for (const Split& s : sp) ends.push_back({s.tau, s.vertex});
        ends.push_back({1.0, e.v1});
        for (std::size_t k = 0; k + 1 < ends.size(); ++k) {
            const double t0 = ends[k].first, t1 = ends[k + 1].first;
            BEdge ne;
            ne.v0 = ends[k].second;
            ne.v1 = ends[k + 1].second;
            if (e.proc) {
                if (t0 == 0.0 && t1 == 1.0) ne.proc = e.proc;
                else if (!(ne.proc = e.proc->part(t0, t1))) return fail("a procedural edge cannot be split");
            } else {
                ne.on[0] = e.on[0];
                ne.on[1] = e.on[1];
                const double a = num::mid(e.t0), b = num::mid(e.t1);
                ne.t0 = t0 == 0.0 ? e.t0 : Interval(a + t0 * (b - a));
                ne.t1 = t1 == 1.0 ? e.t1 : Interval(a + t1 * (b - a));
            }
            pieces[ei].push_back({t0, t1, static_cast<int>(E.size())});
            E.push_back(ne);
        }
    }

    // 2b. Pieces of the curves inside both faces.
    std::vector<std::vector<int>> face_curves(static_cast<std::size_t>(nF)); // output edges lying inside each face
    struct CurveEdgeInfo {
        int fa, fp;
    };
    std::map<int, CurveEdgeInfo> curve_edges;
    for (std::size_t ci = 0; ci < curves.size(); ++ci) {
        const CurveRec& cr = curves[ci];
        const TracedEdge& te = *cr.te;
        const auto sys = te.sys();
        auto& on = curve_on[ci];
        std::sort(on.begin(), on.end());
        if (cr.closed && on.empty()) {
            solve::Pt4 x;
            te.path.point(sys, 0.5, x);
            if (inside_face(cr.fa, x[0], x[1]) && inside_face(cr.fp, x[2], x[3]))
                return fail("an intersection curve is closed inside both faces; not supported yet");
            continue;
        }
        std::vector<double> cuts{0.0};
        for (const auto& o : on) cuts.push_back(o.first);
        cuts.push_back(1.0);
        for (std::size_t k = 0; k + 1 < cuts.size(); ++k) {
            const double s0 = cuts[k], s1 = cuts[k + 1];
            if ((k == 0 || k + 2 == cuts.size()) && s1 - s0 < 1e-9) continue;
            solve::Pt4 x;
            te.path.point(sys, 0.5 * (s0 + s1), x);
            if (!(inside_face(cr.fa, x[0], x[1]) && inside_face(cr.fp, x[2], x[3]))) continue;
            if (k == 0 || k + 2 == cuts.size()) return fail("an intersection curve leaves a face away from its edges");
            BEdge ne;
            ne.v0 = on[k - 1].second;
            ne.v1 = on[k].second;
            ne.proc = te.part(s0, s1);
            const int id = static_cast<int>(E.size());
            E.push_back(ne);
            face_curves[static_cast<std::size_t>(cr.fa)].push_back(id);
            face_curves[static_cast<std::size_t>(cr.fp)].push_back(id);
            curve_edges[id] = {cr.fa, cr.fp};
        }
    }
    res.curves = static_cast<int>(curve_edges.size());

    // 3. Regions of every face: the arrangement of its boundary pieces and
    // its curves, in its parameters.
    struct Half {
        int edge;
        bool rev;        ///< traversed from v1 to v0
        bool boundary;   ///< a piece of the face's boundary
        bool outward;    ///< a boundary piece traversed against the loop
        int from, to;
        double angle;    ///< direction at `from`, in the face's parameters
    };
    struct Region {
        int face;
        std::vector<std::vector<BCoedge>> loops; ///< outer loop first, then holes
        std::vector<int> halves;
        int cls = -1; ///< 1 inside the other solid, 0 outside, -1 unknown
    };
    std::vector<Region> regions;
    std::vector<std::vector<int>> regions_of_edge(E.size()); // regions that use an output edge on their boundary
    for (int f = 0; f < nF; ++f) {
        const BFace& F = M.faces[static_cast<std::size_t>(f)];
        std::vector<Half> H;
        auto add_half = [&](int edge, bool rev, bool boundary, bool outward) {
            const BEdge& e = E[static_cast<std::size_t>(edge)];
            Half h{edge, rev, boundary, outward, rev ? e.v1 : e.v0, rev ? e.v0 : e.v1, 0.0};
            // Direction at `from`: the edge's tangent at that end, in the face's parameters.
            double u, v, du, dv;
            bool ok;
            if (e.proc) {
                if (F.proc >= 0) {
                    auto st = [&](double t, double& s, double& w) { return e.proc->on_proc(F.proc, t, s, w); };
                    const double tau = rev ? 1.0 : 0.0;
                    ok = st(tau, u, v) && piece_derivative(e.proc->breaks(), tau, st, du, dv);
                } else {
                    ok = e.proc->on_carrier(F.carrier, rev ? 1.0 : 0.0, u, v, du, dv);
                }
            } else {
                ok = edge_uv(F, e, rev ? 1.0 : 0.0, u, v, du, dv);
            }
            if (!ok) return false;
            if (rev) { du = -du; dv = -dv; }
            h.angle = std::atan2(dv, du);
            H.push_back(h);
            return true;
        };
        for (const BCoedge& ce : F.loops[0])
            for (const EPiece& p : pieces[static_cast<std::size_t>(ce.edge)]) {
                if (!add_half(p.edge, ce.reversed, true, false) || !add_half(p.edge, !ce.reversed, true, true))
                    return fail("a boundary piece cannot be read in its face's parameters");
            }
        for (int ce : face_curves[static_cast<std::size_t>(f)])
            if (!add_half(ce, false, false, false) || !add_half(ce, true, false, false))
                return fail("an intersection curve cannot be read in its face's parameters");
        // Outgoing halves at each vertex, by angle.
        std::map<int, std::vector<int>> out;
        for (std::size_t i = 0; i < H.size(); ++i) out[H[i].from].push_back(static_cast<int>(i));
        for (auto& [vx, list] : out) {
            std::sort(list.begin(), list.end(), [&](int a, int b) { return H[static_cast<std::size_t>(a)].angle < H[static_cast<std::size_t>(b)].angle; });
            for (std::size_t k = 0; k + 1 < list.size(); ++k)
                if (H[static_cast<std::size_t>(list[k + 1])].angle - H[static_cast<std::size_t>(list[k])].angle < 1e-9)
                    return fail("two edges leave a vertex in the same direction (tangent contact); not supported yet");
        }
        // Twin of each half: same edge, other direction (boundary halves pair
        // with their outward copy).
        auto twin = [&](int i) {
            for (int j : out[H[static_cast<std::size_t>(i)].to])
                if (H[static_cast<std::size_t>(j)].edge == H[static_cast<std::size_t>(i)].edge && H[static_cast<std::size_t>(j)].rev != H[static_cast<std::size_t>(i)].rev)
                    return j;
            return -1;
        };
        // Next half after arriving by i: the outgoing half just clockwise of
        // the twin's direction (faces on the left).
        std::vector<char> used(H.size(), 0);
        std::vector<std::vector<int>> cycles;
        for (std::size_t i0 = 0; i0 < H.size(); ++i0) {
            if (used[i0]) continue;
            std::vector<int> cyc;
            int i = static_cast<int>(i0);
            bool outside = false;
            for (int guard = 0; guard < 100000; ++guard) {
                if (used[static_cast<std::size_t>(i)]) break;
                used[static_cast<std::size_t>(i)] = 1;
                cyc.push_back(i);
                outside = outside || H[static_cast<std::size_t>(i)].outward;
                const int tw = twin(i);
                if (tw < 0) return fail("an edge in a face's arrangement has no twin");
                const auto& list = out[H[static_cast<std::size_t>(i)].to];
                const auto pos = std::find(list.begin(), list.end(), tw) - list.begin();
                const int n = static_cast<int>(list.size());
                i = list[static_cast<std::size_t>((pos - 1 + n) % n)];
            }
            if (!outside) cycles.push_back(std::move(cyc));
        }
        // A cycle running the same way round as the face's own loop bounds a
        // region; one running the other way is a hole (a closed chain of
        // curves) in the smallest such region around it.
        auto cycle_poly = [&](const std::vector<int>& cyc) {
            std::vector<std::array<double, 2>> poly;
            for (int h : cyc) {
                const Half& hh = H[static_cast<std::size_t>(h)];
                const BEdge& e = E[static_cast<std::size_t>(hh.edge)];
                for (int k = 0; k < 16; ++k) {
                    const double tau = hh.rev ? 1.0 - k / 16.0 : k / 16.0;
                    double u = 0, v = 0, du, dv;
                    if (e.proc && F.proc >= 0) e.proc->on_proc(F.proc, tau, u, v);
                    else if (e.proc) e.proc->on_carrier(F.carrier, tau, u, v, du, dv);
                    else edge_uv(F, e, tau, u, v, du, dv);
                    poly.push_back({u, v});
                }
            }
            return poly;
        };
        const double face_sign = signed_area(outline[static_cast<std::size_t>(f)]) > 0.0 ? 1.0 : -1.0;
        std::vector<std::vector<std::array<double, 2>>> polys;
        std::vector<double> areas;
        for (const auto& cyc : cycles) {
            polys.push_back(cycle_poly(cyc));
            areas.push_back(signed_area(polys.back()) * face_sign);
        }
        // The island inside a hole runs along the same curves the other way:
        // it is never the hole's owner.
        auto share_edge = [&](const std::vector<int>& a, const std::vector<int>& b) {
            for (int x : a)
                for (int y : b)
                    if (H[static_cast<std::size_t>(x)].edge == H[static_cast<std::size_t>(y)].edge) return true;
            return false;
        };
        std::vector<int> owner(cycles.size(), -1);
        for (std::size_t c = 0; c < cycles.size(); ++c) {
            if (areas[c] > 0.0) continue;
            double best = 1e300;
            for (std::size_t o = 0; o < cycles.size(); ++o)
                if (areas[o] > 0.0 && areas[o] < best && !share_edge(cycles[o], cycles[c]) &&
                    winding(polys[o], polys[c][0][0], polys[c][0][1]) != 0) {
                    best = areas[o];
                    owner[c] = static_cast<int>(o);
                }
            if (owner[c] < 0) return fail("a closed chain of curves in a face lies in no region of it");
        }
        for (std::size_t c0 = 0; c0 < cycles.size(); ++c0) {
            if (!(areas[c0] > 0.0)) continue;
            Region R;
            R.face = f;
            std::vector<int> cyc;
            for (std::size_t c = 0; c < cycles.size(); ++c) {
                if (c != c0 && owner[c] != static_cast<int>(c0)) continue;
                std::vector<BCoedge> L;
                for (int h : cycles[c]) {
                    L.push_back({H[static_cast<std::size_t>(h)].edge, H[static_cast<std::size_t>(h)].rev});
                    R.halves.push_back(h);
                    cyc.push_back(h);
                }
                if (c == c0) R.loops.insert(R.loops.begin(), std::move(L));
                else R.loops.push_back(std::move(L));
            }
            // Classification from a curve on its boundary: the region lies on
            // the left of the half; inside the other solid if that side is
            // against the other face's outward normal.
            for (int h : cyc) {
                const Half& hh = H[static_cast<std::size_t>(h)];
                if (hh.boundary) continue;
                const auto& info = curve_edges[hh.edge];
                const int other = f == info.fa ? info.fp : info.fa;
                const BEdge& e = E[static_cast<std::size_t>(hh.edge)];
                const double tau = 0.5;
                // Tangent in space along the half, and this face's parameter normal.
                const Vec3d T = (e.proc->point(tau + 1e-6) - e.proc->point(tau - 1e-6)) * (hh.rev ? -1.0 : 1.0);
                double u, v, du, dv;
                if (F.proc >= 0) e.proc->on_proc(F.proc, tau, u, v);
                else e.proc->on_carrier(F.carrier, tau, u, v, du, dv);
                Vec3d S, Su, Sv;
                fsurf[static_cast<std::size_t>(f)]->eval(u, v, S, Su, Sv);
                const Vec3d left = cross(cross(Su, Sv), T);
                // Outward normal of the other face at X.
                double ou, ov;
                const BFace& O = M.faces[static_cast<std::size_t>(other)];
                if (O.proc >= 0) e.proc->on_proc(O.proc, tau, ou, ov);
                else e.proc->on_carrier(O.carrier, tau, ou, ov, du, dv);
                Vec3d Q, Qu, Qv;
                fsurf[static_cast<std::size_t>(other)]->eval(ou, ov, Q, Qu, Qv);
                Vec3d No = cross(Qu, Qv);
                if (O.flipped) No = No * -1.0;
                const int c = dot(left, No) < 0.0 ? 1 : 0;
                if (R.cls >= 0 && R.cls != c) return fail("a region is both inside and outside the other solid");
                R.cls = c;
            }
            for (int h : cyc)
                if (H[static_cast<std::size_t>(h)].boundary)
                    regions_of_edge[static_cast<std::size_t>(H[static_cast<std::size_t>(h)].edge)].push_back(static_cast<int>(regions.size()));
            regions.push_back(std::move(R));
        }
    }

    // 4. Unclassified regions inherit across shared boundary pieces (same
    // solid). A part of A that no curve touches takes a ray test against P:
    // the parity of certified line-surface roots inside P's faces, from a
    // point on the part's boundary; one test per part, then inheritance again.
    auto propagate = [&] {
        for (bool changed = true; changed;) {
            changed = false;
            for (const auto& rs : regions_of_edge) {
                if (rs.size() != 2) continue;
                Region& a = regions[static_cast<std::size_t>(rs[0])];
                Region& b = regions[static_cast<std::size_t>(rs[1])];
                if (solid_of_face(a.face) != solid_of_face(b.face)) continue;
                if (a.cls < 0 && b.cls >= 0) { a.cls = b.cls; changed = true; }
                else if (b.cls < 0 && a.cls >= 0) { b.cls = a.cls; changed = true; }
            }
        }
    };
    double reach = 0.0;
    bool P_rational = true;
    for (int q = nfA; q < nF; ++q) P_rational = P_rational && M.faces[static_cast<std::size_t>(q)].proc < 0;
    for (int q = nfA; q < nF && P_rational; ++q)
        for (const geom::Vec4d& c : M.carriers[static_cast<std::size_t>(M.faces[static_cast<std::size_t>(q)].carrier)].c)
            reach = std::max(reach, std::fabs(c.x / c.w) + std::fabs(c.y / c.w) + std::fabs(c.z / c.w));
    auto inside_P = [&](const Vec3d& X, bool& inside) {
        static constexpr double dirs[4][3] = {{0.5377, 0.8339, 0.1258}, {-0.2769, 0.3319, 0.9018},
                                              {0.7254, -0.6131, 0.3127}, {-0.6490, -0.3035, 0.6976}};
        const double Lr = 2.0 * (reach + std::fabs(X.x) + std::fabs(X.y) + std::fabs(X.z)) + 1.0;
        for (const auto& d : dirs) {
            const PolyPatchD ray = solve::make_line(X, X + Vec3d{d[0], d[1], d[2]} * Lr);
            int count = 0;
            bool valid = true;
            for (int q = nfA; q < nF && valid; ++q) {
                const PolyPatchD& C = M.carriers[static_cast<std::size_t>(M.faces[static_cast<std::size_t>(q)].carrier)];
                const solve::CsiResult cs = solve::intersect_curve_surface(solve::CsiSystem(ray, C));
                if (!cs.complete) { valid = false; break; }
                for (const solve::CsiRoot& r : cs.roots) {
                    const double t = num::mid(r.enc[0]), u = num::mid(r.enc[1]), v = num::mid(r.enc[2]);
                    if (!(r.enc[0].lo > 1e-9 && r.enc[0].hi < 1.0)) { valid = false; break; } // starts on a surface
                    if (!(u > -1e-12 && u < 1.0 + 1e-12 && v > -1e-12 && v < 1.0 + 1e-12)) continue;
                    // Near a face boundary the winding test is a measurement: try another ray.
                    bool near = false;
                    for (double du : {-1e-6, 1e-6})
                        for (double dv : {-1e-6, 1e-6}) near = near || inside_face(q, u + du, v + dv) != inside_face(q, u, v);
                    if (near) { valid = false; break; }
                    if (inside_face(q, u, v)) ++count;
                    (void)t;
                }
            }
            if (valid) { inside = count % 2 == 1; return true; }
        }
        return false;
    };
    propagate();
    for (std::size_t ri = 0; ri < regions.size(); ++ri) {
        Region& R = regions[ri];
        if (R.cls >= 0) continue;
        if (solid_of_face(R.face) == 1) return fail("a part of the second operand never meets the first's boundary; not supported yet");
        const BEdge& e0 = E[static_cast<std::size_t>(R.loops[0][0].edge)];
        const Vec3d X = e0.proc ? e0.proc->point(0.5) : edge_point(e0, 0.5);
        bool inside = false;
        if (!P_rational) return fail("a part of the first operand meets no curve, and the second has procedural faces; not supported yet");
        if (!inside_P(X, inside)) return fail("a part of the first operand could not be placed inside or outside the second");
        R.cls = inside ? 1 : 0;
        propagate();
    }

    // 5. Selection and assembly.
    Brep out;
    out.carriers = M.carriers;
    out.procs = M.procs;
    std::vector<int> vmap(M.vertices.size(), -1), emap(E.size(), -1);
    for (const Region& R : regions) {
        const int s = solid_of_face(R.face);
        bool keep = false, reverse = false;
        if (op == PolyOp::Union) keep = R.cls == 0;
        else if (op == PolyOp::Intersection) keep = R.cls == 1;
        else {
            keep = s == 0 ? R.cls == 0 : R.cls == 1;
            reverse = s == 1;
        }
        if (!keep) continue;
        const BFace& F = M.faces[static_cast<std::size_t>(R.face)];
        BFace nf;
        nf.carrier = F.carrier;
        nf.proc = F.proc;
        nf.flipped = reverse ? !F.flipped : F.flipped;
        for (const auto& loop : R.loops) {
            std::vector<BCoedge> L;
            for (const BCoedge& ce : loop) {
                int& id = emap[static_cast<std::size_t>(ce.edge)];
                if (id < 0) {
                    BEdge ne = E[static_cast<std::size_t>(ce.edge)];
                    for (int* vp : {&ne.v0, &ne.v1}) {
                        int& vm = vmap[static_cast<std::size_t>(*vp)];
                        if (vm < 0) {
                            vm = static_cast<int>(out.vertices.size());
                            out.vertices.push_back(M.vertices[static_cast<std::size_t>(*vp)]);
                        }
                        *vp = vm;
                    }
                    id = static_cast<int>(out.edges.size());
                    out.edges.push_back(ne);
                }
                L.push_back({id, ce.reversed});
            }
            nf.loops.push_back(std::move(L));
        }
        out.faces.push_back(std::move(nf));
    }
    if (out.faces.empty()) return fail("the result is empty");
    res.solid = std::move(out);
    res.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t_start).count();
    return res;
}

} // namespace ak::op
