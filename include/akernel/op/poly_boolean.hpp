// akernel — Booleans between a solid that may have procedural faces and a
// solid whose faces are all flat (docs/FEATURES.md §5.1; decided
// 2026-11-24: cuts by planar solids first).
//
// The rational Boolean engine (op/boolean.hpp) reads polynomial carriers
// only. This one works on any surfaces of the first operand A (polynomial
// or procedural, through geom::Carrier) and needs the second operand P to
// be flat-faced, which keeps every new curve the intersection of one of A's
// surfaces with a plane.
//
// 1. Where A's edges cross P's faces, and where P's edges pierce A's faces
//    (along the curves of the next step), are new vertices.
// 2. Every pair (face of A, face of P) whose enclosures meet is intersected
//    by intersect_carriers (certified, complete). The pieces of each curve
//    inside both faces are traced edges (op/proc_edge.hpp).
// 3. Each face is split by the pieces lying on it: an arrangement in its
//    parameters, traversed with the outgoing edges at each vertex sorted by
//    angle; the cycles that are not outside the face are its regions.
// 4. Regions next to a curve are inside or outside the other solid by the
//    side its outward normal points to; the others inherit it across the
//    uncut edges they share, and A's remaining regions take a ray test
//    against P.
// 5. The result keeps A's regions outside P (union, difference) or inside
//    it (intersection), and P's regions outside A (union), inside A
//    (intersection) or inside A reversed (difference).
//
// Certified: the intersection curves and their completeness. Measured in
// doubles, with refusals near a decision: where edges cross planes, which
// pieces lie inside faces (winding on sampled loops), the order of edges
// around vertices. Refused with a message for now: arc edges, faces of A
// lying in a plane of P, curves closed inside a face, a P component that
// never meets A's boundary.
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
#include <akernel/op/plane_cut.hpp>
#include <akernel/op/proc_edge.hpp>
#include <akernel/solve/carrier_ssi.hpp>

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

/// Parameters of a point on a planar patch (Newton; exact for affine ones).
inline bool plane_uv(const PolyPatchD& P, const Vec3d& X, double& u, double& v) {
    u = v = 0.5;
    const geom::BezierSurface S(P);
    for (int it = 0; it < 60; ++it) {
        const geom::Jet3 J = S.jet(u, v);
        const Vec3d d = X - J[0];
        const Vec3d& Su = J[static_cast<std::size_t>(geom::jet_index(1, 0))];
        const Vec3d& Sv = J[static_cast<std::size_t>(geom::jet_index(0, 1))];
        const double a11 = dot(Su, Su), a12 = dot(Su, Sv), a22 = dot(Sv, Sv);
        const double b1 = dot(Su, d), b2 = dot(Sv, d), det = a11 * a22 - a12 * a12;
        if (!(std::fabs(det) > 0.0)) return false;
        const double du = (a22 * b1 - a12 * b2) / det, dv = (a11 * b2 - a12 * b1) / det;
        u += du;
        v += dv;
        if (std::fabs(du) + std::fabs(dv) < 1e-16) break;
    }
    return std::isfinite(u) && std::isfinite(v);
}

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
    if (!P.procs.empty()) return fail("the second operand must have flat faces only");

    // The merged model: A's carriers then P's, A's procs; vertices, edges and
    // faces of A then of P, with P's carrier indices shifted.
    Brep M;
    const int ncA = static_cast<int>(A.carriers.size());
    M.carriers = A.carriers;
    M.carriers.insert(M.carriers.end(), P.carriers.begin(), P.carriers.end());
    M.procs = A.procs;
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
        for (auto& r : e.on) r.carrier += ncA;
        M.edges.push_back(e);
    }
    M.faces = A.faces;
    for (BFace f : P.faces) {
        f.carrier += ncA;
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
    // P's planes, outward.
    struct Plane {
        Vec3d n;
        double d;
    };
    std::vector<Plane> plane(static_cast<std::size_t>(nF));
    for (int f = nfA; f < nF; ++f) {
        const BFace& F = M.faces[static_cast<std::size_t>(f)];
        const PolyPatchD& C = M.carriers[static_cast<std::size_t>(F.carrier)];
        if (!planar_affine(C) && !same_plane_patches(C, C)) return fail("a face of the second operand is not flat");
        Vec3d Q, Qu, Qv;
        fsurf[static_cast<std::size_t>(f)]->eval(0.5, 0.5, Q, Qu, Qv);
        Vec3d n = cross(Qu, Qv);
        n = n / norm(n);
        if (F.flipped) n = n * -1.0;
        plane[static_cast<std::size_t>(f)] = {n, dot(n, Q)};
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

    // 1a. A's edges crossing P's faces.
    struct Split {
        double tau;
        int vertex;
    };
    std::vector<std::vector<Split>> splits(M.edges.size());
    std::vector<std::pair<int, int>> vertex_on; // (vertex, P face) for crossings of A edges
    for (int ei = 0; ei < neA; ++ei) {
        const BEdge& e = M.edges[static_cast<std::size_t>(ei)];
        const std::vector<double> bt = e.proc ? e.proc->breaks() : std::vector<double>{0.0, 1.0};
        std::vector<double> taus;
        for (std::size_t k = 0; k + 1 < bt.size(); ++k)
            for (int i = 0; i < 32; ++i)
                taus.push_back(bt[k] + (bt[k + 1] - bt[k]) * (i == 0 ? 0.0 : (i + 0.3819660112501051) / 32.5));
        taus.push_back(1.0);
        std::vector<Vec3d> pts(taus.size());
        for (std::size_t i = 0; i < taus.size(); ++i) pts[i] = edge_point(e, taus[i]);
        pts.front() = M.vertices[static_cast<std::size_t>(e.v0)].point;
        pts.back() = M.vertices[static_cast<std::size_t>(e.v1)].point;
        for (int q = nfA; q < nF; ++q) {
            const Plane& H = plane[static_cast<std::size_t>(q)];
            auto g = [&](const Vec3d& X) { return dot(H.n, X) - H.d; };
            for (std::size_t i = 0; i + 1 < taus.size(); ++i) {
                const double ga = g(pts[i]), gb = g(pts[i + 1]);
                if ((ga < 0) == (gb < 0)) continue;
                if (std::fabs(ga) < tol || std::fabs(gb) < tol) return fail("an edge of the first operand passes through a plane of the second at a vertex or sample");
                double a = taus[i], b = taus[i + 1], gl = ga;
                for (int it = 0; it < 80 && b - a > 1e-16; ++it) {
                    const double m = 0.5 * (a + b), gm = g(edge_point(e, m));
                    if ((gm < 0) == (gl < 0)) { a = m; gl = gm; }
                    else b = m;
                }
                const double tau = 0.5 * (a + b);
                const Vec3d X = edge_point(e, tau);
                const int qc = M.faces[static_cast<std::size_t>(q)].carrier;
                double qu, qv;
                if (!plane_uv(M.carriers[static_cast<std::size_t>(qc)], X, qu, qv)) return fail("a point cannot be placed on a plane");
                if (!inside_face(q, qu, qv)) continue;
                const int vx = add_vertex(X);
                add_chart(vx, qc, qu, qv);
                for (int f : edge_faces[static_cast<std::size_t>(ei)]) {
                    const BFace& F = M.faces[static_cast<std::size_t>(f)];
                    if (F.proc >= 0) continue;
                    double u, v, du, dv;
                    if (!edge_uv(F, e, tau, u, v, du, dv)) return fail("an edge cannot be read in its face's parameters");
                    add_chart(vx, F.carrier, u, v);
                }
                splits[static_cast<std::size_t>(ei)].push_back({tau, vx});
                vertex_on.push_back({vx, q});
            }
            for (std::size_t i = 1; i + 1 < taus.size(); ++i)
                if (std::fabs(g(pts[i])) < 1e-7) {
                    // Only a problem where the point is on the face.
                    double qu, qv;
                    if (plane_uv(M.carriers[static_cast<std::size_t>(M.faces[static_cast<std::size_t>(q)].carrier)], pts[i], qu, qv) &&
                        inside_face(q, qu, qv))
                        return fail("an edge of the first operand nearly touches a face of the second");
                }
        }
        for (int q = nfA; q < nF; ++q) {
            const Plane& H = plane[static_cast<std::size_t>(q)];
            for (int vx : {e.v0, e.v1})
                if (std::fabs(dot(H.n, M.vertices[static_cast<std::size_t>(vx)].point) - H.d) < tol) {
                    double qu, qv;
                    if (plane_uv(M.carriers[static_cast<std::size_t>(M.faces[static_cast<std::size_t>(q)].carrier)],
                                 M.vertices[static_cast<std::size_t>(vx)].point, qu, qv) &&
                        inside_face(q, qu, qv))
                        return fail("a vertex of the first operand lies on a face of the second");
                }
        }
    }

    // 2. Curves of (A face, P face) pairs.
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
            const Plane& H = plane[static_cast<std::size_t>(fp)];
            if (fbox_ok[static_cast<std::size_t>(fa)]) {
                const geom::IVec3& b = fbox[static_cast<std::size_t>(fa)];
                const Interval gn = b.x * Interval(H.n.x) + b.y * Interval(H.n.y) + b.z * Interval(H.n.z) - Interval(H.d);
                if (gn.lo > 0.0 || gn.hi < 0.0) continue; // the surface misses the plane
                if (fbox_ok[static_cast<std::size_t>(fp)]) {
                    const geom::IVec3& c = fbox[static_cast<std::size_t>(fp)];
                    if (num::disjoint(b.x, c.x) || num::disjoint(b.y, c.y) || num::disjoint(b.z, c.z)) continue;
                }
                if (gn.hi - gn.lo < 1e-12) return fail("a face of the first operand lies in a plane of the second; not supported yet");
            }
            const solve::SsiResult sr = solve::intersect_carriers(*fsurf[static_cast<std::size_t>(fa)], *fsurf[static_cast<std::size_t>(fp)], opt);
            if (!sr.complete) return fail("an intersection of two faces could not be proven complete");
            for (const solve::SsiCurve& c : sr.curves) {
                auto te = std::make_shared<TracedEdge>();
                te->a = fsurf[static_cast<std::size_t>(fa)];
                te->b = fsurf[static_cast<std::size_t>(fp)];
                const BFace& FA = M.faces[static_cast<std::size_t>(fa)];
                if (FA.proc >= 0) te->proc[0] = FA.proc;
                else te->carrier[0] = FA.carrier;
                te->carrier[1] = M.faces[static_cast<std::size_t>(fp)].carrier;
                te->path.segs = c.segs;
                te->path.make_breaks(te->sys(), *te->a);
                te->path.find_kinks(te->sys(), *te->a, *te->b);
                curves.push_back({fa, fp, te, c.closed});
            }
        }

    // 1b. P's edges piercing A's faces, found along the curves of one of the
    // edge's two faces.
    struct Pierce {
        int vertex, edge, fa;
    };
    std::vector<Pierce> pierces;
    for (int ei = neA; ei < static_cast<int>(M.edges.size()); ++ei) {
        const auto& ef = edge_faces[static_cast<std::size_t>(ei)];
        if (ef.size() != 2) return fail("an edge of the second operand is not used by two faces");
        const int q1 = ef[0], q2 = ef[1];
        const BEdge& e = M.edges[static_cast<std::size_t>(ei)];
        const Plane& H2 = plane[static_cast<std::size_t>(q2)];
        for (const CurveRec& cr : curves) {
            if (cr.fp != q1) continue;
            const TracedEdge& te = *cr.te;
            const auto sys = te.sys();
            const int ns = 512;
            std::vector<double> sg, gv;
            for (int i = 0; i <= ns; ++i) {
                const double s = (i + (i == 0 || i == ns ? 0.0 : 0.3819660112501051)) / ns;
                solve::Pt4 x;
                if (!te.path.point(sys, std::min(1.0, s), x)) return fail("an intersection curve could not be evaluated");
                Vec3d X, Xu, Xv;
                te.a->eval(x[0], x[1], X, Xu, Xv);
                sg.push_back(std::min(1.0, s));
                gv.push_back(dot(H2.n, X) - H2.d);
            }
            // The curve's ends lie on the plane q2 when the face q1 ends at
            // the edge there (the usual case: the curve is traced inside q1's
            // patch); crossings inside the curve are bracketed and bisected.
            const bool end0 = std::fabs(gv.front()) < 1e-9, end1 = std::fabs(gv.back()) < 1e-9;
            std::vector<double> hits;
            if (end0) hits.push_back(0.0);
            if (end1) hits.push_back(1.0);
            for (int i = 0; i < ns; ++i) {
                if ((i == 0 && end0) || (i + 1 == ns && end1)) continue;
                if ((gv[static_cast<std::size_t>(i)] < 0) == (gv[static_cast<std::size_t>(i) + 1] < 0)) continue;
                double a = sg[static_cast<std::size_t>(i)], b = sg[static_cast<std::size_t>(i) + 1], gl = gv[static_cast<std::size_t>(i)];
                for (int it = 0; it < 80 && b - a > 1e-16; ++it) {
                    const double m = 0.5 * (a + b);
                    solve::Pt4 x;
                    te.path.point(sys, m, x);
                    Vec3d X, Xu, Xv;
                    te.a->eval(x[0], x[1], X, Xu, Xv);
                    const double gm = dot(H2.n, X) - H2.d;
                    if ((gm < 0) == (gl < 0)) { a = m; gl = gm; }
                    else b = m;
                }
                hits.push_back(0.5 * (a + b));
            }
            for (double hit : hits) {
                solve::Pt4 x;
                te.path.point(sys, hit, x);
                // On the edge's segment, and inside A's face.
                const BEdge::SideRef* ref = side_on(e, M.faces[static_cast<std::size_t>(q1)].carrier);
                if (!ref) return fail("an edge of the second operand is not a side of its face");
                const double su = x[2], sv = x[3];
                double t = (ref->side == 0 || ref->side == 2) ? su : sv;
                if (ref->reversed) t = 1.0 - t;
                const double t0 = num::mid(e.t0), t1 = num::mid(e.t1);
                const double tau = (t - t0) / (t1 - t0);
                if (!(tau > 1e-9 && tau < 1.0 - 1e-9)) {
                    if (tau > -1e-9 && tau < 1.0 + 1e-9) return fail("a vertex of the second operand lies on a face of the first");
                    continue;
                }
                if (!inside_face(cr.fa, x[0], x[1])) continue;
                Vec3d X, Xu, Xv;
                te.a->eval(x[0], x[1], X, Xu, Xv);
                const int vx = add_vertex(X);
                const BFace& FA = M.faces[static_cast<std::size_t>(cr.fa)];
                if (FA.proc < 0) add_chart(vx, FA.carrier, x[0], x[1]);
                for (int q : {q1, q2}) {
                    const BEdge::SideRef* r = side_on(e, M.faces[static_cast<std::size_t>(q)].carrier);
                    const double tt = t0 + tau * (t1 - t0);
                    double u, v;
                    side_uv(r->side, r->reversed ? 1.0 - tt : tt, u, v);
                    add_chart(vx, M.faces[static_cast<std::size_t>(q)].carrier, u, v);
                }
                splits[static_cast<std::size_t>(ei)].push_back({tau, vx});
                pierces.push_back({vx, ei, cr.fa});
            }
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
    for (const CurveRec& cr : curves) {
        const TracedEdge& te = *cr.te;
        const auto sys = te.sys();
        // Points on this curve: crossings of fa's edges with fp, pierces of
        // fp's edges through fa.
        std::vector<int> pts; // vertices on the curve, located below by distance in space
        const BFace& FA = M.faces[static_cast<std::size_t>(cr.fa)];
        for (const auto& [vx, q] : vertex_on) {
            if (q != cr.fp) continue;
            // Is vx on one of fa's edges?
            bool on_fa = false;
            double u = 0, v = 0;
            for (const BCoedge& ce : FA.loops[0])
                for (const Split& s : splits[static_cast<std::size_t>(ce.edge)])
                    if (s.vertex == vx) {
                        double du, dv;
                        on_fa = edge_uv(FA, M.edges[static_cast<std::size_t>(ce.edge)], s.tau, u, v, du, dv);
                    }
            if (on_fa) pts.push_back(vx);
        }
        for (const Pierce& p : pierces) {
            if (p.fa != cr.fa) continue;
            const auto& ef = edge_faces[static_cast<std::size_t>(p.edge)];
            if (ef[0] != cr.fp && ef[1] != cr.fp) continue;
            pts.push_back(p.vertex);
        }
        // Sample the curve.
        const int ns = 512;
        std::vector<solve::Pt4> smp(ns + 1);
        std::vector<Vec3d> spt(ns + 1);
        for (int i = 0; i <= ns; ++i) {
            if (!te.path.point(sys, static_cast<double>(i) / ns, smp[static_cast<std::size_t>(i)])) return fail("an intersection curve could not be evaluated");
            Vec3d Xu, Xv;
            te.a->eval(smp[static_cast<std::size_t>(i)][0], smp[static_cast<std::size_t>(i)][1], spt[static_cast<std::size_t>(i)], Xu, Xv);
        }
        std::vector<std::pair<double, int>> on;
        for (int vx : pts) {
            const Vec3d X = M.vertices[static_cast<std::size_t>(vx)].point;
            auto dist = [&](double s) {
                solve::Pt4 x;
                te.path.point(sys, s, x);
                Vec3d Y, Yu, Yv;
                te.a->eval(x[0], x[1], Y, Yu, Yv);
                return norm(Y - X);
            };
            int best = 0;
            double bd = 1e300;
            for (int i = 0; i <= ns; ++i) {
                const double dd = norm(spt[static_cast<std::size_t>(i)] - X);
                if (dd < bd) { bd = dd; best = i; }
            }
            double a = std::max(0.0, (best - 1.0) / ns), b = std::min(1.0, (best + 1.0) / ns);
            const double gr = 0.5 * (std::sqrt(5.0) - 1.0);
            double c1 = b - gr * (b - a), c2 = a + gr * (b - a), f1 = dist(c1), f2 = dist(c2);
            for (int it = 0; it < 100 && b - a > 1e-15; ++it) {
                if (f1 < f2) { b = c2; c2 = c1; f2 = f1; c1 = b - gr * (b - a); f1 = dist(c1); }
                else { a = c1; c1 = c2; f1 = f2; c2 = a + gr * (b - a); f2 = dist(c2); }
            }
            const double s = 0.5 * (a + b);
            if (dist(s) < 1e-7) on.push_back({s, vx});
        }
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
        std::vector<BCoedge> loop;
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
            if (outside) continue;
            Region R;
            R.face = f;
            for (int h : cyc) {
                R.loop.push_back({H[static_cast<std::size_t>(h)].edge, H[static_cast<std::size_t>(h)].rev});
                R.halves.push_back(h);
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
    // solid), then A's remaining ones take a ray test against P.
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
    for (Region& R : regions) {
        if (R.cls >= 0) continue;
        if (solid_of_face(R.face) == 1) return fail("a part of the second operand never meets the first's boundary; not supported yet");
        // A point of the region: a vertex of its loop, nudged is not needed
        // for parity against flat faces unless it lies on one.
        const BEdge& e0 = E[static_cast<std::size_t>(R.loop[0].edge)];
        const Vec3d X = e0.proc ? e0.proc->point(0.5) : edge_point(e0, 0.5);
        // Ray in an irrational direction; parity of crossings with P's faces.
        const Vec3d dir = Vec3d{0.5773502691896258, 0.6123724356957945, 0.5400617248673217};
        int crossings = 0;
        for (int q = nfA; q < nF; ++q) {
            const Plane& Hq = plane[static_cast<std::size_t>(q)];
            const double den = dot(Hq.n, dir);
            if (std::fabs(den) < 1e-12) continue;
            const double t = (Hq.d - dot(Hq.n, X)) / den;
            if (!(t > 1e-12)) continue;
            double qu, qv;
            if (!plane_uv(M.carriers[static_cast<std::size_t>(M.faces[static_cast<std::size_t>(q)].carrier)], X + dir * t, qu, qv)) continue;
            if (inside_face(q, qu, qv)) ++crossings;
        }
        R.cls = crossings % 2;
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
        std::vector<BCoedge> L;
        for (const BCoedge& ce : R.loop) {
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
        out.faces.push_back(std::move(nf));
    }
    if (out.faces.empty()) return fail("the result is empty");
    res.solid = std::move(out);
    res.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t_start).count();
    return res;
}

} // namespace ak::op
