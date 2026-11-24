// akernel — cut a solid by a plane, keeping the side n . x <= d. Works on
// solids with procedural faces (docs/FEATURES.md §5.1): the first step of
// Booleans on such solids, decided 2026-11-24 ("cuts by planar solids
// first").
//
// Every face the plane meets is intersected with a planar cap patch by
// intersect_carriers (certified, with a completeness proof), so each
// section curve is a TracedEdge between the face's surface and the cap:
// one definition, read by the face and by the cap in their own parameters.
// The face keeps the region on the kept side; the cap is bounded by the
// section curves.
//
// What is certified and what is not. The section curves are certified
// tubes with completeness. Where the boundary edges cross the plane is
// found in doubles (sign changes along the edge, then bisection), and the
// pieces of a section curve inside a face are told apart by a winding test
// on sampled loops: measurements, refused when close to a decision
// (a vertex or an edge within 1e-9 of the plane, a curve that does not
// cross cleanly).
//
// Limits of this first version, each refused with a message: arc edges
// (intersection edges of rational faces from earlier Booleans); section
// curves closed inside a face; faces with inner loops that the plane
// crosses or splits into more than one region.
#pragma once

#include <algorithm>
#include <chrono>
#include <cmath>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include <akernel/op/brep.hpp>
#include <akernel/op/proc_edge.hpp>
#include <akernel/solve/carrier_ssi.hpp>

namespace ak::op {

struct PlaneCutResult {
    std::optional<Brep> solid;
    std::string error;
    int section_curves = 0; ///< pieces of section curves (edges of the cap)
    double seconds = 0.0;
};

namespace plane_cut_detail {

/// Winding number of (u, v) about a closed polyline.
inline int winding(const std::vector<std::array<double, 2>>& poly, double u, double v) {
    int wn = 0;
    for (std::size_t i = 0; i < poly.size(); ++i) {
        const auto& a = poly[i];
        const auto& b = poly[(i + 1) % poly.size()];
        const double c = (b[0] - a[0]) * (v - a[1]) - (u - a[0]) * (b[1] - a[1]);
        if (a[1] <= v) {
            if (b[1] > v && c > 0) ++wn;
        } else if (b[1] <= v && c < 0) {
            --wn;
        }
    }
    return wn;
}

inline double signed_area(const std::vector<std::array<double, 2>>& poly) {
    double A = 0.0;
    for (std::size_t i = 0; i < poly.size(); ++i) {
        const auto& a = poly[i];
        const auto& b = poly[(i + 1) % poly.size()];
        A += a[0] * b[1] - b[0] * a[1];
    }
    return A / 2;
}

} // namespace plane_cut_detail

inline PlaneCutResult cut_by_plane(const Brep& in, Vec3d n, double d, const solve::SsiOptions& opt = {}) {
    using namespace plane_cut_detail;
    const auto t_start = std::chrono::steady_clock::now();
    PlaneCutResult res;
    auto fail = [&](std::string msg) {
        res.error = std::move(msg);
        res.solid.reset();
        return res;
    };
    const double nl = norm(n);
    if (!(nl > 0.0)) return fail("the plane's normal is zero");
    n = n / nl;
    d /= nl;
    const double tol = 1e-9;
    auto g = [&](const Vec3d& P) { return dot(n, P) - d; };

    for (const BEdge& e : in.edges)
        if (e.arc) return fail("arc edges are not supported by the plane cut yet");

    // Surfaces: polynomial carriers, then procedural ones.
    const int nc = static_cast<int>(in.carriers.size());
    std::vector<std::shared_ptr<const geom::Carrier>> surf;
    for (const PolyPatchD& c : in.carriers) surf.push_back(std::make_shared<geom::PatchCarrier>(c));
    for (const auto& p : in.procs) surf.push_back(p);
    auto sid = [&](const BFace& f) { return f.proc >= 0 ? nc + f.proc : f.carrier; };

    // The cap carrier: a rectangle in the plane, O + u E1 + v E2 with
    // E1 x E2 along n (outward for the kept side), around the solid's
    // vertices and edge samples, with a margin.
    Vec3d e1 = std::fabs(n.x) < 0.9 ? cross(n, Vec3d{1, 0, 0}) : cross(n, Vec3d{0, 1, 0});
    e1 = e1 / norm(e1);
    const Vec3d e2 = cross(n, e1);
    auto edge_point = [&](const BEdge& e, double tau) {
        if (e.proc) return e.proc->point(tau);
        const BEdge::SideRef& ref = e.on[0];
        const double t = num::mid(e.t0) + tau * (num::mid(e.t1) - num::mid(e.t0));
        double u = 0.0, v = 0.0;
        side_uv(ref.side, ref.reversed ? 1.0 - t : t, u, v);
        return geom::BezierSurface(in.carriers[static_cast<std::size_t>(ref.carrier)]).point(u, v);
    };
    double lo1 = 1e300, hi1 = -1e300, lo2 = 1e300, hi2 = -1e300;
    for (const BEdge& e : in.edges)
        for (int k = 0; k <= 32; ++k) {
            const Vec3d P = edge_point(e, k / 32.0);
            lo1 = std::min(lo1, dot(P, e1)); hi1 = std::max(hi1, dot(P, e1));
            lo2 = std::min(lo2, dot(P, e2)); hi2 = std::max(hi2, dot(P, e2));
        }
    // Faces may bulge beyond their edges: a generous margin.
    const double m1 = 0.5 * (hi1 - lo1) + 1e-3, m2 = 0.5 * (hi2 - lo2) + 1e-3;
    lo1 -= m1; hi1 += m1; lo2 -= m2; hi2 += m2;
    const Vec3d O = n * d + e1 * lo1 + e2 * lo2, E1 = e1 * (hi1 - lo1), E2 = e2 * (hi2 - lo2);
    PolyPatchD cap(1, 1);
    auto h = [](const Vec3d& P) { return geom::Vec4d{P.x, P.y, P.z, 1.0}; };
    cap.at(0, 0) = h(O);
    cap.at(1, 0) = h(O + E1);
    cap.at(0, 1) = h(O + E2);
    cap.at(1, 1) = h(O + E1 + E2);
    auto cap_uv = [&](const Vec3d& P, double& u, double& v) {
        u = dot(P - O, E1) / dot(E1, E1);
        v = dot(P - O, E2) / dot(E2, E2);
    };
    auto cap_carrier = std::make_shared<geom::PatchCarrier>(cap);

    Brep out;
    out.carriers = in.carriers;
    out.carriers.push_back(cap);
    const int kCap = nc;
    out.procs = in.procs;

    // Vertices: signs; kept ones are copied.
    std::vector<int> vnew(in.vertices.size(), -1);
    std::vector<double> gv(in.vertices.size());
    for (std::size_t i = 0; i < in.vertices.size(); ++i) {
        gv[i] = g(in.vertices[i].point);
        if (std::fabs(gv[i]) < tol) return fail("a vertex lies on the cutting plane");
        if (gv[i] < 0) {
            vnew[i] = static_cast<int>(out.vertices.size());
            out.vertices.push_back(in.vertices[i]);
        }
    }

    // Faces that use each edge (for charts of new vertices).
    std::vector<std::vector<int>> edge_faces(in.edges.size());
    for (std::size_t f = 0; f < in.faces.size(); ++f)
        for (const auto& loop : in.faces[f].loops)
            for (const BCoedge& ce : loop) edge_faces[static_cast<std::size_t>(ce.edge)].push_back(static_cast<int>(f));

    // Crossings of every edge, and its kept pieces as new edges.
    struct Crossing {
        double tau;
        int vertex; ///< new vertex
    };
    struct Piece {
        double t0, t1;
        int edge; ///< new edge
    };
    std::vector<std::vector<Crossing>> cross_of(in.edges.size());
    std::vector<std::vector<Piece>> pieces_of(in.edges.size());
    for (std::size_t ei = 0; ei < in.edges.size(); ++ei) {
        const BEdge& e = in.edges[ei];
        const std::vector<double> bt = e.proc ? e.proc->breaks() : std::vector<double>{0.0, 1.0};
        std::vector<double> taus;
        for (std::size_t k = 0; k + 1 < bt.size(); ++k)
            for (int i = 0; i < 32; ++i)
                // Offset samples, so that a plane through a simple fraction of
                // the edge does not fall on a sample.
                taus.push_back(bt[k] + (bt[k + 1] - bt[k]) * (i == 0 ? 0.0 : (i + 0.3819660112501051) / 32.5));
        taus.push_back(1.0);
        std::vector<double> gs(taus.size());
        for (std::size_t i = 0; i < taus.size(); ++i) gs[i] = g(edge_point(e, taus[i]));
        gs.front() = gv[static_cast<std::size_t>(e.v0)];
        gs.back() = gv[static_cast<std::size_t>(e.v1)];
        for (std::size_t i = 0; i + 1 < taus.size(); ++i) {
            if (std::fabs(gs[i + 1]) < tol) return fail("an edge passes within 1e-9 of the cutting plane at a sample");
            if ((gs[i] < 0) == (gs[i + 1] < 0)) continue;
            double a = taus[i], b = taus[i + 1], ga = gs[i];
            for (int it = 0; it < 80 && b - a > 1e-16; ++it) {
                const double m = 0.5 * (a + b), gm = g(edge_point(e, m));
                if ((gm < 0) == (ga < 0)) { a = m; ga = gm; }
                else b = m;
            }
            const double tau = 0.5 * (a + b);
            BVertex v;
            v.point = edge_point(e, tau);
            // Charts on the polynomial faces that use the edge, and on the cap.
            for (int f : edge_faces[ei]) {
                const BFace& face = in.faces[static_cast<std::size_t>(f)];
                if (face.proc >= 0 || v.chart(face.carrier)) continue;
                double u, w;
                if (e.proc) {
                    double du, dw;
                    if (!e.proc->on_carrier(face.carrier, tau, u, w, du, dw)) return fail("a procedural edge is not on its face");
                } else {
                    const BEdge::SideRef* ref = side_on(e, face.carrier);
                    if (!ref) return fail("a side edge is not a side of its face's carrier");
                    const double t = num::mid(e.t0) + tau * (num::mid(e.t1) - num::mid(e.t0));
                    side_uv(ref->side, ref->reversed ? 1.0 - t : t, u, w);
                }
                v.charts.push_back({face.carrier, Interval(u), Interval(w)});
            }
            double cu, cv;
            cap_uv(v.point, cu, cv);
            v.charts.push_back({kCap, Interval(cu), Interval(cv)});
            cross_of[ei].push_back({tau, static_cast<int>(out.vertices.size())});
            out.vertices.push_back(v);
        }
        // Near-touches without a sign change.
        for (std::size_t i = 0; i < gs.size(); ++i)
            if (std::fabs(gs[i]) < 1e-7 && i > 0 && i + 1 < gs.size())
                return fail("an edge nearly touches the cutting plane");
        // Kept pieces between consecutive crossings.
        std::vector<std::pair<double, int>> ends; // (tau, vertex)
        ends.push_back({0.0, e.v0});
        for (const Crossing& c : cross_of[ei]) ends.push_back({c.tau, -1 - c.vertex});
        ends.push_back({1.0, e.v1});
        auto newv = [&](int code) { return code >= 0 ? vnew[static_cast<std::size_t>(code)] : -1 - code; };
        for (std::size_t k = 0; k + 1 < ends.size(); ++k) {
            const double t0 = ends[k].first, t1 = ends[k + 1].first;
            if (!(g(edge_point(e, 0.5 * (t0 + t1))) < 0)) continue;
            BEdge ne;
            ne.v0 = newv(ends[k].second);
            ne.v1 = newv(ends[k + 1].second);
            if (ne.v0 < 0 || ne.v1 < 0) return fail("inconsistent signs along an edge");
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
            pieces_of[ei].push_back({t0, t1, static_cast<int>(out.edges.size())});
            out.edges.push_back(ne);
        }
    }

    // A face's parameters of a point of one of its edges.
    auto face_uv = [&](const BFace& face, const BEdge& e, double tau, double& u, double& v) {
        if (face.proc >= 0) return e.proc && e.proc->on_proc(face.proc, tau, u, v);
        if (e.proc) {
            double du, dv;
            return e.proc->on_carrier(face.carrier, tau, u, v, du, dv);
        }
        const BEdge::SideRef* ref = side_on(e, face.carrier);
        if (!ref) return false;
        const double t = num::mid(e.t0) + tau * (num::mid(e.t1) - num::mid(e.t0));
        side_uv(ref->side, ref->reversed ? 1.0 - t : t, u, v);
        return true;
    };

    // Section edges, per crossing vertex: the one leaving it (in its face's
    // loop) and the one arriving.
    std::map<int, std::vector<int>> section_from, section_to; // vertex -> section edges
    std::vector<int> section_edges;

    for (std::size_t fi = 0; fi < in.faces.size(); ++fi) {
        const BFace& face = in.faces[fi];
        bool crossed = false;
        for (const auto& loop : face.loops)
            for (const BCoedge& ce : loop) crossed = crossed || !cross_of[static_cast<std::size_t>(ce.edge)].empty();
        const int S = sid(face);
        if (!crossed) {
            // Whole face kept or dropped; a section curve inside it would be closed.
            const BEdge& e0 = in.edges[static_cast<std::size_t>(face.loops[0][0].edge)];
            if (gv[static_cast<std::size_t>(e0.v0)] > 0) continue;
            BFace nf = face;
            for (auto& loop : nf.loops)
                for (BCoedge& ce : loop) {
                    const auto& ps = pieces_of[static_cast<std::size_t>(ce.edge)];
                    if (ps.size() != 1) return fail("an uncrossed face has an edge that is cut");
                    ce.edge = ps[0].edge;
                }
            out.faces.push_back(nf);
            continue;
        }
        if (face.loops.size() != 1) return fail("faces with inner loops that the plane crosses are not supported yet");

        // The loop as a polyline in the face's parameters, for the winding test.
        std::vector<std::array<double, 2>> poly;
        for (const BCoedge& ce : face.loops[0]) {
            const BEdge& e = in.edges[static_cast<std::size_t>(ce.edge)];
            const int m = 64;
            for (int k = 0; k < m; ++k) {
                const double tau = ce.reversed ? 1.0 - static_cast<double>(k) / m : static_cast<double>(k) / m;
                double u = 0.0, v = 0.0;
                if (!face_uv(face, e, tau, u, v)) return fail("an edge cannot be read in its face's parameters");
                poly.push_back({u, v});
            }
        }

        // Section curves of the face's surface with the cap.
        const solve::SsiResult sr = solve::intersect_carriers(*surf[static_cast<std::size_t>(S)], *cap_carrier, opt);
        if (!sr.complete) return fail("the section of a face with the plane could not be proven complete");

        // Crossing points of this face, in its parameters.
        struct FacePoint {
            int vertex;
            double u = 0.0, v = 0.0;
        };
        std::vector<FacePoint> fps;
        for (const BCoedge& ce : face.loops[0]) {
            const BEdge& e = in.edges[static_cast<std::size_t>(ce.edge)];
            for (const Crossing& c : cross_of[static_cast<std::size_t>(ce.edge)]) {
                double u = 0.0, v = 0.0;
                if (!face_uv(face, e, c.tau, u, v)) return fail("a crossing cannot be read in its face's parameters");
                fps.push_back({c.vertex, u, v});
            }
        }
        std::vector<char> used(fps.size(), 0);
        std::map<int, int> face_section; // crossing vertex -> this face's section edge leaving it

        for (const solve::SsiCurve& curve : sr.curves) {
            auto te = std::make_shared<TracedEdge>();
            te->a = surf[static_cast<std::size_t>(S)];
            te->b = cap_carrier;
            if (face.proc >= 0) te->proc[0] = face.proc;
            else te->carrier[0] = face.carrier;
            te->carrier[1] = kCap;
            te->path.segs = curve.segs;
            te->path.make_breaks(te->sys(), *te->a);
            auto coords = [&](double sg, double& u, double& v) {
                solve::Pt4 x;
                if (!te->path.point(te->sys(), sg, x)) return false;
                u = x[0];
                v = x[1];
                return true;
            };
            // Locate the face's crossing points on this curve: nearest sample,
            // then golden-section refinement of the parameter distance.
            const int ns = 400;
            std::vector<std::array<double, 2>> smp(ns + 1);
            for (int i = 0; i <= ns; ++i)
                if (!coords(static_cast<double>(i) / ns, smp[static_cast<std::size_t>(i)][0], smp[static_cast<std::size_t>(i)][1]))
                    return fail("a section curve could not be evaluated");
            std::vector<std::pair<double, int>> on; // (sigma, index into fps)
            for (std::size_t k = 0; k < fps.size(); ++k) {
                int best = 0;
                double bd = 1e300;
                for (int i = 0; i <= ns; ++i) {
                    const double du = smp[static_cast<std::size_t>(i)][0] - fps[k].u, dv = smp[static_cast<std::size_t>(i)][1] - fps[k].v;
                    if (du * du + dv * dv < bd) { bd = du * du + dv * dv; best = i; }
                }
                double a = std::max(0.0, (best - 1.0) / ns), b = std::min(1.0, (best + 1.0) / ns);
                auto dist = [&](double sg) {
                    double u = 1e300, v = 1e300;
                    coords(sg, u, v);
                    return (u - fps[k].u) * (u - fps[k].u) + (v - fps[k].v) * (v - fps[k].v);
                };
                const double gr = 0.5 * (std::sqrt(5.0) - 1.0);
                double c1 = b - gr * (b - a), c2 = a + gr * (b - a), f1 = dist(c1), f2 = dist(c2);
                for (int it = 0; it < 100 && b - a > 1e-15; ++it) {
                    if (f1 < f2) { b = c2; c2 = c1; f2 = f1; c1 = b - gr * (b - a); f1 = dist(c1); }
                    else { a = c1; c1 = c2; f1 = f2; c2 = a + gr * (b - a); f2 = dist(c2); }
                }
                const double sg = 0.5 * (a + b);
                if (std::sqrt(dist(sg)) < 1e-7) {
                    if (used[k]) return fail("a crossing lies on two section curves");
                    used[k] = 1;
                    on.push_back({sg, static_cast<int>(k)});
                }
            }
            std::sort(on.begin(), on.end());
            if (curve.closed && on.empty()) {
                double u = 0.0, v = 0.0;
                coords(0.5, u, v);
                if (winding(poly, u, v) != 0) return fail("a section curve is closed inside a face; not supported yet");
                continue;
            }
            // Pieces between consecutive crossings; inside ones become edges.
            std::vector<double> cuts{0.0};
            for (const auto& o : on) cuts.push_back(o.first);
            cuts.push_back(1.0);
            for (std::size_t k = 0; k + 1 < cuts.size(); ++k) {
                const double s0 = cuts[k], s1 = cuts[k + 1];
                // A crossing at an end of the curve (the face is its whole
                // patch there) leaves an empty piece.
                if ((k == 0 || k + 2 == cuts.size()) && s1 - s0 < 1e-9) continue;
                double u = 0.0, v = 0.0;
                coords(0.5 * (s0 + s1), u, v);
                const bool inside = winding(poly, u, v) != 0;
                const bool ends_at_crossings = k > 0 && k + 1 < cuts.size() - 1;
                if (!inside) continue;
                if (!ends_at_crossings) return fail("a section curve leaves a face away from its edges");
                const int va = fps[static_cast<std::size_t>(on[k - 1].second)].vertex;
                const int vb = fps[static_cast<std::size_t>(on[k].second)].vertex;
                // Direction: the kept side on the left in the face's
                // parameters, i.e. (Su x Sv) x T towards -n.
                double um, vm;
                coords(0.5 * (s0 + s1), um, vm);
                double u2, v2;
                coords(0.5 * (s0 + s1) + 1e-6 * (s1 - s0), u2, v2);
                Vec3d P, Pu, Pv;
                surf[static_cast<std::size_t>(S)]->eval(um, vm, P, Pu, Pv);
                const Vec3d T = Pu * (u2 - um) + Pv * (v2 - vm);
                const bool forward = dot(cross(cross(Pu, Pv), T), n) < 0.0;
                auto piece = std::static_pointer_cast<const ProcEdge>(te->part(forward ? s0 : s1, forward ? s1 : s0));
                // te spans sigma in [0, 1] (from = 0, to = 1): part maps directly.
                BEdge ne;
                ne.v0 = forward ? va : vb;
                ne.v1 = forward ? vb : va;
                ne.proc = piece;
                const int id = static_cast<int>(out.edges.size());
                out.edges.push_back(ne);
                section_from[ne.v0].push_back(id);
                section_to[ne.v1].push_back(id);
                face_section[ne.v0] = id;
                section_edges.push_back(id);
            }
        }
        for (char u : used)
            if (!u) return fail("a crossing of a face's edge is on no section curve");

        // The face's kept region: walk the loop's kept pieces; at a crossing
        // where the boundary leaves the kept side, follow the section edge
        // that starts there.
        struct Step {
            int edge;
            bool reversed;
            int from, to;
        };
        std::vector<Step> kept;
        for (const BCoedge& ce : face.loops[0]) {
            auto ps = pieces_of[static_cast<std::size_t>(ce.edge)];
            if (ce.reversed) std::reverse(ps.begin(), ps.end());
            for (const Piece& p : ps) {
                const BEdge& ne = out.edges[static_cast<std::size_t>(p.edge)];
                kept.push_back({p.edge, ce.reversed, ce.reversed ? ne.v1 : ne.v0, ce.reversed ? ne.v0 : ne.v1});
            }
        }
        std::map<int, int> starting; // vertex -> index into kept
        for (std::size_t i = 0; i < kept.size(); ++i) starting[kept[i].from] = static_cast<int>(i);
        std::vector<char> done(kept.size(), 0);
        for (std::size_t i0 = 0; i0 < kept.size(); ++i0) {
            if (done[i0]) continue;
            std::vector<BCoedge> L;
            int i = static_cast<int>(i0);
            for (int guard = 0; guard < 100000; ++guard) {
                if (done[static_cast<std::size_t>(i)]) break;
                done[static_cast<std::size_t>(i)] = 1;
                const Step& st = kept[static_cast<std::size_t>(i)];
                L.push_back({st.edge, st.reversed});
                int at = st.to;
                // Through any section edges of this face starting at `at`.
                for (int hop = 0; hop < 1000 && !starting.count(at); ++hop) {
                    const auto it = face_section.find(at);
                    if (it == face_section.end()) return fail("the kept boundary of a face does not close");
                    L.push_back({it->second, false});
                    at = out.edges[static_cast<std::size_t>(it->second)].v1;
                }
                i = starting[at];
            }
            if (L.empty()) continue;
            BFace nf;
            nf.carrier = face.carrier;
            nf.proc = face.proc;
            nf.flipped = face.flipped;
            nf.loops.push_back(std::move(L));
            out.faces.push_back(std::move(nf));
        }
    }
    res.section_curves = static_cast<int>(section_edges.size());
    if (section_edges.empty()) return fail("the plane does not cut the solid");

    // The cap: section edges reversed (each is used forward by its face),
    // chained at the crossing vertices.
    std::vector<char> capped(out.edges.size(), 0);
    std::vector<std::vector<BCoedge>> cap_loops;
    std::vector<double> cap_area;
    std::vector<std::vector<std::array<double, 2>>> cap_poly;
    for (int s0 : section_edges) {
        if (capped[static_cast<std::size_t>(s0)]) continue;
        std::vector<BCoedge> L;
        std::vector<std::array<double, 2>> poly;
        int se = s0;
        for (int guard = 0; guard < 100000; ++guard) {
            capped[static_cast<std::size_t>(se)] = 1;
            L.push_back({se, true});
            const BEdge& E = out.edges[static_cast<std::size_t>(se)];
            for (int k = 0; k < 16; ++k) {
                double u, v, du, dv;
                E.proc->on_carrier(kCap, 1.0 - k / 16.0, u, v, du, dv);
                poly.push_back({u, v});
            }
            const int at = E.v0; // reversed: we arrive at v0
            const auto& nxt = section_to[at];
            if (nxt.size() != 1) return fail("the section curves do not chain into loops");
            se = nxt[0];
            if (se == s0) break;
            if (capped[static_cast<std::size_t>(se)]) return fail("the section curves do not chain into loops");
        }
        cap_area.push_back(signed_area(poly));
        cap_poly.push_back(std::move(poly));
        cap_loops.push_back(std::move(L));
    }
    // Outer loops run counter-clockwise in the cap's parameters (E1 x E2 is
    // along n, the outward normal); holes clockwise, each in the outer loop
    // that contains it.
    std::vector<int> outer_of(cap_loops.size(), -1);
    for (std::size_t i = 0; i < cap_loops.size(); ++i) {
        if (cap_area[i] > 0) continue;
        for (std::size_t j = 0; j < cap_loops.size(); ++j)
            if (cap_area[j] > 0 && winding(cap_poly[j], cap_poly[i][0][0], cap_poly[i][0][1]) != 0) {
                if (outer_of[i] >= 0 && std::fabs(cap_area[j]) > std::fabs(cap_area[static_cast<std::size_t>(outer_of[i])])) continue;
                outer_of[i] = static_cast<int>(j);
            }
        if (outer_of[i] < 0) return fail("a hole of the section is in no outer loop");
    }
    for (std::size_t j = 0; j < cap_loops.size(); ++j) {
        if (!(cap_area[j] > 0)) continue;
        BFace cf;
        cf.carrier = kCap;
        cf.loops.push_back(cap_loops[j]);
        for (std::size_t i = 0; i < cap_loops.size(); ++i)
            if (outer_of[i] == static_cast<int>(j)) cf.loops.push_back(cap_loops[i]);
        out.faces.push_back(std::move(cf));
    }
    res.solid = std::move(out);
    res.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t_start).count();
    return res;
}

} // namespace ak::op
