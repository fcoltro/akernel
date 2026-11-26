// ak_workbench — editing operations: sketches, rigid moves, push/pull of
// flat faces, and fillets/chamfers of picked edges on the bodies the
// workbench built itself (boxes and revolved rings).
//
// Sketches. A closed polyline in a work plane becomes a fillet_detail::Section
// (the cross-section type that op/fillet.hpp extrudes and revolves): the
// polygon is triangulated by ear clipping and each triangle is split into
// three convex quads (centroid and edge midpoints), the section's ruled
// cells; the polygon's edges, each split at its midpoint, are its sides.
// Neighbouring cells take shared points from the same arithmetic (midpoints
// are symmetric in their ends), so shared edges have identical controls and
// the extruded or revolved shell closes by construction.
//
// Push/pull. A flat face normal to a coordinate axis, bounded by straight
// edges, is extruded along its normal as a prism and united with the body
// (pull) or subtracted from it (push): an ordinary Boolean of the kernel.
#pragma once

#include "kernel_ops.hpp"

#include <map>
#include <numbers>

namespace wb {

using P2 = std::array<double, 2>;

// ----------------------------------------------------------------- sketch

inline double signed_area(const std::vector<P2>& p) {
    double a = 0;
    for (std::size_t i = 0; i < p.size(); ++i) {
        const P2& u = p[i];
        const P2& v = p[(i + 1) % p.size()];
        a += u[0] * v[1] - v[0] * u[1];
    }
    return 0.5 * a;
}

inline double cross2(const P2& o, const P2& a, const P2& b) {
    return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0]);
}

inline bool segments_cross(const P2& a, const P2& b, const P2& c, const P2& d) {
    const double d1 = cross2(a, b, c), d2 = cross2(a, b, d), d3 = cross2(c, d, a), d4 = cross2(c, d, b);
    return ((d1 > 0) != (d2 > 0)) && ((d3 > 0) != (d4 > 0)) && d1 != 0 && d2 != 0 && d3 != 0 && d4 != 0;
}

/// Check a sketch polygon: at least 3 distinct points, no repeated points,
/// no crossing edges, non-zero area.
inline std::string check_polygon(const std::vector<P2>& p) {
    if (p.size() < 3) return "a closed sketch needs at least 3 points";
    for (std::size_t i = 0; i < p.size(); ++i)
        for (std::size_t j = i + 1; j < p.size(); ++j)
            if (p[i] == p[j]) return "the sketch has a repeated point";
    const std::size_t n = p.size();
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = i + 1; j < n; ++j) {
            if (j == i + 1 || (i == 0 && j == n - 1)) continue;
            if (segments_cross(p[i], p[(i + 1) % n], p[j], p[(j + 1) % n])) return "the sketch crosses itself";
        }
    if (!(std::fabs(signed_area(p)) > 0)) return "the sketch has no area";
    return {};
}

/// Ear-clipping triangulation of a simple counter-clockwise polygon.
inline std::vector<std::array<int, 3>> triangulate(const std::vector<P2>& p) {
    std::vector<int> idx(p.size());
    for (std::size_t i = 0; i < p.size(); ++i) idx[i] = static_cast<int>(i);
    std::vector<std::array<int, 3>> tris;
    auto P = [&](int i) -> const P2& { return p[static_cast<std::size_t>(i)]; };
    std::size_t guard = 0;
    while (idx.size() > 3 && guard++ < 10000) {
        bool clipped = false;
        const std::size_t m = idx.size();
        for (std::size_t k = 0; k < m; ++k) {
            const int a = idx[(k + m - 1) % m], b = idx[k], c = idx[(k + 1) % m];
            if (!(cross2(P(a), P(b), P(c)) > 0)) continue; // reflex or flat
            bool inside = false;
            for (int q : idx) {
                if (q == a || q == b || q == c) continue;
                if (cross2(P(a), P(b), P(q)) >= 0 && cross2(P(b), P(c), P(q)) >= 0 && cross2(P(c), P(a), P(q)) >= 0) {
                    inside = true;
                    break;
                }
            }
            if (inside) continue;
            tris.push_back({a, b, c});
            idx.erase(idx.begin() + static_cast<std::ptrdiff_t>(k));
            clipped = true;
            break;
        }
        if (!clipped) return {};
    }
    if (idx.size() == 3) tris.push_back({idx[0], idx[1], idx[2]});
    return tris;
}

/// The section of a simple polygon (see the file comment). Points in section
/// coordinates; the polygon is made counter-clockwise first.
inline std::optional<op::fillet_detail::Section> polygon_section(std::vector<P2> p, std::string& error) {
    using namespace op::fillet_detail;
    error = check_polygon(p);
    if (!error.empty()) return std::nullopt;
    if (signed_area(p) < 0) std::reverse(p.begin(), p.end());
    // Drop collinear points (they would give degenerate cells).
    for (std::size_t i = 0; i < p.size() && p.size() > 3;) {
        const std::size_t n = p.size();
        if (cross2(p[(i + n - 1) % n], p[i], p[(i + 1) % n]) == 0) p.erase(p.begin() + static_cast<std::ptrdiff_t>(i));
        else ++i;
    }
    const auto tris = triangulate(p);
    if (tris.empty()) {
        error = "could not triangulate the sketch";
        return std::nullopt;
    }
    auto pt = [](const P2& a) { return spt(a[0], a[1], 1.0); };
    auto mid2 = [](const P2& a, const P2& b) { return P2{(a[0] + b[0]) * 0.5, (a[1] + b[1]) * 0.5}; };
    Section S;
    for (const auto& t : tris) {
        const P2 &a = p[static_cast<std::size_t>(t[0])], &b = p[static_cast<std::size_t>(t[1])], &c = p[static_cast<std::size_t>(t[2])];
        const P2 g{(a[0] + b[0] + c[0]) / 3.0, (a[1] + b[1] + c[1]) / 3.0};
        const P2 mab = mid2(a, b), mbc = mid2(b, c), mca = mid2(c, a);
        // Quads (q0 q1 q2 q3, counter-clockwise) as cells ruled between q0-q1 and q3-q2.
        for (const std::array<P2, 4>& q : {std::array<P2, 4>{a, mab, g, mca}, std::array<P2, 4>{b, mbc, g, mab},
                                            std::array<P2, 4>{c, mca, g, mbc}})
            S.cells.push_back({{pt(q[0]), pt(q[1])}, {pt(q[3]), pt(q[2])}});
    }
    for (std::size_t i = 0; i < p.size(); ++i) {
        const P2& a = p[i];
        const P2& b = p[(i + 1) % p.size()];
        const P2 m = mid2(a, b);
        const double dx = b[0] - a[0], dy = b[1] - a[1], l = std::hypot(dx, dy);
        S.sides.push_back({{pt(a), pt(m)}, dy / l, -dx / l});
        S.sides.push_back({{pt(m), pt(b)}, dy / l, -dx / l});
    }
    return S;
}

/// Work planes for sketches.
enum class Plane { XY, XZ };

/// Extrude a closed sketch. XY: points (x, y) at height `level`, extruded by
/// `distance` along +z (or -z if negative). XZ: points (x, z) in the plane
/// y = level, extruded along -y (or +y if negative).
inline Outcome extrude_sketch(const std::vector<P2>& pts, Plane plane, double level, double distance) {
    std::string err;
    const auto S = polygon_section(pts, err);
    if (!S) return refused("Sketch: " + err);
    if (!(distance != 0)) return refused("Extrude sketch: the distance must not be zero");
    const double c = std::fabs(distance);
    Vec3d o, X, Y, Z;
    if (plane == Plane::XY) {
        X = {1, 0, 0}, Y = {0, 1, 0}, Z = {0, 0, 1};
        o = {0, 0, distance > 0 ? level : level - c};
    } else {
        X = {1, 0, 0}, Y = {0, 0, 1}, Z = {0, -1, 0};
        o = {0, distance > 0 ? level : level + c, 0};
    }
    Outcome r = from_blended(fmt("Extruded sketch (%d points)", static_cast<int>(pts.size())),
                             op::fillet_detail::extrude_section(*S, o, X, Y, Z, c));
    return r;
}

/// Revolve a closed sketch drawn in the XZ plane (x = radius > 0) about the
/// z axis by quarter turns.
inline Outcome revolve_sketch(const std::vector<P2>& pts, int quarters) {
    std::string err;
    const auto S = polygon_section(pts, err);
    if (!S) return refused("Sketch: " + err);
    for (const P2& q : pts)
        if (!(q[0] > 0)) return refused("Revolve sketch: draw it in the XZ plane with x > 0 (the axis is z)");
    if (quarters < 1 || quarters > 4) return refused("Revolve sketch: 1 to 4 quarter turns");
    return from_blended(fmt("Revolved sketch (%d points, %d/4 turn)", static_cast<int>(pts.size()), quarters),
                        op::fillet_detail::revolve_section(*S, 0.0, 0.0, 1.0, 1.0, quarters));
}

// ------------------------------------------------------------- rigid moves

/// Rotation by `deg` degrees about the axis (0, 1, 2 = x, y, z) through
/// `pivot`. Multiples of 90 degrees use the exact signed permutations (any
/// solid); other angles rotate the control points of a solid whose edges are
/// all patch sides (no intersection or procedural edges), keeping every shared
/// side identical so the shell stays closed.
inline Outcome rotate(const Body& A, int axis, double deg, const Vec3d& pivot) {
    if (!A.brep) return refused("Rotate works on solids");
    const double q = deg / 90.0;
    const std::string what = fmt("rotated %g deg about %c", deg, "xyz"[axis]);
    if (q == std::floor(q)) {
        const int turns = static_cast<int>(((static_cast<long long>(q) % 4) + 4) % 4);
        op::SignedPermutation one;
        if (axis == 0) { one.perm = {0, 2, 1}; one.sign = {1, -1, 1}; }
        else if (axis == 1) { one.perm = {2, 1, 0}; one.sign = {1, 1, -1}; }
        else { one.perm = {1, 0, 2}; one.sign = {-1, 1, 1}; }
        one.shift = pivot - one.linear(pivot); // about the pivot (exact for a dyadic pivot)
        op::Brep R = *A.brep;
        for (int k = 0; k < turns; ++k) R = op::transformed(R, one);
        return made(make_solid_body(A.name, std::move(R)), A.name + ", " + what);
    }
    for (const auto& e : A.brep->edges)
        if (e.arc || e.proc)
            return refused("Rotate by " + fmt("%g", deg) +
                           " deg: this solid has intersection or blend edges; use multiples of 90 deg (exact)");
    if (!A.brep->procs.empty()) return refused("Rotate: solids with procedural faces turn by multiples of 90 deg only");
    const double t = deg * std::numbers::pi / 180.0, c = std::cos(t), s = std::sin(t);
    auto rot = [&](double x, double y, double z, double w) {
        // Rotate (x, y, z) - w pivot, then add w pivot back (homogeneous).
        const double px = x - w * pivot.x, py = y - w * pivot.y, pz = z - w * pivot.z;
        double rx = px, ry = py, rz = pz;
        if (axis == 0) { ry = c * py - s * pz; rz = s * py + c * pz; }
        else if (axis == 1) { rz = c * pz - s * px; rx = s * pz + c * px; }
        else { rx = c * px - s * py; ry = s * px + c * py; }
        return Vec4d{rx + w * pivot.x, ry + w * pivot.y, rz + w * pivot.z, w};
    };
    op::Brep R = *A.brep;
    for (auto& C : R.carriers)
        for (Vec4d& h : C.c) h = rot(h.x, h.y, h.z, h.w);
    for (auto& v : R.vertices) {
        const Vec4d p = rot(v.point.x, v.point.y, v.point.z, 1.0);
        v.point = {p.x, p.y, p.z};
        if (v.exact) v.xp = rot(v.xp.x, v.xp.y, v.xp.z, v.xp.w);
        v.exact = false; // rounded now
    }
    R.outward_known = false;
    Outcome o = made(make_solid_body(A.name, std::move(R)), A.name + ", " + what);
    if (o.body) o.body->info.push_back("Rotated by a non-right angle: control points and vertices rounded to doubles");
    return o;
}

// --------------------------------------------------------------- picking info

/// Points of one Brep edge in the kernel mesh.
inline std::vector<Vec3d> edge_points(const Body& b, int edge) {
    std::vector<Vec3d> r;
    if (edge < 0 || edge >= static_cast<int>(b.tess.edges.size())) return r;
    for (int i : b.tess.edges[static_cast<std::size_t>(edge)]) r.push_back(b.tess.points[static_cast<std::size_t>(i)]);
    return r;
}

inline std::string face_info(const Body& b, int f) {
    if (!b.brep || f < 0 || f >= static_cast<int>(b.brep->faces.size())) return "";
    const op::BFace& F = b.brep->faces[static_cast<std::size_t>(f)];
    std::string s = fmt("Face %d: ", f);
    if (F.proc >= 0) s += "procedural (blend) face";
    else {
        const op::PolyPatchD& P = b.brep->carriers[static_cast<std::size_t>(F.carrier)];
        bool rational = false;
        for (const Vec4d& h : P.c) rational = rational || h.w != P.c[0].w;
        s += fmt("patch of degree (%d, %d)%s", P.p, P.q, rational ? ", rational" : "");
    }
    int edges = 0;
    for (const auto& L : F.loops) edges += static_cast<int>(L.size());
    s += fmt(", %d loop(s), %d edges", static_cast<int>(F.loops.size()), edges);
    return s;
}

inline std::string edge_info(const Body& b, int e) {
    if (!b.brep || e < 0 || e >= static_cast<int>(b.brep->edges.size())) return "";
    const op::BEdge& E = b.brep->edges[static_cast<std::size_t>(e)];
    const auto pts = edge_points(b, e);
    double len = 0;
    for (std::size_t i = 1; i < pts.size(); ++i) len += norm(pts[i] - pts[i - 1]);
    const char* kind = E.proc ? "procedural (traced) edge" : E.arc ? "intersection-curve edge" : "patch side";
    std::string s = fmt("Edge %d: %s, length %.6g", e, kind, len);
    if (pts.size() >= 2)
        s += fmt(", from (%.4g, %.4g, %.4g) to (%.4g, %.4g, %.4g)", pts.front().x, pts.front().y, pts.front().z,
                 pts.back().x, pts.back().y, pts.back().z);
    return s;
}

// ---------------------------------------------------------------- push/pull

/// Push (distance < 0) or pull (distance > 0) a flat face that is normal to
/// a coordinate axis and bounded by straight edges.
inline Outcome push_pull(const Body& A, int face, double distance) {
    if (!A.brep || A.tess.tris.empty()) return refused("Push/pull works on solids with a mesh");
    if (face < 0 || face >= static_cast<int>(A.brep->faces.size())) return refused("Push/pull: pick a face");
    if (!(distance != 0)) return refused("Push/pull: the distance must not be zero");
    const op::BFace& F = A.brep->faces[static_cast<std::size_t>(face)];
    if (F.loops.size() != 1) return refused("Push/pull: faces with holes are not supported yet");
    // The face's plane from its mesh: one coordinate constant over all its points.
    const op::Mesh& m = A.tess;
    int axis = -1;
    double level = 0;
    Vec3d nrm{0, 0, 0};
    for (std::size_t t = 0; t < m.tris.size(); ++t) {
        if (m.tri_face[t] != face) continue;
        const auto& T = m.tris[t];
        const Vec3d a = m.points[static_cast<std::size_t>(T[0])], b = m.points[static_cast<std::size_t>(T[1])],
                    c = m.points[static_cast<std::size_t>(T[2])];
        nrm += cross(b - a, c - a);
    }
    const double nl = norm(nrm);
    if (!(nl > 0)) return refused("Push/pull: the face has no area in the mesh");
    nrm = nrm / nl;
    const double comp[3] = {nrm.x, nrm.y, nrm.z};
    for (int k = 0; k < 3; ++k)
        if (std::fabs(comp[k]) > 1 - 1e-12) axis = k;
    if (axis < 0) return refused("Push/pull: only flat faces normal to x, y or z for now");
    // The boundary polygon: loop edges in order, each a straight segment.
    std::vector<Vec3d> poly;
    for (const op::BCoedge& ce : F.loops[0]) {
        auto pts = edge_points(A, ce.edge);
        if (pts.size() < 2) return refused("Push/pull: an edge of the face has no mesh points");
        if (ce.reversed) std::reverse(pts.begin(), pts.end());
        const Vec3d d = pts.back() - pts.front();
        for (const Vec3d& q : pts)
            if (norm(cross(q - pts.front(), d)) > 1e-12 * norm(d) * norm(d) + 1e-15)
                return refused("Push/pull: the face has a curved edge; only straight-edged faces for now");
        poly.push_back(pts.front());
    }
    const double* lv = axis == 0 ? &poly[0].x : axis == 1 ? &poly[0].y : &poly[0].z;
    level = *lv;
    // Exact coordinates: every vertex must lie on the plane exactly.
    for (const Vec3d& q : poly)
        if ((axis == 0 ? q.x : axis == 1 ? q.y : q.z) != level) return refused("Push/pull: the face is not exactly axis-aligned");
    // Face plane coordinates (u, v) and the prism frame.
    const int ua = (axis + 1) % 3, va = (axis + 2) % 3;
    auto coord = [](const Vec3d& q, int k) { return k == 0 ? q.x : k == 1 ? q.y : q.z; };
    auto unit = [](int k) { return Vec3d{k == 0 ? 1.0 : 0.0, k == 1 ? 1.0 : 0.0, k == 2 ? 1.0 : 0.0}; };
    std::vector<P2> pts;
    for (const Vec3d& q : poly) pts.push_back({coord(q, ua), coord(q, va)});
    std::string err;
    const auto S = polygon_section(pts, err);
    if (!S) return refused("Push/pull: " + err);
    const double sgn = comp[axis] > 0 ? 1.0 : -1.0; // outward normal along +axis or -axis
    const bool pull = distance > 0;
    const double c = std::fabs(distance);
    // Prism from the face outward (pull) or inward (push).
    const double start = pull ? level + (sgn > 0 ? 0.0 : -c) : level + (sgn > 0 ? -c : 0.0);
    Vec3d o{0, 0, 0};
    (axis == 0 ? o.x : axis == 1 ? o.y : o.z) = start;
    const op::BlendedSolid prism = op::fillet_detail::extrude_section(*S, o, unit(ua), unit(va), unit(axis), c);
    if (!prism.solid) return refused("Push/pull: prism: " + prism.error);
    Body P = make_solid_body("prism", op::to_brep(*prism.solid));
    Outcome r = boolean(A, P, pull ? 0 : 1);
    if (r.ok && r.body) r.body->name = A.name;
    if (r.ok) r.message = fmt("%s: %s face %d by %g", A.name.c_str(), pull ? "pulled" : "pushed", face, c);
    return r;
}

// ------------------------------------------------- fillet a picked edge

/// Fillet or chamfer one edge of a body the workbench generated (a box or a
/// plain revolved ring), by regenerating it with the edge blended.
inline Outcome blend_picked_edge(const Body& A, int edge, Corner c, double r, double d2) {
    const Recipe& rc = A.recipe;
    if (c == Corner::None) return refused("Choose fillet or chamfer");
    const auto pts = edge_points(A, edge);
    if (pts.size() < 2) return refused("Pick an edge");
    const Vec3d p = pts.front(), q = pts.back();
    if (rc.kind == Recipe::Box) {
        const Vec3d d = q - p;
        const double dc[3] = {std::fabs(d.x), std::fabs(d.y), std::fabs(d.z)};
        const int k = dc[0] > dc[1] && dc[0] > dc[2] ? 0 : dc[1] > dc[2] ? 1 : 2; // edge direction
        const int i = (k + 1) % 3, j = (k + 2) % 3;
        auto C = [](const Vec3d& v, int a) { return a == 0 ? v.x : a == 1 ? v.y : v.z; };
        auto E = [](int a, double len) { return Vec3d{a == 0 ? len : 0, a == 1 ? len : 0, a == 2 ? len : 0}; };
        const bool hi_i = C(p, i) == C(rc.hi, i), hi_j = C(p, j) == C(rc.hi, j);
        const bool lo_i = C(p, i) == C(rc.lo, i), lo_j = C(p, j) == C(rc.lo, j);
        if (!(hi_i || lo_i) || !(hi_j || lo_j)) return refused("That edge is not an edge of the original box");
        // Frame: origin at the corner opposite the edge in the (i, j) section; ex, ey towards the edge.
        const double li = C(rc.hi, i) - C(rc.lo, i), lj = C(rc.hi, j) - C(rc.lo, j), lk = C(rc.hi, k) - C(rc.lo, k);
        Vec3d o = rc.lo;
        Vec3d ex = E(i, li), ey = E(j, lj), ez = E(k, lk);
        if (!hi_i) { (i == 0 ? o.x : i == 1 ? o.y : o.z) += li; ex = ex * -1.0; }
        if (!hi_j) { (j == 0 ? o.x : j == 1 ? o.y : o.z) += lj; ey = ey * -1.0; }
        if (!(dot(cross(ex, ey), ez) > 0)) { // right-handed: swap the section axes
            std::swap(ex, ey);
            std::swap(r, d2);
        }
        Outcome res = from_blended(A.name + (c == Corner::Fillet ? ", edge filleted" : ", edge chamfered"),
                                   op::blend_box_edge(o, ex, ey, ez, edge_cut(c, r, d2)));
        return res;
    }
    if (rc.kind == Recipe::Ring) {
        const double rho = std::hypot(p.x, p.y), z = p.z;
        const double tol = 1e-9 * (1 + rc.r1);
        const bool outer = std::fabs(rho - rc.r1) < tol, inner = std::fabs(rho - rc.r0) < tol;
        const bool top = std::fabs(z - rc.z1) < tol, bottom = std::fabs(z - rc.z0) < tol;
        if (!(outer || inner) || !(top || bottom) || std::fabs(q.z - z) > tol ||
            std::fabs(std::hypot(q.x, q.y) - rho) > tol)
            return refused("That edge is not one of the ring's circular edges");
        const int e = outer ? (top ? 0 : 1) : (top ? 2 : 3);
        return revolve_rect(rc.r0, rc.r1, rc.z0, rc.z1, rc.quarters, e, c, r, d2);
    }
    return refused("Filleting a picked edge works on boxes and plain revolved rings made here; other solids need "
                   "op::round_edge (PR #13)");
}

} // namespace wb
