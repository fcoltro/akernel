// ak_workbench — the kernel calls behind the buttons.
//
// Every tool is a function that returns an Outcome: either a new Body (a
// solid as an op::Brep, or an open NURBS surface) with the kernel's own
// checks already run on it, or the kernel's refusal message. Nothing here
// knows about windows or OpenGL; main.cpp runs these functions on a worker
// thread and draws the RenderMesh they produce.
#pragma once

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <optional>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <akernel/geom/sweep.hpp>
#include <akernel/op/boolean.hpp>
#include <akernel/op/brep_transform.hpp>
#include <akernel/op/certified_volume.hpp>
#include <akernel/op/fillet.hpp>
#include <akernel/op/plane_cut.hpp>
#include <akernel/op/poly_boolean.hpp>
#include <akernel/op/proc_fillet.hpp>
#include <akernel/op/surfaces.hpp>
#include <akernel/op/tessellate.hpp>

namespace wb {

using ak::geom::Vec3d;
using ak::geom::Vec4d;
namespace op = ak::op;
namespace geom = ak::geom;

/// Triangles (position + normal per vertex) and edge segments, in floats for OpenGL.
struct RenderMesh {
    std::vector<float> tri;   ///< x y z nx ny nz, three vertices per triangle
    std::vector<float> lines; ///< x y z, two vertices per segment
    std::vector<int> tri_face;  ///< Brep face of each triangle (-1 for surfaces)
    std::vector<int> line_edge; ///< Brep edge of each segment (-1 for surface boundaries)
    std::array<float, 3> lo{0, 0, 0}, hi{0, 0, 0};
    bool empty() const { return tri.empty() && lines.empty(); }
};

/// What generated a body, so that a picked edge can be blended by
/// regenerating it (op/fillet.hpp blends edges of swept sections).
struct Recipe {
    enum Kind { None, Box, Ring } kind = None;
    Vec3d lo{}, hi{};                  ///< Box
    double r0 = 0, r1 = 0, z0 = 0, z1 = 0; ///< Ring (revolved rectangle)
    int quarters = 4;
};

struct Body {
    int id = 0;                     ///< stable id shown in the UI and used by the MCP tools
    std::string name;
    bool solid = true;              ///< false: an open surface (sweep, loft), not a Brep
    std::optional<op::Brep> brep;   ///< set for solids
    op::Mesh tess;                  ///< the kernel's watertight mesh (solids), for OBJ export
    RenderMesh mesh;
    std::vector<std::string> info;  ///< the inspector's lines: what the kernel certified or measured
    bool visible = true;
    std::array<float, 3> color{0.7f, 0.75f, 0.8f};
    Recipe recipe;                  ///< how it was made, for edits that regenerate it
};

struct Outcome {
    bool ok = false;
    std::string message; ///< a summary, or the kernel's reason for refusing
    std::optional<Body> body;
    double seconds = 0.0;
};

/// printf into a std::string.
template <class... A>
std::string fmt(const char* f, A... a) {
    char buf[512];
    std::snprintf(buf, sizeof buf, f, a...);
    return buf;
}

inline void grow(RenderMesh& m, float x, float y, float z) {
    const float p[3] = {x, y, z};
    const bool first = m.tri.empty() && m.lines.empty();
    for (std::size_t k = 0; k < 3; ++k) {
        m.lo[k] = first ? p[k] : std::min(m.lo[k], p[k]);
        m.hi[k] = first ? p[k] : std::max(m.hi[k], p[k]);
    }
}

inline void push3(std::vector<float>& v, const Vec3d& a) {
    v.push_back(static_cast<float>(a.x));
    v.push_back(static_cast<float>(a.y));
    v.push_back(static_cast<float>(a.z));
}

/// Smooth normals within each Brep face, sharp across edges: a point shared
/// by two faces gets one normal per face.
inline RenderMesh render_mesh(const op::Mesh& m) {
    RenderMesh r;
    std::unordered_map<std::uint64_t, Vec3d> acc;
    auto key = [](int point, int face) {
        return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(point)) << 32) |
               static_cast<std::uint32_t>(face);
    };
    for (std::size_t t = 0; t < m.tris.size(); ++t) {
        const auto& T = m.tris[t];
        const Vec3d& a = m.points[static_cast<std::size_t>(T[0])];
        const Vec3d& b = m.points[static_cast<std::size_t>(T[1])];
        const Vec3d& c = m.points[static_cast<std::size_t>(T[2])];
        const Vec3d n = cross(b - a, c - a); // area weighted
        for (int k = 0; k < 3; ++k) acc[key(T[static_cast<std::size_t>(k)], m.tri_face[t])] += n;
    }
    for (std::size_t t = 0; t < m.tris.size(); ++t)
        for (int k = 0; k < 3; ++k) {
            const int i = m.tris[t][static_cast<std::size_t>(k)];
            const Vec3d& p = m.points[static_cast<std::size_t>(i)];
            Vec3d n = acc[key(i, m.tri_face[t])];
            const double l = norm(n);
            if (l > 0.0) n = n / l;
            grow(r, static_cast<float>(p.x), static_cast<float>(p.y), static_cast<float>(p.z));
            push3(r.tri, p);
            push3(r.tri, n);
            if (k == 0) r.tri_face.push_back(m.tri_face[t]);
        }
    for (std::size_t ei = 0; ei < m.edges.size(); ++ei) {
        const auto& e = m.edges[ei];
        for (std::size_t k = 1; k < e.size(); ++k) {
            push3(r.lines, m.points[static_cast<std::size_t>(e[k - 1])]);
            push3(r.lines, m.points[static_cast<std::size_t>(e[k])]);
            r.line_edge.push_back(static_cast<int>(ei));
        }
    }
    return r;
}

/// Size of a Brep for choosing the tessellation chord: the diagonal of its
/// vertices' box (every closed Brep made here has vertices).
inline double brep_size(const op::Brep& b) {
    Vec3d lo{1e300, 1e300, 1e300}, hi{-1e300, -1e300, -1e300};
    for (const auto& v : b.vertices) {
        lo = {std::min(lo.x, v.point.x), std::min(lo.y, v.point.y), std::min(lo.z, v.point.z)};
        hi = {std::max(hi.x, v.point.x), std::max(hi.y, v.point.y), std::max(hi.z, v.point.z)};
    }
    return b.vertices.empty() ? 1.0 : std::max(norm(hi - lo), 1e-6);
}

/// Run the kernel's checks on a solid and tessellate it. A solid that does
/// not tessellate is still kept (its info says why), but is not drawn.
inline Body make_solid_body(std::string name, op::Brep b) {
    Body body;
    body.name = std::move(name);
    body.solid = true;
    auto& info = body.info;
    info.push_back(fmt("%d faces, %d edges, %d vertices, %d procedural carriers", static_cast<int>(b.faces.size()),
                       static_cast<int>(b.edges.size()), static_cast<int>(b.vertices.size()),
                       static_cast<int>(b.procs.size())));
    const auto rep = op::audit(b);
    int genus = -1;
    for (const auto& s : rep.shells) genus = std::max(genus, s.genus);
    info.push_back(fmt("Topology audit: %s, %d shell(s), Euler characteristic %d%s",
                       rep.closed() ? "closed and oriented" : "NOT closed", static_cast<int>(rep.shells.size()),
                       rep.euler, genus >= 0 ? fmt(", genus %d", genus).c_str() : ""));
    info.push_back(fmt("Volume (quadrature): %.12g", op::volume(b)));
    std::string why;
    if (const auto cv = op::certified_volume(b, &why))
        info.push_back(fmt("Certified volume in [%.12g, %.12g], width %.1e", cv->lo, cv->hi, cv->hi - cv->lo));
    else
        info.push_back("Certified volume: not covered (" + why + ")");
    const double vw = op::vertex_width(b);
    info.push_back(fmt("Certified edge gap bound: %.1e; vertices: %s", op::edge_gap_bound(b),
                       vw < 1e-300 ? "all exact" : fmt("widest enclosure %.1e", vw).c_str()));
    op::TessOptions to;
    to.chord = 1e-3 * brep_size(b);
    const op::TessResult t = op::tessellate(b, to);
    if (t.ok) {
        const op::MeshAudit ma = op::audit(t.mesh);
        info.push_back(fmt("Mesh: %d triangles, %s, chord %.1e, %.0f ms", ma.triangles,
                           ma.watertight() ? "watertight" : "NOT watertight", to.chord, t.stats.ms));
        body.tess = t.mesh;
        body.mesh = render_mesh(t.mesh);
    } else {
        info.push_back("Mesh: not tessellated (" + t.error + ")");
    }
    body.brep = std::move(b);
    return body;
}

/// An open NURBS surface sampled on a grid (sweeps and lofts are surfaces in
/// the kernel today: geom/sweep.hpp).
inline Body make_surface_body(std::string name, const geom::NurbsSurface& S, int n = 64) {
    Body body;
    body.name = std::move(name);
    body.solid = false;
    body.color = {0.85f, 0.7f, 0.45f};
    const double u0 = S.U[static_cast<std::size_t>(S.p)], u1 = S.U[static_cast<std::size_t>(S.nu)];
    const double v0 = S.V[static_cast<std::size_t>(S.q)], v1 = S.V[static_cast<std::size_t>(S.nv)];
    std::vector<Vec3d> P, N;
    for (int i = 0; i <= n; ++i)
        for (int j = 0; j <= n; ++j) {
            const double u = u0 + (u1 - u0) * i / n, v = v0 + (v1 - v0) * j / n;
            const geom::Jet3 J = S.jet(u, v);
            Vec3d nn = cross(J[1], J[2]);
            const double l = norm(nn);
            P.push_back(J[0]);
            N.push_back(l > 0.0 ? nn / l : Vec3d{0, 0, 1});
        }
    auto id = [&](int i, int j) { return static_cast<std::size_t>(i * (n + 1) + j); };
    auto vert = [&](std::size_t k) {
        grow(body.mesh, static_cast<float>(P[k].x), static_cast<float>(P[k].y), static_cast<float>(P[k].z));
        push3(body.mesh.tri, P[k]);
        push3(body.mesh.tri, N[k]);
    };
    auto tri_done = [&] { body.mesh.tri_face.push_back(-1); };
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) {
            for (std::size_t k : {id(i, j), id(i + 1, j), id(i + 1, j + 1)}) vert(k);
            tri_done();
            for (std::size_t k : {id(i, j), id(i + 1, j + 1), id(i, j + 1)}) vert(k);
            tri_done();
        }
    // Boundary curves as edges.
    for (int k = 0; k < n; ++k)
        for (auto [a, b] : {std::pair{id(k, 0), id(k + 1, 0)}, std::pair{id(k, n), id(k + 1, n)},
                            std::pair{id(0, k), id(0, k + 1)}, std::pair{id(n, k), id(n, k + 1)}}) {
            push3(body.mesh.lines, P[a]);
            push3(body.mesh.lines, P[b]);
            body.mesh.line_edge.push_back(-1);
        }
    body.info.push_back(fmt("Open NURBS surface, degree (%d, %d), %d x %d control points", S.p, S.q, S.nu, S.nv));
    auto wr = std::minmax_element(S.cp.begin(), S.cp.end(), [](const Vec4d& a, const Vec4d& b) { return a.w < b.w; });
    body.info.push_back(fmt("Weights in [%g, %g]", wr.first->w, wr.second->w));
    body.info.push_back("Not a solid: sweeps and lofts are surfaces in the kernel today (geom/sweep.hpp)");
    return body;
}

// ------------------------------------------------------------------ timing

template <class F>
Outcome timed(F&& f) {
    const auto t0 = std::chrono::steady_clock::now();
    Outcome o = f();
    o.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    return o;
}

inline Outcome refused(std::string why) {
    Outcome o;
    o.message = std::move(why);
    return o;
}

inline Outcome made(Body b, std::string summary = {}) {
    Outcome o;
    o.ok = true;
    o.message = summary.empty() ? b.name : std::move(summary);
    o.body = std::move(b);
    return o;
}

inline Outcome from_patch_solid(std::string name, const op::PatchSolid& s) {
    if (s.faces.empty()) return refused(name + ": the kernel did not build a closed shell");
    return made(make_solid_body(std::move(name), op::to_brep(s)));
}

inline Outcome from_blended(std::string name, const op::BlendedSolid& s) {
    if (!s.solid) return refused(name + ": " + s.error);
    return from_patch_solid(std::move(name), *s.solid);
}

// ------------------------------------------------------------- primitives

inline Outcome box(const Vec3d& o, const Vec3d& size) {
    if (!(size.x > 0 && size.y > 0 && size.z > 0)) return refused("Box: sizes must be positive");
    Outcome r = from_patch_solid("Box", op::make_box(o, o + size));
    if (r.body) r.body->recipe = {Recipe::Box, o, o + size};
    return r;
}

inline Outcome bump_box(const Vec3d& o, const Vec3d& size, double bump) {
    if (!(size.x > 0 && size.y > 0 && size.z > 0)) return refused("Bump box: sizes must be positive");
    return from_patch_solid("Bump box", op::make_bump_box(o, size.x, size.y, size.z, bump));
}

/// A box whose top is a quadratic B-spline with m x n spans: a wavy freeform top.
inline Outcome spline_box(const Vec3d& o, const Vec3d& size, int m, int n, double amp) {
    if (!(size.x > 0 && size.y > 0 && size.z > 0) || m < 1 || n < 1) return refused("Spline box: bad parameters");
    std::vector<double> h;
    for (int i = 0; i < m + 2; ++i)
        for (int j = 0; j < n + 2; ++j)
            h.push_back(amp * (((i + j) % 2 == 0) ? 1.0 : -1.0) * (i == 0 || j == 0 || i == m + 1 || j == n + 1 ? 0.5 : 1.0));
    return from_patch_solid("Spline box", op::make_spline_box(o, size.x, size.y, size.z, m, n, h));
}

inline Outcome rounded_box(const Vec3d& o, const Vec3d& size, double r) {
    return from_blended("Rounded box", op::rounded_box(o, size.x, size.y, size.z, r));
}

// --------------------------------------------------------- extrude, revolve

enum class Corner { None, Fillet, Chamfer };

inline op::EdgeCut edge_cut(Corner c, double r, double d2) {
    return c == Corner::Fillet ? op::EdgeCut::round(r) : op::EdgeCut::flat(r, d2);
}

/// The plain rectangle [0, a] x [0, b] as a section (one cell, four sides),
/// for extrusion and revolution without a blended corner.
inline op::fillet_detail::Section rectangle_section(double a, double b) {
    using namespace op::fillet_detail;
    const SPt O = spt(0, 0, 1), X = spt(a, 0, 1), XY = spt(a, b, 1), Y = spt(0, b, 1);
    Section S;
    S.cells = {{{O, X}, {Y, XY}}};
    S.sides = {{{O, X}, 0, -1}, {{X, XY}, 1, 0}, {{XY, Y}, 0, 1}, {{Y, O}, -1, 0}};
    return S;
}

/// Extrude the rectangle W x H (in x, y at o) by L along z; its corner at
/// (W, H) optionally rounded or chamfered along the whole length.
inline Outcome extrude_rect(const Vec3d& o, double W, double H, double L, Corner c, double r, double d2) {
    if (!(W > 0 && H > 0 && L > 0)) return refused("Extrude: sizes must be positive");
    if (c == Corner::None) {
        Outcome r = from_blended("Extrusion", op::fillet_detail::extrude_section(rectangle_section(W, H), o, {1, 0, 0},
                                                                                 {0, 1, 0}, {0, 0, 1}, L));
        if (r.body) r.body->recipe = {Recipe::Box, o, o + Vec3d{W, H, L}};
        return r;
    }
    return from_blended(c == Corner::Fillet ? "Extrusion, filleted edge" : "Extrusion, chamfered edge",
                        op::blend_box_edge(o, {W, 0, 0}, {0, H, 0}, {0, 0, L}, edge_cut(c, r, d2)));
}

/// Extrude an L-shaped (step) profile: W x H without the notch x > a1, y > b1;
/// the inside (concave) edge filled by a fillet or chamfer.
inline Outcome extrude_step(const Vec3d& o, double W, double H, double L, double a1, double b1, Corner c, double r,
                            double d2) {
    if (c == Corner::None) return refused("Step: choose a fillet or a chamfer for the inside edge");
    return from_blended(c == Corner::Fillet ? "Step, filleted inside edge" : "Step, chamfered inside edge",
                        op::blend_step_edge(o, {W, 0, 0}, {0, H, 0}, {0, 0, L}, a1, b1, edge_cut(c, r, d2)));
}

/// Revolve the rectangle r0..r1 x z0..z1 about the z axis by `quarters`
/// quarter turns, optionally with one of its four circular edges blended.
inline Outcome revolve_rect(double r0, double r1, double z0, double z1, int quarters, int edge, Corner c, double r,
                            double d2) {
    if (!(r0 > 0 && r1 > r0 && z1 > z0)) return refused("Revolve: needs 0 < r0 < r1 and z0 < z1");
    const char* turn = quarters == 4 ? "full turn" : quarters == 1 ? "quarter turn" : quarters == 2 ? "half turn" : "3/4 turn";
    if (c == Corner::None) {
        Outcome o = from_blended(fmt("Revolution (%s)", turn), op::fillet_detail::revolve_section(
                                                                   rectangle_section(r1 - r0, z1 - z0), r0, z0, 1.0, 1.0, quarters));
        if (o.body) {
            o.body->recipe.kind = Recipe::Ring;
            o.body->recipe.r0 = r0, o.body->recipe.r1 = r1, o.body->recipe.z0 = z0, o.body->recipe.z1 = z1;
            o.body->recipe.quarters = quarters;
        }
        return o;
    }
    return from_blended(fmt("Revolution (%s), %s edge", turn, c == Corner::Fillet ? "filleted" : "chamfered"),
                        op::blend_ring_edge(r0, r1, z0, z1, static_cast<op::RingEdge>(edge), edge_cut(c, r, d2), quarters));
}

// --------------------------------------------------------- sweeps, lofts

/// A circle of radius r (in the plane z = o.z, centred at o) moved along a
/// cubic Bezier path from o by the offsets c1, c2, end, without rotating.
inline Outcome sweep_circle(const Vec3d& o, double r, const Vec3d& c1, const Vec3d& c2, const Vec3d& end) {
    if (!(r > 0)) return refused("Sweep: radius must be positive");
    geom::NurbsCurve profile = geom::exact_circle(r, o.z, 4);
    for (Vec4d& c : profile.cp) c = {c.x + o.x * c.w, c.y + o.y * c.w, c.z, c.w};
    auto at = [&](const Vec3d& d) { return Vec4d{o.x + d.x, o.y + d.y, o.z + d.z, 1}; };
    const geom::NurbsCurve path = geom::bezier_curve({at({0, 0, 0}), at(c1), at(c2), at(end)});
    return made(make_surface_body("Sweep (translational)", geom::translational_sweep(profile, path)));
}

/// Loft (skin) through three circles about the vertical axis through o.
inline Outcome loft_circles(const Vec3d& o, double r0, double r1, double r2, double h) {
    if (!(r0 > 0 && r1 > 0 && r2 > 0 && h > 0)) return refused("Loft: radii and height must be positive");
    auto circle = [&](double r, double z) {
        geom::NurbsCurve c = geom::exact_circle(r, o.z + z, 4);
        for (Vec4d& q : c.cp) q = {q.x + o.x * q.w, q.y + o.y * q.w, q.z, q.w};
        return c;
    };
    const std::vector<geom::NurbsCurve> sec = {circle(r0, 0.0), circle(r1, 0.5 * h), circle(r2, h)};
    const auto s = geom::skin(sec, 2);
    if (!s) return refused("Loft: the kernel refused (singular system or a weight not positive)");
    return made(make_surface_body("Loft (skinned)", s->surface));
}

// ----------------------------------------------------------------- fillet

/// The freeform block of op/proc_fillet.hpp: top z = h + c x s(y) over
/// [0, a] x [0, L] (s the smoothstep), its edge at x = a rounded by a
/// rolling ball of radius r. The blend is a procedural face (BlendCarrier).
inline Outcome freeform_fillet(double a, double L, double h, double c, double r, double cut) {
    op::PolyPatchD T(1, 3);
    static constexpr double s[4] = {0, 0, 1, 1};
    for (int i = 0; i <= 1; ++i)
        for (int j = 0; j <= 3; ++j) T.at(i, j) = {a * i, L * j / 3.0, h + c * a * i * s[j], 1.0};
    op::ProcFilletResult res = op::round_block_edge(T, r, {}, cut);
    if (!res.solid) return refused("Freeform fillet: " + res.error);
    Body b = make_solid_body("Freeform block, filleted edge", std::move(*res.solid));
    b.info.push_back(fmt("Blend spine: %d certified segments", static_cast<int>(res.blend->segments())));
    return made(std::move(b));
}

// ----------------------------------------------------------- operations

inline bool procedural(const op::Brep& b) {
    if (!b.procs.empty()) return true;
    for (const auto& e : b.edges)
        if (e.proc) return true;
    return false;
}

/// A Boolean of two solids: op/boolean.hpp for rational solids, the
/// procedural-face Booleans (op/poly_boolean.hpp) when either has blends.
inline Outcome boolean(const Body& A, const Body& B, int which) {
    static const char* names[3] = {"Union", "Difference", "Intersection"};
    const std::string name = names[which];
    const std::string operands = "Operands: A = " + A.name + ", B = " + B.name;
    if (!A.brep || !B.brep) return refused("Booleans need two solids");
    if (procedural(*A.brep) || procedural(*B.brep)) {
        const op::PolyOp po = which == 0 ? op::PolyOp::Union : which == 1 ? op::PolyOp::Difference : op::PolyOp::Intersection;
        op::PolyBooleanResult r = op::poly_boolean(*A.brep, *B.brep, po);
        if (!r.solid) return refused(name + ": " + r.error);
        if (r.solid->faces.empty()) return refused(name + ": the result is empty");
        Body b = make_solid_body(name, std::move(*r.solid));
        b.info.push_back(operands);
        b.info.push_back(fmt("Procedural-face Boolean: %d certified intersection curve pieces", r.curves));
        return made(std::move(b));
    }
    const op::Op o = which == 0 ? op::Op::Union : which == 1 ? op::Op::Difference : op::Op::Intersection;
    const op::BooleanResult r = op::boolean_op(*A.brep, *B.brep, o);
    if (!r.ok) return refused(name + ": " + r.error);
    if (r.empty) return refused(name + ": the result is empty");
    const bool manifold = r.manifold;
    const int lc = r.line_contacts, pc = r.point_contacts;
    Body b = make_solid_body(name, r.solid);
    b.info.push_back(operands);
    b.info.push_back(manifold ? std::string("Manifold result")
                              : fmt("Non-manifold result: %d line contacts, %d point contacts", lc, pc));
    return made(std::move(b));
}

inline Outcome plane_cut(const Body& A, const Vec3d& n, double d) {
    if (!A.brep) return refused("Plane cut needs a solid");
    op::PlaneCutResult r = op::cut_by_plane(*A.brep, n, d);
    if (!r.solid) return refused("Plane cut: " + r.error);
    if (r.solid->faces.empty()) return refused("Plane cut: nothing is left on the kept side");
    Body b = make_solid_body("Cut " + A.name, std::move(*r.solid));
    b.info.push_back(fmt("Plane cut: %d certified section curve pieces", r.section_curves));
    return made(std::move(b));
}

/// Exact placement: translation, quarter turns and mirrors (op/brep_transform.hpp).
inline Outcome place(const Body& A, const op::SignedPermutation& m, const char* what) {
    if (!A.brep) return refused("Moving works on solids");
    return made(make_solid_body(*what ? A.name + ", " + what : A.name, op::transformed(*A.brep, m)));
}

} // namespace wb
