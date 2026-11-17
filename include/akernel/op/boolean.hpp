// akernel — imprinting and Boolean operations on boundary representations
// (docs/BLUEPRINT.md §3.3).
//
// Operands and results are Breps (op/brep.hpp) with trimmed faces, so a
// result can be an operand. No surface is tessellated and no tolerance is
// used at any step.
//
// Identity by construction (thesis T1). Two operands may share geometry: a
// result combined with one of its own operands shares surfaces, edges and
// vertices with it. Sharing is recognised from how the geometry was made,
// never from distance:
//   - surfaces: identical control points;
//   - sides shared by two surfaces: identical boundary control points;
//   - vertices: identical certified parameters on a surface;
//   - edge pieces: the same curve between the same two vertices.
// Where the two solids share a surface, their faces on it are overlaid in
// the one parameter square, and each region of a face is classified as
// inside, outside, or on the other solid (with the same or the opposite
// outward normal).
//
// Pipeline
//   1. Surface/surface: certified intersection curves of every pair of
//      surfaces that carry faces of different solids, proven complete. A
//      side shared by construction is certified directly as a known curve.
//   2. Edge/surface: certified points where an edge of one solid meets a
//      surface of the other. A point that is an existing vertex (inside the
//      root's uniqueness box) is that vertex. Otherwise it is new if it lies
//      inside a face (the edge pierces the face) or, within a shared
//      surface, on the boundary of one (two edges cross).
//   3. Clipping: each curve of step 1 is cut at every vertex on it. A part
//      that is already an edge of an operand is skipped; a part inside a
//      face of each solid becomes a new edge.
//   4. Splitting: every edge is split at every vertex that lies on it by
//      construction; pieces of the two solids on the same curve between the
//      same vertices are one piece.
//   5. Regions: in each face, its own boundary pieces, the new edges, and the
//      other solid's edge pieces that lie inside it (shared surfaces only)
//      form a planar map. At a vertex the next curve around a region is the
//      first one met turning clockwise, decided with interval cross
//      products.
//   6. Classification: on the other solid's boundary (shared surface), next
//      to a new edge (decided by an interval evaluation on the curve), across
//      uncut edge pieces, or by a certified ray test.
//   7. Selection: regions kept by the operation form the result Brep; the
//      integer audit decides whether it is a closed oriented shell.
//
// Every step either proves what it needs or stops with a message. A curved
// surface touching a plane along sides of its patch is decided exactly
// (op/contact.hpp): the pair contributes no curve, edges touching the plane
// at their ends add no vertex, and the curves and faces that meet
// tangentially there are ordered by the side of the plane the surface lies
// on. Other contact that is not shared by construction (tangency at a point
// or between curved surfaces, two different surfaces that coincide, an edge
// through an edge) is reported as "not in general position".
//
// Not proven: that an operand does not intersect itself.
#pragma once

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <map>
#include <set>
#include <memory>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <akernel/num/expansion.hpp>
#include <akernel/geom/box_index.hpp>
#include <akernel/op/brep.hpp>
#include <akernel/op/contact.hpp>
#include <akernel/op/tsi.hpp>
#include <akernel/solve/contact_point.hpp>
#include <akernel/solve/csi.hpp>
#include <akernel/solve/ssi.hpp>

namespace ak::op {

enum class Op { Intersection, Union, Difference };

struct BooleanOptions {
    solve::SsiOptions ssi;
    solve::CsiOptions csi;
    int refine_depth = 40; ///< bisections allowed when classifying a point against a loop
    bool trace = false;    ///< record what was found in Imprint::trace
    /// Search surface pairs only over the parameter regions that hold faces
    /// (hull of the faces' certified rectangles), not whole carriers. Off by
    /// default: it gave no measurable speed-up on the benchmark chains, and it
    /// changes how the curves are traced (the same curves, other segments).
    bool face_domains = false;
    /// Two flat surfaces in different planes: certify their one line of
    /// intersection directly and prove it complete without subdivision
    /// (falls back to the general search when that proof does not go through).
    bool plane_lines = true;
    /// A surface pair whose intersection could not be proven complete over the
    /// parameter squares is searched again over squares enlarged by this much
    /// (0 switches it off): a curve that only touches both squares at their
    /// boundaries is then traced through the touching point.
    double boundary_margin = 1.0 / 32.0;
    std::size_t boundary_node_limit = 8000; ///< completeness-pass nodes allowed for that second search
    /// Move flat faces on rational carriers (a sector's annular cap) onto
    /// exact parallelograms of their plane where they meet a flat face in the
    /// same plane that is not affinely related to them.
    bool recarry_flat = true;
    /// Prove roots with round coordinates exact and use the exact points of
    /// vertices (to decide that a vertex lies on a surface, and that two
    /// vertices are one). Off: every vertex is an enclosure, as before 0.16.
    bool exact_points = true;
    int boundary_trace_steps = 2000;        ///< tracing steps allowed per branch in that second search
    BooleanOptions() { ssi.verify_node_limit = 400000; }
};

/// Classification of a region of one solid's boundary against the other solid.
enum : int {
    kUnknown = -1,
    kOut = 0,
    kIn = 1,
    kOnSame = 2,     ///< on the other solid's boundary, outward normals agree
    kOnOpposite = 3, ///< on the other solid's boundary, outward normals opposite
};

/// Oriented boundary piece of a region.
struct Piece {
    bool fresh = false; ///< a new intersection edge (true) or a piece of an operand's edge (false)
    int id = -1;        ///< new edge index, or canonical piece index
    bool forward = true; ///< relative to the orientation of `id`
    int src = -1;       ///< operand piece whose geometry is used (old pieces only)
};

/// A cell of a face's planar map.
struct Region {
    int solid = 0, face = -1;
    std::vector<std::vector<Piece>> loops; ///< outer loop first; region on the left of every loop
    std::vector<CycleGeom> geom;           ///< chord-and-box geometry of each loop
    int inside = kUnknown;
};

struct ImprintStats {
    int carrier_pairs = 0, carrier_pairs_solved = 0, edge_carrier_pairs_solved = 0;
    int curves = 0, closed_curves = 0, segments = 0;
    int new_arcs = 0;        ///< edges created by this operation
    int new_vertices = 0;    ///< vertices created by this operation
    int triple_vertices = 0; ///< of those, three-surface points on an arc edge of an operand
    int regions = 0;
    int ray_tests = 0;
    int coarse_arcs = 0;           ///< new arcs whose joins could not be verified for refinement
    double max_vertex_width = 0.0; ///< widest certified enclosure of a new vertex, in parameter units
    // Shared by construction
    int shared_carriers = 0, shared_sides = 0, merged_vertices = 0, unified_pieces = 0, on_regions = 0;
    int coplanar_pairs = 0; ///< pairs of different surfaces proven to lie in one plane
    int plane_lines = 0;    ///< plane pairs whose line was certified directly, where tracing could not finish
    int direct_plane_pairs = 0; ///< plane pairs proven complete directly (one line), without subdivision
    int tangent_sides = 0;      ///< curved surfaces proven to touch a plane only along one of their sides (exact)
    int tangent_ends = 0;       ///< edges proven to touch a plane only at their end vertices (exact)
    int boundary_pairs = 0;     ///< surface pairs proven complete only over enlarged parameter squares
    int recarried_faces = 0;    ///< flat rational carriers whose faces were moved onto parallelograms of their plane
    int curved_coincident_pairs = 0; ///< curved surface pairs proven to be one surface (affine reparameterisation)
    int exact_roots = 0;             ///< new vertices whose position was proven exact (dyadic parameters on a surface)
    int exact_ends = 0;              ///< roots proven to be an end vertex of their edge, by the vertex's exact point
    int iso_edges = 0;               ///< straight edges proven to lie in a curved surface of the other solid, along a line of constant parameter
    int line_contacts = 0;           ///< of their pieces, those along which the other solid only touches the surface inside a face
    int edge_crossings = 0;          ///< edges proven to cross an edge of the other solid in space, at an exact point
    int line_tangencies = 0;         ///< curved surfaces proven to touch a plane along a line inside the patch
    int graph_crossings = 0;         ///< pairs of graphs crossing where they share boundary lines, desingularised exactly
    int point_contacts = 0;          ///< pairs of surfaces touching at one point (a curved one against a plane or a curved one), proven exactly
    int near_misses = 0;             ///< pairs whose facing points are proven apart, with nothing between them in the region the search left
    int edge_touches = 0;            ///< edges in a plane touching a curved surface at one point, meeting it there
    int ssi_pairs_empty = 0;          ///< intersected surface pairs that gave no curve
    double ssi_ms_empty = 0.0, ssi_ms_curves = 0.0; ///< time in surface pairs without and with curves
    std::size_t ssi_nodes_empty = 0, ssi_nodes_curves = 0; ///< completeness-pass nodes, likewise
    // Diagnostics (blueprint §2.2): where the time and the search effort went.
    int stage = 0; ///< last stage entered: 0 operands, 1 surfaces, 2 edges, 3 clipping, 4 splitting, 5 regions, 6 classification
    double ms[7] = {0, 0, 0, 0, 0, 0, 0};
    std::size_t csi_nodes = 0, tsi_nodes = 0, ssi_verify_nodes = 0;
    /// Edge/surface searches: straight edges (curve/surface) and arcs
    /// (three-surface points along a tube), their number and time.
    int edge_searches_straight = 0, edge_searches_arc = 0;
    double edge_ms_straight = 0.0, edge_ms_arc = 0.0;
};

struct Imprint {
    bool ok = false;
    std::string error;
    const Brep* solid[2] = {nullptr, nullptr}; ///< must outlive the imprint
    /// Operands that were changed before the Boolean (flat faces moved onto
    /// affine carriers); `solid` then points here.
    std::shared_ptr<const Brep> owned[2];
    std::vector<PolyPatchD> carriers;          ///< global surfaces; a surface of both operands appears once
    std::vector<int> carrier_of[2];            ///< surface index in solid k -> global surface

    struct Vertex {
        Vec3d point{};
        std::vector<Chart> charts; ///< global surface ids
        bool exact = false;        ///< `xp` is the vertex, exactly (proven when it was made)
        Vec4d xp{};                ///< homogeneous coordinates, each a double
    };
    /// Part of an operand's edge between two consecutive vertices.
    struct OldPiece {
        int solid = 0, edge = -1;
        int v0 = -1, v1 = -1;
        Interval t0, t1;                     ///< side edge: parameter range
        std::shared_ptr<const ArcGeom> geom; ///< arc edge: geometry of the piece
        int canon = -1;                      ///< this piece, or the piece that stands for it on the same curve (of solid 0, or the arc where a side and an arc coincide)
        bool canon_reversed = false;         ///< closed piece: it runs against its canonical piece
        bool coincident = false;             ///< the other solid has a piece on the same curve
    };
    /// A new edge: part of an intersection curve inside a face of each solid.
    struct NewArc {
        int carrier[2] = {-1, -1}; ///< global ids; the curve's coordinates 0,1 and 2,3
        std::shared_ptr<const ArcGeom> geom;
        int v0 = -1, v1 = -1;
        int solid_on[2] = {0, 1}; ///< solid that owns the face on carrier[i]
        int face[2] = {-1, -1};
        /// In the face on carrier[i], is the left of the travel direction
        /// inside the other solid?
        bool left_inside[2] = {false, false};
    };

    std::vector<Vertex> vertices;
    std::vector<OldPiece> pieces;
    std::vector<NewArc> arcs;
    std::vector<std::shared_ptr<const CurveBase>> curves;
    std::vector<Region> regions;
    /// A piece of an edge of one solid that lies inside a face of the other
    /// along a line where its solid only touches that face's surface. The
    /// piece bounds nothing in the face, so the face is not split there, and
    /// a result that keeps both the face's region and the edge has parts
    /// touching along that line without sharing an edge.
    struct LineContact {
        int region = -1; ///< the region of the face that contains the piece
        int piece = -1;  ///< the piece (its canonical index)
    };
    std::vector<LineContact> line_contacts;
    /// Pairs of regions, of the two solids, that touch along a line inside
    /// both (a curved face against a plane tangent to it along a ruling).
    std::vector<std::pair<int, int>> touching_regions;
    /// Pairs of regions, of the two solids, that touch at one point inside
    /// both (a curved face against a plane tangent to it at a point).
    std::vector<std::pair<int, int>> point_touching_regions;
    ImprintStats stats;
    std::vector<std::string> trace; ///< filled when BooleanOptions::trace is set
};

namespace detail {

inline Interval det3i(const Interval m[3][3]) {
    return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) -
           m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
           m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
}

inline bool same_point(const Vec4d& a, const Vec4d& b) {
    return a.x == b.x && a.y == b.y && a.z == b.z && a.w == b.w;
}

/// Identical control points: the same surface by construction.
inline bool same_patch(const PolyPatchD& a, const PolyPatchD& b) {
    if (a.p != b.p || a.q != b.q) return false;
    for (std::size_t i = 0; i < a.c.size(); ++i)
        if (!same_point(a.c[i], b.c[i])) return false;
    return true;
}

/// Do side sa of a and side sb of b have identical control points? If so
/// they are the same curve; `reversed` says the side parameters run opposite.
inline bool shared_side(const PolyPatchD& a, int sa, const PolyPatchD& b, int sb, bool& reversed) {
    const auto ca = side_control(a, sa), cb = side_control(b, sb);
    if (ca.size() != cb.size()) return false;
    bool fwd = true, bwd = true;
    const std::size_t n = ca.size();
    for (std::size_t i = 0; i < n; ++i) {
        fwd = fwd && same_point(ca[i], cb[i]);
        bwd = bwd && same_point(ca[i], cb[n - 1 - i]);
    }
    reversed = !fwd;
    return fwd || bwd;
}

/// Is the chart point on the given side of its square? If so, `tau` is its
/// parameter along the side.
inline Tri on_side(const Chart& c, int side, Interval& tau) {
    const bool along_u = side == 0 || side == 2;
    const Interval& fixed = along_u ? c.v : c.u;
    const double value = (side == 0 || side == 3) ? 0.0 : 1.0;
    tau = along_u ? c.u : c.v;
    if (fixed.lo == value && fixed.hi == value) return Tri::Yes;
    if (fixed.hi < value || fixed.lo > value) return Tri::No;
    return Tri::Unknown;
}

/// Certified tube along a side shared by the two patches of `sys` (side sa
/// of the first, sb of the second). The side's points lie on both surfaces
/// exactly, since the two boundary curves are the same polynomial curve.
/// Certified tube along a known straight curve of the intersection, given by
/// its points `at(tau)`, tau in [0, 1], in the coordinates of `sys`.
inline bool known_curve(const SsiSystem& sys, const std::function<Pt4(double)>& at, const solve::SsiOptions& opt,
                        SsiCurve& out, int pieces = 8) {
    out = SsiCurve{};
    out.known = true;
    std::function<bool(double, double, int)> add = [&](double t0, double t1, int depth) {
        const Pt4 a = at(t0), b = at(t1);
        for (bool generous : {true, false}) {
            const TubeSegment seg = solve::detail::make_segment(a, b, opt, generous);
            if (solve::detail::certify_segment(sys, seg)) {
                out.segs.push_back(seg);
                return true;
            }
        }
        if (depth >= 10) return false;
        const double tm = 0.5 * (t0 + t1);
        return add(t0, tm, depth + 1) && add(tm, t1, depth + 1);
    };
    for (int i = 0; i < pieces; ++i)
        if (!add(static_cast<double>(i) / pieces, static_cast<double>(i + 1) / pieces, 0)) return false;
    return true;
}

inline bool side_curve(const SsiSystem& sys, int sa, int sb, bool reversed, const solve::SsiOptions& opt,
                       SsiCurve& out) {
    return known_curve(sys, [&](double tau) {
        double ua, va, ub, vb;
        side_uv(sa, tau, ua, va);
        side_uv(sb, reversed ? 1.0 - tau : tau, ub, vb);
        return Pt4{ua, va, ub, vb};
    }, opt, out);
}

/// Decide, at a certified point of the curve, on which side of it each
/// surface runs behind the other surface (on the side opposite to Su x Sv).
/// The sign is constant along a regular curve, so one point settles it.
inline bool arc_sides(const SsiSystem& sys, const SsiCurve& curve, bool left_behind[2]) {
    const int n = static_cast<int>(curve.segs.size());
    for (int attempt = 0; attempt < n && attempt < 8; ++attempt) {
        const TubeSegment& seg = curve.segs[static_cast<std::size_t>((n / 2 + attempt) % n)];
        Box4 enc;
        if (!solve::curve_point(sys, seg, 0.5, enc)) continue;
        Interval J[3][4];
        sys.J_over(enc, J);
        Interval t[4];
        for (int i = 0; i < 4; ++i) {
            Interval m[3][3];
            for (int r = 0; r < 3; ++r)
                for (int c = 0, cc = 0; c < 4; ++c)
                    if (c != i) m[r][cc++] = J[r][c];
            const Interval d = det3i(m);
            t[i] = i % 2 == 0 ? d : -d;
        }
        Interval mt(0.0);
        for (int i = 0; i < 4; ++i) mt += Interval(seg.m[i]) * t[i];
        if (num::contains_zero(mt)) continue;
        const bool flip = (mt.hi < 0.0) != seg.reversed;
        Interval T[4];
        for (int i = 0; i < 4; ++i) T[i] = flip ? -t[i] : t[i];

        Vec3<Interval> Su, Sv, Bs, Bt;
        sys.frames(enc, Su, Sv, Bs, Bt);
        const Vec3<Interval> nA = cross(Su, Sv), nB = cross(Bs, Bt);
        // Left of the travel direction (du, dv) in a parameter square is (-dv, du).
        const Vec3<Interval> LA = Sv * T[0] - Su * T[1];
        const Vec3<Interval> LB = Bt * T[2] - Bs * T[3];
        const Interval va = dot(LA, nB), vb = dot(LB, nA);
        if (num::contains_zero(va) || num::contains_zero(vb)) continue;
        left_behind[0] = va.hi < 0.0;
        left_behind[1] = vb.hi < 0.0;
        return true;
    }
    return false;
}

/// Is the position inside the arc (strictly between its ends)?
inline Tri in_arc(const ArcGeom& g, const CurvePos& p) {
    const Tri after_a = before(g.a, p), before_b = before(p, g.b);
    if (!g.wrap) {
        if (after_a == Tri::Yes && before_b == Tri::Yes) return Tri::Yes;
        if (after_a == Tri::No || before_b == Tri::No) return Tri::No;
    } else {
        if (after_a == Tri::Yes || before_b == Tri::Yes) return Tri::Yes;
        if (after_a == Tri::No && before_b == Tri::No) return Tri::No;
    }
    return Tri::Unknown;
}

/// A direction in a parameter plane, as intervals. With `pert`, the exact
/// direction (x, y) turned towards the exact direction (nx, ny) by an
/// infinitesimal angle (symbolic perturbation): products are decided by
/// (x, y) and, on an exact tie, by the first-order term.
struct Dir {
    Interval x, y;
    Interval nx{0.0}, ny{0.0};
    bool pert = false;
    bool exact() const { return x.lo == x.hi && y.lo == y.hi; }
};
inline Interval dir_sign(const Interval& c) {
    if (c.lo > 0.0) return Interval(1.0);
    if (c.hi < 0.0) return Interval(-1.0);
    return c;
}
/// a.x b.y - a.y b.x. Directions along the sides of the square are exact, and
/// interval arithmetic would widen an exact zero by one unit; for those the
/// sign is computed exactly and returned as -1, 0 or 1 (only signs are used).
inline Interval cross2(const Dir& a, const Dir& b) {
    if (a.pert || b.pert) {
        const Dir ea{a.x, a.y}, eb{b.x, b.y};
        const Interval c = cross2(ea, eb);
        if (!(c.lo == 0.0 && c.hi == 0.0)) return dir_sign(c);
        Interval t(0.0);
        if (a.pert) t += cross2(Dir{a.nx, a.ny}, eb);
        if (b.pert) t += cross2(ea, Dir{b.nx, b.ny});
        return dir_sign(t);
    }
    if (a.exact() && b.exact())
        return Interval(static_cast<double>(
            (num::Expansion(a.x.lo) * num::Expansion(b.y.lo) - num::Expansion(a.y.lo) * num::Expansion(b.x.lo)).sign()));
    return a.x * b.y - a.y * b.x;
}
inline Interval dot2(const Dir& a, const Dir& b) {
    if (a.pert || b.pert) {
        const Dir ea{a.x, a.y}, eb{b.x, b.y};
        const Interval c = dot2(ea, eb);
        if (!(c.lo == 0.0 && c.hi == 0.0)) return dir_sign(c);
        Interval t(0.0);
        if (a.pert) t += dot2(Dir{a.nx, a.ny}, eb);
        if (b.pert) t += dot2(ea, Dir{b.nx, b.ny});
        return dir_sign(t);
    }
    if (a.exact() && b.exact())
        return Interval(static_cast<double>(
            (num::Expansion(a.x.lo) * num::Expansion(b.x.lo) + num::Expansion(a.y.lo) * num::Expansion(b.y.lo)).sign()));
    return a.x * b.x + a.y * b.y;
}

/// Where d lies relative to d0: 0 right of it, 1 opposite it, 2 left of it,
/// 3 near opposite it (not decided whether right, opposite or left, but
/// pointing strictly backwards: a curve continuing the way back), -1
/// undecided (including d along d0).
inline int side_of(const Dir& d0, const Dir& d) {
    const Interval c = cross2(d0, d);
    if (c.hi < 0.0) return 0;
    if (c.lo > 0.0) return 2;
    if (c.lo == 0.0 && c.hi == 0.0 && dot2(d0, d).hi < 0.0) return 1;
    if (dot2(d0, d).hi < 0.0) return 3;
    return -1;
}

/// Turning clockwise from d0, is d1 (on side s1 of d0) met before d2 (on side s2)?
inline Tri cw_before(int s1, const Dir& d1, int s2, const Dir& d2) {
    if (s1 < 0 || s2 < 0) return Tri::Unknown;
    if ((s1 == 3) != (s2 == 3)) {
        // One direction (say d1) lies near the opposite of d0: its enclosure
        // points strictly backwards and straddles d0's line, so it contains
        // the exact opposite -d0. There cross2(-d0, d2) = -cross2(d0, d2),
        // whose sign is d2's side: right gives d2 first, left gives d1 first,
        // which is the order for d1 = -d0. If cross2(d1, d2) has one strict
        // sign over the whole enclosure, no direction in it coincides with
        // d2, so the true d1 is on the same side of d2 as -d0 is, and the
        // order is the same. (The other is never near-opposite itself.)
        const Interval c = cross2(d1, d2);
        if (c.hi < 0.0) return Tri::Yes;
        if (c.lo > 0.0) return Tri::No;
        return Tri::Unknown;
    }
    if (s1 == 3) return Tri::Unknown;
    if (s1 != s2) return s1 < s2 ? Tri::Yes : Tri::No;
    if (s1 == 1) return Tri::Unknown;
    const Interval c = cross2(d1, d2);
    if (c.hi < 0.0) return Tri::Yes;
    if (c.lo > 0.0) return Tri::No;
    return Tri::Unknown;
}

/// A certified point strictly inside the region bounded by the loops: the
/// first candidate that is proven inside, or with `skip` the one after
/// `skip` such candidates (for a caller whose first point was unusable).
inline bool interior_point(std::vector<CycleGeom>& loops, int depth, double& u, double& v, int skip = 0) {
    if (loops.empty() || loops[0].links.empty()) return false;
    struct Chord { double mx, my, nx, ny, len; };
    std::vector<Chord> chords;
    double cx = 0.0, cy = 0.0;
    for (const Link& L : loops[0].links) {
        const double dx = L.p1[0] - L.p0[0], dy = L.p1[1] - L.p0[1];
        const double len = std::sqrt(dx * dx + dy * dy);
        cx += L.p0[0];
        cy += L.p0[1];
        if (len > 0.0) chords.push_back({0.5 * (L.p0[0] + L.p1[0]), 0.5 * (L.p0[1] + L.p1[1]), -dy / len, dx / len, len});
    }
    cx /= static_cast<double>(loops[0].links.size());
    cy /= static_cast<double>(loops[0].links.size());
    auto try_at = [&](double x, double y) {
        if (locate(loops, Interval(x), Interval(y), depth) != Where::Inside) return false;
        if (skip-- > 0) return false;
        u = x;
        v = y;
        return true;
    };
    if (try_at(cx, cy)) return true;
    // The region is on the left of its outer loop.
    const std::size_t stride = std::max<std::size_t>(1, chords.size() / 8);
    for (double frac : {0.3, 0.1, 0.03, 0.01, 0.003, 0.001})
        for (std::size_t k = 0; k < chords.size(); k += stride) {
            const Chord& c = chords[k];
            const double d = frac * c.len;
            if (try_at(c.mx + d * c.nx, c.my + d * c.ny)) return true;
        }
    return false;
}

class Builder {
    Imprint& im;
    const BooleanOptions& opt;

    std::vector<FaceGeom> fgeom[2];            ///< per solid and face: its loops, for point location
    std::vector<std::vector<int>> faces_on[2]; ///< per solid and global surface: its faces
    std::vector<Vec3<Interval>> face_box[2];   ///< per solid and face: a box that contains the face
    std::vector<PolyPatchI> icarrier;          ///< per global surface: its interval form
    std::vector<Vec3<Interval>> surface_box;   ///< per global surface: hull of the boxes of its faces and edges
    std::vector<char> surface_used;            ///< per global surface: faces or edges of either solid lie on it
    std::vector<Vec3<Interval>> edge_box[2];   ///< per solid and edge: a box that contains the edge
    std::vector<PolyPatchI> face_patch[2];     ///< per solid and face: its carrier over a rectangle holding the face
    std::vector<int> piece_ids[2];             ///< per solid: its edge pieces (not shared), in increasing order
    geom::BoxIndex piece_index[2];             ///< per solid: those pieces by the boxes of their edges
    std::vector<std::array<double, 4>> surface_rect; ///< per global surface: u0 u1 v0 v1 hull of its faces' parameter rectangles
    std::vector<int> vmap[2];                  ///< per solid: vertex -> global vertex

    enum { kOld = 0, kPierce = 1, kCross = 2, kSeam = 3 };
    std::vector<int> vkind, arc_ends;
    std::map<int, int> expected_ends; ///< pierce vertex -> faces of the pierced edge: the new curves that must end there
    std::vector<int> edge_faces[2];   ///< per solid and edge: how many faces use it

    struct Curve {
        int x = -1, y = -1; ///< global surfaces, x < y; the curve's coordinates 0,1 and 2,3
        std::shared_ptr<const CurveBase> base;
        bool left_behind[2] = {false, false};
    };
    std::vector<Curve> curves;

    struct Split {
        int vertex = -1;
        Interval t;   ///< on a side edge
        CurvePos pos; ///< on an arc edge
    };
    std::vector<std::vector<int>> edge_pieces[2]; ///< per solid and edge: pieces from v0 to v1

    struct Node {
        double uv[2] = {0.0, 0.0};
        Interval box[2];
        std::vector<int> out;
    };
    struct HalfEdge {
        int from = -1, to = -1;
        Piece piece;
        bool fwd = true;  ///< relative to the geometry used (piece `src`, or the new arc)
        bool own = false; ///< a piece of the face's own boundary
        int twin = -1;
    };
    enum class Hit { Found, None, Undecided };

    // Surfaces that are the same plane without being the same patch (decided
    // exactly, op/planar.hpp) form a class; maps translate parameters.
    std::vector<int> rep;        ///< per global surface: its class representative
    std::vector<char> is_planar; ///< per global surface: an affine parallelogram
    std::vector<std::vector<int>> members; ///< per representative: the surfaces of its class
    std::map<std::pair<int, int>, Affine> maps;

    /// Proven contacts (op/contact.hpp): curved surface `a` meets the plane of
    /// surface `b` exactly in its side `side`, and nowhere else.
    struct Touch {
        int a = -1, side = -1, b = -1;
    };
    /// Curved surface a meets plane b exactly along a line of constant
    /// parameter inside a's patch, touching it (op/contact.hpp, line_contact).
    struct LineTouch {
        int a = -1, b = -1;
        bool fixed_u = true;
        Interval r;
    };
    std::vector<LineTouch> line_touches;
    /// A curved surface `a` touching the surface `b` at one point: its exact
    /// parameters on `a`, the point, and, where `b` is curved too, its exact
    /// parameters on `b` (on a plane they follow from the point).
    struct PointTouch {
        int a = -1, b = -1;
        double u = 0.0, v = 0.0;
        geom::Vec3<Interval> P;
        bool b_curved = false;
        double s = 0.0, t = 0.0;
        bool alone = false;  ///< the surfaces meet nowhere else (no curve between them)
        bool has_xp = false; ///< xp is the point, exactly
        Vec4d xp{};
    };
    std::vector<PointTouch> point_touches;
    std::vector<Touch> touches;

    const Brep& S(int k) const { return *im.solid[k]; }
    int gid(int k, int carrier) const { return im.carrier_of[k][static_cast<std::size_t>(carrier)]; }
    bool equiv(int a, int b) const { return rep[static_cast<std::size_t>(a)] == rep[static_cast<std::size_t>(b)]; }
    const std::vector<int>& same_plane(int g) const { return members[static_cast<std::size_t>(rep[static_cast<std::size_t>(g)])]; }
    /// Parameters of surface `from` in those of surface `to` (same class).
    const Affine& map(int from, int to) {
        const auto key = std::make_pair(from, to);
        auto it = maps.find(key);
        if (it != maps.end()) return it->second;
        Affine m = Affine::identity();
        if (from != to) same_surface_map(im.carriers[static_cast<std::size_t>(from)], im.carriers[static_cast<std::size_t>(to)], m);
        return maps[key] = m;
    }
    /// Certified tube along the line where planar surfaces x and y meet, over
    /// the part of it near both parameter squares.
    bool plane_line(int x, int y, const SsiSystem& sys, SsiCurve& out) const {
        const PolyPatchD& X = im.carriers[static_cast<std::size_t>(x)];
        const PolyPatchD& Y = im.carriers[static_cast<std::size_t>(y)];
        auto sub = [](const Vec4d& a, const Vec4d& b) { return Vec3d{a.x - b.x, a.y - b.y, a.z - b.z}; };
        const Vec4d Y00 = corner_point(Y, 0, 0), X00 = corner_point(X, 0, 0);
        const Vec3d n = cross(sub(corner_point(Y, 1, 0), Y00), sub(corner_point(Y, 0, 1), Y00));
        const Vec3d e1 = sub(corner_point(X, 1, 0), X00), e2 = sub(corner_point(X, 0, 1), X00);
        // The line in X's parameters: a u + b v + c = 0.
        const double a = dot(n, e1), b = dot(n, e2), c = dot(n, sub(X00, Y00));
        if (!(std::fabs(a) + std::fabs(b) > 0.0)) return false;
        // A point on it and its direction, then clip to the squares, enlarged.
        const double nn = a * a + b * b;
        const double pu = -a * c / nn, pv = -b * c / nn, du = -b, dv = a;
        const double m = 0.1;
        double lo = -1e300, hi = 1e300;
        auto clip = [&](double p, double d, double l, double h) {
            if (d == 0.0) {
                if (p < l || p > h) { lo = 1.0; hi = 0.0; }
                return;
            }
            double s0 = (l - p) / d, s1 = (h - p) / d;
            if (s0 > s1) std::swap(s0, s1);
            lo = std::max(lo, s0);
            hi = std::min(hi, s1);
        };
        clip(pu, du, -m, 1 + m);
        clip(pv, dv, -m, 1 + m);
        const PolyPatchI Xi = geom::to_interval(X);
        auto y_params = [&](double u, double v, Interval& s, Interval& t) {
            const auto h = Xi.eval(Interval(u), Interval(v));
            return plane_params(Y, {h.x / h.w, h.y / h.w, h.z / h.w}, s, t);
        };
        // Y's parameters are affine along the line: clip by them too.
        Interval s0, t0, s1, t1;
        if (!(lo < hi) || !y_params(pu, pv, s0, t0) || !y_params(pu + du, pv + dv, s1, t1)) return false;
        clip(num::mid(s0), num::mid(s1) - num::mid(s0), -m, 1 + m);
        clip(num::mid(t0), num::mid(t1) - num::mid(t0), -m, 1 + m);
        if (!(lo < hi)) return false;
        const double l0 = lo, l1 = hi;
        return known_curve(sys, [&](double tau) {
            const double s = l0 + tau * (l1 - l0), u = pu + s * du, v = pv + s * dv;
            Interval ys(0.5), yt(0.5);
            y_params(u, v, ys, yt);
            return Pt4{u, v, num::mid(ys), num::mid(yt)};
        }, opt.ssi, out);
    }
    /// Two flat surfaces x and y in different planes, intersected without
    /// subdivision. Both are affine parallelograms, X(u, v) = X0 + u E1 + v E2
    /// and Y(s, t) = Y0 + s F1 + t F2, so when they are not degenerate and their
    /// planes are not parallel, the whole solution set of X(u, v) = Y(s, t) in
    /// (u, v, s, t) is one straight line l = P + sigma d (computed here in
    /// intervals: P and d enclose a true point and the true direction).
    ///
    /// A certified tube segment has exactly one solution on each of its
    /// slices; every solution is on l, so that solution is the point of l on
    /// the slice. The segment therefore covers the sigma-range between its two
    /// end slices, which intervals bound from inside. If those ranges overlap
    /// and contain an outer bound of the sigma-range where l is inside the
    /// domain, the intersection is complete. When any step cannot be
    /// certified, return false: the caller then uses the general search.
    bool plane_pair(int x, int y, const SsiSystem& sys, const solve::SsiOptions& sopt, solve::SsiResult& res) const {
        using planar_detail::corner;
        using V3 = Vec3<Interval>;
        auto point = [](const Vec4d& c) {
            const Interval w(c.w);
            return V3{Interval(c.x) / w, Interval(c.y) / w, Interval(c.z) / w};
        };
        auto positive = [](const Interval& i) { return i.lo > 0.0; };
        const PolyPatchD& X = im.carriers[static_cast<std::size_t>(x)];
        const PolyPatchD& Y = im.carriers[static_cast<std::size_t>(y)];
        const V3 X0 = point(corner(X, 0, 0)), E1 = point(corner(X, 1, 0)) - X0, E2 = point(corner(X, 0, 1)) - X0;
        const V3 Y0 = point(corner(Y, 0, 0)), F1 = point(corner(Y, 1, 0)) - Y0, F2 = point(corner(Y, 0, 1)) - Y0;
        const V3 nX = cross(E1, E2), nY = cross(F1, F2);
        if (!positive(dot(nX, nX)) || !positive(dot(nY, nY))) return false; // degenerate parallelogram
        // In X's parameters the line is a u + b v + c = 0.
        const Interval a = dot(nY, E1), b = dot(nY, E2), c = dot(nY, X0 - Y0);
        const Interval nn = a * a + b * b;
        if (!positive(nn)) return false; // parallel planes (or too close to tell)
        // Y's parameters of a point of Y's plane: the normal equations.
        const Interval g11 = dot(F1, F1), g12 = dot(F1, F2), g22 = dot(F2, F2);
        const Interval det = g11 * g22 - g12 * g12;
        if (!positive(det)) return false;
        auto st = [&](const V3& r, Interval& s, Interval& t) {
            const Interval r1 = dot(F1, r), r2 = dot(F2, r);
            s = (g22 * r1 - g12 * r2) / det;
            t = (g11 * r2 - g12 * r1) / det;
        };
        Box4 P, d;
        P[0] = -(a * c) / nn;
        P[1] = -(b * c) / nn;
        d[0] = -b;
        d[1] = a;
        st(X0 + E1 * P[0] + E2 * P[1] - Y0, P[2], P[3]); // a point of l
        st(E1 * d[0] + E2 * d[1], d[2], d[3]);           // its direction
        // Outer bound of the sigma-range where l is inside the domain.
        double lo = -std::numeric_limits<double>::infinity(), hi = std::numeric_limits<double>::infinity();
        for (std::size_t k = 0; k < 4; ++k) {
            if (num::contains_zero(d[k])) continue;
            const Interval q0 = (Interval(sopt.dom_lo[k]) - P[k]) / d[k], q1 = (Interval(sopt.dom_hi[k]) - P[k]) / d[k];
            lo = std::max(lo, std::min(q0.lo, q1.lo));
            hi = std::min(hi, std::max(q0.hi, q1.hi));
        }
        res = solve::SsiResult{};
        if (lo > hi) { res.complete = true; return true; } // l misses the domain
        if (!std::isfinite(lo) || !std::isfinite(hi)) return false;
        // A certified tube over a little more than that range.
        Pt4 pm, dm;
        double speed = 0.0;
        for (std::size_t k = 0; k < 4; ++k) {
            pm[k] = num::mid(P[k]);
            dm[k] = num::mid(d[k]);
            speed = std::max(speed, std::fabs(dm[k]));
        }
        if (!(speed > 0.0)) return false;
        const double pad = 0.05 * (hi - lo) + 0.02 / speed;
        const double s0 = lo - pad, s1 = hi + pad;
        // Segments about as long as traced ones (h_max in parameter space),
        // so that later certified points on parts of the curve find them
        // as they would on a traced curve.
        double extent = 0.0;
        for (std::size_t k = 0; k < 4; ++k) extent += (s1 - s0) * (s1 - s0) * dm[k] * dm[k];
        const int pieces = static_cast<int>(std::min(4096.0, std::max(1.0, std::ceil(std::sqrt(extent) / opt.ssi.h_max))));
        SsiCurve line;
        if (!known_curve(sys, [&](double tau) {
                const double sg = s0 + tau * (s1 - s0);
                return Pt4{pm[0] + sg * dm[0], pm[1] + sg * dm[1], pm[2] + sg * dm[2], pm[3] + sg * dm[3]};
            }, opt.ssi, line, pieces))
            return false;
        // The sigma-range each segment covers, bounded from inside.
        std::vector<std::pair<double, double>> covered;
        for (const solve::TubeSegment& seg : line.segs) {
            Interval md(0.0), mp(0.0);
            for (std::size_t k = 0; k < 4; ++k) {
                md += Interval(seg.m[k]) * d[k];
                mp += Interval(seg.m[k]) * (P[k] - Interval(seg.a[k]));
            }
            if (num::contains_zero(md)) return false;
            const Interval at_lo = (Interval(seg.c_lo) - mp) / md, at_hi = (Interval(seg.c_hi) - mp) / md;
            if (md.lo > 0.0) covered.push_back({at_lo.hi, at_hi.lo});
            else covered.push_back({at_hi.hi, at_lo.lo});
        }
        std::sort(covered.begin(), covered.end());
        double reach = lo;
        for (const auto& [from, to] : covered) {
            if (from > reach) break;
            reach = std::max(reach, to);
        }
        if (!(reach >= hi)) return false;
        line.known = false; // a curve to clip, not a shared side
        res.curves.push_back(std::move(line));
        res.complete = true;
        return true;
    }

    /// Do faces of both solids lie on this surface or one in the same plane?
    bool shared(int g) const {
        bool a = false, b = false;
        for (int h : same_plane(g)) {
            a = a || !faces_on[0][h].empty();
            b = b || !faces_on[1][h].empty();
        }
        return a && b;
    }
    /// Do two vertices both lie on three planes that meet in exactly one point?
    /// A chart on a planar surface says the vertex lies in its plane. Planes
    /// are compared by class (exactly the same plane); the three normals are
    /// shown independent by an interval determinant.
    bool three_planes(const std::vector<Chart>& a, const std::vector<Chart>& b) const {
        std::vector<int> common;
        for (const Chart& c : a) {
            if (!is_planar[static_cast<std::size_t>(c.carrier)]) continue;
            const int r = rep[static_cast<std::size_t>(c.carrier)];
            if (std::find(common.begin(), common.end(), r) != common.end()) continue;
            for (const Chart& d : b)
                if (is_planar[static_cast<std::size_t>(d.carrier)] && rep[static_cast<std::size_t>(d.carrier)] == r) {
                    common.push_back(r);
                    break;
                }
        }
        if (common.size() < 3) return false;
        auto normal = [&](int g) {
            const PolyPatchD& P = im.carriers[static_cast<std::size_t>(g)];
            auto e = [](const Vec4d& p, const Vec4d& q) {
                return Vec3<Interval>{Interval(p.x) - Interval(q.x), Interval(p.y) - Interval(q.y), Interval(p.z) - Interval(q.z)};
            };
            using planar_detail::corner;
            return cross(e(corner(P, 1, 0), corner(P, 0, 0)), e(corner(P, 0, 1), corner(P, 0, 0)));
        };
        for (std::size_t i = 0; i < common.size(); ++i)
            for (std::size_t j = i + 1; j < common.size(); ++j)
                for (std::size_t l = j + 1; l < common.size(); ++l) {
                    const Interval det = dot(normal(common[i]), cross(normal(common[j]), normal(common[l])));
                    if (!num::contains_zero(det)) return true;
                }
        return false;
    }
    /// Does solid k have faces on g or on a surface in g's plane?
    bool faces_in_plane(int k, int g) const {
        for (int h : same_plane(g))
            if (!faces_on[k][h].empty()) return true;
        return false;
    }
    /// A certified box around a vertex, from one of its charts.
    Vec3<Interval> vertex_box(int v) const {
        const Chart& c = im.vertices[static_cast<std::size_t>(v)].charts.front();
        const auto h = icarrier[static_cast<std::size_t>(c.carrier)].eval(c.u, c.v);
        return {h.x / h.w, h.y / h.w, h.z / h.w};
    }
    /// Give a vertex parameters on every surface in the same plane as one it
    /// has, where that surface has faces or edges near it. (A surface's
    /// charts serve its faces and edges; a vertex whose box misses all of
    /// them cannot lie on any. Without this, a vertex on a plane tiled into n
    /// faces would carry n charts.) Unused surfaces get a chart always.
    void spread(int v) {
        const std::vector<Chart> own = im.vertices[static_cast<std::size_t>(v)].charts;
        bool boxed = false;
        Vec3<Interval> box;
        for (const Chart& c : own)
            for (int h : same_plane(c.carrier)) {
                if (chart(v, h)) continue;
                if (surface_used[static_cast<std::size_t>(h)]) {
                    if (!boxed) { box = vertex_box(v); boxed = true; }
                    if (solve::detail::boxes_disjoint(box, surface_box[static_cast<std::size_t>(h)])) continue;
                }
                Chart m;
                m.carrier = h;
                map(c.carrier, h).apply(c.u, c.v, m.u, m.v);
                im.vertices[static_cast<std::size_t>(v)].charts.push_back(m);
            }
    }
    /// Side `side` of surface a lies exactly in the plane of b, and a touches
    /// that plane only there (proven by the caller). Every vertex on that side
    /// (an exact side coordinate on a) lies in the plane: give it parameters
    /// on b, enclosures of the exact ones, and on the surfaces of b's plane.
    void on_plane_side(int a, int side, int b) {
        touches.push_back({a, side, b});
        for (int v = 0; v < static_cast<int>(im.vertices.size()); ++v) {
            touch_charts(v);
            exact_charts(v);
        }
    }
    /// For every proven contact (a along its side s with b's plane): a vertex
    /// with an exact coordinate on that side of a lies in the plane; give it
    /// parameters on b (enclosures of the exact ones) and on b's class.
    void touch_charts(int v) {
        for (const Touch& t : touches) {
            const Chart* c = chart(v, t.a);
            if (!c || chart(v, t.b)) continue;
            Interval tau;
            if (on_side(*c, t.side, tau) != Tri::Yes) continue;
            const auto h = icarrier[static_cast<std::size_t>(t.a)].eval(c->u, c->v);
            const Vec3<Interval> P{h.x / h.w, h.y / h.w, h.z / h.w};
            if (surface_used[static_cast<std::size_t>(t.b)] && solve::detail::boxes_disjoint(P, surface_box[static_cast<std::size_t>(t.b)]))
                continue;
            Chart m;
            m.carrier = t.b;
            if (!plane_params(im.carriers[static_cast<std::size_t>(t.b)], P, m.u, m.v)) continue;
            im.vertices[static_cast<std::size_t>(v)].charts.push_back(m);
            spread(v);
        }
    }
    /// The side of curved surface g along which it meets the plane of b, if
    /// it meets that plane nowhere else, even slightly beyond its square
    /// (proven exactly, op/contact.hpp; cached). -1 if there is none.
    std::map<std::pair<int, int>, int> only_side_cache;
    static constexpr double kSnapMargin = 1.0 / 64.0;
    int only_side(int g, int b) {
        if (!is_planar[static_cast<std::size_t>(b)] || equiv(g, b)) return -1;
        const auto key = std::make_pair(g, b);
        const auto it = only_side_cache.find(key);
        if (it != only_side_cache.end()) return it->second;
        int found = -1;
        const PolyPatchD& G = im.carriers[static_cast<std::size_t>(g)];
        const PolyPatchD& B = im.carriers[static_cast<std::size_t>(b)];
        // Two different planes meet in one line: a side of g in b's plane is
        // all of it, over the whole plane of g.
        const bool flat = is_planar[static_cast<std::size_t>(g)];
        for (int s = 0; s < 4 && found < 0; ++s)
            if (side_in_plane(G, s, B) && (flat || side_contact(G, s, B, 12, kSnapMargin).only_side)) found = s;
        only_side_cache[key] = found;
        return found;
    }
    /// Make exact what is exact. A vertex on a plane b and on a curved surface
    /// g that meets b's plane only along its side s (within the margin of
    /// g's square) is on that side: its coordinate across s is the side's
    /// value. The side is then a side of every surface that shares it by
    /// construction (identical control points), and the vertex's parameters
    /// there follow exactly. Each chart keeps its intersection with the
    /// exact form, which contains the vertex's true parameters.
    /// A new vertex that is exactly on a side of a surface (its parameter
    /// across the side is the side's value) lies on every surface that shares
    /// that side by construction: it gets parameters there. (A box's top edge
    /// is a side of its top and of a side face; a vertex found on the top's
    /// side is on the side face too.)
    void share_side_charts(int v) {
        std::vector<Chart> add;
        const std::vector<Chart> ch = im.vertices[static_cast<std::size_t>(v)].charts;
        for (const Chart& c : ch)
            for (int sd = 0; sd < 4; ++sd) {
                Interval tau;
                if (detail::on_side(c, sd, tau) != Tri::Yes || tau.lo < 0.0 || tau.hi > 1.0) continue;
                for (const auto& [h, t, rev] : sharing_side(c.carrier, sd)) {
                    bool has = false;
                    for (const Chart& x : ch) has = has || x.carrier == h;
                    for (const Chart& x : add) has = has || x.carrier == h;
                    if (has) continue;
                    Chart n;
                    n.carrier = h;
                    side_uv(t, rev ? Interval(1.0) - tau : tau, n.u, n.v);
                    add.push_back(n);
                }
            }
        if (add.empty()) return;
        for (const Chart& x : add) im.vertices[static_cast<std::size_t>(v)].charts.push_back(x);
        spread(v);
    }
    void exact_charts(int v) {
        if (touches.empty()) return; // only after a contact has been proven
        std::vector<Chart>& ch = im.vertices[static_cast<std::size_t>(v)].charts;
        std::vector<Chart> added; // charts on surfaces that share a side, appended at the end
        // x := x intersected with y; false (and x unchanged) if they are disjoint.
        auto narrow = [](Interval& x, const Interval& y) {
            const double lo = std::max(x.lo, y.lo), hi = std::min(x.hi, y.hi);
            if (lo > hi) return false;
            x = Interval(lo, hi);
            return true;
        };
        for (std::size_t i = 0; i < ch.size(); ++i) {
            const int g = ch[i].carrier;
            const bool flat = is_planar[static_cast<std::size_t>(g)];
            for (std::size_t j = 0; j < ch.size(); ++j) {
                const int b = ch[j].carrier;
                if (j == i || !is_planar[static_cast<std::size_t>(b)]) continue;
                const int s = only_side(g, b);
                if (s < 0) continue;
                Chart& c = ch[i];
                if (!flat && (c.u.lo < -kSnapMargin || c.u.hi > 1.0 + kSnapMargin || c.v.lo < -kSnapMargin ||
                              c.v.hi > 1.0 + kSnapMargin))
                    continue; // beyond the square the proof covers
                const double value = (s == 0 || s == 3) ? 0.0 : 1.0;
                Interval& across = (s == 0 || s == 2) ? c.v : c.u;
                if (!num::contains(across, value)) continue;
                across = Interval(value);
                // The side's points are those of every side identical to it.
                const Interval tau = (s == 0 || s == 2) ? c.u : c.v;
                for (std::size_t l = 0; l < ch.size(); ++l) {
                    if (l == i) continue;
                    for (int t = 0; t < 4; ++t) {
                        bool rev = false;
                        if (!shared_side(im.carriers[static_cast<std::size_t>(g)], s, im.carriers[static_cast<std::size_t>(ch[l].carrier)], t, rev))
                            continue;
                        Interval u, w;
                        side_uv(t, rev ? Interval(1.0) - tau : tau, u, w);
                        if (narrow(u, ch[l].u) && narrow(w, ch[l].v)) { ch[l].u = u; ch[l].v = w; }
                        break;
                    }
                }
                // A surface that shares the side and has no chart yet: the vertex
                // is on it when its position along the side is within the side.
                if (tau.lo >= 0.0 && tau.hi <= 1.0)
                    for (const auto& [h, t, rev] : sharing_side(g, s)) {
                        bool has = false;
                        for (const Chart& x : ch) has = has || x.carrier == h;
                        for (const Chart& x : added) has = has || x.carrier == h;
                        if (has) continue;
                        Chart n;
                        n.carrier = h;
                        side_uv(t, rev ? Interval(1.0) - tau : tau, n.u, n.v);
                        added.push_back(n);
                    }
                break;
            }
        }
        if (!added.empty()) {
            for (const Chart& x : added) ch.push_back(x);
            spread(v);
        }
    }
    /// Sides of other global surfaces identical to side s of g (by
    /// construction: identical control points), with their direction.
    std::map<std::pair<int, int>, std::vector<std::tuple<int, int, bool>>> sharing_cache;
    const std::vector<std::tuple<int, int, bool>>& sharing_side(int g, int s) {
        const auto key = std::make_pair(g, s);
        const auto it = sharing_cache.find(key);
        if (it != sharing_cache.end()) return it->second;
        std::vector<std::tuple<int, int, bool>> out;
        for (int h = 0; h < static_cast<int>(im.carriers.size()); ++h) {
            if (h == g) continue;
            for (int t = 0; t < 4; ++t) {
                bool rev = false;
                if (shared_side(im.carriers[static_cast<std::size_t>(g)], s, im.carriers[static_cast<std::size_t>(h)], t, rev)) {
                    out.emplace_back(h, t, rev);
                    break;
                }
            }
        }
        return sharing_cache[key] = std::move(out);
    }
    bool fail(const std::string& msg) { im.error = msg; return false; }
    static std::string str(int i) { return std::to_string(i); }
    void note(const std::string& line) { if (opt.trace) im.trace.push_back(line); }
    /// A surface for the trace: its id and p (an exact parallelogram), f
    /// (flat, on a rational carrier) or c (curved).
    std::string sk(int s) const {
        const PolyPatchD& P = im.carriers[static_cast<std::size_t>(s)];
        return str(s) + (is_planar[static_cast<std::size_t>(s)] ? "p" : (same_plane_patches(P, P) ? "f" : "c"));
    }

    const Chart* chart(int v, int g) const {
        for (const Chart& c : im.vertices[static_cast<std::size_t>(v)].charts)
            if (c.carrier == g) return &c;
        return nullptr;
    }
    void edge_carriers(int k, int e, int& g0, int& g1) const {
        const BEdge& ed = S(k).edges[static_cast<std::size_t>(e)];
        g0 = gid(k, ed.arc ? ed.carrier[0] : ed.on[0].carrier);
        g1 = gid(k, ed.arc ? ed.carrier[1] : ed.on[1].carrier);
    }

    /// The face of solid k on global surface g that contains a certified parameter box.
    Hit find_face(int k, int g, const Interval& u, const Interval& v, int& face) {
        face = -1;
        for (int f : faces_on[k][g]) {
            const Where w = locate(fgeom[k][f], u, v, opt.refine_depth);
            if (w == Where::Inside) { face = f; return Hit::Found; }
            if (w != Where::Outside) return Hit::Undecided;
        }
        return Hit::None;
    }
    /// The same over every surface in the plane of g: the face of solid k that
    /// contains the point, on whichever of those surfaces it lies. `reversed`
    /// says the face's surface is parameterised with the opposite orientation.
    Hit find_face_in_plane(int k, int g, const Interval& u, const Interval& v, int& face, bool& reversed) {
        face = -1;
        reversed = false;
        for (int h : same_plane(g)) {
            if (faces_on[k][h].empty()) continue;
            Interval uh = u, vh = v;
            if (h != g) map(g, h).apply(u, v, uh, vh);
            const Hit hit = find_face(k, h, uh, vh, face);
            if (hit == Hit::Undecided) return hit;
            if (hit == Hit::Found) {
                if (h != g) {
                    const Interval d = map(g, h).det();
                    if (num::contains_zero(d)) return Hit::Undecided;
                    reversed = d.hi < 0.0;
                }
                return hit;
            }
        }
        return Hit::None;
    }
    /// Does face f of solid k lie strictly on one side of the plane of the
    /// planar surface h? The face lies in the convex hull of the control
    /// points of its carrier over its rectangle (positive weights), so it
    /// does when all of them do.
    bool off_plane(int k, int f, int h) const {
        using planar_detail::corner;
        const PolyPatchD& P = im.carriers[static_cast<std::size_t>(h)];
        auto e = [](const Vec4d& p, const Vec4d& q) {
            return Vec3<Interval>{Interval(p.x) - Interval(q.x), Interval(p.y) - Interval(q.y), Interval(p.z) - Interval(q.z)};
        };
        const Vec4d c00 = corner(P, 0, 0);
        const Vec3<Interval> n = cross(e(corner(P, 1, 0), c00), e(corner(P, 0, 1), c00));
        const Interval w0(c00.w);
        const Interval d = dot(n, Vec3<Interval>{Interval(c00.x) / w0, Interval(c00.y) / w0, Interval(c00.z) / w0});
        int sign = 0;
        for (const auto& c : face_patch[k][static_cast<std::size_t>(f)].c) {
            if (!(c.w.lo > 0.0)) return false;
            const Interval s = dot(n, Vec3<Interval>{c.x, c.y, c.z}) - d * c.w; // w (n . x - d)
            const int si = s.lo > 0.0 ? 1 : s.hi < 0.0 ? -1 : 0;
            if (si == 0 || (sign != 0 && si != sign)) return false;
            sign = si;
        }
        return sign != 0;
    }

public:
    Builder(Imprint& imprint, const BooleanOptions& options) : im(imprint), opt(options) {}

    bool run() {
        using Stage = bool (Builder::*)();
        const Stage stages[7] = {&Builder::prepare,     &Builder::intersect_carriers, &Builder::intersect_edges,
                                 &Builder::clip_curves, &Builder::split_edges,        &Builder::build_regions,
                                 &Builder::classify};
        for (int i = 0; i < 7; ++i) {
            im.stats.stage = i;
            const auto t0 = std::chrono::steady_clock::now();
            const bool ok = (this->*stages[i])();
            im.stats.ms[i] = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
            if (!ok) return false;
        }
        return true;
    }

private:
    // --- 0. operands, shared surfaces and shared vertices -------------------------
    bool prepare() {

        for (int k = 0; k < 2; ++k) {
            const Brep& s = S(k);
            if (s.faces.empty()) return fail("solid " + str(k) + " has no faces");
            const topo::AuditReport a = audit(s);
            if (!a.closed()) return fail("solid " + str(k) + " is not a closed oriented boundary");
            if (!a.manifold()) {
                std::string why;
                if (!check_edge_fans(s, why)) return fail("solid " + str(k) + ": " + why);
            }
            if (!s.outward_known && !(volume(s, true) > 0.0)) return fail("solid " + str(k) + " is inside out");
            for (const BEdge& e : s.edges)
                if (!e.arc && e.on[0].reversed)
                    return fail("solid " + str(k) + ": a side edge runs against its first side (unsupported)");
        }
        // Surfaces: a surface of solid 1 with the same control points as one of
        // solid 0 is that surface.
        for (int k = 0; k < 2; ++k) {
            im.carrier_of[k].clear();
            for (const PolyPatchD& c : S(k).carriers) {
                int g = -1;
                if (k == 1)
                    for (std::size_t j = 0; j < S(0).carriers.size() && g < 0; ++j)
                        if (same_patch(c, S(0).carriers[j])) g = im.carrier_of[0][j];
                if (g < 0) {
                    g = static_cast<int>(im.carriers.size());
                    im.carriers.push_back(c);
                }
                im.carrier_of[k].push_back(g);
            }
        }
        const int nG = static_cast<int>(im.carriers.size());
        surface_rect.assign(static_cast<std::size_t>(nG), {1.0, 0.0, 1.0, 0.0});
        for (int k = 0; k < 2; ++k) {
            const Brep& s = S(k);
            const int nF = static_cast<int>(s.faces.size());
            faces_on[k].assign(static_cast<std::size_t>(nG), {});
            fgeom[k].resize(static_cast<std::size_t>(nF));
            for (int f = 0; f < nF; ++f) {
                faces_on[k][gid(k, s.faces[f].carrier)].push_back(f);
                std::string err;
                if (!face_geometry(s, f, fgeom[k][f], err)) return fail("solid " + str(k) + ": " + err);
                // The face lies in the parameter rectangle spanned by its
                // loops, so the carrier over that rectangle encloses it.
                Interval u(1.0, 0.0), v(1.0, 0.0);
                for (const CycleGeom& cyc : fgeom[k][f].loops)
                    for (const Link& L : cyc.links) {
                        u = Interval(std::min(u.lo, L.bu.lo), std::max(u.hi, L.bu.hi));
                        v = Interval(std::min(v.lo, L.bv.lo), std::max(v.hi, L.bv.hi));
                    }
                const double u0 = std::max(u.lo, 0.0), u1 = std::min(u.hi, 1.0);
                const double v0 = std::max(v.lo, 0.0), v1 = std::min(v.hi, 1.0);
                // The surface's region that holds faces: the hull of their
                // parameter rectangles (each encloses its face's loops).
                auto& rect = surface_rect[static_cast<std::size_t>(gid(k, s.faces[f].carrier))];
                rect = {std::min(rect[0], u0), std::max(rect[1], u1), std::min(rect[2], v0), std::max(rect[3], v1)};
                const PolyPatchI patch = geom::to_interval(s.carriers[s.faces[f].carrier]);
                face_patch[k].push_back(u0 < u1 && v0 < v1 ? patch.restricted(u0, u1, v0, v1) : patch);
                face_box[k].push_back(geom::projected_bounds(face_patch[k].back()));
            }
        }
        icarrier.clear();
        for (const PolyPatchD& c : im.carriers) icarrier.push_back(geom::to_interval(c));
        // A box around each edge.
        for (int k = 0; k < 2; ++k) {
            const Brep& s = S(k);
            edge_box[k].clear();
            for (const BEdge& edge : s.edges) {
                Vec3<Interval> eb;
                if (!edge.arc) {
                    const PolyPatchD c = solve::make_curve(side_control(s.carriers[static_cast<std::size_t>(edge.on[0].carrier)], edge.on[0].side));
                    eb = geom::projected_bounds(geom::to_interval(c).restricted(std::max(edge.t0.lo, 0.0), std::min(edge.t1.hi, 1.0), 0.0, 1.0));
                } else {
                    const CurveBase& base = *edge.geom->base;
                    bool first = true;
                    for (const Span& sp : edge.geom->spans) {
                        const auto p = base.sys->point_a(base.curve->segs[static_cast<std::size_t>(sp.seg)].X);
                        eb = first ? p : Vec3<Interval>{num::hull(eb.x, p.x), num::hull(eb.y, p.y), num::hull(eb.z, p.z)};
                        first = false;
                    }
                }
                edge_box[k].push_back(eb);
            }
        }
        // Where each surface is used: the hull of the boxes of the faces and
        // the edges on it. (An edge can lie on a side of a surface whose own
        // faces are elsewhere, its face being on another surface in the same
        // plane; and on a surface without faces, a seam.)
        surface_box.assign(static_cast<std::size_t>(nG), {});
        surface_used.assign(static_cast<std::size_t>(nG), 0);
        auto widen = [&](int g, const Vec3<Interval>& eb) {
            auto& b = surface_box[static_cast<std::size_t>(g)];
            if (!surface_used[static_cast<std::size_t>(g)]) { b = eb; surface_used[static_cast<std::size_t>(g)] = 1; return; }
            b = {num::hull(b.x, eb.x), num::hull(b.y, eb.y), num::hull(b.z, eb.z)};
        };
        for (int k = 0; k < 2; ++k) {
            for (int g = 0; g < nG; ++g)
                for (int f : faces_on[k][static_cast<std::size_t>(g)]) widen(g, face_box[k][static_cast<std::size_t>(f)]);
            for (int e = 0; e < static_cast<int>(S(k).edges.size()); ++e) {
                int g0 = -1, g1 = -1;
                edge_carriers(k, e, g0, g1);
                widen(g0, edge_box[k][static_cast<std::size_t>(e)]);
                widen(g1, edge_box[k][static_cast<std::size_t>(e)]);
            }
        }

        // Surfaces in the same plane (exact), and curved surfaces that are
        // affine reparameterisations of one another: union-find. Surfaces
        // without faces count too: an edge can lie on one (a seam).
        rep.resize(static_cast<std::size_t>(nG));
        is_planar.assign(static_cast<std::size_t>(nG), 0);
        for (int g = 0; g < nG; ++g) {
            rep[static_cast<std::size_t>(g)] = g;
            is_planar[static_cast<std::size_t>(g)] = planar_affine(im.carriers[static_cast<std::size_t>(g)]) ? 1 : 0;
        }
        auto root = [&](int g) {
            while (rep[static_cast<std::size_t>(g)] != g) g = rep[static_cast<std::size_t>(g)];
            return g;
        };
        // Planes. Coplanarity is an equivalence, so each surface is tested
        // against one member of each plane found so far, not against every
        // surface; and a cheap enclosure of the plane (unit normal and
        // distance from the origin, both up to sign) rules most planes out
        // before the exact test. Exactly coplanar surfaces have overlapping
        // enclosures, so nothing is missed.
        struct PlaneKey { Interval k[4]; };
        auto plane_key = [&](int g) {
            using planar_detail::corner;
            const PolyPatchD& P = im.carriers[static_cast<std::size_t>(g)];
            auto e = [](const Vec4d& p, const Vec4d& q) {
                return Vec3<Interval>{Interval(p.x) - Interval(q.x), Interval(p.y) - Interval(q.y), Interval(p.z) - Interval(q.z)};
            };
            const Vec4d c00 = corner(P, 0, 0);
            const Vec3<Interval> n = cross(e(corner(P, 1, 0), c00), e(corner(P, 0, 1), c00));
            const Interval len = num::sqrt(dot(n, n)), w(c00.w);
            const Vec3<Interval> o{Interval(c00.x) / w, Interval(c00.y) / w, Interval(c00.z) / w};
            return PlaneKey{{num::abs(n.x / len), num::abs(n.y / len), num::abs(n.z / len), num::abs(dot(n, o) / len)}};
        };
        auto may_coincide = [](const PlaneKey& a, const PlaneKey& b) {
            for (int i = 0; i < 4; ++i)
                if (num::disjoint(a.k[i], b.k[i])) return false;
            return true;
        };
        std::vector<int> plane_rep;
        std::vector<PlaneKey> plane_rep_key;
        std::vector<int> curved;
        for (int g = 0; g < nG; ++g) {
            if (!is_planar[static_cast<std::size_t>(g)]) { curved.push_back(g); continue; }
            const PlaneKey key = plane_key(g);
            bool joined = false;
            for (std::size_t r = 0; r < plane_rep.size() && !joined; ++r) {
                if (!may_coincide(key, plane_rep_key[r])) continue;
                if (!coplanar(im.carriers[static_cast<std::size_t>(plane_rep[r])], im.carriers[static_cast<std::size_t>(g)])) continue;
                rep[static_cast<std::size_t>(g)] = plane_rep[r];
                ++im.stats.coplanar_pairs;
                joined = true;
            }
            if (!joined) { plane_rep.push_back(g); plane_rep_key.push_back(key); }
        }
        // Curved (or rational flat) surfaces: one patch an affine
        // reparameterisation of the other, proven exactly. Only worth trying
        // where the two can meet at all. Bilinear patches are left out: flat
        // ones are planes above, and others (a rotated box, rounded) do not
        // coincide with anything but themselves.
        std::vector<Vec3<Interval>> curved_box;
        for (int g : curved) curved_box.push_back(geom::projected_bounds(icarrier[static_cast<std::size_t>(g)]));
        for (std::size_t i = 0; i < curved.size(); ++i)
            for (std::size_t j = i + 1; j < curved.size(); ++j) {
                const int a = curved[i], b = curved[j];
                if (root(a) == root(b)) continue;
                const PolyPatchD& A = im.carriers[static_cast<std::size_t>(a)];
                const PolyPatchD& B = im.carriers[static_cast<std::size_t>(b)];
                const bool same_degree = A.p == B.p && A.q == B.q, swapped_degree = A.p == B.q && A.q == B.p;
                if ((!same_degree && !swapped_degree) || (A.p == 1 && A.q == 1)) continue;
                if (solve::detail::boxes_disjoint(curved_box[i], curved_box[j])) continue;
                Affine m;
                // Also with the parameters exchanged (an inner and an outer
                // face of one cylinder, made by different solids).
                if (!affine_reparam(A, B, m) && !affine_reparam(B, A, m) && !affine_reparam_any(A, B, m) &&
                    !affine_reparam_any(B, A, m))
                    continue;
                ++im.stats.curved_coincident_pairs;
                rep[static_cast<std::size_t>(root(b))] = root(a);
            }
        members.assign(static_cast<std::size_t>(nG), {});
        for (int g = 0; g < nG; ++g) {
            rep[static_cast<std::size_t>(g)] = root(g);
            members[static_cast<std::size_t>(rep[static_cast<std::size_t>(g)])].push_back(g);
        }
        for (int g = 0; g < nG; ++g)
            if (rep[static_cast<std::size_t>(g)] == g && shared(g)) ++im.stats.shared_carriers;

        // Vertices: a vertex of solid 1 with the same certified parameters on a
        // shared surface as a vertex of solid 0, or at exactly the same corner
        // point, is that vertex.
        auto corner_of = [&](const std::vector<Chart>& charts, Vec4d& p) {
            for (const Chart& c : charts)
                if (exact_corner(im.carriers[static_cast<std::size_t>(c.carrier)], c.u, c.v, p)) return true;
            return false;
        };
        auto chart_box = [&](const std::vector<Chart>& charts) {
            const Chart& c = charts.front();
            const auto h = icarrier[static_cast<std::size_t>(c.carrier)].eval(c.u, c.v);
            return Vec3<Interval>{h.x / h.w, h.y / h.w, h.z / h.w};
        };
        for (int k = 0; k < 2; ++k) {
            const Brep& s = S(k);
            vmap[k].assign(s.vertices.size(), -1);
            const int n0 = static_cast<int>(im.vertices.size());
            // The vertices of solid 0 by their boxes: a vertex of solid 1 can
            // only be one whose box meets its own (both contain the point).
            geom::BoxIndex index;
            if (k == 1) {
                std::vector<Vec3<Interval>> boxes;
                for (int w = 0; w < n0; ++w) boxes.push_back(vertex_box(w));
                index = geom::BoxIndex(std::move(boxes));
            }
            for (std::size_t i = 0; i < s.vertices.size(); ++i) {
                std::vector<Chart> charts;
                for (const Chart& c : s.vertices[i].charts) charts.push_back({gid(k, c.carrier), c.u, c.v});
                int match = -1;
                std::vector<int> near;
                if (k == 1) near = charts.empty() ? index.query({Interval::entire(), Interval::entire(), Interval::entire()})
                                                  : index.query(chart_box(charts));
                for (int w : near) {
                    if (match >= 0) break;
                    for (const Chart& c : charts) {
                        const Chart* d = chart(w, c.carrier);
                        if (d && d->u.lo == c.u.lo && d->u.hi == c.u.hi && d->v.lo == c.v.lo && d->v.hi == c.v.hi) {
                            match = w;
                            break;
                        }
                    }
                }
                Vec4d p;
                if (match < 0 && k == 1 && corner_of(charts, p))
                    for (int w : near) {
                        if (match >= 0) break;
                        Vec4d q;
                        if (corner_of(im.vertices[static_cast<std::size_t>(w)].charts, q) && planar_detail::same_point_h(q, p))
                            match = w;
                    }
                // Both lie, by construction, on the same three planes, and those
                // planes meet in a single point: they are that point.
                for (int w : near) {
                    if (match >= 0) break;
                    if (three_planes(charts, im.vertices[static_cast<std::size_t>(w)].charts)) match = w;
                }
                // Both have an exact point (recorded, or that of a surface on
                // which the parameters are exact, a corner say), and it is the
                // same point. Looked at only where nothing else matched.
                Vec4d xi = s.vertices[i].xp;
                bool exact_i = s.vertices[i].exact;
                if (opt.exact_points && match < 0 && !near.empty()) {
                    for (std::size_t c = 0; c < charts.size() && !exact_i; ++c)
                        exact_i = is_point(charts[c].u) && is_point(charts[c].v) &&
                                  exact_point_at(im.carriers[static_cast<std::size_t>(charts[c].carrier)], charts[c].u.lo, charts[c].v.lo, xi);
                    if (exact_i)
                        for (int w : near) {
                            if (match >= 0) break;
                            Vec4d xw;
                            if (exact_point(w, xw) && planar_detail::same_point_h(xw, xi)) match = w;
                        }
                }
                if (match >= 0) {
                    Imprint::Vertex& o = im.vertices[static_cast<std::size_t>(match)];
                    for (const Chart& c : charts)
                        if (!chart(match, c.carrier)) o.charts.push_back(c);
                    if (!o.exact && exact_i) { o.exact = true; o.xp = xi; }
                    vmap[k][i] = match;
                    ++im.stats.merged_vertices;
                    continue;
                }
                vmap[k][i] = static_cast<int>(im.vertices.size());
                im.vertices.push_back({s.vertices[i].point, std::move(charts)});
                im.vertices.back().exact = exact_i;
                im.vertices.back().xp = xi;
            }
        }
        for (int v = 0; v < static_cast<int>(im.vertices.size()); ++v) spread(v);
        vkind.assign(im.vertices.size(), kOld);
        arc_ends.assign(im.vertices.size(), 0);
        return true;
    }

    /// Does a straight side lie in the plane of g? (Its surfaces are in that
    /// plane, or the side itself lies in it exactly.)
    bool side_lies_in(int g0, int side0, int g1, int g) const {
        if (equiv(g, g0) || equiv(g, g1)) return true;
        if (!is_planar[static_cast<std::size_t>(g)]) return false;
        return side_in_plane(im.carriers[static_cast<std::size_t>(g0)], side0, im.carriers[static_cast<std::size_t>(g)]);
    }
    /// Is a vertex on side `side0` of g0 (which g1 shares), and where? Exact
    /// parameters decide directly. Parameters carried over from a surface in
    /// the same plane are intervals; then the vertex lies on both planes g0
    /// and g1, whose intersection is exactly the side's line.
    Tri vertex_on_side(int v, int g0, int side0, int g1, Interval& tau) const {
        const Chart* a = chart(v, g0);
        if (!a || !chart(v, g1)) return Tri::No;
        const Tri st = on_side(*a, side0, tau);
        if (st != Tri::Unknown) return st;
        if (is_planar[static_cast<std::size_t>(g0)] && is_planar[static_cast<std::size_t>(g1)] && !equiv(g0, g1)) return Tri::Yes;
        return Tri::Unknown;
    }

    // --- 1. surface/surface --------------------------------------------------------
    bool intersect_carriers() {
        const int nG = static_cast<int>(im.carriers.size());
        auto interacts = [&](int a, int b) { // faces of solid 0 on a against faces of solid 1 on b
            for (int fa : faces_on[0][a])
                for (int fb : faces_on[1][b]) {
                    if (solve::detail::boxes_disjoint(face_box[0][fa], face_box[1][fb])) continue;
                    // A face entirely on one side of the other's plane cannot meet it.
                    if (is_planar[static_cast<std::size_t>(a)] && off_plane(1, fb, a)) continue;
                    if (is_planar[static_cast<std::size_t>(b)] && off_plane(0, fa, b)) continue;
                    return true;
                }
            return false;
        };
        // Pairs of surfaces that carry faces of the two solids, not in one
        // plane (those are overlaid in step 5, not intersected): counted, and
        // the ones whose faces have meeting boxes found through an index of
        // the face boxes of solid 1. Processed in increasing (x, y) order.
        {
            auto pairs = [](long long a, long long b, long long m) { return a * b - m - m * (m - 1) / 2; };
            long long total = 0;
            std::map<int, std::array<long long, 3>> per_class; // |A|, |B|, |A and B| per class
            long long na = 0, nb = 0, nab = 0;
            for (int g = 0; g < nG; ++g) {
                const bool a = !faces_on[0][g].empty(), b = !faces_on[1][g].empty();
                na += a; nb += b; nab += a && b;
                auto& c = per_class[rep[static_cast<std::size_t>(g)]];
                c[0] += a; c[1] += b; c[2] += a && b;
            }
            total = pairs(na, nb, nab);
            for (const auto& kv : per_class) total -= pairs(kv.second[0], kv.second[1], kv.second[2]);
            im.stats.carrier_pairs = static_cast<int>(total);
        }
        std::set<std::pair<int, int>> candidates;
        {
            const geom::BoxIndex index(face_box[1]);
            for (int fa = 0; fa < static_cast<int>(face_box[0].size()); ++fa) {
                const int x = gid(0, S(0).faces[static_cast<std::size_t>(fa)].carrier);
                for (int fb : index.query(face_box[0][static_cast<std::size_t>(fa)])) {
                    const int y = gid(1, S(1).faces[static_cast<std::size_t>(fb)].carrier);
                    if (x == y || equiv(x, y)) continue;
                    candidates.insert({std::min(x, y), std::max(x, y)});
                }
            }
        }
        for (const auto& [x, y] : candidates) {
            {
                const bool xy = !faces_on[0][x].empty() && !faces_on[1][y].empty();
                const bool yx = !faces_on[0][y].empty() && !faces_on[1][x].empty();
                if (!((xy && interacts(x, y)) || (yx && interacts(y, x)))) continue;
                ++im.stats.carrier_pairs_solved;
                const std::string where = "surfaces " + str(x) + " and " + str(y);
                // Two flat surfaces in one plane that are not one class: their
                // parameterisations are not related affinely, which the overlay
                // needs. Named, so that it is not mistaken for tangency.
                auto one_plane = [&]() {
                    return same_plane_patches(im.carriers[static_cast<std::size_t>(x)], im.carriers[static_cast<std::size_t>(y)]);
                };
                const char* kOnePlane = ": the two surfaces lie in one plane, but their parameterisations are not "
                                        "related affinely (a rational flat face, such as a sector's annular cap, "
                                        "against another flat face); unsupported";
                auto sys = std::make_shared<SsiSystem>(im.carriers[x], im.carriers[y]);
                // Search only where faces are: the hull of the face rectangles on
                // each surface. Curve parts outside it lie outside every face and
                // would be clipped away; completeness is proven over this box.
                solve::SsiOptions sopt = opt.ssi;
                if (opt.face_domains) {
                    const auto& rx = surface_rect[static_cast<std::size_t>(x)];
                    const auto& ry = surface_rect[static_cast<std::size_t>(y)];
                    if (rx[0] <= rx[1] && rx[2] <= rx[3] && ry[0] <= ry[1] && ry[2] <= ry[3]) {
                        sopt.dom_lo = {rx[0], rx[2], ry[0], ry[2]};
                        sopt.dom_hi = {rx[1], rx[3], ry[1], ry[3]};
                    }
                }
                std::vector<SsiCurve> known;
                bool shared_a[4] = {false, false, false, false}, shared_b[4] = {false, false, false, false};
                for (int sa = 0; sa < 4; ++sa)
                    for (int sb = 0; sb < 4; ++sb) {
                        bool rev = false;
                        if (!shared_side(im.carriers[x], sa, im.carriers[y], sb, rev)) continue;
                        SsiCurve c;
                        if (!side_curve(*sys, sa, sb, rev, opt.ssi, c))
                            return fail(where + ": they share a side but do not meet transversally along it "
                                                "(tangent along the side; not in general position)");
                        known.push_back(std::move(c));
                        shared_a[sa] = shared_b[sb] = true;
                        ++im.stats.shared_sides;
                    }
                // A straight side of one lying exactly in the plane of the other
                // is a known line of the intersection, on the boundary of one
                // domain at least, where subdivision alone cannot finish.
                bool touch_only = false; // proven: the surfaces meet only along that side
                for (int which = 0; which < 2 && !touch_only; ++which) {
                    const int a = which == 0 ? x : y, b = which == 0 ? y : x;
                    if (!is_planar[static_cast<std::size_t>(b)]) continue;
                    for (int s = 0; s < 4; ++s) {
                        if ((which == 0 ? shared_a : shared_b)[s]) continue;
                        const PolyPatchD& A = im.carriers[static_cast<std::size_t>(a)];
                        const PolyPatchD& B = im.carriers[static_cast<std::size_t>(b)];
                        if (!side_in_plane(A, s, B)) continue;
                        const PolyPatchI Ai = geom::to_interval(A);
                        SsiCurve c;
                        const bool ok = known_curve(*sys, [&](double tau) {
                            double ua, va;
                            side_uv(s, tau, ua, va);
                            const auto h = Ai.eval(Interval(ua), Interval(va));
                            const geom::Vec3<Interval> P{h.x / h.w, h.y / h.w, h.z / h.w};
                            Interval ub(0.5), vb(0.5);
                            plane_params(B, P, ub, vb);
                            return which == 0 ? Pt4{ua, va, num::mid(ub), num::mid(vb)} : Pt4{num::mid(ub), num::mid(vb), ua, va};
                        }, opt.ssi, c);
                        if (!ok) {
                            // No certified tube: the surfaces touch along the side.
                            // Decided exactly instead: if a meets b's plane in that
                            // side and nowhere else, the side is the whole
                            // intersection of the two surfaces. It is an edge of a's
                            // faces already (or outside them), so the pair adds no
                            // curve; the side's vertices lie in b's plane.
                            const SidesContact sc = one_plane() ? SidesContact{} : sides_contact(A, B);
                            if (sc.only_sides && (sc.mask >> s & 1)) {
                                touch_only = true;
                                ++im.stats.tangent_sides;
                                note(where + ": surface " + str(a) + " meets the plane of " + str(b) + " only along its side" +
                                     (sc.mask == (1 << s) ? " " + str(s) : "s of mask " + str(sc.mask)) + " (exact)");
                                for (int t = 0; t < 4; ++t)
                                    if (sc.mask >> t & 1) on_plane_side(a, t, b);
                                break;
                            }
                            return fail(where + (one_plane() ? std::string(kOnePlane)
                                                             : std::string(": a side of one lies in the other but they do not meet "
                                                                           "transversally along it (tangent along it; not in general position)")));
                        }
                        known.push_back(std::move(c));
                        ++im.stats.shared_sides;
                    }
                }
                if (touch_only) continue;
                // A curved surface touching the plane along a line inside its
                // patch (a plane tangent to a cylinder along a ruling), proven
                // exactly: the pair meets only on that line and does not cross
                // there, so it contributes no curve.
                if (known.empty() && is_planar[static_cast<std::size_t>(x)] != is_planar[static_cast<std::size_t>(y)]) {
                    const int a = is_planar[static_cast<std::size_t>(x)] ? y : x, b = a == x ? y : x;
                    const PlaneLineContact lc = line_contact(im.carriers[static_cast<std::size_t>(a)], im.carriers[static_cast<std::size_t>(b)]);
                    if (lc.found) {
                        line_touches.push_back({a, b, lc.fixed_u, lc.r});
                        ++im.stats.line_tangencies;
                        note(where + ": surface " + str(a) + " touches the plane of " + str(b) + " along its line " +
                             (lc.fixed_u ? "u" : "v") + " = " + std::to_string(num::mid(lc.r)) + " (exact)");
                        continue;
                    }
                }
                bool iso_side_tangent = false;
                // A straight side of a lying in a curved surface b along one of
                // b's lines of constant parameter (proven exactly,
                // find_iso_line): a known line of the intersection, which may run
                // along the boundary of a's square only.
                for (int which = 0; which < 2; ++which) {
                    const int a = which == 0 ? x : y, b = which == 0 ? y : x;
                    if (is_planar[static_cast<std::size_t>(b)]) continue;
                    const PolyPatchD& A = im.carriers[static_cast<std::size_t>(a)];
                    for (int s = 0; s < 4; ++s) {
                        if ((which == 0 ? shared_a : shared_b)[s]) continue;
                        const auto cp = side_control(A, s);
                        IsoLine line;
                        if (!find_iso_line(cp, b, line)) continue;
                        const double la[3] = {line.A.x, line.A.y, line.A.z}, lb[3] = {line.B.x, line.B.y, line.B.z};
                        const double p0[3] = {cp.front().x, cp.front().y, cp.front().z}, p1[3] = {cp.back().x, cp.back().y, cp.back().z};
                        int c = 0;
                        for (int i = 1; i < 3; ++i)
                            if (std::fabs(lb[i] - la[i]) > std::fabs(lb[c] - la[c])) c = i;
                        const double fixed = line.a / line.b;
                        SsiCurve curve;
                        const bool ok = known_curve(*sys, [&](double tau) {
                            double ua, va;
                            side_uv(s, tau, ua, va);
                            const double xc = p0[c] + tau * (p1[c] - p0[c]);
                            const double f = (line.A.w * xc - la[c]) / (lb[c] - la[c]);
                            const double ub = line.fixed_u ? fixed : f, vb = line.fixed_u ? f : fixed;
                            return which == 0 ? Pt4{ua, va, ub, vb} : Pt4{ub, vb, ua, va};
                        }, opt.ssi, curve);
                        if (!ok) {
                            // Tangent somewhere along it: two graphs crossing where they
                            // share this line, decided below by exact factorisation
                            // (graph_crossing), or reported then.
                            iso_side_tangent = true;
                            continue;
                        }
                        known.push_back(std::move(curve));
                        (which == 0 ? shared_a : shared_b)[s] = true;
                        ++im.stats.shared_sides;
                        note(where + ": side " + str(s) + " of surface " + str(a) + " lies in surface " + str(b) +
                             " along its line " + (line.fixed_u ? "u" : "v") + " = " + std::to_string(fixed) + " (exact)");
                    }
                }
                // The same for a curved surface b: a side of a shared by
                // construction with a surface z that is the same surface as b
                // (an affine reparameterisation, proven exactly) lies on b; b's
                // parameters along it come from z's through the map.
                for (int which = 0; which < 2; ++which) {
                    const int a = which == 0 ? x : y, b = which == 0 ? y : x;
                    if (is_planar[static_cast<std::size_t>(b)]) continue;
                    const PolyPatchD& A = im.carriers[static_cast<std::size_t>(a)];
                    for (int s = 0; s < 4; ++s) {
                        if ((which == 0 ? shared_a : shared_b)[s]) continue;
                        int z = -1, sz = -1;
                        bool rev = false;
                        for (int h : same_plane(b)) {
                            if (h == b || equiv(h, a)) continue;
                            for (int t = 0; t < 4 && z < 0; ++t)
                                if (shared_side(A, s, im.carriers[static_cast<std::size_t>(h)], t, rev)) { z = h; sz = t; }
                            if (z >= 0) break;
                        }
                        if (z < 0) continue;
                        const Affine& zb = map(z, b);
                        SsiCurve c;
                        const bool ok = known_curve(*sys, [&](double tau) {
                            double ua, va, uz, vz;
                            side_uv(s, tau, ua, va);
                            side_uv(sz, rev ? 1.0 - tau : tau, uz, vz);
                            Interval ub, vb;
                            zb.apply(Interval(uz), Interval(vz), ub, vb);
                            return which == 0 ? Pt4{ua, va, num::mid(ub), num::mid(vb)} : Pt4{num::mid(ub), num::mid(vb), ua, va};
                        }, opt.ssi, c);
                        if (!ok)
                            return fail(where + ": a side of one lies on the other but they do not meet "
                                                "transversally along it (not in general position)");
                        known.push_back(std::move(c));
                        (which == 0 ? shared_a : shared_b)[s] = true;
                        ++im.stats.shared_sides;
                    }
                }
                const auto ssi_t0 = std::chrono::steady_clock::now();
                solve::SsiResult res;
                const bool direct = opt.plane_lines && known.empty() && is_planar[static_cast<std::size_t>(x)] &&
                                    is_planar[static_cast<std::size_t>(y)] && plane_pair(x, y, *sys, sopt, res);
                if (direct) ++im.stats.direct_plane_pairs;
                else res = solve::intersect(*sys, sopt, known);
                im.stats.ssi_verify_nodes += res.stats.verify_nodes;
                {
                    const double ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - ssi_t0).count();
                    bool any = false;
                    for (const auto& c : res.curves) any = any || !c.known;
                    if (any) { im.stats.ssi_ms_curves += ms; im.stats.ssi_nodes_curves += res.stats.verify_nodes; }
                    else { im.stats.ssi_ms_empty += ms; im.stats.ssi_nodes_empty += res.stats.verify_nodes; ++im.stats.ssi_pairs_empty; }
                }
                if (!res.complete && is_planar[static_cast<std::size_t>(x)] && is_planar[static_cast<std::size_t>(y)]) {
                    // Two planes meet in exactly one line. Tracing cannot follow it
                    // when it only touches the domains (at a corner, say); the line
                    // can be certified directly instead, and the completeness proof
                    // then covers the touching point.
                    SsiCurve line;
                    if (plane_line(x, y, *sys, line)) {
                        std::vector<SsiCurve> with_line = known;
                        with_line.push_back(std::move(line));
                        res = solve::intersect(*sys, sopt, with_line);
                        im.stats.ssi_verify_nodes += res.stats.verify_nodes;
                        if (res.complete) res.curves[known.size()].known = false; // a curve to clip, not a shared side
                        ++im.stats.plane_lines;
                    }
                }
                if (!res.complete && known.empty() && opt.boundary_margin > 0.0 && !one_plane()) {
                    // A curve that meets both parameter squares only on their
                    // boundaries (through a common corner of two faces, say) cannot
                    // be traced inside them, and the completeness proof cannot
                    // close around that point. Over a slightly larger box the curve
                    // runs through it, its tube covers it, and completeness there
                    // implies completeness over the original box. The parts outside
                    // the squares lie outside every face and are clipped away.
                    // (Not with shared sides: their continuations beyond the
                    // squares would be traced as curves along the known ones.)
                    solve::SsiOptions wide = sopt;
                    // Where it works it needs few nodes (under a thousand on the
                    // benchmarks); a pair that is not of this kind is cut short.
                    wide.verify_node_limit = std::min<std::size_t>(wide.verify_node_limit, opt.boundary_node_limit);
                    // Tracing towards a true tangency creeps on with ever shorter
                    // steps (tens of thousands of segments, each then checked by
                    // the completeness pass); the curves this search is for take
                    // a few dozen.
                    wide.trace_steps = std::min(wide.trace_steps, opt.boundary_trace_steps);
                    for (std::size_t i = 0; i < 4; ++i) {
                        wide.dom_lo[i] -= opt.boundary_margin;
                        wide.dom_hi[i] += opt.boundary_margin;
                    }
                    const auto wide_t0 = std::chrono::steady_clock::now();
                    solve::SsiResult r2 = solve::intersect(*sys, wide);
                    const double wide_ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - wide_t0).count();
                    im.stats.ssi_verify_nodes += r2.stats.verify_nodes;
                    std::size_t wide_segs = 0;
                    for (const auto& c : r2.curves) wide_segs += c.segs.size();
                    note(where + ": search over a box enlarged by " + std::to_string(opt.boundary_margin) +
                         (r2.complete ? " proved it complete" : " did not finish") + " (" + str(static_cast<int>(r2.curves.size())) +
                         " curves, " + str(static_cast<int>(wide_segs)) + " segments, " + str(static_cast<int>(r2.stats.steps_rejected)) +
                         " steps rejected, " +
                         str(static_cast<int>(r2.stats.verify_nodes)) + " nodes, " + str(static_cast<int>(wide_ms)) + " ms; first search " +
                         str(static_cast<int>(res.stats.verify_nodes)) + " nodes)");
                    if (r2.complete) {
                        res = std::move(r2);
                        ++im.stats.boundary_pairs;
                    }
                }
                // Two graphs crossing where they share boundary lines: decided by
                // exact factorisation (graph_crossing), with the crossing curve
                // traced against a shifted surface on which it is transversal.
                bool swap_xy = false; // the curves' coordinates are y's, then x's
                if (!res.complete && known.empty()) {
                    graph_crossing_refusal.clear();
                    if (!graph_crossing(x, y, sys, sopt, res, where) && graph_crossing_refusal.empty() &&
                        graph_crossing(y, x, sys, sopt, res, where))
                        swap_xy = true;
                    if (!res.complete && !graph_crossing_refusal.empty()) return fail(where + ": " + graph_crossing_refusal);
                }
                // A curved surface touching a plane at one point, inside the
                // region the search could not finish: proven exactly.
                if (!res.complete && known.empty() &&
                    (point_touch(x, y, res, where) || curved_touch(x, y, res, where) || near_miss(x, y, res, where)))
                    res.complete = true;
                if (!res.complete && iso_side_tangent)
                    return fail(where + ": a side of one lies in the other along one of its lines but they do not "
                                        "meet transversally along it (tangent along it; not in general position)");
                if (!res.complete && opt.trace) {
                    // Where the search could not finish: the hull of the boxes
                    // left, its position in space, and the surfaces' kinds.
                    const Box4& H = res.unverified_hull;
                    const auto pa = sys->point_a(H), pb = sys->point_b(H);
                    std::string s = where + ": not complete; surfaces " + sk(x) + " (" + str(im.carriers[x].p) + "x" +
                                    str(im.carriers[x].q) + "), " + sk(y) + " (" + str(im.carriers[y].p) + "x" +
                                    str(im.carriers[y].q) + "); " + str(static_cast<int>(res.unverified_boxes)) + " boxes in";
                    for (int i = 0; i < 4; ++i)
                        s += " [" + std::to_string(H[static_cast<std::size_t>(i)].lo) + ", " + std::to_string(H[static_cast<std::size_t>(i)].hi) + "]";
                    s += "; at (" + std::to_string(num::mid(pa.x)) + ", " + std::to_string(num::mid(pa.y)) + ", " +
                         std::to_string(num::mid(pa.z)) + ") / (" + std::to_string(num::mid(pb.x)) + ", " +
                         std::to_string(num::mid(pb.y)) + ", " + std::to_string(num::mid(pb.z)) + "); " +
                         str(static_cast<int>(res.curves.size())) + " curves found";
                    note(s);
                }
                if (!res.complete)
                    return fail(where + (one_plane() ? std::string(kOnePlane)
                                                     : std::string(": the intersection could not be proven complete "
                                                                   "(tangential contact, overlapping surfaces, or contact on the "
                                                                   "boundary of both; not in general position)")));
                for (SsiCurve& c : res.curves) {
                    if (c.known) continue; // a shared side: its pieces are operand edges already
                    if (c.end_front == solve::EndKind::Singular || c.end_back == solve::EndKind::Singular)
                        return fail(where + ": an intersection curve ends at a suspected tangency");
                    im.stats.segments += static_cast<int>(c.segs.size());
                    if (c.closed) ++im.stats.closed_curves;
                    Curve cu;
                    cu.x = swap_xy ? y : x;
                    cu.y = swap_xy ? x : y;
                    auto shared_c = std::make_shared<const SsiCurve>(std::move(c));
                    if (!arc_sides(*sys, *shared_c, cu.left_behind))
                        return fail(where + ": could not prove on which side of the intersection curve "
                                            "each surface lies");
                    cu.base = make_curve_base(sys, shared_c);
                    im.curves.push_back(cu.base);
                    curves.push_back(std::move(cu));
                }
            }
        }
        im.stats.curves = static_cast<int>(curves.size());
        return true;
    }

    // --- 2. edge/surface ------------------------------------------------------------
    /// An existing vertex that is this root of "side of g0 against g": it lies
    /// on both by construction and inside the root's uniqueness box.
    bool identify_side_root(int g0, int side0, int g1, int g, const solve::Box3& unique) const {
        for (int v = 0; v < static_cast<int>(im.vertices.size()); ++v) {
            const Chart* c = chart(v, g);
            if (!c) continue;
            Interval tau;
            if (vertex_on_side(v, g0, side0, g1, tau) != Tri::Yes) continue;
            if (num::subset(tau, unique[0]) && num::subset(c->u, unique[1]) && num::subset(c->v, unique[2]))
                return true;
        }
        return false;
    }
    /// The same for a three-surface root on (g0, g1, g).
    bool identify_root6(int g0, int g1, int g, const Box6& unique) const {
        for (int v = 0; v < static_cast<int>(im.vertices.size()); ++v) {
            const Chart* a = chart(v, g0);
            const Chart* b = chart(v, g1);
            const Chart* c = chart(v, g);
            if (!a || !b || !c) continue;
            const Interval box[6] = {a->u, a->v, b->u, b->v, c->u, c->v};
            bool in = true;
            for (int i = 0; i < 6 && in; ++i) in = num::subset(box[i], unique[static_cast<std::size_t>(i)]);
            if (in) return true;
        }
        return false;
    }
    // --- exact points ---------------------------------------------------------------
    // A root whose parameters are dyadic numbers can be known exactly: the
    // candidate (the enclosure's midpoint snapped to a dyadic grid) is put
    // into every surface with exact arithmetic, and if the surfaces give one
    // and the same point there, that point is a common point inside the
    // root's uniqueness box, hence the root. Its parameters are then points,
    // not enclosures, and its position is an exact homogeneous point. With
    // that, "does this vertex lie on that surface" is decided exactly: the
    // coincidences that models with round dimensions are full of.

    /// The exact point of a vertex: recorded when the vertex was made, or the
    /// point of a surface on which it has exact parameters, if representable.
    bool exact_point(int v, Vec4d& E) const {
        if (!opt.exact_points) return false;
        if (im.vertices[static_cast<std::size_t>(v)].exact) {
            E = im.vertices[static_cast<std::size_t>(v)].xp;
            return true;
        }
        for (const Chart& c : im.vertices[static_cast<std::size_t>(v)].charts)
            if (is_point(c.u) && is_point(c.v) &&
                exact_point_at(im.carriers[static_cast<std::size_t>(c.carrier)], c.u.lo, c.v.lo, E))
                return true;
        return false;
    }
    /// Parameters on surface h of the exact point E, if E lies exactly on h.
    /// A plane: an exact test, then the parameters as enclosures, made exact
    /// where they are dyadic. Any other surface: the candidate is the hint's
    /// midpoint snapped to the dyadic grid, and h's exact point there must be E.
    bool exact_on_surface(const Vec4d& E, int h, const Interval& hint_u, const Interval& hint_v, Chart& out) const {
        const PolyPatchD& P = im.carriers[static_cast<std::size_t>(h)];
        out.carrier = h;
        Vec4d X;
        if (is_planar[static_cast<std::size_t>(h)]) {
            if (!planar_detail::on_plane(P, E)) return false;
            const Interval w(E.w);
            if (!plane_params(P, {Interval(E.x) / w, Interval(E.y) / w, Interval(E.z) / w}, out.u, out.v)) return false;
            const double us = snap_dyadic(num::mid(out.u)), vs = snap_dyadic(num::mid(out.v));
            if (num::contains(out.u, us) && num::contains(out.v, vs) && exact_point_at(P, us, vs, X) &&
                planar_detail::same_point_h(X, E)) {
                out.u = Interval(us);
                out.v = Interval(vs);
            } else {
                // One parameter alone may be dyadic (v = 1 on a side, u a
                // third): E lies on the parallelogram's line at that value,
                // between its exact points there at the other parameter 0 and 1.
                for (int c = 0; c < 2; ++c) {
                    const double val = c == 0 ? us : vs;
                    if (!num::contains(c == 0 ? out.u : out.v, val)) continue;
                    Vec4d A, B;
                    if (!(c == 0 ? exact_point_at(P, val, 0.0, A) && exact_point_at(P, val, 1.0, B)
                                 : exact_point_at(P, 0.0, val, A) && exact_point_at(P, 1.0, val, B)))
                        continue;
                    if (A.w != 1.0 || B.w != 1.0 || !on_line(A, B, E)) continue;
                    (c == 0 ? out.u : out.v) = Interval(val);
                }
            }
            return true;
        }
        const double us = snap_dyadic(num::mid(hint_u)), vs = snap_dyadic(num::mid(hint_v));
        if (exact_point_at(P, us, vs, X) && planar_detail::same_point_h(X, E)) {
            out.u = Interval(us);
            out.v = Interval(vs);
            return true;
        }
        // Parameters that are rational but not dyadic (a point of a circle
        // at the parameter 1/3): proven the same way; they stay enclosures.
        double au = 0.0, bu = 1.0, av = 0.0, bv = 1.0;
        if (!rational_near(num::mid(hint_u), au, bu) || !rational_near(num::mid(hint_v), av, bv)) return false;
        if (!exact_point_at(P, au, bu, av, bv, X) || !planar_detail::same_point_h(X, E)) return false;
        out.u = Interval(au) / Interval(bu);
        out.v = Interval(av) / Interval(bv);
        return true;
    }
    /// Parameters on surface h of the exact point E, if E lies exactly on h,
    /// without a hint: on a curved surface the candidate comes from locating
    /// E on h in doubles.
    bool exact_chart(const Vec4d& E, int h, Chart& out) const {
        if (is_planar[static_cast<std::size_t>(h)]) return exact_on_surface(E, h, Interval(0.0), Interval(0.0), out);
        const geom::BezierSurface Sh(im.carriers[static_cast<std::size_t>(h)]);
        const Vec3d x{E.x / E.w, E.y / E.w, E.z / E.w};
        double u = 0.0, v = 0.0;
        if (!surfaces_detail::locate_on_patch(Sh, x, 1.0 + std::fabs(x.x) + std::fabs(x.y) + std::fabs(x.z), u, v)) return false;
        return exact_on_surface(E, h, Interval(u), Interval(v), out);
    }
    /// Does the exact point E lie on edge e of solid k, strictly between its
    /// ends? It must lie exactly on both of the edge's surfaces; the charts
    /// there are returned. A side: its position along the side, from the
    /// exact parameter across it. An arc: its position on the certified curve.
    bool exact_on_edge(const Vec4d& E, int k, int e, Chart out[2]) {
        const BEdge& edge = S(k).edges[static_cast<std::size_t>(e)];
        int g0 = -1, g1 = -1;
        edge_carriers(k, e, g0, g1);
        if (!exact_chart(E, g0, out[0]) || !exact_chart(E, g1, out[1])) return false;
        if (!edge.arc) {
            Interval tau;
            if (on_side(out[0], edge.on[0].side, tau) != Tri::Yes) return false;
            const Interval te = edge.on[0].reversed ? Interval(1.0) - tau : tau;
            return edge.t0.hi < te.lo && te.hi < edge.t1.lo;
        }
        bool on = false;
        CurvePos pos;
        if (!locate_on_curve(*edge.geom->base, {out[0].u, out[0].v, out[1].u, out[1].v}, on, pos) || !on) return false;
        return in_arc(*edge.geom, pos) == Tri::Yes;
    }

    /// Is vertex v the root that is the only common point of the surfaces
    /// gs[0..n-1] inside `box` (two parameters per surface)? Proven when the
    /// vertex has parameters inside the box on every one of them: its own
    /// enclosures, or, on a surface where it has none, those of its exact
    /// point, which must lie exactly on that surface (`enc`, the root's
    /// enclosure, gives the candidate). The vertex is then a common point
    /// inside the box, so it is the root, and it gets the parameters it lacked.
    bool exact_vertex_is_root(int v, const int* gs, int n, const Interval* enc, const Interval* box) {
        std::vector<Chart> add;
        Vec4d E;
        bool have = false, tried = false;
        for (int i = 0; i < n; ++i) {
            Chart c;
            if (const Chart* mine = chart(v, gs[i])) {
                c = *mine;
            } else {
                if (!tried) { have = exact_point(v, E); tried = true; }
                if (!have || !exact_on_surface(E, gs[i], enc[2 * i], enc[2 * i + 1], c)) return false;
                add.push_back(c);
            }
            if (!num::subset(c.u, box[2 * i]) || !num::subset(c.v, box[2 * i + 1])) return false;
        }
        if (add.empty()) return true;
        for (const Chart& c : add) im.vertices[static_cast<std::size_t>(v)].charts.push_back(c);
        spread(v);
        return true;
    }
    /// The same for a root (t, u, v) of a straight edge's side against g: an
    /// end vertex of the edge whose side parameter is exact and inside the
    /// box, and which lies on g inside the box, is the root.
    bool exact_end_is_root(int k, int e, int g, const solve::CsiRoot& root) {
        const BEdge& edge = S(k).edges[static_cast<std::size_t>(e)];
        for (int end = 0; end < 2; ++end) {
            const Interval& te = end == 0 ? edge.t0 : edge.t1;
            // The end's true parameter lies in te (exact, or an enclosure such
            // as 1/3): the vertex is a common point inside the box if te is.
            if (!num::subset(te, root.unique[0])) continue;
            const int v = vmap[k][static_cast<std::size_t>(end == 0 ? edge.v0 : edge.v1)];
            const Interval enc[2] = {root.enc[1], root.enc[2]}, box[2] = {root.unique[1], root.unique[2]};
            if (exact_vertex_is_root(v, &g, 1, enc, box)) {
                note("edge " + str(e) + " of solid " + str(k) + " ends on surface " + str(g) + " at vertex " + str(v) + " (exact point)");
                ++im.stats.exact_ends;
                return true;
            }
        }
        return false;
    }
    /// An exact root of side `ref` of g0 against g. The side parameter must be
    /// dyadic; the side's exact point there must lie exactly on g, at
    /// parameters inside the uniqueness box (dyadic, or rational, or any on a
    /// plane). That point is then a solution inside the box, hence the root.
    /// `ts` is the side parameter and `on_g` the parameters on g: points
    /// where they are dyadic, enclosures otherwise.
    bool exact_side_root(int g0, const BEdge::SideRef& ref, int g, const solve::CsiRoot& root, double& ts, Chart& on_g) const {
        if (!opt.exact_points) return false;
        ts = snap_dyadic(num::mid(root.enc[0]));
        if (!num::contains(root.enc[0], ts) || !num::contains(root.unique[0], ts) || ts < 0.0 || ts > 1.0) return false;
        Interval cu, cv;
        side_uv(ref.side, Interval(ref.reversed ? 1.0 - ts : ts), cu, cv);
        Vec4d A;
        if (!is_point(cu) || !is_point(cv) || !exact_point_at(im.carriers[static_cast<std::size_t>(g0)], cu.lo, cv.lo, A)) return false;
        return exact_on_surface(A, g, root.enc[1], root.enc[2], on_g) && num::subset(on_g.u, root.unique[1]) &&
               num::subset(on_g.v, root.unique[2]);
    }
    /// A three-surface root that is a point with dyadic coordinates although
    /// its parameters are dyadic on none of the surfaces (the point
    /// (1, 0.75, 0) at a third of an edge and a third of a circle). The
    /// candidate X is the root's position `at` rounded to the dyadic grid; it
    /// must lie exactly on each surface at parameters inside the uniqueness
    /// box. `out` are those parameters, as tight enclosures.
    bool exact_space_root(const int gs[3], const TsiRoot& root, const Vec3<Interval>& at, Vec4d& X, Chart out[3]) const {
        if (!opt.exact_points) return false;
        X = {snap_dyadic(num::mid(at.x)), snap_dyadic(num::mid(at.y)), snap_dyadic(num::mid(at.z)), 1.0};
        if (!num::contains(at.x, X.x) || !num::contains(at.y, X.y) || !num::contains(at.z, X.z)) return false;
        for (int j = 0; j < 3; ++j) {
            const std::size_t ij = static_cast<std::size_t>(2 * j);
            if (!exact_on_surface(X, gs[j], root.enc[ij], root.enc[ij + 1], out[j]) || !num::subset(out[j].u, root.unique[ij]) ||
                !num::subset(out[j].v, root.unique[ij + 1]))
                return false;
        }
        return true;
    }
    /// Is a side (its control points) a straight segment in its linear
    /// parameterisation: unit weights, and every control point where degree
    /// elevation of the two ends puts it (exact)?
    static bool linear_side(const std::vector<Vec4d>& cp) {
        const int n = static_cast<int>(cp.size()) - 1;
        if (n < 1 || detail::same_point(cp.front(), cp.back())) return false;
        for (const Vec4d& c : cp)
            if (c.w != 1.0) return false;
        using E = num::Expansion;
        for (int i = 1; i < n; ++i) {
            const double a[3] = {cp[static_cast<std::size_t>(i)].x, cp[static_cast<std::size_t>(i)].y, cp[static_cast<std::size_t>(i)].z};
            const double p0[3] = {cp.front().x, cp.front().y, cp.front().z}, p1[3] = {cp.back().x, cp.back().y, cp.back().z};
            for (int c = 0; c < 3; ++c)
                if ((E(static_cast<double>(n)) * E(a[c]) - E(static_cast<double>(n - i)) * E(p0[c]) - E(static_cast<double>(i)) * E(p1[c])).sign() != 0)
                    return false;
        }
        return true;
    }
    /// The same for a root of a straight edge's side against g: X must lie
    /// exactly on the side's line, at a side parameter `s` inside the box,
    /// and exactly on g at parameters inside the box.
    bool exact_space_side_root(int g0, const BEdge::SideRef& ref, int g, const solve::CsiRoot& root, const Vec3<Interval>& at,
                               Vec4d& X, Interval& s, Chart& on_g) const {
        if (!opt.exact_points) return false;
        X = {snap_dyadic(num::mid(at.x)), snap_dyadic(num::mid(at.y)), snap_dyadic(num::mid(at.z)), 1.0};
        if (!num::contains(at.x, X.x) || !num::contains(at.y, X.y) || !num::contains(at.z, X.z)) return false;
        const auto cp = side_control(im.carriers[static_cast<std::size_t>(g0)], ref.side);
        if (!linear_side(cp) || !collinear(cp.front(), cp.back(), X)) return false;
        const double p0[3] = {cp.front().x, cp.front().y, cp.front().z}, p1[3] = {cp.back().x, cp.back().y, cp.back().z};
        const double x[3] = {X.x, X.y, X.z};
        int c = 0;
        for (int i = 1; i < 3; ++i)
            if (std::fabs(p1[i] - p0[i]) > std::fabs(p1[c] - p0[c])) c = i;
        s = (Interval(x[c]) - Interval(p0[c])) / (Interval(p1[c]) - Interval(p0[c]));
        if (!num::subset(s, root.unique[0])) return false;
        return exact_on_surface(X, g, root.enc[1], root.enc[2], on_g) && num::subset(on_g.u, root.unique[1]) &&
               num::subset(on_g.v, root.unique[2]);
    }
    // --- two graphs crossing where they share a boundary line ------------------------
    // Two polynomial patches that are graphs over the xy plane, with x affine
    // in one parameter and y in the other (the tops of boxes with a hill), can
    // share a boundary line: both are at the same height there. Then the
    // difference of their heights, D, vanishes on that line, and the curve
    // where they cross elsewhere meets the line at points where the two
    // surfaces are tangent. No interval certificate covers such a point.
    // Exact algebra does: in A's parameters D is a polynomial with exact
    // coefficients (B's heights restricted to A's square by rational
    // blossoming); every side of A's square on which D vanishes is divided
    // out exactly, leaving k with D = (side factors) k, the factors positive
    // inside the square. The crossing curve is k = 0, and it is the
    // transversal intersection of A with C = B + kappa e_z, where kappa is
    // chosen so that A's height minus C's is a positive multiple of k (C is
    // kept as an interval enclosure of its exact control points). Tracing A
    // against C certifies the curve, ends included; its coordinates are A's
    // and B's parameters (C has B's x and y), its points lie on A and on B
    // (where k = 0, C = B), and the side on which each surface lies agrees
    // with A against B (both are signs of the same height difference times a
    // positive factor). The shared lines are edges of the operands or lie
    // outside their faces; their crossings with the curve become vertices.

    struct Rat { double a = 0.0, b = 1.0; }; ///< a / b, b a positive integer
    bool graph_crossing_seen = false; ///< a pair was decided by graph_crossing
    std::string graph_crossing_refusal; ///< why graph_crossing refused the last pair, if it decided to

    /// x = x0 + x1 u and y = y0 + y1 v exactly, weights 1? (A graph patch.)
    static bool graph_axes(const PolyPatchD& P, double x[2], double y[2]) {
        using E = num::Expansion;
        for (const Vec4d& c : P.c)
            if (c.w != 1.0) return false;
        x[0] = P.at(0, 0).x;
        y[0] = P.at(0, 0).y;
        const E dx = E(P.at(P.p, 0).x) - E(x[0]), dy = E(P.at(0, P.q).y) - E(y[0]);
        if (dx.size() > 1 || dy.size() > 1) return false;
        x[1] = dx.estimate();
        y[1] = dy.estimate();
        if (x[1] == 0.0 || y[1] == 0.0) return false;
        for (int i = 0; i <= P.p; ++i)
            for (int j = 0; j <= P.q; ++j) {
                const Vec4d& c = P.at(i, j);
                if ((E(static_cast<double>(P.p)) * E(c.x) - E(static_cast<double>(P.p)) * E(x[0]) - E(static_cast<double>(i)) * E(x[1])).sign() != 0)
                    return false;
                if ((E(static_cast<double>(P.q)) * E(c.y) - E(static_cast<double>(P.q)) * E(y[0]) - E(static_cast<double>(j)) * E(y[1])).sign() != 0)
                    return false;
            }
        return true;
    }
    /// (p - p0) / d as an exact rational with a small denominator, if it is one.
    static bool ratio(double p, double p0, double d, Rat& r) {
        using E = num::Expansion;
        double a = 0.0, b = 1.0;
        if (!rational_near((p - p0) / d, a, b) || b > 64.0) return false;
        if ((E(a) * E(d) - E(b) * (E(p) - E(p0))).sign() != 0) return false;
        r = {a, b};
        return true;
    }
    /// Coefficients c (degree n) of a polynomial on [0, 1], re-expressed on
    /// [x0, x1] for rationals x0, x1: coefficient k is the blossom at (n - k)
    /// copies of x0 and k of x1, computed in homogeneous form and brought to
    /// the common positive factor (b0 b1)^n.
    static std::vector<num::Expansion> restrict_rat(const std::vector<num::Expansion>& c, Rat x0, Rat x1) {
        using E = num::Expansion;
        const int n = static_cast<int>(c.size()) - 1;
        std::vector<E> out(c.size());
        for (int k = 0; k <= n; ++k) {
            std::vector<E> pts = c;
            for (int r = 1; r <= n; ++r) {
                const Rat x = r <= k ? x1 : x0;
                for (int i = 0; i + r <= n; ++i)
                    pts[static_cast<std::size_t>(i)] = E(x.b - x.a) * pts[static_cast<std::size_t>(i)] + E(x.a) * pts[static_cast<std::size_t>(i + 1)];
            }
            out[static_cast<std::size_t>(k)] = pts[0] * E(std::pow(x0.b, k) * std::pow(x1.b, n - k));
        }
        return out;
    }
    static ExactBernstein restrict_rat2(const ExactBernstein& f, Rat u0, Rat u1, Rat v0, Rat v1) {
        ExactBernstein r = f;
        for (int j = 0; j <= f.q; ++j) {
            std::vector<num::Expansion> row;
            for (int i = 0; i <= f.p; ++i) row.push_back(r.at(i, j));
            row = restrict_rat(row, u0, u1);
            for (int i = 0; i <= f.p; ++i) r.at(i, j) = row[static_cast<std::size_t>(i)];
        }
        for (int i = 0; i <= f.p; ++i) {
            std::vector<num::Expansion> col;
            for (int j = 0; j <= f.q; ++j) col.push_back(r.at(i, j));
            col = restrict_rat(col, v0, v1);
            for (int j = 0; j <= f.q; ++j) r.at(i, j) = col[static_cast<std::size_t>(j)];
        }
        return r;
    }
    /// f raised to degree (p, q) (each step times its new degree, a positive factor).
    static ExactBernstein elevate_to(ExactBernstein f, int p, int q) {
        using E = num::Expansion;
        while (f.p < p) {
            ExactBernstein g(f.p + 1, f.q);
            const int m = f.p + 1;
            for (int j = 0; j <= f.q; ++j)
                for (int i = 0; i <= m; ++i) {
                    E x;
                    if (i > 0) x = x + E(static_cast<double>(i)) * f.at(i - 1, j);
                    if (i < m) x = x + E(static_cast<double>(m - i)) * f.at(i, j);
                    g.at(i, j) = x;
                }
            f = g;
        }
        while (f.q < q) {
            ExactBernstein g(f.p, f.q + 1);
            const int m = f.q + 1;
            for (int i = 0; i <= f.p; ++i)
                for (int j = 0; j <= m; ++j) {
                    E x;
                    if (j > 0) x = x + E(static_cast<double>(j)) * f.at(i, j - 1);
                    if (j < m) x = x + E(static_cast<double>(m - j)) * f.at(i, j);
                    g.at(i, j) = x;
                }
            f = g;
        }
        return f;
    }
    static Interval enclose(const num::Expansion& e) {
        Interval s(0.0);
        for (double x : e.components()) s = s + Interval(x);
        return s;
    }
    /// Roots of a polynomial (Bernstein coefficients, enclosed) on [lo, hi],
    /// each enclosed and proven unique by an interval Newton step. False if
    /// they could not all be isolated.
    static bool roots_1d(const std::vector<Interval>& c, double lo, double hi, int depth, std::vector<Interval>& out) {
        bool pos = true, neg = true;
        for (const Interval& x : c) { pos = pos && x.lo > 0.0; neg = neg && x.hi < 0.0; }
        if (pos || neg) return true;
        const int n = static_cast<int>(c.size()) - 1;
        // Derivative coefficients: n (c_{i+1} - c_i), scaled by 1 / (hi - lo).
        bool dpos = true, dneg = true;
        for (int i = 0; i < n; ++i) {
            const Interval d = (c[static_cast<std::size_t>(i + 1)] - c[static_cast<std::size_t>(i)]);
            dpos = dpos && d.lo > 0.0;
            dneg = dneg && d.hi < 0.0;
        }
        if (dpos || dneg) {
            // Monotone: at most one root; it exists if the ends differ in sign.
            const Interval f0 = c.front(), f1 = c.back();
            if ((f0.hi < 0.0 && f1.lo > 0.0) || (f0.lo > 0.0 && f1.hi < 0.0)) {
                // Bisect the enclosure down, keeping the sign change.
                double a = lo, b = hi;
                std::vector<Interval> cc = c;
                for (int it = 0; it < 60 && b - a > 1e-15 * std::max(1.0, std::fabs(a)); ++it) {
                    // value at the midpoint, by de Casteljau
                    std::vector<Interval> t = cc;
                    for (int r = 1; r <= n; ++r)
                        for (int i = 0; i + r <= n; ++i) t[static_cast<std::size_t>(i)] = (t[static_cast<std::size_t>(i)] + t[static_cast<std::size_t>(i + 1)]) * Interval(0.5);
                    const Interval fm = t[0];
                    if (num::contains_zero(fm)) break;
                    // split at the middle, keep the half with the sign change
                    std::vector<Interval> left(cc.size()), right(cc.size()), w = cc;
                    for (int r = 0; r <= n; ++r) {
                        left[static_cast<std::size_t>(r)] = w[0];
                        right[static_cast<std::size_t>(n - r)] = w[static_cast<std::size_t>(n - r)];
                        for (int i = 0; i + r < n; ++i) w[static_cast<std::size_t>(i)] = (w[static_cast<std::size_t>(i)] + w[static_cast<std::size_t>(i + 1)]) * Interval(0.5);
                    }
                    const double m = 0.5 * (a + b);
                    const bool left_change = (left.front().hi < 0.0) != (fm.hi < 0.0);
                    if (left_change) { cc = left; b = m; } else { cc = right; a = m; }
                }
                out.push_back(Interval(a, b));
                return true;
            }
            return !(num::contains_zero(c.front()) || num::contains_zero(c.back()));
        }
        if (depth <= 0) return false;
        std::vector<Interval> left(c.size()), right(c.size()), w = c;
        for (int r = 0; r <= n; ++r) {
            left[static_cast<std::size_t>(r)] = w[0];
            right[static_cast<std::size_t>(n - r)] = w[static_cast<std::size_t>(n - r)];
            for (int i = 0; i + r < n; ++i) w[static_cast<std::size_t>(i)] = (w[static_cast<std::size_t>(i)] + w[static_cast<std::size_t>(i + 1)]) * Interval(0.5);
        }
        // The midpoint itself must not be a root (it would be counted twice).
        if (num::contains_zero(left.back())) return false;
        const double m = 0.5 * (lo + hi);
        return roots_1d(left, lo, m, depth - 1, out) && roots_1d(right, m, hi, depth - 1, out);
    }

    /// The pair (x, y) as two graphs crossing where they share boundary
    /// lines. On success `res` holds the crossing curves, proven complete,
    /// `sys` is their system (A against C), and vertices where the curves
    /// meet the shared lines have been added.
    bool graph_crossing(int x, int y, std::shared_ptr<SsiSystem>& sys, const solve::SsiOptions& sopt, solve::SsiResult& res,
                        const std::string& where) {
        using E = num::Expansion;
        const PolyPatchD& A = im.carriers[static_cast<std::size_t>(x)];
        const PolyPatchD& B = im.carriers[static_cast<std::size_t>(y)];
        if (A.p != B.p || A.q != B.q) return false;
        double ax[2], ay[2], bx[2], by[2];
        if (!graph_axes(A, ax, ay) || !graph_axes(B, bx, by)) return false;
        // A's square in B's parameters, and B's square in A's.
        Rat s0, s1, t0, t1, u0, u1, v0, v1;
        if (!ratio(ax[0], bx[0], bx[1], s0) || !ratio(ax[0] + ax[1], bx[0], bx[1], s1) || !ratio(ay[0], by[0], by[1], t0) ||
            !ratio(ay[0] + ay[1], by[0], by[1], t1) || !ratio(bx[0], ax[0], ax[1], u0) || !ratio(bx[0] + bx[1], ax[0], ax[1], u1) ||
            !ratio(by[0], ay[0], ay[1], v0) || !ratio(by[0] + by[1], ay[0], ay[1], v1))
            return false;
        if ((E(ax[0]) + E(ax[1]) - E(ax[0] + ax[1])).sign() != 0 || (E(bx[0]) + E(bx[1]) - E(bx[0] + bx[1])).sign() != 0 ||
            (E(ay[0]) + E(ay[1]) - E(ay[0] + ay[1])).sign() != 0 || (E(by[0]) + E(by[1]) - E(by[0] + by[1])).sign() != 0)
            return false;
        // Heights as exact polynomials; D = S zA - (zB restricted to A's square).
        ExactBernstein zA(A.p, A.q), zB(B.p, B.q);
        for (int i = 0; i <= A.p; ++i)
            for (int j = 0; j <= A.q; ++j) { zA.at(i, j) = E(A.at(i, j).z); zB.at(i, j) = E(B.at(i, j).z); }
        const ExactBernstein zBA = restrict_rat2(zB, s0, s1, t0, t1);
        const double S = std::pow(s0.b * s1.b, A.p) * std::pow(t0.b * t1.b, A.q);
        ExactBernstein D(A.p, A.q);
        for (std::size_t i = 0; i < D.c.size(); ++i) D.c[i] = E(S) * zA.c[i] - zBA.c[i];
        // Divide out the sides of A's square on which D vanishes.
        ExactBernstein k = D;
        int mask = 0;
        for (bool again = true; again;) {
            again = false;
            for (int sd = 0; sd < 4; ++sd) {
                const bool across_v = sd == 0 || sd == 2;
                if ((across_v ? k.q : k.p) == 0 || !contact_detail::side_zero(k, sd)) continue;
                k = contact_detail::divide_side(k, sd);
                mask |= 1 << sd;
                again = true;
            }
        }
        if (mask == 0) return false;
        note(where + ": two graphs sharing boundary lines (sides of mask " + str(mask) + " of surface " + str(x) + ")");
        // A double root of k along a shared side is a point where the crossing
        // curve touches the shared line instead of leaving through it: both
        // tops have the same slope there. For a quadratic this is decided
        // exactly (zero discriminant, root in the closed side). A discriminant
        // that is not zero but tiny against its terms (relative 2^-40) is the
        // same picture blurred by rounding of the data (two crossings, or none,
        // closer to each other than any enclosure can separate); both refused.
        for (int sd = 0; sd < 4; ++sd) {
            if (!(mask >> sd & 1)) continue;
            const bool along_u = sd == 0 || sd == 2;
            if ((along_u ? k.p : k.q) != 2) continue;
            E c[3];
            for (int i = 0; i <= 2; ++i) c[i] = sd == 0 ? k.at(i, 0) : sd == 1 ? k.at(k.p, i) : sd == 2 ? k.at(i, k.q) : k.at(0, i);
            const int s0c = c[0].sign(), s1c = c[1].sign(), s2c = c[2].sign();
            if (s0c == 0 && s1c == 0 && s2c == 0) continue;
            // A (near) double root lies in the side iff the end values do not
            // have opposite signs and the middle coefficient opposes them.
            if (!(s0c == 0 || s2c == 0 || (s0c == s2c && s1c == -s0c))) continue;
            const E disc = c[1] * c[1] - c[0] * c[2];
            const Interval dm = enclose(disc), sz = enclose(c[1] * c[1]) + num::abs(enclose(c[0] * c[2]));
            const int ds = disc.sign();
            if (ds != 0 && !(num::abs(dm).hi <= std::ldexp(sz.lo, -40))) continue;
            graph_crossing_refusal =
                ds == 0 ? "the two tops cross along a curve that touches the boundary line they share (tangent to it at a point where both "
                          "have the same slope; not in general position)"
                : ds > 0 ? "the two tops cross along a curve that meets the boundary line they share twice, closer together than the data's "
                           "rounding can separate (a touch blurred by rounding; not in general position)"
                         : "the two tops cross along a curve that passes the boundary line they share closer than the data's rounding can "
                           "separate, without meeting it (a touch blurred by rounding; not in general position)";
            note(where + ": " + graph_crossing_refusal);
            return false;
        }
        solve::SsiResult out;
        out.complete = true;
        if (contact_detail::strict_sign(k, 12) == 0) {
            // C = B + kappa e_z with kappa, on A's square, (D - k_e) / S.
            const ExactBernstein ke = elevate_to(k, A.p, A.q);
            ExactBernstein K(A.p, A.q);
            for (std::size_t i = 0; i < K.c.size(); ++i) K.c[i] = D.c[i] - ke.c[i];
            const ExactBernstein KB = restrict_rat2(K, u0, u1, v0, v1);
            const double S2 = std::pow(u0.b * u1.b, A.p) * std::pow(v0.b * v1.b, A.q);
            PolyPatchI Ci = geom::to_interval(B);
            PolyPatchD Cm = B;
            for (int i = 0; i <= B.p; ++i)
                for (int j = 0; j <= B.q; ++j) {
                    const Interval kz = enclose(KB.at(i, j)) / (Interval(S) * Interval(S2));
                    Ci.at(i, j).z = Ci.at(i, j).z + kz;
                    Cm.at(i, j).z = num::mid(Ci.at(i, j).z);
                }
            auto sys2 = std::make_shared<SsiSystem>(A, Cm, Ci);
            out = solve::intersect(*sys2, sopt);
            im.stats.ssi_verify_nodes += out.stats.verify_nodes;
            {
                std::string s2 = where + ": the crossing curve, A against the shifted surface: " + (out.complete ? "complete" : "not complete") +
                                 ", " + str(static_cast<int>(out.curves.size())) + " curves";
                if (!out.complete) {
                    s2 += "; " + str(static_cast<int>(out.unverified_boxes)) + " boxes in";
                    for (int i = 0; i < 4; ++i)
                        s2 += " [" + std::to_string(out.unverified_hull[static_cast<std::size_t>(i)].lo) + ", " +
                              std::to_string(out.unverified_hull[static_cast<std::size_t>(i)].hi) + "]";
                    s2 += "; side mask " + str(mask) + ", k of degree " + str(k.p) + "x" + str(k.q);
                }
                note(s2);
            }
            if (!out.complete) return false;
            for (const SsiCurve& c : out.curves)
                if (c.end_front == solve::EndKind::Singular || c.end_back == solve::EndKind::Singular) return false;
            // Where the curve meets a shared side: k restricted to that side.
            for (int sd = 0; sd < 4; ++sd) {
                if (!(mask >> sd & 1)) continue;
                std::vector<Interval> c;
                const bool along_u = sd == 0 || sd == 2;
                const int n = along_u ? k.p : k.q;
                for (int i = 0; i <= n; ++i)
                    c.push_back(enclose(sd == 0 ? k.at(i, 0) : sd == 1 ? k.at(k.p, i) : sd == 2 ? k.at(i, k.q) : k.at(0, i)));
                std::vector<Interval> roots;
                if (!roots_1d(c, 0.0, 1.0, 40, roots)) return false;
                for (const Interval& tau : roots) {
                    Interval ua, va;
                    side_uv(sd, tau, ua, va);
                    // B's parameters: s = s0 + (s1 - s0) u, exactly where u is.
                    auto lin = [](Rat r0, Rat r1, const Interval& w) {
                        // At an end of A's range the value is that rational exactly:
                        // a point where it is 0 or 1 (a side of B's square too).
                        if (w.lo == w.hi && (w.lo == 0.0 || w.lo == 1.0)) {
                            const Rat r = w.lo == 0.0 ? r0 : r1;
                            if (r.a == 0.0) return Interval(0.0);
                            if (r.a == r.b) return Interval(1.0);
                            return Interval(r.a) / Interval(r.b);
                        }
                        return Interval(r0.a) / Interval(r0.b) + w * (Interval(r1.a) / Interval(r1.b) - Interval(r0.a) / Interval(r0.b));
                    };
                    const Interval sb = lin(s0, s1, ua), tb = lin(t0, t1, va);
                    if (sb.hi < 0.0 || sb.lo > 1.0 || tb.hi < 0.0 || tb.lo > 1.0) continue; // outside B
                    if (!(sb.lo >= 0.0 && sb.hi <= 1.0 && tb.lo >= 0.0 && tb.hi <= 1.0)) return false; // on B's boundary
                    // Already a vertex?
                    bool known_v = false;
                    for (int v = 0; v < static_cast<int>(im.vertices.size()) && !known_v; ++v) {
                        const Chart* ca = chart(v, x);
                        const Chart* cb = chart(v, y);
                        known_v = ca && cb && !num::disjoint(ca->u, ua) && !num::disjoint(ca->v, va);
                    }
                    if (known_v) continue;
                    const auto h = icarrier[static_cast<std::size_t>(x)].eval(ua, va);
                    std::vector<Chart> charts = {{x, ua, va}, {y, sb, tb}};
                    // On the surfaces that share this side of x (the side faces of
                    // a box below its top), and the side of y it lies on, if any.
                    auto add_sharing = [&](int g, int side, const Interval& t) {
                        for (const auto& [h2, t2, rev] : sharing_side(g, side)) {
                            bool has = false;
                            for (const Chart& c2 : charts) has = has || c2.carrier == h2;
                            if (has) continue;
                            Chart n2;
                            n2.carrier = h2;
                            side_uv(t2, rev ? Interval(1.0) - t : t, n2.u, n2.v);
                            charts.push_back(n2);
                        }
                    };
                    add_sharing(x, sd, tau);
                    for (int sb_side = 0; sb_side < 4; ++sb_side) {
                        Interval t_on;
                        if (detail::on_side(Chart{y, sb, tb}, sb_side, t_on) == Tri::Yes) add_sharing(y, sb_side, t_on);
                    }
                    im.vertices.push_back({{num::mid(h.x / h.w), num::mid(h.y / h.w), num::mid(h.z / h.w)}, std::move(charts)});
                    vkind.push_back(kCross);
                    arc_ends.push_back(0);
                    const int nv = static_cast<int>(im.vertices.size()) - 1;
                    spread(nv);
                    exact_charts(nv);
                    share_side_charts(nv);
                    ++im.stats.new_vertices;
                    {
                        const Vec3d& P = im.vertices.back().point;
                        char buf[96];
                        std::snprintf(buf, sizeof buf, " (%.9g %.9g %.9g)", P.x, P.y, P.z);
                        note(where + ": the crossing curve meets side " + str(sd) + " of surface " + str(x) + " at vertex " + str(nv) + buf);
                    }
                }
            }
            sys = sys2;
        }
        res = std::move(out);
        ++im.stats.graph_crossings;
        graph_crossing_seen = true;
        return true;
    }

    /// A curved patch that touches a plane at one point, inside the region the
    /// surface/surface search could not finish (roadmap phase 2; RESULTS.md
    /// §3, "Tangent points, certified"). The plane function g of the patch
    /// (exact Bernstein form, homogeneous, so rational patches are covered)
    /// vanishes together with its gradient at dyadic parameters p*, and g's
    /// Hessian is definite over a dyadic box U containing every box the
    /// search left: g is then strictly convex or concave on U, so p* is its
    /// only zero there. Since the search proved everything outside those
    /// boxes, the pair meets in the curves found and in that one point.
    bool point_touch(int x, int y, const solve::SsiResult& res, const std::string& where) {
        using E = num::Expansion;
        using contact_detail::restrict_dir;
        const bool px = is_planar[static_cast<std::size_t>(x)], py = is_planar[static_cast<std::size_t>(y)];
        if (px == py) return false;
        const int c = px ? y : x, pl = px ? x : y;
        const PolyPatchD& X = im.carriers[static_cast<std::size_t>(c)];
        const PolyPatchD& Y = im.carriers[static_cast<std::size_t>(pl)];
        if (!planar_affine(Y)) return false;
        const Interval hu = res.unverified_hull[px ? 2u : 0u], hv = res.unverified_hull[px ? 3u : 1u];
        if (!(hu.hi - hu.lo < 1e-3 && hv.hi - hv.lo < 1e-3)) return false; // a small region only
        const ExactBernstein g = plane_function(X, Y);
        auto deriv = [](const ExactBernstein& f, bool in_u) {
            const int n = in_u ? f.p : f.q;
            ExactBernstein d(in_u ? std::max(f.p - 1, 0) : f.p, in_u ? f.q : std::max(f.q - 1, 0));
            if (n == 0) {
                for (auto& e : d.c) e = E(0.0);
                return d;
            }
            for (int i = 0; i <= d.p; ++i)
                for (int j = 0; j <= d.q; ++j)
                    d.at(i, j) = (in_u ? f.at(i + 1, j) - f.at(i, j) : f.at(i, j + 1) - f.at(i, j)) * E(static_cast<double>(n));
            return d;
        };
        auto value = [&](const ExactBernstein& f, double u, double v) {
            return restrict_dir(restrict_dir(f, true, u, u), false, v, v).c[0];
        };
        // The box U: the hull of the boxes left, widened, with dyadic ends.
        auto down = [](double t) { return std::ldexp(std::floor(std::ldexp(t, 24)), -24); };
        auto up = [](double t) { return std::ldexp(std::ceil(std::ldexp(t, 24)), -24); };
        const double mu = std::max(hu.hi - hu.lo, std::ldexp(1.0, -20)), mv = std::max(hv.hi - hv.lo, std::ldexp(1.0, -20));
        const double u0 = down(hu.lo - mu), u1 = up(hu.hi + mu), v0 = down(hv.lo - mv), v1 = up(hv.hi + mv);
        // p*: the dyadic point of fewest bits near the middle where g and its
        // gradient vanish exactly.
        const ExactBernstein gu = deriv(g, true), gv = deriv(g, false);
        double us = 0.0, vs = 0.0;
        bool found = false;
        for (int bits = 1; bits <= 40 && !found; ++bits) {
            us = std::ldexp(std::nearbyint(std::ldexp(num::mid(hu), bits)), -bits);
            vs = std::ldexp(std::nearbyint(std::ldexp(num::mid(hv), bits)), -bits);
            if (!(us >= u0 && us <= u1 && vs >= v0 && vs <= v1)) continue;
            found = value(g, us, vs).sign() == 0 && value(gu, us, vs).sign() == 0 && value(gv, us, vs).sign() == 0;
        }
        // Definite Hessian over U, from the coefficient hulls, with a lower
        // bound m on its eigenvalues' magnitude.
        auto hull_over = [&](const ExactBernstein& f) {
            const ExactBernstein r = restrict_dir(restrict_dir(f, true, u0, u1), false, v0, v1);
            Interval h = enclose(r.c[0]);
            for (const E& e : r.c) h = num::hull(h, enclose(e));
            return h;
        };
        const ExactBernstein guu = deriv(gu, true), guv = deriv(gu, false), gvv = deriv(gv, false);
        const Interval huu = hull_over(guu), huv = hull_over(guv), hvv = hull_over(gvv);
        const Interval Hs[2][2] = {{huu, huv}, {huv, hvv}};
        double m = 0.0;
        const int sg = solve::contact_detail::definite(Hs, m);
        if (sg == 0) {
            note(where + ": the surface's curvature against the plane is not proven definite over the region left (not an isolated "
                         "contact or near miss)");
            return false;
        }
        if (!found) {
            // A near miss: Newton in doubles towards g's critical point p (a
            // double, hence dyadic); g and its gradient there exactly; and by
            // Taylor's theorem with the eigenvalue bound, g(p + d) >= g(p) +
            // grad g(p).d + m |d|^2 / 2 >= g(p) - |grad g(p)|^2 / (2 m) on U
            // (mirrored for a concave g). If that keeps g(p)'s sign, g has no
            // zero on U: the surface does not reach the plane there.
            auto dval = [&](const ExactBernstein& f, double u, double v) { return num::mid(enclose(value(f, u, v))); };
            double pu = num::mid(hu), pv = num::mid(hv);
            for (int it = 0; it < 8; ++it) {
                const double a = dval(guu, pu, pv), b = dval(guv, pu, pv), cc = dval(gvv, pu, pv);
                const double fu = dval(gu, pu, pv), fv = dval(gv, pu, pv), d = a * cc - b * b;
                if (!(std::fabs(d) > 0.0)) break;
                const double du = (cc * fu - b * fv) / d, dv = (a * fv - b * fu) / d;
                pu -= du;
                pv -= dv;
                if (std::fabs(du) + std::fabs(dv) < 1e-17) break;
            }
            if (!(pu >= u0 && pu <= u1 && pv >= v0 && pv <= v1)) return false;
            const Interval g0 = enclose(value(g, pu, pv)), g1 = enclose(value(gu, pu, pv)), g2 = enclose(value(gv, pu, pv));
            const Interval drop = (num::sqr(g1) + num::sqr(g2)) / (Interval(2.0) * Interval(m));
            const bool apart = sg > 0 ? (g0 - drop).lo > 0.0 : (g0 + drop).hi < 0.0;
            if (!apart) {
                note(where + ": the surface comes close to the plane, but is not proven to stay off it over the region left");
                return false;
            }
            ++im.stats.near_misses;
            char buf[200];
            std::snprintf(buf, sizeof buf, ": a near miss, proven exactly: surface %d stays off the plane %d over the region left "
                                           "(plane function %.3e at its closest point, against at most %.3e of drop)",
                          c, pl, num::mid(g0), drop.hi);
            note(where + buf);
            return true;
        }
        const auto h = icarrier[static_cast<std::size_t>(c)].eval(Interval(us), Interval(vs));
        PointTouch pt{c, pl, us, vs, {h.x / h.w, h.y / h.w, h.z / h.w}};
        pt.alone = res.curves.empty();
        pt.has_xp = exact_point_at(X, us, vs, pt.xp);
        point_touches.push_back(pt);
        ++im.stats.point_contacts;
        note(where + ": surface " + str(c) + " touches the plane " + str(pl) + " at one point, at parameters (" + std::to_string(us) +
             ", " + std::to_string(vs) + "), proven exactly, and nowhere else in the region left");
        return true;
    }

    /// Two curved patches touching at one point, inside the region the search
    /// could not finish: the contact is proven exact by solve::classify_contact
    /// (a regular root of the deflated system, checked at dyadic parameters in
    /// homogeneous expansion arithmetic) and the only intersection in boxes
    /// around every box the search left by solve::unique_touch_in (the normal
    /// offset from one to the other strictly convex there).
    bool curved_touch(int x, int y, const solve::SsiResult& res, const std::string& where) {
        if (is_planar[static_cast<std::size_t>(x)] || is_planar[static_cast<std::size_t>(y)]) return false;
        const Box4& H = res.unverified_hull;
        for (const Interval& h : H)
            if (!(h.hi - h.lo < 1e-3)) return false; // a small region only
        const PolyPatchD& A = im.carriers[static_cast<std::size_t>(x)];
        const PolyPatchD& B = im.carriers[static_cast<std::size_t>(y)];
        const solve::ContactResult r = solve::classify_contact(A, B, {num::mid(H[0]), num::mid(H[1]), num::mid(H[2]), num::mid(H[3]), 0.0});
        if (r.kind != solve::ContactKind::Touch) {
            note(where + ": the contact of the two surfaces there is classified as " + solve::to_string(r.kind) +
                 (r.why.empty() ? std::string() : " (" + r.why + ")"));
            return false;
        }
        auto widen = [](const Interval& h) {
            const double m = std::max(h.hi - h.lo, std::ldexp(1.0, -20));
            return Interval(h.lo - m, h.hi + m);
        };
        std::string why;
        if (!solve::unique_touch_in(A, B, r.point, widen(H[0]), widen(H[1]), widen(H[2]), widen(H[3]), why)) {
            note(where + ": an exact contact of the two surfaces, but not proven the only intersection there: " + why);
            return false;
        }
        const auto h = icarrier[static_cast<std::size_t>(x)].eval(Interval(r.point[0]), Interval(r.point[1]));
        point_touches.push_back({x, y, r.point[0], r.point[1], {h.x / h.w, h.y / h.w, h.z / h.w}, true, r.point[2], r.point[3]});
        ++im.stats.point_contacts;
        note(where + ": surfaces " + str(x) + " and " + str(y) + " touch at one point, at parameters (" + std::to_string(r.point[0]) + ", " +
             std::to_string(r.point[1]) + ") and (" + std::to_string(r.point[2]) + ", " + std::to_string(r.point[3]) +
             "), proven exactly, and nowhere else in the region left");
        return true;
    }

    /// A near miss: a curved patch and another patch (curved, or a plane)
    /// that come close without meeting, inside the region the search could
    /// not finish. solve::classify_contact proves a facing pair with a nonzero
    /// gap there, and solve::separated_in proves that the normal offset from
    /// the curved patch keeps that gap's sign over boxes around every box the
    /// search left: the pair has no intersection there.
    bool near_miss(int x, int y, const solve::SsiResult& res, const std::string& where) {
        const bool px = is_planar[static_cast<std::size_t>(x)], py = is_planar[static_cast<std::size_t>(y)];
        if (px && py) return false;
        const Box4& H = res.unverified_hull;
        for (const Interval& h : H)
            if (!(h.hi - h.lo < 1e-3)) return false; // a small region only
        // A is the curved one (the offset is taken along its normal).
        const int ia = px ? y : x, ib = px ? x : y;
        const Interval hu = H[px ? 2u : 0u], hv = H[px ? 3u : 1u], hs = H[px ? 0u : 2u], ht = H[px ? 1u : 3u];
        const PolyPatchD& A = im.carriers[static_cast<std::size_t>(ia)];
        const PolyPatchD& B = im.carriers[static_cast<std::size_t>(ib)];
        const solve::ContactResult r = solve::classify_contact(A, B, {num::mid(hu), num::mid(hv), num::mid(hs), num::mid(ht), 0.0});
        if (r.kind != solve::ContactKind::Separated) return false;
        auto widen = [](const Interval& h) {
            const double m = std::max(h.hi - h.lo, std::ldexp(1.0, -20));
            return Interval(h.lo - m, h.hi + m);
        };
        std::array<double, 5> c;
        for (int i = 0; i < 5; ++i) c[static_cast<std::size_t>(i)] = num::mid(r.enc[static_cast<std::size_t>(i)]);
        std::string why;
        const int side = solve::separated_in(A, B, c, widen(hu), widen(hv), widen(hs), widen(ht), why);
        if (side == 0) {
            note(where + ": a facing pair with a gap, but not proven apart over the region left: " + why);
            return false;
        }
        ++im.stats.near_misses;
        char buf[160];
        std::snprintf(buf, sizeof buf, ": a near miss, gap in [%.3e, %.3e] along the normal of surface %d, proven not to close over the region left",
                      r.enc[4].lo, r.enc[4].hi, ia);
        note(where + buf);
        return true;
    }

    // --- an edge in a curved surface -------------------------------------------------
    // A straight edge of one solid can lie in a curved surface of the other:
    // along a straight line of a cylinder, say. The edge/surface search then
    // finds no isolated root and cannot finish. Where the edge runs along a
    // line of constant parameter of the surface, this is decided exactly, and
    // the edge is treated as an edge lying in a plane is: it is not searched
    // against that surface, its vertices get parameters there, it enters the
    // planar maps of the surface's faces as a straight piece, the surface's
    // curves with the edge's own two surfaces are the edge itself along it,
    // and classification does not cross it.
    struct IsoLine {
        int g = -1;          ///< the surface the edge lies in
        bool fixed_u = true; ///< along u = a / b (else along v = a / b)
        double a = 0.0, b = 1.0;
        Vec4d A{}, B{};      ///< the line's points at the free parameter 0 and 1 (homogeneous, one weight)
    };
    std::map<std::pair<int, int>, std::vector<IsoLine>> iso_lines; ///< (solid, edge) -> the surfaces it lies in

    /// Is the homogeneous point X on the line through P0 and P1 (unit
    /// weights)? Exact: (X - w P0) x (P1 - P0) = 0.
    static bool on_line(const Vec4d& P0, const Vec4d& P1, const Vec4d& X) {
        using E = num::Expansion;
        const E w(X.w);
        const E r[3] = {E(X.x) - w * E(P0.x), E(X.y) - w * E(P0.y), E(X.z) - w * E(P0.z)};
        const E d[3] = {E(P1.x) - E(P0.x), E(P1.y) - E(P0.y), E(P1.z) - E(P0.z)};
        return (r[1] * d[2] - r[2] * d[1]).sign() == 0 && (r[2] * d[0] - r[0] * d[2]).sign() == 0 &&
               (r[0] * d[1] - r[1] * d[0]).sign() == 0;
    }
    /// Does the straight edge e of solid k lie in the curved surface g along
    /// a line of constant parameter? The candidate comes from the parameters
    /// of the edge's two end points on g (in doubles): one of them the same
    /// at both ends, close to a fraction a / b. Proof: g has degree 1 in the
    /// other parameter, so its curve at a / b is the segment between its
    /// exact points A and B there, linearly parameterised when their weights
    /// are equal; and A and B lie exactly on the edge's line.
    bool find_iso_line(int k, int e, int g, IsoLine& L) const {
        const BEdge& edge = S(k).edges[static_cast<std::size_t>(e)];
        if (edge.arc) return false;
        return find_iso_line(side_control(S(k).carriers[static_cast<std::size_t>(edge.on[0].carrier)], edge.on[0].side), g, L);
    }
    /// The same for a segment given by the control points of a side.
    bool find_iso_line(const std::vector<Vec4d>& cp, int g, IsoLine& L) const {
        if (is_planar[static_cast<std::size_t>(g)]) return false;
        const PolyPatchD& G = im.carriers[static_cast<std::size_t>(g)];
        if (!linear_side(cp)) return false;
        const Vec4d P0 = cp.front(), P1 = cp.back();
        const geom::BezierSurface Sg(G);
        const double scale = 1.0 + std::fabs(P0.x) + std::fabs(P0.y) + std::fabs(P0.z) + std::fabs(P1.x) + std::fabs(P1.y) + std::fabs(P1.z);
        double u0 = 0.0, v0 = 0.0, u1 = 0.0, v1 = 0.0;
        if (!surfaces_detail::locate_on_patch(Sg, {P0.x, P0.y, P0.z}, scale, u0, v0) ||
            !surfaces_detail::locate_on_patch(Sg, {P1.x, P1.y, P1.z}, scale, u1, v1))
            return false;
        for (int fixed = 0; fixed < 2; ++fixed) { // 0: u is constant along the edge; 1: v is
            const double f0 = fixed == 0 ? u0 : v0, f1 = fixed == 0 ? u1 : v1;
            if (std::fabs(f0 - f1) > 1e-8) continue;
            double a = 0.0, b = 1.0;
            if (!rational_near(0.5 * (f0 + f1), a, b)) continue;
            Vec4d A, B;
            if ((fixed == 0 ? G.q : G.p) == 1) {
                const bool ok = fixed == 0 ? exact_point_at(G, a, b, 0.0, 1.0, A) && exact_point_at(G, a, b, 1.0, 1.0, B)
                                           : exact_point_at(G, 0.0, 1.0, a, b, A) && exact_point_at(G, 1.0, 1.0, a, b, B);
                if (!ok || A.w != B.w || planar_detail::same_point_h(A, B)) continue;
            } else {
                // Of higher degree along the line: only a side of the square,
                // whose control points must be a straight segment, linearly
                // parameterised (the side of a box's top with a hill on it).
                if (b != 1.0 || (a != 0.0 && a != 1.0)) continue;
                const int side = fixed == 0 ? (a == 0.0 ? 3 : 1) : (a == 0.0 ? 0 : 2);
                const auto sc = side_control(G, side);
                if (!linear_side(sc)) continue;
                A = sc.front();
                B = sc.back();
            }
            if (!on_line(P0, P1, A) || !on_line(P0, P1, B)) continue;
            L.g = g;
            L.fixed_u = fixed == 0;
            L.a = a;
            L.b = b;
            L.A = A;
            L.B = B;
            return true;
        }
        return false;
    }
    /// The record of edge e of solid k lying in a surface of g's class, if any.
    const IsoLine* iso_line_on(int k, int e, int g) const {
        const auto it = iso_lines.find({k, e});
        if (it == iso_lines.end()) return nullptr;
        for (const IsoLine& L : it->second)
            if (equiv(L.g, g)) return &L;
        return nullptr;
    }
    /// Parameters on L.g of the point of edge e (solid k) at the edge
    /// parameter t: the line's constant parameter and, along it, the
    /// position between A and B. Exact where the position is dyadic.
    Chart iso_chart(int k, int e, const IsoLine& L, const Interval& t) const {
        const BEdge& edge = S(k).edges[static_cast<std::size_t>(e)];
        const auto cp = side_control(S(k).carriers[static_cast<std::size_t>(edge.on[0].carrier)], edge.on[0].side);
        const double p0[3] = {cp.front().x, cp.front().y, cp.front().z}, p1[3] = {cp.back().x, cp.back().y, cp.back().z};
        const double la[3] = {L.A.x, L.A.y, L.A.z}, lb[3] = {L.B.x, L.B.y, L.B.z};
        int c = 0;
        for (int i = 1; i < 3; ++i)
            if (std::fabs(lb[i] - la[i]) > std::fabs(lb[c] - la[c])) c = i;
        const Interval s = edge.on[0].reversed ? Interval(1.0) - t : t; // the side's own parameter
        const Interval x = Interval(p0[c]) + s * (Interval(p1[c]) - Interval(p0[c]));
        Interval free = (Interval(L.A.w) * x - Interval(la[c])) / (Interval(lb[c]) - Interval(la[c]));
        if (is_point(s)) {
            // A + f (B - A) = w (P0 + s (P1 - P0)) in every coordinate, exactly?
            using E = num::Expansion;
            const double f = snap_dyadic(num::mid(free));
            bool same = num::contains(free, f);
            for (int i = 0; i < 3 && same; ++i)
                same = (E(la[i]) + E(f) * (E(lb[i]) - E(la[i])) - E(L.A.w) * (E(p0[i]) + E(s.lo) * (E(p1[i]) - E(p0[i])))).sign() == 0;
            if (same) free = Interval(f);
        }
        const Interval fixed = L.b == 1.0 ? Interval(L.a) : Interval(L.a) / Interval(L.b);
        return L.fixed_u ? Chart{L.g, fixed, free} : Chart{L.g, free, fixed};
    }
    /// Give vertex v, on edge e of solid k at the edge parameter t, its
    /// parameters on the surface the edge lies in, if it has none there.
    void iso_vertex_chart(int k, int e, const IsoLine& L, int v, const Interval& t) {
        if (chart(v, L.g)) return;
        im.vertices[static_cast<std::size_t>(v)].charts.push_back(iso_chart(k, e, L, t));
        spread(v);
    }

    /// Does the solid of edge e (solid k), which lies in L.g, stay on one side
    /// of that surface along the edge (Yes: it touches it), or does the
    /// surface pass between the edge's faces (No)? The faces along the edge
    /// lie on its two surfaces, each inside its parameter square, so each
    /// leaves the edge across its side, into the square. Their signs against
    /// the surface's normal are taken at one point of the edge, with
    /// intervals. One point is enough: a sign can change along the edge only
    /// where L.g is tangent to one of the two surfaces, and their curves
    /// with L.g along the edge are certified regular.
    std::map<std::tuple<int, int, int>, Tri> touching_cache;
    Tri iso_touching(int k, int e, const IsoLine& L) {
        const auto key = std::make_tuple(k, e, L.g);
        const auto it = touching_cache.find(key);
        if (it != touching_cache.end()) return it->second;
        Tri result = Tri::Unknown;
        const BEdge& edge = S(k).edges[static_cast<std::size_t>(e)];
        int g0 = -1, g1 = -1;
        edge_carriers(k, e, g0, g1);
        bool on[2] = {false, false}, other = false;
        for (const BFace& f : S(k).faces)
            for (const auto& loop : f.loops)
                for (const BCoedge& ce : loop) {
                    if (ce.edge != e) continue;
                    const int h = gid(k, f.carrier);
                    if (h == g0) on[0] = true;
                    else if (h == g1) on[1] = true;
                    else other = true;
                }
        if (!other && (on[0] || on[1])) {
            const Interval tm(0.5 * (num::mid(edge.t0) + num::mid(edge.t1)));
            const Chart lc = iso_chart(k, e, L, tm);
            int sign = 0;
            bool decided = true, crossing = false;
            for (int i = 0; i < 2 && decided; ++i) {
                if (!on[i]) continue;
                const SsiSystem sys(im.carriers[static_cast<std::size_t>(L.g)], im.carriers[static_cast<std::size_t>(i == 0 ? g0 : g1)]);
                Interval su, sv;
                side_uv(edge.on[i].side, edge.on[i].reversed ? Interval(1.0) - tm : tm, su, sv);
                const Box4 X = {lc.u, lc.v, su, sv};
                Interval J[3][4];
                sys.J_over(X, J);
                // Columns 0 and 1: positive multiples of L.g's tangents; 2 and 3:
                // negative multiples of the other surface's.
                const Interval n[3] = {J[1][0] * J[2][1] - J[2][0] * J[1][1], J[2][0] * J[0][1] - J[0][0] * J[2][1],
                                       J[0][0] * J[1][1] - J[1][0] * J[0][1]};
                const int side = edge.on[i].side, col = side == 0 || side == 2 ? 3 : 2;
                const Interval dot = n[0] * J[0][col] + n[1] * J[1][col] + n[2] * J[2][col];
                if (num::contains_zero(dot)) { decided = false; break; }
                // Into the square: +v at side 0, -u at 1, -v at 2, +u at 3.
                const int inward = side == 0 || side == 3 ? -1 : 1; // against the column's negative factor
                const int si = inward * (dot.lo > 0.0 ? 1 : -1);
                if (sign == 0) sign = si;
                else if (si != sign) crossing = true;
            }
            if (decided) result = crossing ? Tri::No : Tri::Yes;
        }
        return touching_cache[key] = result;
    }

    /// An exact three-surface root on gs[0..2]. One surface with dyadic
    /// parameters anchors the exact point, which must lie exactly on the
    /// other two inside the uniqueness box. `out` are the parameters on the
    /// three surfaces (points where dyadic).
    bool exact_triple_root(const int gs[3], const TsiRoot& root, Chart out[3]) const {
        if (!opt.exact_points) return false;
        for (int a = 0; a < 3; ++a) {
            const std::size_t ia = static_cast<std::size_t>(2 * a);
            const double us = snap_dyadic(num::mid(root.enc[ia])), vs = snap_dyadic(num::mid(root.enc[ia + 1]));
            if (!num::contains(root.enc[ia], us) || !num::contains(root.unique[ia], us) || !num::contains(root.enc[ia + 1], vs) ||
                !num::contains(root.unique[ia + 1], vs))
                continue;
            Vec4d A;
            if (!exact_point_at(im.carriers[static_cast<std::size_t>(gs[a])], us, vs, A)) continue;
            bool ok = true;
            for (int j = 0; j < 3 && ok; ++j) {
                const std::size_t ij = static_cast<std::size_t>(2 * j);
                if (j == a) { out[j] = {gs[j], Interval(us), Interval(vs)}; continue; }
                ok = exact_on_surface(A, gs[j], root.enc[ij], root.enc[ij + 1], out[j]) && num::subset(out[j].u, root.unique[ij]) &&
                     num::subset(out[j].v, root.unique[ij + 1]);
            }
            if (ok) return true;
        }
        return false;
    }

    /// Is a new root (charts on g0, g1, g) an existing vertex, by a proof that
    /// a box containing both holds exactly one common point of the three
    /// surfaces? Used where the cheaper tests cannot decide: a crossing found
    /// from two edges of different carriers of one surface, whose parameters
    /// are intervals straddling a side.
    bool identify_point3(int g0, int g1, int g, const std::vector<Chart>& root) {
        // Once proven the same point, the root's charts enclose the vertex as
        // well as its own do: keep the intersection (a side coordinate that is
        // exact on the root becomes exact on the vertex).
        auto merge = [&](int v) {
            for (const Chart& rc : root) {
                Chart* mine = nullptr;
                for (Chart& c : im.vertices[static_cast<std::size_t>(v)].charts)
                    if (c.carrier == rc.carrier) mine = &c;
                if (!mine) { im.vertices[static_cast<std::size_t>(v)].charts.push_back(rc); continue; }
                const double ul = std::max(mine->u.lo, rc.u.lo), uh = std::min(mine->u.hi, rc.u.hi);
                const double vl = std::max(mine->v.lo, rc.v.lo), vh = std::min(mine->v.hi, rc.v.hi);
                if (ul <= uh && vl <= vh) { mine->u = Interval(ul, uh); mine->v = Interval(vl, vh); }
            }
            spread(v);
        };
        for (int v = 0; v < static_cast<int>(im.vertices.size()); ++v) {
            const Chart* a = chart(v, g0);
            const Chart* b = chart(v, g1);
            const Chart* c = chart(v, g);
            if (!a || !b || !c) continue;
            const Interval mine[6] = {a->u, a->v, b->u, b->v, c->u, c->v};
            const Interval theirs[6] = {root[0].u, root[0].v, root[1].u, root[1].v, root[2].u, root[2].v};
            Box6 Y;
            Pt6 y{};
            bool near = true;
            for (int i = 0; i < 6 && near; ++i) {
                const double lo = std::min(mine[i].lo, theirs[i].lo), hi = std::max(mine[i].hi, theirs[i].hi);
                near = hi - lo < 1e-6;
                const double pad = 1e-9 + 0.5 * (hi - lo);
                Y[static_cast<std::size_t>(i)] = Interval(lo - pad, hi + pad);
                y[static_cast<std::size_t>(i)] = 0.5 * (lo + hi);
            }
            if (!near) continue;
            const SsiSystem edge_sys(im.carriers[static_cast<std::size_t>(g0)], im.carriers[static_cast<std::size_t>(g1)]);
            const TsiSystem sys(edge_sys, im.carriers[static_cast<std::size_t>(g)]);
            Box6 K;
            if (tsi_newton(sys, y)) {
                for (int i = 0; i < 6; ++i)
                    y[static_cast<std::size_t>(i)] = std::min(std::max(y[static_cast<std::size_t>(i)], Y[static_cast<std::size_t>(i)].lo),
                                                              Y[static_cast<std::size_t>(i)].hi);
            }
            if (tsi_krawczyk(sys, y, Y, K)) {
                merge(v);
                return true;
            }
        }
        return false;
    }
    /// A root that is not an existing vertex. `charts` are on g0, g1 and g.
    /// Does the straight side `side0` of g0 (which g1 shares) lie exactly in
    /// the plane of a surface that carries faces of solid k, other than
    /// through g0's or g1's own plane? (A line along which a curved surface of
    /// one solid touches a face of the other.)
    bool side_in_face_plane(int k, int g0, int side0, int g1) const {
        for (std::size_t r = 0; r < members.size(); ++r) {
            if (members[r].empty() || !is_planar[static_cast<std::size_t>(members[r].front())]) continue;
            const int h = members[r].front();
            if (equiv(h, g0) || equiv(h, g1) || !faces_in_plane(k, h)) continue;
            if (side_in_plane(im.carriers[static_cast<std::size_t>(g0)], side0, im.carriers[static_cast<std::size_t>(h)])) return true;
        }
        return false;
    }

    /// `xp`: the root's exact point, where it was proven in space.
    bool add_root(int k, int g, int g0, int g1, std::vector<Chart> charts, const Vec3d& point, bool triple,
                  const std::string& where, int faces_of_edge = 2, int side0 = -1, const Vec4d* xp = nullptr) {
        // A root whose exact point is that of an existing vertex is that
        // vertex (found from another edge, or an operand's own vertex): it
        // gets the root's parameters where it has none, and nothing is added.
        Vec4d E{};
        bool exact = false;
        if (opt.exact_points) {
            exact = xp != nullptr;
            if (exact) E = *xp;
            for (std::size_t c = 0; c < charts.size() && !exact; ++c)
                exact = is_point(charts[c].u) && is_point(charts[c].v) &&
                        exact_point_at(im.carriers[static_cast<std::size_t>(charts[c].carrier)], charts[c].u.lo, charts[c].v.lo, E);
            if (exact)
                for (int w = 0; w < static_cast<int>(im.vertices.size()); ++w) {
                    Vec4d W;
                    if (!exact_point(w, W) || !planar_detail::same_point_h(W, E)) continue;
                    Imprint::Vertex& vx = im.vertices[static_cast<std::size_t>(w)];
                    for (const Chart& c : charts)
                        if (!chart(w, c.carrier)) vx.charts.push_back(c);
                    if (!vx.exact) { vx.exact = true; vx.xp = E; }
                    spread(w);
                    ++im.stats.merged_vertices;
                    note(where + ": the root is vertex " + str(w) + " (same exact point)");
                    return true;
                }
        }
        int face = -1;
        const Hit hit = find_face(1 - k, g, charts[2].u, charts[2].v, face);
        if (hit == Hit::None) return true;
        int kind = kPierce;
        if (hit == Hit::Undecided) {
            // On the boundary of a face: when the edge lies in a surface that
            // carries faces of the other solid, that is where it crosses an edge
            // of the other solid inside that surface. (The edge's own solid may
            // have no face there: the edge can be a seam.) The same where the
            // edge, a straight side, lies exactly in such a plane without being
            // on a surface of it (a curved surface touching that plane along it).
            // Or the edge lies in a plane that holds a side of the surface hit:
            // the face's boundary there runs in that plane too (a box's top
            // edge in the plane of a cap the edge lies on), and the two edges
            // cross inside it.
            auto plane_holds_side = [&]() {
                for (int gp : {g0, g1}) {
                    if (!is_planar[static_cast<std::size_t>(gp)]) continue;
                    for (int s = 0; s < 4; ++s)
                        if (side_in_plane(im.carriers[static_cast<std::size_t>(g)], s, im.carriers[static_cast<std::size_t>(gp)])) return true;
                }
                return false;
            };
            // Or, in space: the root's exact point lies on an edge of the other
            // solid that bounds a face on g, strictly inside that edge (two
            // edges crossing at a point with round coordinates). The vertex
            // then lies on both edges and on all four surfaces.
            auto through_edge = [&]() {
                if (!exact) return false;
                for (int e2 = 0; e2 < static_cast<int>(S(1 - k).edges.size()); ++e2) {
                    int h0 = -1, h1 = -1;
                    edge_carriers(1 - k, e2, h0, h1);
                    if (!equiv(h0, g) && !equiv(h1, g)) continue;
                    Chart on2[2];
                    if (!exact_on_edge(E, 1 - k, e2, on2)) continue;
                    for (const Chart& c : on2) {
                        bool has = false;
                        for (const Chart& d : charts) has = has || d.carrier == c.carrier;
                        if (!has) charts.push_back(c);
                    }
                    ++im.stats.edge_crossings;
                    note(where + ": the edge crosses edge " + str(e2) + " of solid " + str(1 - k) + " (exact point)");
                    return true;
                }
                return false;
            };
            if (!faces_in_plane(1 - k, g0) && !faces_in_plane(1 - k, g1) &&
                !(side0 >= 0 && side_in_face_plane(1 - k, g0, side0, g1)) && !plane_holds_side() && !through_edge()) {
                if (opt.trace) {
                    std::string s = where + ": root at (" + std::to_string(point.x) + ", " + std::to_string(point.y) + ", " +
                                    std::to_string(point.z) + "); edge surfaces " + str(g0) + (is_planar[static_cast<std::size_t>(g0)] ? "p" : "") +
                                    " class " + str(rep[static_cast<std::size_t>(g0)]) + ", " + str(g1) +
                                    (is_planar[static_cast<std::size_t>(g1)] ? "p" : "") + " class " + str(rep[static_cast<std::size_t>(g1)]) +
                                    "; surface " + str(g) + " class " + str(rep[static_cast<std::size_t>(g)]) + "; faces of each solid in those planes:";
                    for (int c : {g0, g1, g}) s += " [" + str(faces_in_plane(0, c) ? 1 : 0) + str(faces_in_plane(1, c) ? 1 : 0) + "]";
                    note(s);
                }
                return fail(where + ": the edge meets the boundary of a face (edge through edge; not in "
                                    "general position)");
            }
            kind = kCross;
        }
        for (const Chart& c : charts)
            im.stats.max_vertex_width =
                std::max(im.stats.max_vertex_width, std::max(num::width(c.u), num::width(c.v)));
        im.vertices.push_back({point, std::move(charts)});
        {
            // Its exact point: proven in space, or the point of a surface on
            // which its parameters were proven exact.
            const int nv = static_cast<int>(im.vertices.size()) - 1;
            Vec4d X = E;
            if (xp) X = *xp;
            if (xp || exact || exact_point(nv, X)) {
                im.vertices.back().exact = true;
                im.vertices.back().xp = X;
                ++im.stats.exact_roots;
            }
        }
        vkind.push_back(kind);
        arc_ends.push_back(0);
        if (kind == kPierce) expected_ends[static_cast<int>(im.vertices.size()) - 1] = faces_of_edge;
        spread(static_cast<int>(im.vertices.size()) - 1);
        touch_charts(static_cast<int>(im.vertices.size()) - 1);
        exact_charts(static_cast<int>(im.vertices.size()) - 1);
        if (!line_touches.empty() || graph_crossing_seen) share_side_charts(static_cast<int>(im.vertices.size()) - 1);
        ++im.stats.new_vertices;
        if (triple) ++im.stats.triple_vertices;
        return true;
    }

    /// Does the straight edge e of solid k (the curve C, its side) meet the
    /// plane g only at its end vertices? Decided exactly (op/contact.hpp): the
    /// curve's plane function vanishes at most at the ends of the side, and an
    /// end where it does is an end of the edge, with exact parameter. That
    /// vertex then lies in the plane: it gets parameters on g.
    bool touches_at_ends(int k, int e, int g, const PolyPatchD& C) {
        const BEdge& edge = S(k).edges[static_cast<std::size_t>(e)];
        const EndContact ec = end_contact(C, im.carriers[static_cast<std::size_t>(g)]);
        if (!ec.decided) return false;
        std::vector<int> on_plane;
        for (int end = 0; end < 2; ++end) {
            if ((end == 0 ? ec.order0 : ec.order1) == 0) continue;
            const double t = end;
            if (edge.t1.hi < t || edge.t0.lo > t) continue; // the touching end of the side is not on this edge
            int v = -1;
            if (edge.t0.lo == t && edge.t0.hi == t) v = vmap[k][static_cast<std::size_t>(edge.v0)];
            else if (edge.t1.lo == t && edge.t1.hi == t) v = vmap[k][static_cast<std::size_t>(edge.v1)];
            if (v < 0) return false; // the side touches the plane inside the edge
            on_plane.push_back(v);
            if (chart(v, g)) continue;
            const Vec4d& h = C.c[end == 0 ? 0 : static_cast<std::size_t>(C.p)];
            const Interval w(h.w);
            const Vec3<Interval> P{Interval(h.x) / w, Interval(h.y) / w, Interval(h.z) / w};
            Chart m;
            m.carrier = g;
            if (!plane_params(im.carriers[static_cast<std::size_t>(g)], P, m.u, m.v)) return false;
            im.vertices[static_cast<std::size_t>(v)].charts.push_back(m);
            spread(v);
        }
        ++im.stats.tangent_ends;
        note("edge " + str(e) + " of solid " + str(k) + " meets surface " + str(g) + " only at " +
             str(static_cast<int>(on_plane.size())) + " of its end vertices (exact)");
        return true;
    }

    bool intersect_edges() {
        for (int k = 0; k < 2; ++k) {
            edge_faces[k].assign(S(k).edges.size(), 0);
            for (const BFace& f : S(k).faces)
                for (const auto& loop : f.loops)
                    for (const BCoedge& ce : loop) ++edge_faces[k][static_cast<std::size_t>(ce.edge)];
        }
        for (int k = 0; k < 2; ++k) {
            const Brep& own = S(k);
            // Surfaces of the other solid with a face near an edge, from an
            // index of its face boxes, in increasing order.
            const geom::BoxIndex findex(face_box[1 - k]);
            for (int e = 0; e < static_cast<int>(own.edges.size()); ++e) {
                const BEdge& edge = own.edges[static_cast<std::size_t>(e)];
                int g0 = -1, g1 = -1;
                edge_carriers(k, e, g0, g1);
                std::vector<int> near_surfaces;
                for (int fo : findex.query(edge_box[k][static_cast<std::size_t>(e)]))
                    near_surfaces.push_back(gid(1 - k, S(1 - k).faces[static_cast<std::size_t>(fo)].carrier));
                std::sort(near_surfaces.begin(), near_surfaces.end());
                near_surfaces.erase(std::unique(near_surfaces.begin(), near_surfaces.end()), near_surfaces.end());                ArcEnclosureCache enclosures; // shared by this edge's searches against all surfaces

                const std::string name = "edge " + str(e) + " of solid " + str(k);
                PolyPatchD side_curve_patch;
                double tlo = 0.0, thi = 1.0;
                if (!edge.arc) {
                    side_curve_patch = solve::make_curve(
                        side_control(own.carriers[static_cast<std::size_t>(edge.on[0].carrier)], edge.on[0].side));
                    tlo = std::max(edge.t0.lo, 0.0);
                    thi = std::min(edge.t1.hi, 1.0);
                }
                for (int g : near_surfaces) {
                    // The edge lies in g: g is (in the plane of) one of its surfaces.
                    if (equiv(g, g0) || equiv(g, g1)) continue;
                    if (!edge.arc) {
                        // The edge is a side of g as well, or lies exactly in g's plane.
                        bool lies = side_lies_in(g0, edge.on[0].side, g1, g);
                        for (int s = 0; s < 4 && !lies; ++s) {
                            bool rev = false;
                            lies = shared_side(im.carriers[static_cast<std::size_t>(g)], s,
                                               im.carriers[static_cast<std::size_t>(g0)], edge.on[0].side, rev);
                        }
                        if (lies) continue;
                    }
                    // The edge lies in the plane of b, and g meets that plane only
                    // along its side s (proven exactly): every point the edge
                    // shares with g is on that side. Within g's faces that side is
                    // covered by edges of the other solid, which lie in the plane
                    // too and are searched against this edge's surfaces; those
                    // searches give the crossings. The edge only touches g there.
                    // Conversely, an arc on a surface a that meets the plane g
                    // only along its side s meets g only where it reaches that
                    // side: on the boundary of a's face, so at its own vertices,
                    // which have parameters on g already (touch_charts).
                    {
                        bool in_touch_plane = false;
                        for (const Touch& t : touches) {
                            if (t.a == g && (equiv(t.b, g0) || equiv(t.b, g1))) in_touch_plane = true;
                            // Or a side lying exactly in that plane without being on a
                            // surface of its class (a circle of a cap left on the cap's
                            // old rational carrier, say).
                            if (t.a == g && !edge.arc && is_planar[static_cast<std::size_t>(t.b)] &&
                                side_in_plane(im.carriers[static_cast<std::size_t>(g0)], edge.on[0].side,
                                              im.carriers[static_cast<std::size_t>(t.b)]))
                                in_touch_plane = true;
                            if (edge.arc && t.b == g && (t.a == g0 || t.a == g1)) in_touch_plane = true;
                        }
                        // A surface touching a plane along a line inside it: an
                        // edge in the plane meets the surface only on that line,
                        // and an edge on the surface meets the plane only there;
                        // neither crosses, so neither needs a vertex.
                        for (const LineTouch& t : line_touches) {
                            if (t.a == g && (equiv(t.b, g0) || equiv(t.b, g1))) in_touch_plane = true;
                            if (equiv(t.b, g) && (t.a == g0 || t.a == g1)) in_touch_plane = true;
                        }
                        if (in_touch_plane) continue;
                    }
                    // The edge lies in a plane that the curved surface g touches at
                    // one point p and meets nowhere else (proven exactly): every
                    // point the edge shares with g is p. It is a root exactly when
                    // p lies on the edge (decided exactly), and then a vertex where
                    // the edge touches g; otherwise there is none.
                    if (!edge.arc) {
                        const PointTouch* pt = nullptr;
                        for (const PointTouch& t : point_touches)
                            if (t.a == g && !t.b_curved && t.alone && (equiv(t.b, g0) || equiv(t.b, g1))) pt = &t;
                        if (pt) {
                            const std::string where = name + " and surface " + str(g);
                            if (!pt->has_xp) return fail(where + ": the edge lies in a plane touching the surface at a point that is not representable");
                            Chart ce[2];
                            if (exact_on_edge(pt->xp, k, e, ce)) {
                                std::vector<Chart> charts = {ce[0], ce[1], Chart{g, Interval(pt->u), Interval(pt->v)}};
                                const Vec4d& X = pt->xp;
                                if (!identify_point3(g0, g1, g, charts)) {
                                    const std::size_t before = im.vertices.size();
                                    if (!add_root(k, g, g0, g1, std::move(charts), {X.x / X.w, X.y / X.w, X.z / X.w}, false, where,
                                                  edge_faces[k][static_cast<std::size_t>(e)], edge.on[0].side, &X))
                                        return false;
                                    // Not a pierce: the edge touches g, so no new curve
                                    // has to end here (the pierce count does not apply).
                                    if (im.vertices.size() > before) {
                                        const int nv = static_cast<int>(im.vertices.size()) - 1;
                                        vkind[static_cast<std::size_t>(nv)] = kCross;
                                        expected_ends.erase(nv);
                                    }
                                    ++im.stats.edge_touches;
                                    note(where + ": the edge touches the surface at the point where its plane does (exact)");
                                }
                            } else {
                                for (int ev : {edge.v0, edge.v1}) {
                                    const auto& V = im.vertices[static_cast<std::size_t>(vmap[k][static_cast<std::size_t>(ev)])];
                                    using E = num::Expansion;
                                    auto same = [&](double a, double aw, double b, double bw) { return (E(a) * E(bw) - E(b) * E(aw)).sign() == 0; };
                                    if (V.exact && same(V.xp.x, V.xp.w, pt->xp.x, pt->xp.w) && same(V.xp.y, V.xp.w, pt->xp.y, pt->xp.w) &&
                                        same(V.xp.z, V.xp.w, pt->xp.z, pt->xp.w))
                                        return fail(where + ": the surface touches a plane exactly at a vertex of the edge (not in general position)");
                                }
                            }
                            continue;
                        }
                    }
                    ++im.stats.edge_carrier_pairs_solved;
                    const std::string where = name + " and surface " + str(g);
                    if (!edge.arc) {
                        const solve::CsiSystem cs(side_curve_patch, im.carriers[static_cast<std::size_t>(g)]);
                        const auto t0 = std::chrono::steady_clock::now();
                        const solve::CsiResult res = solve::intersect_curve_surface(cs, opt.csi, tlo, thi);
                        im.stats.edge_ms_straight += std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
                        ++im.stats.edge_searches_straight;
                        im.stats.csi_nodes += res.nodes;
                        if (!res.complete) {
                            // A tangential root cannot be isolated by intervals. For a
                            // plane, decide exactly whether the edge meets it only at
                            // its end vertices (then there is nothing new here).
                            if (is_planar[static_cast<std::size_t>(g)] && touches_at_ends(k, e, g, side_curve_patch)) continue;
                            // Or the edge lies in the surface, along a line of
                            // constant parameter (exact): nothing to find here.
                            // Its end vertices lie on g.
                            IsoLine line;
                            if (find_iso_line(k, e, g, line)) {
                                iso_lines[{k, e}].push_back(line);
                                iso_vertex_chart(k, e, line, vmap[k][static_cast<std::size_t>(edge.v0)], edge.t0);
                                iso_vertex_chart(k, e, line, vmap[k][static_cast<std::size_t>(edge.v1)], edge.t1);
                                ++im.stats.iso_edges;
                                note(where + ": the edge lies in the surface, along its line " + (line.fixed_u ? "u" : "v") + " = " +
                                     std::to_string(line.a / line.b) + " (exact)");
                                continue;
                            }
                            if (opt.trace) {
                                const auto a = cs.curve_point(Interval(tlo)), b = cs.curve_point(Interval(thi));
                                note(where + ": the search was not complete; side on surfaces " + sk(g0) + ", " + sk(g1) +
                                     (straight_piece_edge(k, e) ? ", straight" : ", curved") + ", from (" +
                                     std::to_string(num::mid(a.x)) + ", " + std::to_string(num::mid(a.y)) + ", " +
                                     std::to_string(num::mid(a.z)) + ") to (" + std::to_string(num::mid(b.x)) + ", " +
                                     std::to_string(num::mid(b.y)) + ", " + std::to_string(num::mid(b.z)) + "); surface " + sk(g) +
                                     "; " + str(static_cast<int>(res.roots.size())) + " roots found");
                            }
                            return fail(where + ": could not be resolved (the edge touches the surface or lies "
                                                "in it; not in general position)");
                        }
                        for (const solve::CsiRoot& root : res.roots) {
                            const Interval& t = root.enc[0];
                            if (t.hi < edge.t0.lo || t.lo > edge.t1.hi) continue;
                            if (identify_side_root(g0, edge.on[0].side, g1, g, root.unique)) continue;
                            if (!(edge.t0.hi < t.lo && t.hi < edge.t1.lo)) {
                                // A root that cannot be separated from an end of the
                                // edge. For a plane, decide exactly whether the side
                                // meets it only at its ends: then the root is that end
                                // vertex, which lies in the plane and gets parameters
                                // on it, and the edge has no other point there.
                                if (is_planar[static_cast<std::size_t>(g)] && touches_at_ends(k, e, g, side_curve_patch)) break;
                                // Or the end vertex has an exact point that lies
                                // exactly on g inside the root's uniqueness box.
                                if (exact_end_is_root(k, e, g, root)) continue;
                                if (opt.trace) {
                                    const auto p = cs.curve_point(t);
                                    std::string s = where + ": root at (" + std::to_string(num::mid(p.x)) + ", " +
                                                    std::to_string(num::mid(p.y)) + ", " + std::to_string(num::mid(p.z)) +
                                                    "), edge parameter [" + std::to_string(t.lo) + ", " + std::to_string(t.hi) +
                                                    "] of [" + std::to_string(edge.t0.lo) + ", " + std::to_string(edge.t1.hi) +
                                                    "]; edge surfaces " + str(g0) + ", " + str(g1) + "; end vertices";
                                    for (int ev : {edge.v0, edge.v1}) {
                                        const int v = vmap[k][static_cast<std::size_t>(ev)];
                                        s += " " + str(v) + " (on";
                                        for (const Chart& c : im.vertices[static_cast<std::size_t>(v)].charts) s += " " + str(c.carrier);
                                        s += ")";
                                    }
                                    note(s);
                                }
                                return fail(where + ": the surface passes through a vertex of the edge "
                                                    "(not in general position)");
                            }
                            std::vector<Chart> charts(3);
                            // A root proven exact: its dyadic parameters are points,
                            // not enclosures.
                            // Or a point with dyadic coordinates whose parameters
                            // are not dyadic: exact in space, enclosures here.
                            const auto p = cs.curve_point(t);
                            double ts = 0.0;
                            Chart exact_g;
                            Vec4d X{};
                            Interval s_exact;
                            const bool exact = exact_side_root(g0, edge.on[0], g, root, ts, exact_g);
                            const bool in_space = !exact && exact_space_side_root(g0, edge.on[0], g, root, p, X, s_exact, exact_g);
                            const Interval tc = in_space ? s_exact : t;
                            for (int i = 0; i < 2; ++i) {
                                charts[static_cast<std::size_t>(i)].carrier = i == 0 ? g0 : g1;
                                const Interval tau = exact ? Interval(edge.on[i].reversed ? 1.0 - ts : ts)
                                                           : (edge.on[i].reversed ? Interval(1.0) - tc : tc);
                                side_uv(edge.on[i].side, tau, charts[static_cast<std::size_t>(i)].u, charts[static_cast<std::size_t>(i)].v);
                            }
                            charts[2] = exact || in_space ? exact_g : Chart{g, root.enc[1], root.enc[2]};
                            if (identify_point3(g0, g1, g, charts)) continue;
                            if (!add_root(k, g, g0, g1, std::move(charts), {num::mid(p.x), num::mid(p.y), num::mid(p.z)},
                                          false, where, edge_faces[k][static_cast<std::size_t>(e)], edge.on[0].side, in_space ? &X : nullptr))
                                return false;
                        }
                    } else {
                        const CurveBase& base = *edge.geom->base;
                        const auto t0 = std::chrono::steady_clock::now();
                        const TsiResult res = intersect_arc_surface(*edge.geom, im.carriers[static_cast<std::size_t>(g)], opt.csi, enclosures);
                        im.stats.edge_ms_arc += std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
                        ++im.stats.edge_searches_arc;
                        im.stats.tsi_nodes += res.nodes;
                        if (!res.complete && opt.trace) {
                            const auto& V0 = own.vertices[static_cast<std::size_t>(edge.v0)].point;
                            const auto& V1 = own.vertices[static_cast<std::size_t>(edge.v1)].point;
                            char buf[200];
                            std::snprintf(buf, sizeof buf, ": the arc of surfaces %d and %d from (%.6g %.6g %.6g) to (%.6g %.6g %.6g)", g0, g1,
                                          V0.x, V0.y, V0.z, V1.x, V1.y, V1.z);
                            note(where + buf);
                        }
                        if (!res.complete)
                            return fail(where + ": the three-surface search could not be resolved (" +
                                        str(static_cast<int>(res.unverified_boxes)) + " boxes left after " +
                                        str(static_cast<int>(res.nodes)) + " nodes; the edge touches the surface, "
                                        "lies in it, or the search limit was reached)");
                        for (const TsiRoot& root : res.roots) {
                            const Box4 P = {root.enc[0], root.enc[1], root.enc[2], root.enc[3]};
                            bool on = false;
                            CurvePos pos;
                            if (!locate_on_curve(base, P, on, pos))
                                return fail(where + ": could not place a crossing on the edge");
                            if (!on) continue;
                            const Tri in = in_arc(*edge.geom, pos);
                            if (in == Tri::No) continue;
                            if (identify_root6(g0, g1, g, root.unique)) continue;
                            if (in == Tri::Unknown) {
                                // A root that cannot be separated from an end of the
                                // edge: it is that end vertex if the vertex has an
                                // exact point that lies exactly on g, inside the
                                // root's uniqueness box.
                                bool is_end = false;
                                const int gs[3] = {g0, g1, g};
                                for (int ev : {edge.v0, edge.v1}) {
                                    const int v = vmap[k][static_cast<std::size_t>(ev)];
                                    if (!exact_vertex_is_root(v, gs, 3, root.enc.data(), root.unique.data())) continue;
                                    note(name + " ends on surface " + str(g) + " at vertex " + str(v) + " (exact point)");
                                    ++im.stats.exact_ends;
                                    is_end = true;
                                    break;
                                }
                                if (is_end) continue;
                                if (opt.trace) {
                                    const auto p = base.sys->point_a(P);
                                    std::string s = where + ": root at (" + std::to_string(num::mid(p.x)) + ", " +
                                                    std::to_string(num::mid(p.y)) + ", " + std::to_string(num::mid(p.z)) +
                                                    "); arc on surfaces " + sk(g0) + ", " + sk(g1) + "; surface " + sk(g) + "; end vertices";
                                    for (int ev : {edge.v0, edge.v1}) {
                                        const int v = vmap[k][static_cast<std::size_t>(ev)];
                                        s += " " + str(v) + " (on";
                                        for (const Chart& c : im.vertices[static_cast<std::size_t>(v)].charts)
                                            s += " " + sk(c.carrier) + (is_point(c.u) && is_point(c.v) ? "!" : "");
                                        Vec4d E;
                                        s += exact_point(v, E) ? "; exact point)" : "; no exact point)";
                                    }
                                    note(s);
                                }
                                return fail(where + ": the surface passes through a vertex of the edge "
                                                    "(not in general position)");
                            }
                            std::vector<Chart> charts = {{g0, root.enc[0], root.enc[1]},
                                                         {g1, root.enc[2], root.enc[3]},
                                                         {g, root.enc[4], root.enc[5]}};
                            const auto p = base.sys->point_a(P);
                            const int gs3[3] = {g0, g1, g};
                            Chart ex3[3];
                            Vec4d X{};
                            bool in_space = false;
                            if (exact_triple_root(gs3, root, ex3) || (in_space = exact_space_root(gs3, root, p, X, ex3)))
                                for (std::size_t i = 0; i < 3; ++i) charts[i] = ex3[i];
                            if (identify_point3(g0, g1, g, charts)) continue;
                            if (!add_root(k, g, g0, g1, std::move(charts), {num::mid(p.x), num::mid(p.y), num::mid(p.z)},
                                          true, where, edge_faces[k][static_cast<std::size_t>(e)], -1, in_space ? &X : nullptr))
                                return false;
                        }
                    }
                }
            }
        }
        return true;
    }

    // --- 3. clip the curves to the faces ---------------------------------------------
    /// Does a certified point of the curve on surfaces (x, y) lie on an
    /// operand's arc edge on the same two surfaces?
    Tri on_operand_edge(int x, int y, const Box4& mid, const CurveBase& curve) {
        // A straight edge that lies in a curved surface of the other solid
        // (iso_lines): that surface's curve with either of the edge's own
        // surfaces runs along the edge. Proof that this curve is that line:
        // a point of the line, in both surfaces' parameters, lies on the
        // certified curve, which has one solution on each slice of its tube.
        for (const auto& entry : iso_lines) {
            const int k = entry.first.first, e = entry.first.second;
            const BEdge& edge = S(k).edges[static_cast<std::size_t>(e)];
            int g0 = -1, g1 = -1;
            edge_carriers(k, e, g0, g1);
            for (const IsoLine& L : entry.second)
                for (int i = 0; i < 2; ++i) {
                    const int gi = i == 0 ? g0 : g1;
                    const bool same = equiv(gi, x) && equiv(L.g, y), swapped = equiv(gi, y) && equiv(L.g, x);
                    if (!same && !swapped) continue;
                    // The curve point's parameters on the edge's own surface.
                    Interval a0, a1;
                    if (same) map(x, gi).apply(mid[0], mid[1], a0, a1);
                    else map(y, gi).apply(mid[2], mid[3], a0, a1);
                    Interval tau;
                    if (on_side(Chart{gi, a0, a1}, edge.on[i].side, tau) == Tri::No) continue;
                    const Interval te = edge.on[i].reversed ? Interval(1.0) - tau : tau;
                    if (te.hi < edge.t0.lo || te.lo > edge.t1.hi) continue;
                    // The line's point at that edge parameter, on both surfaces.
                    const Interval tm(std::min(1.0, std::max(0.0, num::mid(te))));
                    Interval su, sv, lu, lv;
                    side_uv(edge.on[i].side, edge.on[i].reversed ? Interval(1.0) - tm : tm, su, sv);
                    const Chart lc = iso_chart(k, e, L, tm);
                    Box4 Q;
                    if (same) {
                        map(gi, x).apply(Interval(su), Interval(sv), Q[0], Q[1]);
                        map(L.g, y).apply(lc.u, lc.v, Q[2], Q[3]);
                    } else {
                        map(L.g, x).apply(lc.u, lc.v, Q[0], Q[1]);
                        map(gi, y).apply(Interval(su), Interval(sv), Q[2], Q[3]);
                    }
                    bool on = false;
                    CurvePos pos;
                    if (!locate_on_curve(curve, Q, on, pos) || !on) return Tri::Unknown;
                    if (edge.t0.hi < te.lo && te.hi < edge.t1.lo) return Tri::Yes;
                    return Tri::Unknown;
                }
        }
        // A curved surface touching a plane b along a line L (line_touches),
        // where L also lies in the other plane of this curve: the curve runs
        // along L there, and L is the line of the operands' straight edges on
        // that plane and b. Proof that this part of the curve is L: a point of
        // L, in both surfaces' parameters, lies on the certified curve, which
        // has one solution per slice. Then the part is an existing edge if
        // that point lies within an edge on (that plane, b).
        for (const LineTouch& t : line_touches) {
            const int plane = y == t.a ? x : (x == t.a ? y : -1);
            if (plane < 0 || !is_planar[static_cast<std::size_t>(plane)] || equiv(plane, t.b)) continue;
            const PolyPatchD& A = im.carriers[static_cast<std::size_t>(t.a)];
            double ra = 0.0, rb = 1.0;
            if (!rational_near(num::mid(t.r), ra, rb)) continue;
            Vec4d L0, L1;
            if (!(t.fixed_u ? exact_point_at(A, ra, rb, 0.0, 1.0, L0) && exact_point_at(A, ra, rb, 1.0, 1.0, L1)
                            : exact_point_at(A, 0.0, 1.0, ra, rb, L0) && exact_point_at(A, 1.0, 1.0, ra, rb, L1)))
                continue;
            const PolyPatchD& Pl = im.carriers[static_cast<std::size_t>(plane)];
            if (!planar_detail::on_plane(Pl, L0) || !planar_detail::on_plane(Pl, L1)) continue; // L in this plane
            // The point of L at the curve point's position along L.
            const std::size_t ca = y == t.a ? 2 : 0, cp = y == t.a ? 0 : 2;
            const Interval along = t.fixed_u ? mid[ca + 1] : mid[ca];
            const Interval au = t.fixed_u ? t.r : along, av = t.fixed_u ? along : t.r;
            const auto h = icarrier[static_cast<std::size_t>(t.a)].eval(au, av);
            const Vec3<Interval> P3{h.x / h.w, h.y / h.w, h.z / h.w};
            Box4 Q;
            Q[ca] = au;
            Q[ca + 1] = av;
            if (!plane_params(Pl, P3, Q[cp], Q[cp + 1])) continue;
            bool on = false;
            CurvePos pos;
            if (!locate_on_curve(curve, Q, on, pos) || !on) continue;
            for (int k = 0; k < 2; ++k)
                for (int e = 0; e < static_cast<int>(S(k).edges.size()); ++e) {
                    const BEdge& edge = S(k).edges[static_cast<std::size_t>(e)];
                    int g0 = -1, g1 = -1;
                    edge_carriers(k, e, g0, g1);
                    const bool here = (equiv(g0, plane) && equiv(g1, t.b)) || (equiv(g0, t.b) && equiv(g1, plane));
                    if (!here) continue;
                    if (edge.arc) {
                        // The curve of two planes, made by an earlier operation:
                        // the point is on it (both planes hold L); where along it?
                        Box4 E4;
                        if (!plane_params(im.carriers[static_cast<std::size_t>(g0)], P3, E4[0], E4[1]) ||
                            !plane_params(im.carriers[static_cast<std::size_t>(g1)], P3, E4[2], E4[3]))
                            continue;
                        bool on_e = false;
                        CurvePos epos;
                        if (!locate_on_curve(*edge.geom->base, E4, on_e, epos) || !on_e) continue;
                        const Tri in = in_arc(*edge.geom, epos);
                        if (in != Tri::No) return in;
                        continue;
                    }
                    // Its position along the edge, from the point's parameters on
                    // the edge's first carrier (a plane, holding the line).
                    Interval eu, ev;
                    if (!plane_params(im.carriers[static_cast<std::size_t>(g0)], P3, eu, ev)) continue;
                    const int side = edge.on[0].side;
                    Interval tau = side == 0 || side == 2 ? eu : ev;
                    if (edge.on[0].reversed) tau = Interval(1.0) - tau;
                    if (tau.hi < edge.t0.lo || tau.lo > edge.t1.hi) continue;
                    if (edge.t0.hi < tau.lo && tau.hi < edge.t1.lo) return Tri::Yes;
                    return Tri::Unknown;
                }
        }
        for (int k = 0; k < 2; ++k)
            for (int e = 0; e < static_cast<int>(S(k).edges.size()); ++e) {
                const BEdge& edge = S(k).edges[static_cast<std::size_t>(e)];
                int g0 = -1, g1 = -1;
                edge_carriers(k, e, g0, g1);
                // The edge's surfaces are x and y, or in their planes.
                const bool same = equiv(g0, x) && equiv(g1, y), swapped = equiv(g0, y) && equiv(g1, x);
                if (!same && !swapped) continue;
                Box4 P;
                const int cx = same ? g0 : g1, cy = same ? g1 : g0; // the edge's surfaces matching x and y
                Interval a0, a1, b0, b1;
                map(x, cx).apply(mid[0], mid[1], a0, a1);
                map(y, cy).apply(mid[2], mid[3], b0, b1);
                P = same ? Box4{a0, a1, b0, b1} : Box4{b0, b1, a0, a1};
                if (!edge.arc) {
                    // Two distinct planes meet in one line, so the curve runs along
                    // the edge's line; is the point within the edge?
                    Interval tau;
                    const Chart c{g0, P[0], P[1]};
                    if (on_side(c, edge.on[0].side, tau) == Tri::No) continue;
                    if (tau.hi < edge.t0.lo || tau.lo > edge.t1.hi) continue;
                    if (edge.t0.hi < tau.lo && tau.hi < edge.t1.lo) return Tri::Yes;
                    return Tri::Unknown;
                }
                bool on = false;
                CurvePos pos;
                if (!locate_on_curve(*edge.geom->base, P, on, pos)) return Tri::Unknown;
                if (!on) continue;
                const Tri in = in_arc(*edge.geom, pos);
                if (in != Tri::No) return in;
            }
        return Tri::No;
    }

    bool clip_curves() {
        if (opt.trace) {
            std::string s = "shared surfaces:";
            for (int g = 0; g < static_cast<int>(im.carriers.size()); ++g)
                if (shared(g)) s += " " + str(g);
            note(s);
            for (int v = 0; v < static_cast<int>(im.vertices.size()); ++v) {
                if (vkind[static_cast<std::size_t>(v)] == kOld) continue;
                std::string line = "vertex " + str(v) + (vkind[static_cast<std::size_t>(v)] == kPierce ? " (pierce)" : " (cross)") + " on surfaces";
                for (const Chart& c : im.vertices[static_cast<std::size_t>(v)].charts) line += " " + str(c.carrier);
                note(line);
            }
        }
        for (const Curve& cu : curves) {
            const int x = cu.x, y = cu.y;
            const std::string name = "the curve between surfaces " + str(x) + " and " + str(y);
            struct OnCurve {
                int vertex = -1;
                CurvePos pos;
            };
            std::vector<OnCurve> hits;
            for (int v = 0; v < static_cast<int>(im.vertices.size()); ++v) {
                const Chart* a = chart(v, x);
                const Chart* b = chart(v, y);
                if (!a || !b) continue;
                OnCurve h;
                h.vertex = v;
                bool on = false;
                if (!locate_on_curve(*cu.base, {a->u, a->v, b->u, b->v}, on, h.pos))
                    return fail(name + ": could not decide whether a vertex lies on it");
                if (on) hits.push_back(h);
            }
            for (std::size_t i = 1; i < hits.size(); ++i)
                for (std::size_t j = i; j > 0; --j) {
                    const Tri t = before(hits[j].pos, hits[j - 1].pos);
                    if (t == Tri::Unknown) {
                        if (opt.trace) {
                            std::string s = name + ": vertices " + str(hits[j - 1].vertex) + " and " + str(hits[j].vertex) + " not ordered;";
                            for (int v : {hits[j - 1].vertex, hits[j].vertex}) {
                                const Imprint::Vertex& vx = im.vertices[static_cast<std::size_t>(v)];
                                Vec4d E;
                                s += " " + str(v) + " at (" + std::to_string(vx.point.x) + ", " + std::to_string(vx.point.y) + ", " +
                                     std::to_string(vx.point.z) + ") kind " + str(vkind[static_cast<std::size_t>(v)]) + " on";
                                for (const Chart& c : vx.charts)
                                    s += " " + sk(c.carrier) + "[" + std::to_string(num::width(c.u)) + "," + std::to_string(num::width(c.v)) + "]";
                                s += exact_point(v, E) ? " exact;" : " inexact;";
                            }
                            note(s);
                        }
                        return fail(name + ": two vertices on it could not be ordered");
                    }
                    if (t == Tri::No) break;
                    std::swap(hits[j], hits[j - 1]);
                }

            const int m = static_cast<int>(hits.size());
            if (opt.trace) {
                std::string s = "  " + name + (cu.base->closed ? " (closed)" : "") + ": vertices on it:";
                for (const OnCurve& h : hits) s += " " + str(h.vertex);
                note(s);
            }
            struct Part {
                int i = -1, j = -1; ///< indices into hits, or -1 for the whole closed curve
                bool wrap = false;
            };
            std::vector<Part> parts;
            if (cu.base->closed && m == 0) parts.push_back({-1, -1, false});
            if (cu.base->closed && m == 1) parts.push_back({0, 0, true}); // the whole loop, from the vertex round to it
            for (int i = 0; i + 1 < m; ++i) parts.push_back({i, i + 1, false});
            if (cu.base->closed && m >= 2) parts.push_back({m - 1, 0, true});

            for (const Part& part : parts) {
                const std::shared_ptr<const ArcGeom> arc =
                    part.i < 0 ? make_full_arc(cu.base) : make_arc(cu.base, hits[part.i].pos, hits[part.j].pos, part.wrap);
                Box4 mid;
                if (!arc_point(*arc, mid)) return fail(name + ": no certified point could be placed on a part of it");
                const Tri on_edge = on_operand_edge(x, y, mid, *cu.base);
                if (on_edge == Tri::Unknown) return fail(name + ": could not decide whether a part of it is an existing edge");
                const std::string part_name = "  " + name + ", part " +
                                              (part.i < 0 ? std::string("(whole loop)")
                                                          : str(hits[part.i].vertex) + " -> " + str(hits[part.j].vertex));
                if (on_edge == Tri::Yes) {
                    note(part_name + ": an existing edge");
                    continue; // an edge of an operand already
                }

                Imprint::NewArc na;
                int found = 0;
                for (int sx = 0; sx < 2; ++sx) { // the solid that owns the face on x
                    if (faces_on[sx][x].empty() || faces_on[1 - sx][y].empty()) continue;
                    int fx = -1, fy = -1;
                    const Hit hx = find_face(sx, x, mid[0], mid[1], fx);
                    const Hit hy = find_face(1 - sx, y, mid[2], mid[3], fy);
                    if (hx == Hit::None || hy == Hit::None) continue;
                    if (hx == Hit::Undecided || hy == Hit::Undecided) {
                        if (opt.trace) {
                            const auto p = cu.base->sys->point_a(mid);
                            note(part_name + ": midpoint (" + std::to_string(num::mid(p.x)) + ", " + std::to_string(num::mid(p.y)) + ", " +
                                 std::to_string(num::mid(p.z)) + ") on " + sk(x) + " [" + std::to_string(num::mid(mid[0])) + ", " +
                                 std::to_string(num::mid(mid[1])) + "] " + (hx == Hit::Undecided ? "undecided" : "decided") + ", on " + sk(y) +
                                 " [" + std::to_string(num::mid(mid[2])) + ", " + std::to_string(num::mid(mid[3])) + "] " +
                                 (hy == Hit::Undecided ? "undecided" : "decided"));
                        }
                        return fail(name + ": could not decide whether a part of it lies inside a face");
                    }
                    ++found;
                    na.solid_on[0] = sx;
                    na.solid_on[1] = 1 - sx;
                    na.face[0] = fx;
                    na.face[1] = fy;
                    // "Behind a surface" is inside its solid unless its face is flipped.
                    na.left_inside[0] = cu.left_behind[0] != S(1 - sx).faces[static_cast<std::size_t>(fy)].flipped;
                    na.left_inside[1] = cu.left_behind[1] != S(sx).faces[static_cast<std::size_t>(fx)].flipped;
                }
                note(part_name + (found ? ": new edge" : ": not inside a face of each solid"));
                if (found == 0) continue;
                if (found > 1) return fail(name + ": a part of it lies in faces of both solids on both surfaces");
                na.carrier[0] = x;
                na.carrier[1] = y;
                na.geom = arc;
                if (part.i < 0) {
                    // A closed curve that meets no vertex: mark it with a seam vertex.
                    const TubeSegment& seg = cu.base->curve->segs.front();
                    Pt4 xs = seg_start_point(seg);
                    Box4 enc;
                    if (!solve::certify_point(*cu.base->sys, seg.m, seg.a, seg_start_slice(seg), xs, enc))
                        return fail(name + ": its seam point could not be certified");
                    const auto p = cu.base->sys->point_a(enc);
                    na.v0 = na.v1 = static_cast<int>(im.vertices.size());
                    im.vertices.push_back({{num::mid(p.x), num::mid(p.y), num::mid(p.z)}, {{x, enc[0], enc[1]}, {y, enc[2], enc[3]}}});
                    vkind.push_back(kSeam);
                    arc_ends.push_back(0);
                    spread(na.v0);
                } else {
                    na.v0 = hits[part.i].vertex;
                    na.v1 = hits[part.j].vertex;
                    ++arc_ends[static_cast<std::size_t>(na.v0)];
                    ++arc_ends[static_cast<std::size_t>(na.v1)];
                }
                if (!arc->refinable) ++im.stats.coarse_arcs;
                im.arcs.push_back(std::move(na));
            }
        }
        im.stats.new_arcs = static_cast<int>(im.arcs.size());
        return true;
    }

    // --- 4. split the operands' edges at every vertex on them ----------------------
    /// Is the piece a straight segment? A side of degree 1, a side of higher
    /// degree whose control points lie exactly on one line (a straight edge
    /// of a degree-elevated flat face, such as the base of a box with a hill
    /// on top), or the intersection of two planes.
    bool straight_piece(const Imprint::OldPiece& p) const {
        const BEdge& e = S(p.solid).edges[static_cast<std::size_t>(p.edge)];
        int g0 = -1, g1 = -1;
        edge_carriers(p.solid, p.edge, g0, g1);
        if (!e.arc) {
            const PolyPatchD& c = im.carriers[static_cast<std::size_t>(g0)];
            if ((e.on[0].side == 0 || e.on[0].side == 2 ? c.p : c.q) == 1) return true;
            const auto cp = side_control(c, e.on[0].side);
            for (const Vec4d& x : cp)
                if (x.w != cp.front().w || !(x.w > 0.0) || !collinear(cp.front(), cp.back(), x)) return false;
            return !detail::same_point(cp.front(), cp.back());
        }
        if (is_planar[static_cast<std::size_t>(g0)] && is_planar[static_cast<std::size_t>(g1)]) return !equiv(g0, g1);
        // Two flat surfaces in different planes, one or both on a rational
        // carrier (a cap that was not moved): their curve is a line as well.
        return is_flat(g0) && is_flat(g1) &&
               !same_plane_patches(im.carriers[static_cast<std::size_t>(g0)], im.carriers[static_cast<std::size_t>(g1)]);
    }
    /// Is the side edge e of solid k a straight segment (see straight_piece)?
    bool straight_piece_edge(int k, int e) const {
        Imprint::OldPiece p;
        p.solid = k;
        p.edge = e;
        return straight_piece(p);
    }
    /// Do all control points of surface g lie in one plane (exact)? True for
    /// parallelograms and for flat faces on rational carriers.
    mutable std::vector<signed char> flat_known; ///< per surface: -1 not tested yet
    bool is_flat(int g) const {
        if (is_planar[static_cast<std::size_t>(g)]) return true;
        if (flat_known.size() != im.carriers.size()) flat_known.assign(im.carriers.size(), -1);
        signed char& f = flat_known[static_cast<std::size_t>(g)];
        if (f < 0) f = same_plane_patches(im.carriers[static_cast<std::size_t>(g)], im.carriers[static_cast<std::size_t>(g)]) ? 1 : 0;
        return f == 1;
    }

    bool split_edges() {
        // Where an edge pierces a face away from any shared surface, one new
        // curve ends per face of the edge: two, or more where parts of the
        // solid touch along it. (Around shared surfaces the counts vary; the
        // audit of the result is the check there.)
        for (std::size_t v = 0; v < im.vertices.size(); ++v) {
            if (vkind[v] != kPierce) continue;
            bool plain = true;
            for (const Chart& c : im.vertices[v].charts) plain = plain && !shared(c.carrier);
            const auto want = expected_ends.find(static_cast<int>(v));
            const int expect = want == expected_ends.end() ? 2 : want->second;
            if (plain && arc_ends[v] != expect)
                return fail("a vertex where an edge pierces a face has " + str(arc_ends[v]) +
                            " new intersection curves ending at it; expected " + str(expect));
        }

        // A vertex on an edge has its box meet the edge's box.
        std::vector<Vec3<Interval>> vboxes;
        for (int v = 0; v < static_cast<int>(im.vertices.size()); ++v) vboxes.push_back(vertex_box(v));
        const geom::BoxIndex vindex(std::move(vboxes));
        for (int k = 0; k < 2; ++k) {
            const Brep& own = S(k);
            const int nE = static_cast<int>(own.edges.size());
            edge_pieces[k].assign(static_cast<std::size_t>(nE), {});
            for (int e = 0; e < nE; ++e) {
                const BEdge& edge = own.edges[static_cast<std::size_t>(e)];
                const std::string name = "edge " + str(e) + " of solid " + str(k);
                int g0 = -1, g1 = -1;
                edge_carriers(k, e, g0, g1);
                const int ev0 = vmap[k][static_cast<std::size_t>(edge.v0)], ev1 = vmap[k][static_cast<std::size_t>(edge.v1)];
                std::vector<Split> list;
                for (int v : vindex.query(edge_box[k][static_cast<std::size_t>(e)])) {
                    if (v == ev0 || v == ev1) continue;
                    const Chart* a = chart(v, g0);
                    const Chart* b = chart(v, g1);
                    if (!a || !b) continue;
                    Split sp;
                    sp.vertex = v;
                    if (!edge.arc) {
                        Tri st = vertex_on_side(v, g0, edge.on[0].side, g1, sp.t);
                        if (st == Tri::Unknown) {
                            // A vertex with an exact point: its parameters on g0,
                            // proven exact where they can be, decide it.
                            Vec4d E;
                            Chart c0;
                            if (exact_point(v, E) && exact_chart(E, g0, c0)) {
                                Interval tau;
                                const Tri sx = on_side(c0, edge.on[0].side, tau);
                                if (sx != Tri::Unknown) { st = sx; sp.t = tau; }
                            }
                        }
                        if (st == Tri::No) continue;
                        // Beyond the ends of the edge it does not matter whether
                        // the vertex is on the side's line (a collinear edge
                        // further along it, such as the next tile's).
                        if (sp.t.hi < edge.t0.lo || sp.t.lo > edge.t1.hi) continue;
                        if (st == Tri::Unknown) return fail(name + ": could not decide whether a vertex lies on it");
                        if (!(edge.t0.hi < sp.t.lo && sp.t.hi < edge.t1.lo))
                            return fail(name + ": a vertex on it could not be separated from its ends");
                    } else {
                        bool on = false;
                        if (!locate_on_curve(*edge.geom->base, {a->u, a->v, b->u, b->v}, on, sp.pos))
                            return fail(name + ": could not decide whether a vertex lies on it");
                        if (!on) continue;
                        const Tri in = in_arc(*edge.geom, sp.pos);
                        if (in == Tri::No) continue;
                        if (in == Tri::Unknown) return fail(name + ": a vertex on it could not be separated from its ends");
                    }
                    list.push_back(sp);
                }
                // Certified order along the edge, from v0 to v1.
                auto precedes = [&](const Split& x, const Split& y) -> Tri {
                    if (!edge.arc) {
                        if (x.t.hi < y.t.lo) return Tri::Yes;
                        if (y.t.hi < x.t.lo) return Tri::No;
                        return Tri::Unknown;
                    }
                    if (edge.geom->wrap) {
                        const Tri wx = before(x.pos, edge.geom->a), wy = before(y.pos, edge.geom->a);
                        if (wx == Tri::Unknown || wy == Tri::Unknown) return Tri::Unknown;
                        if (wx != wy) return wx == Tri::No ? Tri::Yes : Tri::No;
                    }
                    return before(x.pos, y.pos);
                };
                for (std::size_t i = 1; i < list.size(); ++i)
                    for (std::size_t j = i; j > 0; --j) {
                        const Tri t = precedes(list[j], list[j - 1]);
                        if (t == Tri::Unknown) return fail(name + ": two vertices could not be ordered along it");
                        if (t == Tri::No) break;
                        std::swap(list[j], list[j - 1]);
                    }
                int prev = ev0;
                Interval tprev = edge.t0;
                CurvePos pprev = edge.arc ? edge.geom->a : CurvePos{};
                for (std::size_t i = 0; i <= list.size(); ++i) {
                    const bool last = i == list.size();
                    Imprint::OldPiece piece;
                    piece.solid = k;
                    piece.edge = e;
                    piece.v0 = prev;
                    piece.v1 = last ? ev1 : list[i].vertex;
                    if (!edge.arc) {
                        piece.t0 = tprev;
                        piece.t1 = last ? edge.t1 : list[i].t;
                        tprev = piece.t1;
                    } else if (list.empty()) {
                        piece.geom = edge.geom;
                    } else {
                        const CurvePos pnext = last ? edge.geom->b : list[i].pos;
                        piece.geom = make_arc(edge.geom->base, pprev, pnext, before(pprev, pnext) != Tri::Yes);
                        pprev = pnext;
                    }
                    prev = piece.v1;
                    piece.canon = static_cast<int>(im.pieces.size());
                    edge_pieces[k][static_cast<std::size_t>(e)].push_back(piece.canon);
                    im.pieces.push_back(std::move(piece));
                }
            }
        }

        // Pieces of the two solids on the same curve between the same vertices
        // are one piece.
        auto side_key = [&](int k, const BEdge& e) {
            const std::array<int, 2> a = {gid(k, e.on[0].carrier), e.on[0].side};
            const std::array<int, 2> b = {gid(k, e.on[1].carrier), e.on[1].side};
            return std::min(a, b);
        };
        std::map<std::array<int, 6>, std::vector<int>> by_key;
        auto key_of = [&](const Imprint::OldPiece& p) {
            const BEdge& e = S(p.solid).edges[static_cast<std::size_t>(p.edge)];
            std::array<int, 6> key{};
            if (straight_piece(p)) {
                // A straight segment is fixed by its two end points, whatever
                // surfaces it was made on.
                key = {2, 0, 0, 0, std::min(p.v0, p.v1), std::max(p.v0, p.v1)};
            } else if (!e.arc) {
                // By surface class, so that sides of different carriers of one
                // surface meet; whether they are the same curve is checked below.
                const int r0 = rep[static_cast<std::size_t>(gid(p.solid, e.on[0].carrier))];
                const int r1 = rep[static_cast<std::size_t>(gid(p.solid, e.on[1].carrier))];
                key = {0, std::min(r0, r1), std::max(r0, r1), 0, std::min(p.v0, p.v1), std::max(p.v0, p.v1)};
            } else {
                int g0 = -1, g1 = -1;
                edge_carriers(p.solid, p.edge, g0, g1);
                key = {1, std::min(g0, g1), std::max(g0, g1), 0, std::min(p.v0, p.v1), std::max(p.v0, p.v1)};
            }
            return key;
        };
        // Arcs on different but equivalent surfaces (one plane, one cylinder
        // under an exact reparameterisation) by the classes of their surfaces,
        // used only where the surfaces themselves give no match.
        auto class_key = [&](const Imprint::OldPiece& p) {
            int g0 = -1, g1 = -1;
            edge_carriers(p.solid, p.edge, g0, g1);
            const int r0 = rep[static_cast<std::size_t>(g0)], r1 = rep[static_cast<std::size_t>(g1)];
            return std::array<int, 6>{3, std::min(r0, r1), std::max(r0, r1), 0, std::min(p.v0, p.v1), std::max(p.v0, p.v1)};
        };
        std::map<std::array<int, 6>, std::vector<int>> by_class;
        std::map<std::pair<int, int>, std::vector<int>> mixed; // curved pieces of solid 0 by their end vertices
        for (std::size_t i = 0; i < im.pieces.size(); ++i)
            if (im.pieces[i].solid == 0) {
                by_key[key_of(im.pieces[i])].push_back(static_cast<int>(i));
                const BEdge& e = S(0).edges[static_cast<std::size_t>(im.pieces[i].edge)];
                if (e.arc && !straight_piece(im.pieces[i])) by_class[class_key(im.pieces[i])].push_back(static_cast<int>(i));
                if ((!e.arc || !straight_piece(im.pieces[i])) && im.pieces[i].v0 != im.pieces[i].v1)
                    mixed[{std::min(im.pieces[i].v0, im.pieces[i].v1), std::max(im.pieces[i].v0, im.pieces[i].v1)}].push_back(
                        static_cast<int>(i));
            }
        for (std::size_t i = 0; i < im.pieces.size(); ++i) {
            Imprint::OldPiece& p = im.pieces[i];
            if (p.solid != 1) continue;
            const BEdge& e = S(1).edges[static_cast<std::size_t>(p.edge)];
            // The candidates: by the surfaces themselves, else (arcs only) by
            // their classes. (Each map's iterator is compared only with its own end.)
            const std::vector<int>* candidates = nullptr;
            bool by_classes = false;
            if (const auto it = by_key.find(key_of(p)); it != by_key.end()) {
                candidates = &it->second;
            } else if (e.arc && !straight_piece(p)) {
                if (const auto jt = by_class.find(class_key(p)); jt != by_class.end()) {
                    candidates = &jt->second;
                    by_classes = true;
                }
            }
            static const std::vector<int> none;
            for (int q : candidates ? *candidates : none) {
                Imprint::OldPiece& o = im.pieces[static_cast<std::size_t>(q)];
                if (by_classes) {
                    // p's midpoint, carried into o's surfaces through the exact
                    // maps between equivalent surfaces, on o's arc?
                    Box4 mid;
                    if (!arc_point(*p.geom, mid)) return fail("no certified point could be placed on an edge piece");
                    int g0 = -1, g1 = -1, h0 = -1, h1 = -1;
                    edge_carriers(1, p.edge, g0, g1);
                    edge_carriers(0, o.edge, h0, h1);
                    const bool swap = !(equiv(g0, h0) && equiv(g1, h1));
                    if (swap && !(equiv(g0, h1) && equiv(g1, h0))) continue;
                    Box4 P;
                    map(g0, swap ? h1 : h0).apply(mid[0], mid[1], P[swap ? 2 : 0], P[swap ? 3 : 1]);
                    map(g1, swap ? h0 : h1).apply(mid[2], mid[3], P[swap ? 0 : 2], P[swap ? 1 : 3]);
                    bool on = false;
                    CurvePos pos;
                    if (!locate_on_curve(*o.geom->base, P, on, pos)) return fail("could not compare two edge pieces");
                    if (!on) continue;
                    const Tri in = in_arc(*o.geom, pos);
                    if (in == Tri::Unknown) return fail("could not compare two edge pieces");
                    if (in == Tri::No) continue;
                    if (p.v0 == p.v1) return fail("could not decide the direction of two coinciding closed edges");
                } else if (e.arc && !straight_piece(p)) {
                    // Two arcs between the same vertices may be different branches.
                    Box4 mid;
                    if (!arc_point(*p.geom, mid)) return fail("no certified point could be placed on an edge piece");
                    int g0 = -1, g1 = -1, h0 = -1, h1 = -1;
                    edge_carriers(1, p.edge, g0, g1);
                    edge_carriers(0, o.edge, h0, h1);
                    const Box4 P = g0 == h0 ? mid : Box4{mid[2], mid[3], mid[0], mid[1]};
                    bool on = false;
                    CurvePos pos;
                    if (!locate_on_curve(*o.geom->base, P, on, pos)) return fail("could not compare two edge pieces");
                    if (!on) continue;
                    const Tri in = in_arc(*o.geom, pos);
                    if (in == Tri::Unknown) return fail("could not compare two edge pieces");
                    if (in == Tri::No) continue;
                } else if (!e.arc && !straight_piece(p)) {
                    const BEdge& oe = S(0).edges[static_cast<std::size_t>(o.edge)];
                    if (side_key(1, e) != side_key(0, oe)) {
                        // Sides of different carriers of the same surfaces: the
                        // same curve if p's midpoint, carried into o's surface,
                        // lies exactly on o's side, inside o's range.
                        const Tri same = same_side_curve(p, e, o, oe);
                        if (same == Tri::Unknown) return fail("could not compare two edge pieces");
                        if (same == Tri::No) continue;
                    }
                }
                if (p.v0 == p.v1) {
                    // Two closed pieces of one curve through one vertex: their
                    // end vertices cannot say whether they run the same way.
                    const Tri rev = closed_relation(p, o);
                    if (rev == Tri::Unknown) return fail("could not decide the direction of two coinciding closed edges");
                    p.canon_reversed = rev == Tri::Yes;
                }
                p.canon = q;
                p.coincident = o.coincident = true;
                ++im.stats.unified_pieces;
                break;
            }
            // A curved side against an arc between the same two vertices: the
            // circle of a cap that was moved onto a plane is the intersection
            // curve of its cylinder with that plane, and the other solid may
            // have the same circle as a side of a cylinder patch.
            // (A straight side counts too: the curve of a plane and a cylinder
            // along one of its rulings is straight, though it is an arc.)
            if (p.coincident || (e.arc && straight_piece(p)) || p.v0 == p.v1) continue;
            const auto mt = mixed.find({std::min(p.v0, p.v1), std::max(p.v0, p.v1)});
            if (mt == mixed.end()) continue;
            for (int q : mt->second) {
                Imprint::OldPiece& o = im.pieces[static_cast<std::size_t>(q)];
                const BEdge& oe = S(0).edges[static_cast<std::size_t>(o.edge)];
                if (o.coincident || oe.arc == e.arc) continue;
                const Tri same = e.arc ? side_on_arc(o, oe, p) : side_on_arc(p, e, o);
                if (same == Tri::Unknown) return fail("could not compare two edge pieces");
                if (same == Tri::No) continue;
                // The arc stands for both: faces on either of its surfaces, or
                // on surfaces equivalent to them, can carry it, while the side
                // may belong to a carrier that was left behind (the cap's old
                // rational carrier).
                if (e.arc) o.canon = static_cast<int>(i);
                else p.canon = q;
                p.coincident = o.coincident = true;
                ++im.stats.unified_pieces;
                break;
            }
        }
        return true;
    }

    /// Is the curved side piece s (of edge se) the same curve as the arc piece
    /// a of the other solid? The side's midpoint is a point of both of the
    /// arc's surfaces: of one because the side's surface is that surface
    /// (same class), of the other for the same reason or because the side
    /// lies exactly in that plane. Carried into the arc's parameters, it must
    /// lie on the arc's certified curve, inside the piece.
    Tri side_on_arc(const Imprint::OldPiece& s, const BEdge& se, const Imprint::OldPiece& a) {
        int a0 = -1, a1 = -1;
        edge_carriers(a.solid, a.edge, a0, a1);
        const Interval t(0.5 * (num::mid(s.t0) + num::mid(s.t1)));
        Box4 P;
        for (int i = 0; i < 2; ++i) {
            const int ai = i == 0 ? a0 : a1;
            bool have = false;
            for (int ks = 0; ks < 2 && !have; ++ks) {
                const int gs = gid(s.solid, se.on[ks].carrier);
                if (!equiv(gs, ai)) continue;
                Interval u, v;
                side_uv(se.on[ks].side, se.on[ks].reversed ? Interval(1.0) - t : t, u, v);
                if (gs != ai) map(gs, ai).apply(Interval(u), Interval(v), u, v);
                P[static_cast<std::size_t>(2 * i)] = u;
                P[static_cast<std::size_t>(2 * i + 1)] = v;
                have = true;
            }
            if (!have && is_planar[static_cast<std::size_t>(ai)]) {
                const int gs = gid(s.solid, se.on[0].carrier);
                if (!side_in_plane(im.carriers[static_cast<std::size_t>(gs)], se.on[0].side, im.carriers[static_cast<std::size_t>(ai)]))
                    return Tri::No;
                Interval u, v;
                side_uv(se.on[0].side, se.on[0].reversed ? Interval(1.0) - t : t, u, v);
                const auto h = icarrier[static_cast<std::size_t>(gs)].eval(u, v);
                if (!plane_params(im.carriers[static_cast<std::size_t>(ai)], {h.x / h.w, h.y / h.w, h.z / h.w},
                                  P[static_cast<std::size_t>(2 * i)], P[static_cast<std::size_t>(2 * i + 1)]))
                    return Tri::Unknown;
                have = true;
            }
            if (!have) return Tri::No;
        }
        bool on = false;
        CurvePos pos;
        if (!locate_on_curve(*a.geom->base, P, on, pos)) return Tri::Unknown;
        if (!on) return Tri::No;
        return in_arc(*a.geom, pos);
    }

    /// Is side piece p (of solid 1, edge e) the same curve as side piece o (of
    /// solid 0, edge oe), where their sides belong to different carriers of
    /// the same surfaces? p's midpoint is carried into each of o's carriers of
    /// the same class and must lie exactly on o's side there, strictly inside
    /// o's parameter range.
    Tri same_side_curve(const Imprint::OldPiece& p, const BEdge& e, const Imprint::OldPiece& o, const BEdge& oe) {
        const double t = 0.5 * (num::mid(p.t0) + num::mid(p.t1));
        for (int kp = 0; kp < 2; ++kp) {
            const int gp = gid(1, e.on[kp].carrier);
            Interval u, v;
            side_uv(e.on[kp].side, e.on[kp].reversed ? Interval(1.0) - Interval(t) : Interval(t), u, v);
            for (int ko = 0; ko < 2; ++ko) {
                const int go = gid(0, oe.on[ko].carrier);
                if (!equiv(gp, go)) continue;
                Chart c{go, u, v};
                if (gp != go) map(gp, go).apply(u, v, c.u, c.v);
                Interval tau;
                const Tri on = on_side(c, oe.on[ko].side, tau);
                if (on != Tri::Yes) return on;
                const Interval te = oe.on[ko].reversed ? Interval(1.0) - tau : tau;
                if (o.t0.hi < te.lo && te.hi < o.t1.lo) return Tri::Yes;
                if (te.hi < o.t0.lo || te.lo > o.t1.hi) return Tri::No;
                return Tri::Unknown;
            }
        }
        return Tri::No; // no carrier of p is of the same surface as one of o
    }

    /// Do two closed pieces on one curve (p of solid 1, o of solid 0) run in
    /// opposite directions? Pieces of one certified curve both run in its
    /// travel order. Otherwise the curve's tangent at a certified point of p is
    /// enclosed, and its sign against the slicing direction of each tube,
    /// with each span's direction, gives each piece's direction of travel.
    Tri closed_relation(const Imprint::OldPiece& p, const Imprint::OldPiece& o) {
        if (p.geom->base == o.geom->base) return Tri::No;
        int g0 = -1, g1 = -1, h0 = -1, h1 = -1;
        edge_carriers(1, p.edge, g0, g1);
        edge_carriers(0, o.edge, h0, h1);
        const bool swap = g0 != h0;
        if (swap ? (g0 != h1 || g1 != h0) : g1 != h1) return Tri::Unknown;
        const ArcGeom& pa = *p.geom;
        const Span& sp = pa.spans[pa.spans.size() / 2];
        const TubeSegment& seg = pa.base->curve->segs[static_cast<std::size_t>(sp.seg)];
        Box4 mid;
        if (!arc_point(pa, mid)) return Tri::Unknown;
        Interval J[3][4];
        pa.base->sys->J_over(mid, J);
        Interval t[4];
        for (int i = 0; i < 4; ++i) {
            Interval m[3][3];
            for (int r = 0; r < 3; ++r)
                for (int c = 0, cc = 0; c < 4; ++c)
                    if (c != i) m[r][cc++] = J[r][c];
            const Interval det = det3i(m);
            t[i] = i % 2 == 0 ? det : -det;
        }
        Interval mt(0.0);
        for (int i = 0; i < 4; ++i) mt += Interval(seg.m[i]) * t[i];
        if (num::contains_zero(mt) || sp.cb == sp.ca) return Tri::Unknown;
        const int dir_p = ((sp.cb > sp.ca) == (mt.lo > 0.0)) ? 1 : -1;

        const Box4 P = swap ? Box4{mid[2], mid[3], mid[0], mid[1]} : mid;
        bool on = false;
        CurvePos pos;
        if (!locate_on_curve(*o.geom->base, P, on, pos) || !on) return Tri::Unknown;
        const Span& so = o.geom->base->spans[static_cast<std::size_t>(pos.span)];
        const TubeSegment& sego = o.geom->base->curve->segs[static_cast<std::size_t>(so.seg)];
        const Interval to[4] = {swap ? t[2] : t[0], swap ? t[3] : t[1], swap ? t[0] : t[2], swap ? t[1] : t[3]};
        Interval mo(0.0);
        for (int i = 0; i < 4; ++i) mo += Interval(sego.m[i]) * to[i];
        if (num::contains_zero(mo) || so.cb == so.ca) return Tri::Unknown;
        const int dir_o = ((so.cb > so.ca) == (mo.lo > 0.0)) ? 1 : -1;
        return dir_p == dir_o ? Tri::No : Tri::Yes;
    }

    // --- 5. planar map of each face ----------------------------------------------------
    /// Direction in which an arc leaves a vertex, in the parameter plane of surface g.
    bool arc_dir(const ArcGeom& arc, int c0, int c1, int g, bool forward, int vertex, Dir& d) {
        int co = c0 == g ? 0 : (c1 == g ? 2 : -1);
        int via = -1; // a surface of the arc in g's plane, when it is not g itself
        if (co < 0) {
            if (equiv(c0, g)) { co = 0; via = c0; }
            else if (equiv(c1, g)) { co = 2; via = c1; }
            else return false;
        }
        const Chart* a = chart(vertex, c0);
        const Chart* b = chart(vertex, c1);
        if (!a || !b) return false;
        const Box4 P = {a->u, a->v, b->u, b->v};
        Interval J[3][4];
        arc.base->sys->J_over(P, J);
        Interval t[4];
        for (int i = 0; i < 4; ++i) {
            Interval m[3][3];
            for (int r = 0; r < 3; ++r)
                for (int c = 0, cc = 0; c < 4; ++c)
                    if (c != i) m[r][cc++] = J[r][c];
            const Interval det = det3i(m);
            t[i] = i % 2 == 0 ? det : -det;
        }
        const Span& sp = forward ? arc.spans.front() : arc.spans.back();
        const TubeSegment& seg = arc.base->curve->segs[static_cast<std::size_t>(sp.seg)];
        Interval mt(0.0);
        for (int i = 0; i < 4; ++i) mt += Interval(seg.m[i]) * t[i];
        if (num::contains_zero(mt)) return false;
        const bool flip = (mt.hi < 0.0) != seg.reversed; // t points against the travel direction
        const bool negate = flip == forward;              // leaving backward runs against the travel
        d = {negate ? -t[co] : t[co], negate ? -t[co + 1] : t[co + 1]};
        if (via >= 0) {
            Dir m;
            map(via, g).apply_linear(d.x, d.y, m.x, m.y);
            d = m;
        }
        return true;
    }

    /// Geometry of an arc in the parameters of face surface g: the curve
    /// coordinate offset, and a map if the arc lies on a surface in g's plane.
    bool arc_frame(int c0, int c1, int g, int& co, const Affine*& m) {
        m = nullptr;
        if (c0 == g) { co = 0; return true; }
        if (c1 == g) { co = 2; return true; }
        if (equiv(c0, g)) { co = 0; m = &map(c0, g); return true; }
        if (equiv(c1, g)) { co = 2; m = &map(c1, g); return true; }
        return false;
    }

    /// Direction in which a half-edge leaves its start node.
    bool leave_dir(int g, const HalfEdge& h, Dir& d) {
        if (h.piece.fresh) {
            const Imprint::NewArc& a = im.arcs[static_cast<std::size_t>(h.piece.id)];
            return arc_dir(*a.geom, a.carrier[0], a.carrier[1], g, h.fwd, h.from, d);
        }
        const Imprint::OldPiece& p = im.pieces[static_cast<std::size_t>(h.piece.src)];
        const BEdge& e = S(p.solid).edges[static_cast<std::size_t>(p.edge)];
        int g0 = -1, g1 = -1;
        edge_carriers(p.solid, p.edge, g0, g1);
        if (e.arc) return arc_dir(*p.geom, g0, g1, g, h.fwd, h.from, d);
        // A straight piece along a side of the square: an exact direction.
        const BEdge::SideRef* ref = g0 == g ? &e.on[0] : (g1 == g ? &e.on[1] : nullptr);
        if (!ref) {
            // A straight piece elsewhere in g's plane: the difference of its end
            // points' parameters encloses its direction.
            const Chart* a = chart(h.from, g);
            const Chart* b = chart(h.to, g);
            if (!a || !b) return false;
            d = {b->u - a->u, b->v - a->v};
            // The piece lies on its edge's line: where the edge's own end
            // vertices have exact parameters here, its chord is the exact
            // direction, signed by the piece's.
            Chart c0{}, c1{};
            if (chart_in_class(vmap[p.solid][static_cast<std::size_t>(e.v0)], g, c0) &&
                chart_in_class(vmap[p.solid][static_cast<std::size_t>(e.v1)], g, c1) && is_point(c0.u) && is_point(c0.v) &&
                is_point(c1.u) && is_point(c1.v)) {
                const Chart* e0 = &c0;
                const Chart* e1 = &c1;
                double dx = 0.0, dy = 0.0;
                if (flat_detail_sum(e1->u.lo, -e0->u.lo, dx) && flat_detail_sum(e1->v.lo, -e0->v.lo, dy) && (dx != 0.0 || dy != 0.0)) {
                    const Dir de{Interval(dx), Interval(dy)};
                    const Interval s = dot2(de, d);
                    if (s.lo > 0.0) d = de;
                    else if (s.hi < 0.0) d = Dir{Interval(-dx), Interval(-dy)};
                }
            }
            return true;
        }
        const bool along_u = ref->side == 0 || ref->side == 2;
        double sgn = ref->reversed ? -1.0 : 1.0; // edge parameter against side parameter
        if (!h.fwd) sgn = -sgn;
        d = {Interval(along_u ? sgn : 0.0), Interval(along_u ? 0.0 : sgn)};
        return true;
    }

    /// At a vertex on a side of g's parameter square, a curve of the planar
    /// map may leave tangent to that side: a plane's trace along the circle
    /// of a cylinder it touches, say. Its direction enclosure then straddles
    /// the side's direction e and cannot be ordered against the side by
    /// interval signs. But the curve lies in the face, hence in the closed
    /// square, so its true direction is e turned towards the inward normal n
    /// by an angle somewhere in [0, theta], theta bounded by the enclosure
    /// (0 for a true tangency, where the curve still runs on n's side of the
    /// side to second order). If every other direction at the vertex lies
    /// outside that cone, provably, then e turned by an infinitesimal angle
    /// orders the same way: the direction is replaced by that symbolic
    /// perturbation of e (Dir::pert), compared exactly.
    /// Two curves tangent at one vertex are not resolved (returns false).
    bool tangent_proxies(int g, int vid, Dir& d0, std::vector<Dir>& cd) const {
        const Chart* c = chart(vid, g);
        if (!c) return false;
        struct SideDir {
            Dir e, n;
        };
        std::vector<SideDir> sides;
        auto exact_at = [](const Interval& x, double v) { return x.lo == v && x.hi == v; };
        const Interval z(0.0), p1(1.0), m1(-1.0);
        if (exact_at(c->u, 0.0)) { sides.push_back({{z, p1}, {p1, z}}); sides.push_back({{z, m1}, {p1, z}}); }
        if (exact_at(c->u, 1.0)) { sides.push_back({{z, p1}, {m1, z}}); sides.push_back({{z, m1}, {m1, z}}); }
        if (exact_at(c->v, 0.0)) { sides.push_back({{p1, z}, {z, p1}}); sides.push_back({{m1, z}, {z, p1}}); }
        if (exact_at(c->v, 1.0)) { sides.push_back({{p1, z}, {z, m1}}); sides.push_back({{m1, z}, {z, m1}}); }
        if (sides.empty()) return false;
        std::vector<Dir*> all = {&d0};
        for (Dir& d : cd) all.push_back(&d);
        std::vector<int> base(all.size(), -1);
        int tangent = 0;
        for (std::size_t i = 0; i < all.size(); ++i) {
            const Dir& d = *all[i];
            if (d.exact()) continue;
            for (std::size_t k = 0; k < sides.size(); ++k) {
                const Interval cr = cross2(sides[k].e, d), dt = dot2(sides[k].e, d);
                if (cr.lo <= 0.0 && cr.hi >= 0.0 && dt.lo > 0.0) { base[i] = static_cast<int>(k); ++tangent; break; }
            }
        }
        if (tangent == 0) return false;
        std::vector<std::array<Dir, 2>> en;
        for (const SideDir& sd : sides) en.push_back({sd.e, sd.n});
        return apply_proxies(all, base, en);
    }
    /// The cone argument of tangent_proxies: each direction i with base[i] >= 0
    /// is the exact direction e = en[base[i]][0] turned towards the exact
    /// normal n = en[base[i]][1] by an angle in [0, theta]. If every other
    /// direction is provably outside that cone, replace it by e + eps n.
    static bool apply_proxies(std::vector<Dir*>& all, const std::vector<int>& base, const std::vector<std::array<Dir, 2>>& en) {
        for (std::size_t i = 0; i < all.size(); ++i) {
            if (base[i] < 0) continue;
            const Dir& e = en[static_cast<std::size_t>(base[i])][0];
            const Dir& n = en[static_cast<std::size_t>(base[i])][1];
            const Interval en_cross = cross2(e, n); // exact: both are exact directions
            if (en_cross.lo == 0.0 && en_cross.hi == 0.0) return false;
            const int nside = en_cross.lo > 0.0 ? 1 : -1;
            for (std::size_t j = 0; j < all.size(); ++j) {
                if (j == i) continue;
                if (base[j] >= 0) return false; // two curves tangent at one vertex
                const Dir& o = *all[j];
                if (o.exact() && cross2(e, o).lo == 0.0 && cross2(e, o).hi == 0.0 && dot2(e, o).lo > 0.0) continue; // along e itself
                const Interval ce = cross2(e, o), ci = cross2(*all[i], o);
                const bool behind = nside > 0 ? ce.hi < 0.0 : ce.lo > 0.0;                       // far side of e from n
                const bool beyond = (nside > 0 ? ce.lo > 0.0 && ci.lo > 0.0 : ce.hi < 0.0 && ci.hi < 0.0); // past the whole cone
                const bool opposite = dot2(e, o).hi < 0.0;                                          // the cone is within 90 degrees of e
                if (!behind && !beyond && !opposite) return false;
            }
        }
        // e turned towards n by an infinitesimal angle (symbolic, see Dir).
        for (std::size_t i = 0; i < all.size(); ++i) {
            if (base[i] < 0) continue;
            const Dir& e = en[static_cast<std::size_t>(base[i])][0];
            const Dir& n = en[static_cast<std::size_t>(base[i])][1];
            *all[i] = Dir{e.x, e.y, n.x, n.y, true};
        }
        return true;
    }

    /// The vertex's parameters on g: its chart there, or one carried over
    /// from a surface of g's class (a vertex away from g's faces has none).
    bool chart_in_class(int v, int g, Chart& out) {
        if (const Chart* c = chart(v, g)) {
            out = *c;
            // An enclosure where the exact point is known: prefer the exact
            // parameters (they lie in the enclosure, which is checked).
            Chart x;
            if ((!is_point(c->u) || !is_point(c->v)) && exact_plane_chart(v, g, x) && num::contains(c->u, x.u.lo) &&
                num::contains(c->v, x.v.lo))
                out = x;
            return true;
        }
        for (const Chart& c : im.vertices[static_cast<std::size_t>(v)].charts)
            if (c.carrier != g && equiv(c.carrier, g)) {
                out.carrier = g;
                map(c.carrier, g).apply(c.u, c.v, out.u, out.v);
                return true;
            }
        return exact_plane_chart(v, g, out);
    }
    /// Exact parameters on the plane g of a vertex that is an exact point
    /// (a corner control point with a power-of-two weight) lying in it:
    /// computed in doubles and accepted only if g evaluates there to the point
    /// exactly. False otherwise (not an error: the caller has nothing exact).
    bool exact_plane_chart(int v, int g, Chart& out) {
        if (!is_planar[static_cast<std::size_t>(g)]) return false;
        const PolyPatchD& P = im.carriers[static_cast<std::size_t>(g)];
        using planar_detail::corner;
        const Vec4d c00 = corner(P, 0, 0), c10 = corner(P, 1, 0), c01 = corner(P, 0, 1);
        auto pow2 = [](double w) {
            int e = 0;
            return w > 0.0 && std::frexp(w, &e) == 0.5;
        };
        if (!pow2(c00.w) || c00.w != c10.w || c00.w != c01.w) return false;
        for (const Chart& c : im.vertices[static_cast<std::size_t>(v)].charts) {
            Vec4d p;
            if (!exact_corner(im.carriers[static_cast<std::size_t>(c.carrier)], c.u, c.v, p) || !pow2(p.w)) continue;
            const double q[3] = {p.x / p.w, p.y / p.w, p.z / p.w}; // exact: power-of-two weight
            const geom::Vec3<Interval> X{Interval(q[0]), Interval(q[1]), Interval(q[2])};
            Interval s, t;
            if (!plane_params(P, X, s, t)) return false;
            // The exact parameters of a dyadic point on a dyadic parallelogram
            // are dyadic: snap the estimate to a fine dyadic grid, then check.
            auto snap40 = [](double x) { return std::ldexp(std::nearbyint(std::ldexp(x, 40)), -40); };
            const double sm = snap40(num::mid(s)), tm = snap40(num::mid(t));
            const double a[3] = {c00.x / c00.w, c00.y / c00.w, c00.z / c00.w}, b[3] = {c10.x / c10.w, c10.y / c10.w, c10.z / c10.w},
                         d[3] = {c01.x / c01.w, c01.y / c01.w, c01.z / c01.w};
            bool exact = true;
            for (int k = 0; k < 3 && exact; ++k) {
                using E = num::Expansion;
                const E val = E(a[k]) + E(sm) * (E(b[k]) - E(a[k])) + E(tm) * (E(d[k]) - E(a[k]));
                exact = (val - E(q[k])).sign() == 0;
            }
            if (!exact) return false;
            out = {g, Interval(sm), Interval(tm)};
            return true;
        }
        return false;
    }
    /// a + b if it is exact.
    static bool flat_detail_sum(double a, double b, double& s) {
        s = a + b;
        return (num::Expansion(s) - num::Expansion(a) - num::Expansion(b)).sign() == 0;
    }

    /// Surfaces of the curve a half-edge runs on, and its arc geometry if any.
    void piece_surfaces(const HalfEdge& x, int& c0, int& c1, const ArcGeom*& geom) const {
        geom = nullptr;
        if (x.piece.fresh) {
            const Imprint::NewArc& a = im.arcs[static_cast<std::size_t>(x.piece.id)];
            c0 = a.carrier[0];
            c1 = a.carrier[1];
            geom = a.geom.get();
            return;
        }
        const Imprint::OldPiece& p = im.pieces[static_cast<std::size_t>(x.piece.src)];
        edge_carriers(p.solid, p.edge, c0, c1);
        if (S(p.solid).edges[static_cast<std::size_t>(p.edge)].arc) geom = p.geom.get();
    }

    /// The same for a curve tangent at a vertex to a straight piece that lies
    /// in a plane b, where the curve lies on a surface a that touches b only
    /// along a side (a proven contact, `touches`). Then the curve never
    /// crosses b's trace in this face's plane, so its side of that line is
    /// read off any certified point of it, and the cone argument applies with
    /// e the line's exact direction and n the normal towards that side.
    bool line_tangent_proxies(int g, int vid, Dir& d0, std::vector<Dir>& cd, const HalfEdge& back,
                              const std::vector<int>& cand, const std::vector<HalfEdge>& he) {
        if (touches.empty()) return false;
        const Chart* vc = chart(vid, g);
        if (!vc) return false;
        std::vector<Dir*> all = {&d0};
        std::vector<const HalfEdge*> hs = {&back};
        for (std::size_t c = 0; c < cd.size(); ++c) {
            all.push_back(&cd[c]);
            hs.push_back(&he[static_cast<std::size_t>(cand[c])]);
        }
        std::vector<int> base(all.size(), -1);
        std::vector<std::array<Dir, 2>> en;
        int tangent = 0;
        for (std::size_t i = 0; i < all.size(); ++i) {
            // A curve (its direction an enclosure, or exactly the line's when
            // the tangent was computed exactly).
            int a0 = -1, a1 = -1;
            const ArcGeom* geom = nullptr;
            piece_surfaces(*hs[i], a0, a1, geom);
            if (!geom) continue;
            for (std::size_t j = 0; j < all.size() && base[i] < 0; ++j) {
                if (j == i || !all[j]->exact()) continue;
                int b0 = -1, b1 = -1;
                const ArcGeom* lg = nullptr;
                piece_surfaces(*hs[j], b0, b1, lg);
                if (lg) continue; // a straight piece
                bool proven = false;
                // The straight piece's side, if it is one: it may lie exactly in
                // the touched plane without being on a surface of its class.
                const BEdge::SideRef* side_ref = nullptr;
                if (!hs[j]->piece.fresh) {
                    const Imprint::OldPiece& lp = im.pieces[static_cast<std::size_t>(hs[j]->piece.src)];
                    const BEdge& le = S(lp.solid).edges[static_cast<std::size_t>(lp.edge)];
                    if (!le.arc) side_ref = &le.on[0];
                }
                for (const Touch& t : touches) {
                    if (t.a != a0 && t.a != a1) continue;
                    if (equiv(t.b, b0) || equiv(t.b, b1)) proven = true;
                    else if (side_ref && is_planar[static_cast<std::size_t>(t.b)]) {
                        const Imprint::OldPiece& lp = im.pieces[static_cast<std::size_t>(hs[j]->piece.src)];
                        const int sc = gid(lp.solid, side_ref->carrier);
                        proven = proven || side_in_plane(im.carriers[static_cast<std::size_t>(sc)], side_ref->side,
                                                         im.carriers[static_cast<std::size_t>(t.b)]);
                    }
                }
                note("  line tangent: curve on " + str(a0) + "," + str(a1) + " against line on " + str(b0) + "," + str(b1) +
                     (proven ? ", contact proven" : ", no proven contact"));
                if (!proven) continue;
                // Along the line, one way or the other.
                for (double sgn : {1.0, -1.0}) {
                    const Dir e{Interval(sgn * all[j]->x.lo), Interval(sgn * all[j]->y.lo)};
                    const Interval cr = cross2(e, *all[i]), dt = dot2(e, *all[i]);
                    if (opt.trace) {
                        char buf[256];
                        std::snprintf(buf, sizeof buf, "  line tangent: line (%.17g, %.17g), curve [%.17g, %.17g] x [%.17g, %.17g]%s, cross [%.3g, %.3g], dot [%.3g, %.3g]",
                                      e.x.lo, e.y.lo, all[i]->x.lo, all[i]->x.hi, all[i]->y.lo, all[i]->y.hi, all[i]->pert ? " (perturbed)" : "",
                                      cr.lo, cr.hi, dt.lo, dt.hi);
                        note(buf);
                    }
                    if (!(cr.lo <= 0.0 && cr.hi >= 0.0 && dt.lo > 0.0)) continue;
                    // Which side of the line the curve runs on: a certified point of it.
                    Box4 mid;
                    int co = 0;
                    const Affine* m = nullptr;
                    if (!arc_point(*geom, mid) || !arc_frame(a0, a1, g, co, m)) return false;
                    Interval mu = mid[static_cast<std::size_t>(co)], mv = mid[static_cast<std::size_t>(co + 1)];
                    if (m) m->apply(Interval(mu), Interval(mv), mu, mv);
                    const Interval side = cross2(e, Dir{mu - vc->u, mv - vc->v});
                    if (num::contains_zero(side)) { note("  line tangent: side of the curve undecided"); return false; }
                    const Dir n = side.lo > 0.0 ? Dir{-e.y, e.x} : Dir{e.y, -e.x};
                    base[i] = static_cast<int>(en.size());
                    en.push_back({e, n});
                    ++tangent;
                    break;
                }
            }
        }
        if (tangent == 0) { note("  line tangent: no tangent pair"); return false; }
        const bool ok = apply_proxies(all, base, en);
        if (!ok) note("  line tangent: another direction lies in the cone");
        return ok;
    }

    bool build_regions() {
        // Edge pieces by the boxes of their edges, per solid: a face looks up
        // the other solid's pieces near it.
        for (int k = 0; k < 2; ++k) {
            std::vector<Vec3<Interval>> boxes;
            piece_ids[k].clear();
            for (int pid = 0; pid < static_cast<int>(im.pieces.size()); ++pid) {
                const Imprint::OldPiece& p = im.pieces[static_cast<std::size_t>(pid)];
                if (p.solid != k || p.coincident) continue;
                piece_ids[k].push_back(pid);
                boxes.push_back(edge_box[k][static_cast<std::size_t>(p.edge)]);
            }
            piece_index[k] = geom::BoxIndex(std::move(boxes));
        }
        for (int k = 0; k < 2; ++k)
            for (int f = 0; f < static_cast<int>(S(k).faces.size()); ++f)
                if (!build_face(k, f)) return false;

        im.stats.regions = static_cast<int>(im.regions.size());
        return true;
    }

    bool build_face(int k, int f) {
        const Brep& own = S(k);
        const BFace& face = own.faces[static_cast<std::size_t>(f)];
        const int g = gid(k, face.carrier);
        const std::string name = "face " + str(f) + " of solid " + str(k);
        std::map<int, Node> nodes;
        std::vector<HalfEdge> he;

        auto node = [&](int vid) -> Node* {
            const auto it = nodes.find(vid);
            if (it != nodes.end()) return &it->second;
            const Chart* c = chart(vid, g);
            if (!c) return nullptr;
            Node nd;
            nd.box[0] = c->u;
            nd.box[1] = c->v;
            nd.uv[0] = num::mid(c->u);
            nd.uv[1] = num::mid(c->v);
            return &(nodes[vid] = nd);
        };
        auto add = [&](int from, int to, const Piece& piece, bool fwd, bool own_piece) -> int {
            Node* a = node(from);
            if (!a || !node(to)) return -1;
            const int id = static_cast<int>(he.size());
            nodes[from].out.push_back(id);
            HalfEdge h;
            h.from = from;
            h.to = to;
            h.piece = piece;
            h.fwd = fwd;
            h.own = own_piece;
            he.push_back(h);
            return id;
        };
        auto add_pair = [&](int v0, int v1, Piece piece) {
            const int a = add(v0, v1, piece, true, false);
            Piece back = piece;
            back.forward = !piece.forward;
            const int b = add(v1, v0, back, false, false);
            if (a < 0 || b < 0) return false;
            he[static_cast<std::size_t>(a)].twin = b;
            he[static_cast<std::size_t>(b)].twin = a;
            return true;
        };


        // A vertex on the face's own boundary lies on the face's surface. Where
        // that is a plane and the vertex has no parameters on it (it was found
        // on an edge stored as a side of other surfaces through the same
        // line), they are computed from its certified position.
        auto ensure_plane_chart = [&](int vid) {
            if (chart(vid, g) || !is_planar[static_cast<std::size_t>(g)]) return;
            const Vec3<Interval> P = vertex_box(vid);
            Chart m;
            m.carrier = g;
            if (!plane_params(im.carriers[static_cast<std::size_t>(g)], P, m.u, m.v)) return;
            im.vertices[static_cast<std::size_t>(vid)].charts.push_back(m);
        };
        // The face's own loops, split at the vertices on them.
        for (const auto& loop : face.loops)
            for (const BCoedge& ce : loop)
                for (int id : edge_pieces[k][static_cast<std::size_t>(ce.edge)]) {
                    ensure_plane_chart(im.pieces[static_cast<std::size_t>(id)].v0);
                    ensure_plane_chart(im.pieces[static_cast<std::size_t>(id)].v1);
                }
        for (const auto& loop : face.loops)
            for (const BCoedge& ce : loop) {
                const std::vector<int>& ids = edge_pieces[k][static_cast<std::size_t>(ce.edge)];
                const int m = static_cast<int>(ids.size());
                for (int i = 0; i < m; ++i) {
                    const int id = ids[static_cast<std::size_t>(ce.reversed ? m - 1 - i : i)];
                    const Imprint::OldPiece& p = im.pieces[static_cast<std::size_t>(id)];
                    const Imprint::OldPiece& c = im.pieces[static_cast<std::size_t>(p.canon)];
                    const bool fwd = !ce.reversed;
                    // Does the canonical piece run the same way? Read from the end
                    // vertices, or, for a closed piece, decided when it was unified.
                    const bool same = p.v0 != p.v1 ? c.v0 == p.v0 : !p.canon_reversed;
                    const Piece piece{false, p.canon, same ? fwd : !fwd, id};
                    if (add(fwd ? p.v0 : p.v1, fwd ? p.v1 : p.v0, piece, fwd, true) < 0) {
                        if (opt.trace) {
                            const int bad = !chart(p.v0, g) ? p.v0 : p.v1;
                            int e0 = -1, e1 = -1;
                            edge_carriers(p.solid, p.edge, e0, e1);
                            std::string s = name + " (surface " + sk(g) + ", class " + str(rep[static_cast<std::size_t>(g)]) + "): vertex " +
                                            str(bad) + " of its edge " + str(p.edge) + " (on " + sk(e0) + " class " +
                                            str(rep[static_cast<std::size_t>(e0)]) + ", " + sk(e1) + " class " + str(rep[static_cast<std::size_t>(e1)]) +
                                            ") has no parameters here (face plane and edge planes coplanar: " +
                                            str(coplanar(im.carriers[static_cast<std::size_t>(g)], im.carriers[static_cast<std::size_t>(e0)]) ? 1 : 0) +
                                            str(coplanar(im.carriers[static_cast<std::size_t>(g)], im.carriers[static_cast<std::size_t>(e1)]) ? 1 : 0) +
                                            "); it has them on";
                            for (const Chart& c2 : im.vertices[static_cast<std::size_t>(bad)].charts)
                                s += " " + sk(c2.carrier) + "[" + std::to_string(num::mid(c2.u)) + "," + std::to_string(num::mid(c2.v)) + "]";
                            const Vec3d pt = im.vertices[static_cast<std::size_t>(bad)].point;
                            s += "; at (" + std::to_string(pt.x) + ", " + std::to_string(pt.y) + ", " + std::to_string(pt.z) + ")";
                            note(s);
                        }
                        return fail(name + ": a vertex has no parameters on this face");
                    }
                }
            }
        // The new edges in this face.
        for (int ai = 0; ai < static_cast<int>(im.arcs.size()); ++ai) {
            const Imprint::NewArc& arc = im.arcs[static_cast<std::size_t>(ai)];
            for (int i = 0; i < 2; ++i) {
                if (arc.solid_on[i] != k || arc.face[i] != f) continue;
                if (!add_pair(arc.v0, arc.v1, Piece{true, ai, true, -1}))
                    return fail(name + ": a new vertex has no parameters on this face");
            }
        }
        // Edge pieces of the other solid that lie on this face's surface and
        // inside the face. The surface is usually shared, but an edge can lie
        // on a surface that carries no face of its solid: a seam between two
        // faces on one surface, left where an earlier operation cut along it.
        for (int hit : piece_index[1 - k].query(face_box[k][static_cast<std::size_t>(f)])) {
                const int pid = piece_ids[1 - k][static_cast<std::size_t>(hit)];
                const Imprint::OldPiece& p = im.pieces[static_cast<std::size_t>(pid)];
                const BEdge& e = S(p.solid).edges[static_cast<std::size_t>(p.edge)];
                int g0 = -1, g1 = -1;
                edge_carriers(p.solid, p.edge, g0, g1);
                Interval mu, mv; // a certified point of the piece, in g's parameters
                if (!e.arc) {
                    // A straight side lying in g's plane: of a surface in that plane,
                    // or itself in it. (A side of g's own square is on this face's
                    // boundary or outside it, and then the test below says so.)
                    if (!side_lies_in(g0, e.on[0].side, g1, g)) {
                        // Or a straight edge lying in this (curved) surface along a
                        // line of constant parameter: a straight piece here too.
                        // Every vertex on it lies on the surface.
                        const IsoLine* line = iso_line_on(p.solid, p.edge, g);
                        if (!line) continue;
                        iso_vertex_chart(p.solid, p.edge, *line, p.v0, p.t0);
                        iso_vertex_chart(p.solid, p.edge, *line, p.v1, p.t1);
                    }
                    const Chart* a = chart(p.v0, g);
                    const Chart* b = chart(p.v1, g);
                    if (!a || !b) {
                        // Charts are spread only to surfaces whose faces are near
                        // the vertex. An end far from this face puts the piece
                        // outside it: a piece does not cross a face's boundary
                        // (it was split where it does).
                        const int far = !a ? p.v0 : p.v1;
                        if (solve::detail::boxes_disjoint(vertex_box(far), face_box[k][static_cast<std::size_t>(f)])) continue;
                        return fail(name + ": a vertex of the other solid has no parameters on this face");
                    }
                    mu = (a->u + b->u) * Interval(0.5);
                    mv = (a->v + b->v) * Interval(0.5);
                } else {
                    int co = -1;
                    const Affine* m = nullptr;
                    if (!arc_frame(g0, g1, g, co, m)) continue;
                    Box4 mid;
                    if (!arc_point(*p.geom, mid)) return fail(name + ": no certified point could be placed on an edge piece");
                    mu = mid[static_cast<std::size_t>(co)];
                    mv = mid[static_cast<std::size_t>(co + 1)];
                    if (m) m->apply(Interval(mu), Interval(mv), mu, mv);
                }
                Where w = locate(fgeom[k][f], mu, mv, opt.refine_depth);
                if (!e.arc && w == Where::Undecided) {
                    // A piece does not cross this face's boundary (it was split
                    // where it does), so any point of it decides: try others
                    // where the midpoint lies in the box of a boundary curve.
                    const Chart* a = chart(p.v0, g);
                    const Chart* b = chart(p.v1, g);
                    for (double s : {0.25, 0.75, 0.125, 0.875, 0.375, 0.625}) {
                        const Interval is(s);
                        const Where ws = locate(fgeom[k][f], a->u + is * (b->u - a->u), a->v + is * (b->v - a->v), opt.refine_depth);
                        if (ws == Where::Undecided) continue;
                        w = ws;
                        break;
                    }
                }
                if (opt.trace && !e.arc && !equiv(g, g0) && !equiv(g, g1))
                    note(name + ": piece " + str(pid) + " of the other solid (" + str(p.v0) + " -> " + str(p.v1) +
                         ") lies in its plane; its midpoint is " +
                         (w == Where::Inside ? "inside" : w == Where::Outside ? "outside" : w == Where::Boundary ? "on the boundary" : "undecided"));
                // A straight piece along this face's boundary, proven on a
                // boundary piece. (Pieces that are one curve with a piece of this
                // solid are not looked at here at all.) An undecided piece is
                // reported below, never taken for a boundary piece.
                if (!e.arc && w == Where::Boundary) continue;
                if (w == Where::Outside) continue;
                if (w != Where::Inside) {
                    if (opt.trace) {
                        std::string s = name + " (surface " + sk(g) + "): piece " + str(pid) + " of the other solid's edge " +
                                        str(p.edge) + " (" + (e.arc ? "arc" : "side") + (straight_piece(p) ? ", straight," : "") +
                                        " on surfaces " + sk(g0) + ", " + sk(g1) + "), vertices " + str(p.v0) + " -> " + str(p.v1) +
                                        ", is undecided; own boundary:";
                        for (const HalfEdge& x : he) {
                            if (!x.own) continue;
                            const Imprint::OldPiece& q = im.pieces[static_cast<std::size_t>(x.piece.src)];
                            int h0 = -1, h1 = -1;
                            edge_carriers(q.solid, q.edge, h0, h1);
                            s += " " + str(x.from) + ">" + str(x.to) + (S(q.solid).edges[static_cast<std::size_t>(q.edge)].arc ? " arc " : " side ") +
                                 sk(h0) + "," + sk(h1) + ";";
                        }
                        note(s);
                    }
                    return fail(name + ": could not decide whether an edge of the other solid lies inside it");
                }
                if (!add_pair(p.v0, p.v1, Piece{false, pid, true, pid}))
                    return fail(name + ": a vertex of the other solid has no parameters on this face");
            }

        // A straight piece of the other solid that lies in this face along a
        // line of its surface (iso_lines), where that solid only touches the
        // surface (proven), and that ends inside the face with nothing
        // continuing it: it bounds nothing here. Such pieces are taken out of
        // the map and recorded as line contacts. A piece that is continued, or
        // one where the surface passes between the edge's faces, stays.
        std::vector<int> contact_pieces;
        {
            std::vector<char> removed(he.size(), 0);
            for (bool again = true; again;) {
                again = false;
                for (std::size_t hi = 0; hi < he.size(); ++hi) {
                    const HalfEdge& h = he[hi];
                    if (removed[hi] || h.own || h.piece.fresh || h.twin < 0) continue;
                    bool continued = false;
                    for (int o : nodes[h.to].out)
                        continued = continued || (o != h.twin && !removed[static_cast<std::size_t>(o)]);
                    if (continued) continue;
                    const Imprint::OldPiece& p = im.pieces[static_cast<std::size_t>(h.piece.src)];
                    const IsoLine* line = iso_line_on(p.solid, p.edge, g);
                    if (!line || iso_touching(p.solid, p.edge, *line) != Tri::Yes) continue;
                    removed[hi] = removed[static_cast<std::size_t>(h.twin)] = 1;
                    contact_pieces.push_back(h.piece.src);
                    again = true;
                }
            }
            if (!contact_pieces.empty()) {
                std::vector<int> index(he.size(), -1);
                std::vector<HalfEdge> kept;
                for (std::size_t hi = 0; hi < he.size(); ++hi)
                    if (!removed[hi]) { index[hi] = static_cast<int>(kept.size()); kept.push_back(he[hi]); }
                for (HalfEdge& h : kept)
                    if (h.twin >= 0) h.twin = index[static_cast<std::size_t>(h.twin)];
                he = std::move(kept);
                for (auto& kv : nodes) kv.second.out.clear();
                for (std::size_t hi = 0; hi < he.size(); ++hi) nodes[he[hi].from].out.push_back(static_cast<int>(hi));
            }
        }

        // Successor of each half-edge around the region on its left: at the end
        // node, the first outgoing curve met turning clockwise from the way back.
        const int nH = static_cast<int>(he.size());
        std::vector<Dir> dirs(static_cast<std::size_t>(nH));
        for (int hi = 0; hi < nH; ++hi)
            if (!leave_dir(g, he[static_cast<std::size_t>(hi)], dirs[static_cast<std::size_t>(hi)]))
                return fail(name + ": the direction of a curve at a vertex could not be determined");
        // The intersection curve a half-edge runs on, as its pair of surfaces
        // (-1 for a straight piece). Two half-edges on the same curve that meet
        // at a vertex continue each other: their directions are opposite by
        // construction, which intervals alone could not prove.
        auto curve_of = [&](const HalfEdge& h) -> std::pair<int, int> {
            int g0 = -1, g1 = -1;
            if (h.piece.fresh) {
                const Imprint::NewArc& a = im.arcs[static_cast<std::size_t>(h.piece.id)];
                g0 = a.carrier[0];
                g1 = a.carrier[1];
            } else {
                const Imprint::OldPiece* p = &im.pieces[static_cast<std::size_t>(h.piece.src)];
                if (!S(p->solid).edges[static_cast<std::size_t>(p->edge)].arc) {
                    // A curved side that is one curve with an arc of the other
                    // solid (which then stands for it) runs on that arc's curve.
                    p = &im.pieces[static_cast<std::size_t>(h.piece.id)];
                    if (!S(p->solid).edges[static_cast<std::size_t>(p->edge)].arc) return {-1, -1};
                }
                edge_carriers(p->solid, p->edge, g0, g1);
            }
            g0 = rep[static_cast<std::size_t>(g0)];
            g1 = rep[static_cast<std::size_t>(g1)];
            return {std::min(g0, g1), std::max(g0, g1)};
        };
        // Straight pieces on one line also continue each other. A side has two
        // exact points of its line (its end control points); an edge on two
        // planes lies on their line. Decided exactly.
        struct Line {
            bool has_points = false;
            Vec4d p, q;
            int plane[2] = {-1, -1};
        };
        auto line_of = [&](const HalfEdge& h, Line& L) {
            int g0 = -1, g1 = -1;
            if (h.piece.fresh) {
                const Imprint::NewArc& a = im.arcs[static_cast<std::size_t>(h.piece.id)];
                g0 = a.carrier[0];
                g1 = a.carrier[1];
            } else {
                const Imprint::OldPiece& p = im.pieces[static_cast<std::size_t>(h.piece.src)];
                const BEdge& e = S(p.solid).edges[static_cast<std::size_t>(p.edge)];
                edge_carriers(p.solid, p.edge, g0, g1);
                if (!e.arc) {
                    if (!straight_piece(p)) return false;
                    const auto cp = side_control(im.carriers[static_cast<std::size_t>(g0)], e.on[0].side);
                    L.has_points = true;
                    L.p = cp.front();
                    L.q = cp.back();
                    return L.p.w == 1.0 && L.q.w == 1.0;
                }
            }
            // Two planes, or two flat surfaces in different planes of which
            // one or both are on a rational carrier.
            if (is_planar[static_cast<std::size_t>(g0)] && is_planar[static_cast<std::size_t>(g1)]) {
                if (equiv(g0, g1)) return false;
            } else if (!is_flat(g0) || !is_flat(g1) ||
                       same_plane_patches(im.carriers[static_cast<std::size_t>(g0)], im.carriers[static_cast<std::size_t>(g1)])) {
                return false;
            }
            L.plane[0] = g0;
            L.plane[1] = g1;
            return true;
        };
        auto same_line = [&](const HalfEdge& a, const HalfEdge& b) {
            Line la, lb;
            if (!line_of(a, la) || !line_of(b, lb)) return false;
            if (la.has_points && lb.has_points) return collinear(la.p, la.q, lb.p) && collinear(la.p, la.q, lb.q);
            if (!la.has_points && !lb.has_points)
                return std::min(rep[static_cast<std::size_t>(la.plane[0])], rep[static_cast<std::size_t>(la.plane[1])]) ==
                           std::min(rep[static_cast<std::size_t>(lb.plane[0])], rep[static_cast<std::size_t>(lb.plane[1])]) &&
                       std::max(rep[static_cast<std::size_t>(la.plane[0])], rep[static_cast<std::size_t>(la.plane[1])]) ==
                           std::max(rep[static_cast<std::size_t>(lb.plane[0])], rep[static_cast<std::size_t>(lb.plane[1])]);
            const Line& pts = la.has_points ? la : lb;
            const Line& pl = la.has_points ? lb : la;
            for (int i = 0; i < 2; ++i) {
                const PolyPatchD& P = im.carriers[static_cast<std::size_t>(pl.plane[i])];
                if (!planar_detail::on_plane(P, pts.p) || !planar_detail::on_plane(P, pts.q)) return false;
            }
            return true;
        };
        std::vector<int> next(static_cast<std::size_t>(nH), -1);
        for (int hi = 0; hi < nH; ++hi) {
            const HalfEdge& h = he[static_cast<std::size_t>(hi)];
            std::vector<int> cand;
            for (int o : nodes[h.to].out)
                if (o != h.twin) cand.push_back(o);
            if (cand.empty()) {
                if (opt.trace) {
                    note(name + " (surface " + str(g) + "): nothing leaves vertex " + str(h.to) + " after the half-edge from " +
                         str(h.from));
                    for (const HalfEdge& x : he) {
                        std::string kind = x.piece.fresh ? "new arc " + str(x.piece.id) : "";
                        if (!x.piece.fresh) {
                            const Imprint::OldPiece& p = im.pieces[static_cast<std::size_t>(x.piece.src)];
                            int g0 = -1, g1 = -1;
                            edge_carriers(p.solid, p.edge, g0, g1);
                            kind = std::string(x.own ? "own" : "other solid's") + " piece " + str(x.piece.src) + " of edge " +
                                   str(p.edge) + " (" + (S(p.solid).edges[static_cast<std::size_t>(p.edge)].arc ? "arc" : "side") +
                                   " on surfaces " + str(g0) + "," + str(g1) + (p.coincident ? ", shared" : "") + ")";
                        }
                        note("  " + str(x.from) + " -> " + str(x.to) + ": " + kind);
                    }
                    std::string vlist = "  vertex " + str(h.to) + " charts on surfaces";
                    for (const Chart& c : im.vertices[static_cast<std::size_t>(h.to)].charts) vlist += " " + str(c.carrier);
                    note(vlist);
                }
                return fail(name + ": the planar map is not closed");
            }
            if (cand.size() == 1) { next[static_cast<std::size_t>(hi)] = cand[0]; continue; }
            HalfEdge back = h;
            back.from = h.to;
            back.to = h.from;
            back.fwd = !h.fwd;
            Dir d0;
            if (!leave_dir(g, back, d0)) return fail(name + ": the direction of a curve at a vertex could not be determined");
            const auto on_curve = curve_of(h);
            std::vector<Dir> cd(cand.size());
            for (std::size_t c = 0; c < cand.size(); ++c) cd[c] = dirs[static_cast<std::size_t>(cand[c])];
            std::vector<char> is_continuation(cand.size(), 0);
            for (std::size_t c = 0; c < cand.size(); ++c) {
                const HalfEdge& o = he[static_cast<std::size_t>(cand[c])];
                // Two pieces of one operand edge, split at this vertex, continue
                // each other: the edge's curve is smooth there.
                bool same_edge = false;
                if (!h.piece.fresh && !o.piece.fresh) {
                    const Imprint::OldPiece& ph = im.pieces[static_cast<std::size_t>(h.piece.src)];
                    const Imprint::OldPiece& po = im.pieces[static_cast<std::size_t>(o.piece.src)];
                    same_edge = ph.solid == po.solid && ph.edge == po.edge && h.piece.src != o.piece.src;
                }
                is_continuation[c] = same_edge || (on_curve.first >= 0 && curve_of(o) == on_curve) || same_line(h, o);
            }
            auto order = [&](std::size_t& best_out) {
                std::vector<int> side(cand.size());
                for (std::size_t c = 0; c < cand.size(); ++c) side[c] = is_continuation[c] ? 1 : side_of(d0, cd[c]);
                std::size_t best = 0;
                for (std::size_t c = 1; c < cand.size(); ++c) {
                    const Tri t = cw_before(side[c], cd[c], side[best], cd[best]);
                    if (t == Tri::Unknown) return false;
                    if (t == Tri::Yes) best = c;
                }
                best_out = best;
                return true;
            };
            std::size_t best = 0;
            if (!order(best)) {
                // A curve tangent to a side of the square at this vertex.
                const Dir d0_kept = d0;
                const std::vector<Dir> cd_kept = cd;
                bool ordered = tangent_proxies(g, h.to, d0, cd) && order(best);
                if (!ordered) {
                    d0 = d0_kept;
                    cd = cd_kept;
                    ordered = line_tangent_proxies(g, h.to, d0, cd, back, cand, he) && order(best);
                }
                if (!ordered) {
                    if (opt.trace) {
                        const Chart* c = chart(h.to, g);
                        auto iv = [](const Interval& x) { return "[" + std::to_string(x.lo) + ", " + std::to_string(x.hi) + "]"; };
                        std::string s = name + " (surface " + str(g) + "): vertex " + str(h.to);
                        if (c) s += " at " + iv(c->u) + " x " + iv(c->v);
                        s += "; back " + iv(d0.x) + "," + iv(d0.y) + "; out:";
                        for (std::size_t q = 0; q < cand.size(); ++q) {
                            const HalfEdge& o = he[static_cast<std::size_t>(cand[q])];
                            s += " {" + std::string(o.piece.fresh ? "new" : (o.own ? "own" : "other")) + " " + iv(cd[q].x) + "," + iv(cd[q].y);
                            if (!o.piece.fresh) {
                                const Imprint::OldPiece& op_ = im.pieces[static_cast<std::size_t>(o.piece.src)];
                                const BEdge& oe = S(op_.solid).edges[static_cast<std::size_t>(op_.edge)];
                                for (int ev : {oe.v0, oe.v1}) {
                                    Chart ecc{};
                                    const bool has = chart_in_class(vmap[op_.solid][static_cast<std::size_t>(ev)], g, ecc);
                                    s += " end " + str(vmap[op_.solid][static_cast<std::size_t>(ev)]) + (has ? " " + iv(ecc.u) + "x" + iv(ecc.v) : " no chart");
                                }
                            }
                            s += "}";
                        }
                        note(s);
                    }
                    return fail(name + ": the curves at a vertex could not be ordered (tangential?)");
                }
                note(name + ": curves at vertex " + str(h.to) + " ordered with a curve tangent to a side of the square");
            }
            next[static_cast<std::size_t>(hi)] = cand[best];
        }

        // Cycles of the map. Each becomes a loop of some region.
        struct Cycle {
            std::vector<Piece> pieces;
            CycleGeom geom;
            bool straight = false;          ///< has a piece on a side of the parameter square
            const ArcGeom* probe = nullptr; ///< an arc of the cycle
            int probe_coord = 0;
            const Affine* probe_map = nullptr;
            int probe_he = -1;
            bool have_point = false;
            Interval pu, pv; ///< certified point of the cycle
        };
        std::vector<Cycle> cycles;
        std::vector<int> cycle_of(static_cast<std::size_t>(nH), -1);
        for (int start = 0; start < nH; ++start) {
            if (cycle_of[static_cast<std::size_t>(start)] >= 0) continue;
            Cycle cyc;
            const int index = static_cast<int>(cycles.size());
            int hi = start;
            for (int guard = 0; guard <= nH; ++guard) {
                if (cycle_of[static_cast<std::size_t>(hi)] >= 0) break;
                cycle_of[static_cast<std::size_t>(hi)] = index;
                const HalfEdge& h = he[static_cast<std::size_t>(hi)];
                cyc.pieces.push_back(h.piece);
                const ArcGeom* arc = nullptr;
                int coord = 0, v0 = -1, v1 = -1;
                const Affine* amap = nullptr;
                bool square_side = false;
                if (h.piece.fresh) {
                    const Imprint::NewArc& na = im.arcs[static_cast<std::size_t>(h.piece.id)];
                    arc = na.geom.get();
                    if (!arc_frame(na.carrier[0], na.carrier[1], g, coord, amap))
                        return fail(name + ": a new edge does not lie on this face's surface");
                    v0 = na.v0;
                    v1 = na.v1;
                } else {
                    const Imprint::OldPiece& p = im.pieces[static_cast<std::size_t>(h.piece.src)];
                    const BEdge& e = S(p.solid).edges[static_cast<std::size_t>(p.edge)];
                    int g0 = -1, g1 = -1;
                    edge_carriers(p.solid, p.edge, g0, g1);
                    if (e.arc) {
                        arc = p.geom.get();
                        if (!arc_frame(g0, g1, g, coord, amap))
                            return fail(name + ": an edge does not lie on this face's surface");
                        v0 = p.v0;
                        v1 = p.v1;
                    } else {
                        square_side = h.own && (g0 == g || g1 == g);
                    }
                }
                if (arc) {
                    append_arc_links(cyc.geom.links, *arc, coord, h.fwd, nodes[v0].uv, nodes[v1].uv, amap);
                    if (!cyc.probe) {
                        cyc.probe = arc;
                        cyc.probe_coord = coord;
                        cyc.probe_map = amap;
                        cyc.probe_he = hi;
                        cyc.have_point = false; // an arc point is preferred to a straight piece's
                    }
                } else {
                    if (square_side) cyc.straight = true;
                    const Node& a = nodes[h.from];
                    const Node& b = nodes[h.to];
                    append_straight_link(cyc.geom.links, a.uv, a.box, b.uv, b.box);
                    if (!cyc.probe && cyc.probe_he < 0) {
                        // A certified point of the cycle: the middle of a straight piece.
                        cyc.probe_he = hi;
                        cyc.pu = (a.box[0] + b.box[0]) * Interval(0.5);
                        cyc.pv = (a.box[1] + b.box[1]) * Interval(0.5);
                        cyc.have_point = true;
                    }
                }
                hi = next[static_cast<std::size_t>(hi)];
            }
            if (hi != start) return fail(name + ": a boundary cycle of the planar map does not close");
            cycles.push_back(std::move(cyc));
        }

        // A cycle with a piece on a side of the square is an outer loop: the
        // face cannot extend beyond the square. Any other cycle is a simple
        // closed curve whose orientation says whether it is an outer loop
        // (counter-clockwise) or a hole (clockwise).
        std::vector<int> outers, holes;
        for (std::size_t ci = 0; ci < cycles.size(); ++ci) {
            Cycle& cyc = cycles[ci];
            bool outer = cyc.straight;
            if (!outer) {
                const int o = cycle_orientation(cyc.geom, opt.refine_depth);
                if (o == 0) return fail(name + ": the orientation of a closed loop could not be decided");
                outer = o > 0;
            }
            (outer ? outers : holes).push_back(static_cast<int>(ci));
        }
        auto cycle_point = [&](Cycle& cyc) {
            if (cyc.have_point) return true;
            Box4 enc;
            if (!cyc.probe || !arc_point(*cyc.probe, enc)) return false;
            cyc.pu = enc[static_cast<std::size_t>(cyc.probe_coord)];
            cyc.pv = enc[static_cast<std::size_t>(cyc.probe_coord + 1)];
            if (cyc.probe_map) cyc.probe_map->apply(Interval(cyc.pu), Interval(cyc.pv), cyc.pu, cyc.pv);
            cyc.have_point = true;
            return true;
        };
        // Does outer loop `o` enclose the certified point of cycle `c`?
        auto encloses = [&](int o, Cycle& c, bool& yes) {
            const WindingResult w = winding(cycles[static_cast<std::size_t>(o)].geom, c.pu, c.pv, opt.refine_depth);
            yes = w.decided && w.wn == 1;
            return w.decided;
        };
        std::vector<std::vector<int>> holes_of(cycles.size());
        for (int hc : holes) {
            Cycle& hole = cycles[static_cast<std::size_t>(hc)];
            if (!cycle_point(hole)) return fail(name + ": no certified point could be placed on a hole");
            // The cell on the other side of the probed arc lies inside the hole.
            const int tw = he[static_cast<std::size_t>(hole.probe_he)].twin;
            const int twin = tw >= 0 ? cycle_of[static_cast<std::size_t>(tw)] : -1;
            std::vector<int> containers;
            for (int o : outers) {
                if (o == twin) continue;
                bool in = false;
                if (!encloses(o, hole, in)) return fail(name + ": could not decide which region contains a closed loop");
                if (in) containers.push_back(o);
            }
            // Outer loops that enclose the same point are nested and disjoint;
            // the innermost is enclosed by all the others.
            int best = -1, best_depth = -1;
            for (int o : containers) {
                int depth = 0;
                Cycle& oc = cycles[static_cast<std::size_t>(o)];
                if (oc.probe || oc.have_point) {
                    if (!cycle_point(oc)) return fail(name + ": no certified point could be placed on a loop");
                    for (int p : containers) {
                        if (p == o) continue;
                        bool nested = false;
                        if (!encloses(p, oc, nested)) return fail(name + ": could not decide how closed loops are nested");
                        if (nested) ++depth;
                    }
                }
                if (depth > best_depth) { best_depth = depth; best = o; }
            }
            if (best < 0) return fail(name + ": a closed loop lies in no region");
            holes_of[static_cast<std::size_t>(best)].push_back(hc);
        }
        for (int o : outers) {
            Region reg;
            reg.solid = k;
            reg.face = f;
            reg.loops.push_back(cycles[static_cast<std::size_t>(o)].pieces);
            reg.geom.push_back(cycles[static_cast<std::size_t>(o)].geom);
            for (int hc : holes_of[static_cast<std::size_t>(o)]) {
                reg.loops.push_back(cycles[static_cast<std::size_t>(hc)].pieces);
                reg.geom.push_back(cycles[static_cast<std::size_t>(hc)].geom);
            }
            im.regions.push_back(std::move(reg));
        }
        // The region that contains each line contact.
        for (int pid : contact_pieces) {
            const Imprint::OldPiece& p = im.pieces[static_cast<std::size_t>(pid)];
            const Chart* a = chart(p.v0, g);
            const Chart* b = chart(p.v1, g);
            if (!a || !b) return fail(name + ": a vertex of a line contact has no parameters on this face");
            int in = -1;
            for (std::size_t r = im.regions.size() - outers.size(); r < im.regions.size() && in < 0; ++r)
                for (double s : {0.5, 0.25, 0.75}) {
                    const Interval is(s);
                    if (locate(im.regions[r].geom, a->u + is * (b->u - a->u), a->v + is * (b->v - a->v), opt.refine_depth) != Where::Inside)
                        continue;
                    in = static_cast<int>(r);
                    break;
                }
            if (in < 0) return fail(name + ": could not place a line contact in a region of the face");
            im.line_contacts.push_back({in, p.canon});
            ++im.stats.line_contacts;
        }
        return true;
    }

    // --- 6. inside / outside / on ----------------------------------------------------
    /// Parity of the crossings of a ray from P with solid k.
    bool point_in_solid(int k, const Vec3d& P, bool& inside) {
        double reach = 0.0;
        for (const PolyPatchD& f : S(k).carriers)
            for (const Vec4d& c : f.c) reach = std::max(reach, geom::norm(Vec3d{c.x / c.w, c.y / c.w, c.z / c.w} - P));
        const double L = 2.0 * reach + 1.0;
        static constexpr double dirs[5][3] = {{0.5377, 0.8339, 0.1258},
                                              {-0.2769, 0.3319, 0.9018},
                                              {0.7254, -0.6131, 0.3127},
                                              {-0.6490, -0.3035, 0.6976},
                                              {0.1143, -0.5571, -0.8225}};
        for (const auto& d : dirs) {
            const PolyPatchD ray = solve::make_line(P, P + Vec3d{d[0], d[1], d[2]} * L);
            int count = 0;
            bool valid = true;
            for (int g = 0; g < static_cast<int>(im.carriers.size()) && valid; ++g) {
                if (faces_on[k][g].empty()) continue;
                const solve::CsiSystem cs(ray, im.carriers[static_cast<std::size_t>(g)]);
                const solve::CsiResult res = solve::intersect_curve_surface(cs, opt.csi);
                if (!res.complete) valid = false;
                for (const solve::CsiRoot& r : res.roots) {
                    int face = -1;
                    const Hit hit = find_face(k, g, r.enc[1], r.enc[2], face);
                    if (hit == Hit::None) continue; // e.g. the start point, on a surface but outside its faces
                    if (hit == Hit::Undecided || !(r.enc[0].lo > 0.0 && r.enc[0].hi < 1.0)) valid = false;
                    else ++count;
                }
            }
            if (valid) { inside = count % 2 == 1; return true; }
        }
        return false;
    }

    bool classify() {
        const int nR = static_cast<int>(im.regions.size());
        // Regions of the two solids that touch along a line inside both: a
        // point of the line proven inside each of them. (Recorded for the
        // result; classification itself needs nothing here, since the curved
        // face lies on one side of the plane.)
        for (const LineTouch& t : line_touches)
            for (int r1 = 0; r1 < nR; ++r1) {
                Region& R1 = im.regions[static_cast<std::size_t>(r1)];
                if (gid(R1.solid, S(R1.solid).faces[static_cast<std::size_t>(R1.face)].carrier) != t.a) continue;
                for (int r2 = 0; r2 < nR; ++r2) {
                    Region& R2 = im.regions[static_cast<std::size_t>(r2)];
                    if (R2.solid == R1.solid) continue;
                    const int gb = gid(R2.solid, S(R2.solid).faces[static_cast<std::size_t>(R2.face)].carrier);
                    if (!equiv(gb, t.b)) continue;
                    bool found = false;
                    for (double sv : {0.5, 0.25, 0.75, 0.125, 0.875, 0.375, 0.625}) {
                        const Interval ua = t.fixed_u ? t.r : Interval(sv), va = t.fixed_u ? Interval(sv) : t.r;
                        if (locate(R1.geom, ua, va, opt.refine_depth) != Where::Inside) continue;
                        const auto h = icarrier[static_cast<std::size_t>(t.a)].eval(ua, va);
                        Interval pu, pv;
                        if (!plane_params(im.carriers[static_cast<std::size_t>(gb)], {h.x / h.w, h.y / h.w, h.z / h.w}, pu, pv)) continue;
                        if (locate(R2.geom, pu, pv, opt.refine_depth) != Where::Inside) continue;
                        found = true;
                        break;
                    }
                    if (found) im.touching_regions.push_back({r1, r2});
                }
            }
        // Regions of the two solids that touch at one point inside both. A
        // point on a face's boundary would need a vertex there: refused.
        for (const PointTouch& t : point_touches)
            for (int r1 = 0; r1 < nR; ++r1) {
                Region& R1 = im.regions[static_cast<std::size_t>(r1)];
                if (gid(R1.solid, S(R1.solid).faces[static_cast<std::size_t>(R1.face)].carrier) != t.a) continue;
                const Where w1 = locate(R1.geom, Interval(t.u), Interval(t.v), opt.refine_depth);
                if (w1 == Where::Outside) continue;
                if (w1 != Where::Inside)
                    return fail(std::string("a curved face touches a surface at one point ") +
                                (w1 == Where::Boundary ? "on the boundary of the face" : "not located in the face"));
                for (int r2 = 0; r2 < nR; ++r2) {
                    Region& R2 = im.regions[static_cast<std::size_t>(r2)];
                    if (R2.solid == R1.solid) continue;
                    const int gb = gid(R2.solid, S(R2.solid).faces[static_cast<std::size_t>(R2.face)].carrier);
                    if (!equiv(gb, t.b)) continue;
                    Interval pu(t.s), pv(t.t);
                    if (!t.b_curved && !plane_params(im.carriers[static_cast<std::size_t>(gb)], t.P, pu, pv))
                        return fail("a point where a curved face touches a plane could not be placed in the plane's face");
                    Where w2 = locate(R2.geom, pu, pv, opt.refine_depth);
                    // Exactly on the face's boundary: interval point location cannot
                    // say so, but the structure does when the point is a vertex of
                    // the imprint with parameters on the plane (made where an edge
                    // of the plane's solid touches the curved surface).
                    if (w2 == Where::Undecided && t.has_xp && !t.b_curved) {
                        using E = num::Expansion;
                        auto same = [&](double a, double aw, double b, double bw) { return (E(a) * E(bw) - E(b) * E(aw)).sign() == 0; };
                        for (int vi = 0; vi < static_cast<int>(im.vertices.size()); ++vi) {
                            const auto& V = im.vertices[static_cast<std::size_t>(vi)];
                            if (!V.exact || !chart(vi, gb)) continue;
                            if (same(V.xp.x, V.xp.w, t.xp.x, t.xp.w) && same(V.xp.y, V.xp.w, t.xp.y, t.xp.w) &&
                                same(V.xp.z, V.xp.w, t.xp.z, t.xp.w))
                                w2 = Where::Boundary;
                        }
                    }
                    if (w2 == Where::Outside) continue;
                    // On the boundary of a flat face: the point lies on an edge of
                    // the plane's solid (that edge has a vertex there, where it
                    // touches the curved face); the regions still touch there.
                    if (w2 == Where::Boundary && !t.b_curved && t.alone) {
                        im.point_touching_regions.push_back({r1, r2});
                        continue;
                    }
                    if (w2 != Where::Inside)
                        return fail(std::string("two faces touch at one point ") +
                                    (w2 == Where::Boundary ? "on the boundary of the second" : "not located in the second"));
                    im.point_touching_regions.push_back({r1, r2});
                }
            }
        // On the other solid's boundary: a region of a shared surface inside a
        // face of the other solid there.
        for (Region& reg : im.regions) {
            const BFace& face = S(reg.solid).faces[static_cast<std::size_t>(reg.face)];
            const int g = gid(reg.solid, face.carrier);
            if (!faces_in_plane(1 - reg.solid, g)) continue;
            double u = 0.0, v = 0.0;
            if (!interior_point(reg.geom, opt.refine_depth, u, v))
                return fail("no interior point could be found for a region of a shared surface");
            int other = -1;
            bool reversed = false;
            const Hit hit = find_face_in_plane(1 - reg.solid, g, Interval(u), Interval(v), other, reversed);
            if (hit == Hit::Undecided) {
                if (opt.trace) {
                    note("region of face " + str(reg.face) + " of solid " + str(reg.solid) + " on surface " + str(g) +
                         ": its point (" + std::to_string(u) + ", " + std::to_string(v) + ") is undecided against the other solid");
                    for (const auto& loop : reg.loops) {
                        std::string s = "  loop:";
                        for (const Piece& p : loop) {
                            if (p.fresh) { s += " new" + str(p.id); continue; }
                            const Imprint::OldPiece& op_ = im.pieces[static_cast<std::size_t>(p.src)];
                            s += std::string(" ") + (op_.solid == reg.solid ? "own" : "other") + str(p.src) + "(" +
                                 str(op_.v0) + "-" + str(op_.v1) + ")";
                        }
                        note(s);
                    }
                }
                return fail("could not decide whether a region lies on the other solid's boundary");
            }
            if (hit == Hit::Found) {
                // Outward normals agree unless exactly one of: this face flipped, the
                // other face flipped, the other surface parameterised the other way.
                const bool other_flipped = S(1 - reg.solid).faces[static_cast<std::size_t>(other)].flipped;
                reg.inside = (other_flipped != face.flipped) != reversed ? kOnOpposite : kOnSame;
                ++im.stats.on_regions;
            }
        }
        // Next to a new edge: decided on the curve.
        std::vector<int> queue;
        for (int ri = 0; ri < nR; ++ri) {
            Region& reg = im.regions[static_cast<std::size_t>(ri)];
            for (const auto& loop : reg.loops)
                for (const Piece& p : loop) {
                    if (!p.fresh) continue;
                    const Imprint::NewArc& a = im.arcs[static_cast<std::size_t>(p.id)];
                    const int i = a.solid_on[0] == reg.solid ? 0 : 1;
                    const int in = a.left_inside[i] == p.forward ? kIn : kOut;
                    if (reg.inside != kUnknown && reg.inside != in)
                        return fail("the curves around one region disagree about inside and outside");
                    reg.inside = in;
                }
            if (reg.inside == kIn || reg.inside == kOut) queue.push_back(ri);
        }
        // Across uncut pieces of the region's own edges: a piece that does not
        // lie on the other solid's boundary does not cross it. A piece can lie
        // on it without being shared with it: it may run inside a face of the
        // other solid on a shared surface, with both its neighbouring regions
        // on a third surface (a seam between two faces on one surface).
        std::map<int, std::vector<int>> uses; // own piece -> regions of its solid that use it
        for (int ri = 0; ri < nR; ++ri)
            for (const auto& loop : im.regions[static_cast<std::size_t>(ri)].loops)
                for (const Piece& p : loop)
                    if (!p.fresh && im.pieces[static_cast<std::size_t>(p.src)].solid == im.regions[static_cast<std::size_t>(ri)].solid)
                        uses[p.src].push_back(ri);
        // Two regions per piece, or an even number at an edge where parts of
        // the solid touch. A piece that is not on the other solid's boundary
        // has every region around it on the same side of that solid, however
        // many there are.
        for (const auto& kv : uses)
            if (kv.second.size() < 2 || kv.second.size() % 2 != 0)
                return fail("an edge piece is not shared by an even number of regions of its solid");
        std::map<int, bool> crossable;
        for (const auto& kv : uses) {
            const Imprint::OldPiece& piece = im.pieces[static_cast<std::size_t>(kv.first)];
            bool ok = !piece.coincident;
            int g0 = -1, g1 = -1;
            edge_carriers(piece.solid, piece.edge, g0, g1);
            const BEdge& e = S(piece.solid).edges[static_cast<std::size_t>(piece.edge)];
            // A straight edge lying exactly in a plane of the other solid's
            // faces may be on its boundary: not crossed; its regions are
            // decided otherwise. (Where a curved surface touches that plane, or
            // where the edge lies in a face of the other solid and its own two
            // faces leave it on the two sides of that face, the edge leaving a
            // box through one of its edges, say.)
            if (ok && !e.arc && side_in_face_plane(1 - piece.solid, g0, e.on[0].side, g1)) ok = false;
            // The same for an edge lying in a curved surface of the other solid.
            if (ok && iso_lines.count({piece.solid, piece.edge})) ok = false;
            for (int g : {g0, g1}) {
                if (!ok || !faces_in_plane(1 - piece.solid, g)) continue;
                Interval u, v;
                if (e.arc) {
                    Box4 mid;
                    if (!arc_point(*piece.geom, mid)) return fail("no certified point could be placed on an edge piece");
                    const std::size_t co = g == g0 ? 0 : 2;
                    u = mid[co];
                    v = mid[co + 1];
                } else {
                    const BEdge::SideRef& ref = g == g0 ? e.on[0] : e.on[1];
                    const double t = 0.5 * (num::mid(piece.t0) + num::mid(piece.t1));
                    side_uv(ref.side, ref.reversed ? Interval(1.0) - Interval(t) : Interval(t), u, v);
                }
                int face = -1;
                bool rev = false;
                if (find_face_in_plane(1 - piece.solid, g, u, v, face, rev) != Hit::None) ok = false;
            }
            crossable[kv.first] = ok;
        }
        auto flood = [&]() {
            while (!queue.empty()) {
                const int ri = queue.back();
                queue.pop_back();
                const Region& reg = im.regions[static_cast<std::size_t>(ri)];
                for (const auto& loop : reg.loops)
                    for (const Piece& p : loop) {
                        if (p.fresh) continue;
                        const Imprint::OldPiece& piece = im.pieces[static_cast<std::size_t>(p.src)];
                        if (piece.solid != reg.solid || !crossable[p.src]) continue;
                        for (int o : uses[p.src]) {
                            Region& other = im.regions[static_cast<std::size_t>(o)];
                            if (o == ri || other.inside >= kOnSame) continue;
                            if (other.inside == kUnknown) {
                                other.inside = reg.inside;
                                queue.push_back(o);
                            } else if (other.inside != reg.inside) {
                                if (opt.trace) {
                                    int g0 = -1, g1 = -1;
                                    edge_carriers(piece.solid, piece.edge, g0, g1);
                                    note("regions " + str(ri) + " (solid " + str(reg.solid) + " face " + str(reg.face) + ", " +
                                         (reg.inside == kIn ? "in" : "out") + ") and " + str(o) + " (face " + str(other.face) + ", " +
                                         (other.inside == kIn ? "in" : "out") + ") across piece " + str(p.src) + " of edge " +
                                         str(piece.edge) + " (" + (S(piece.solid).edges[static_cast<std::size_t>(piece.edge)].arc ? "arc" : "side") +
                                         " on " + sk(g0) + ", " + sk(g1) + "), vertices " + str(piece.v0) + " -> " + str(piece.v1));
                                }
                                return false;
                            }
                        }
                    }
            }
            return true;
        };
        const char* clash = "two regions that share an uncut edge piece disagree about inside and outside";
        if (!flood()) return fail(clash);
        // What is left: a certified ray test from a point inside the region.
        for (int ri = 0; ri < nR; ++ri) {
            Region& reg = im.regions[static_cast<std::size_t>(ri)];
            if (reg.inside != kUnknown) continue;
            // A point of the region that lies on the other solid's boundary
            // (on a line where the solids touch inside the face, say) cannot be
            // classified by a ray; another point of the region is tried.
            bool in = false, classified = false;
            for (int attempt = 0; attempt < 6 && !classified; ++attempt) {
                double u = 0.0, v = 0.0;
                if (!interior_point(reg.geom, opt.refine_depth, u, v, attempt)) {
                    if (attempt == 0) return fail("no interior point could be found for a region");
                    break;
                }
                const Vec4d h = im.carriers[static_cast<std::size_t>(gid(reg.solid, S(reg.solid).faces[static_cast<std::size_t>(reg.face)].carrier))].eval(u, v);
                classified = point_in_solid(1 - reg.solid, {h.x / h.w, h.y / h.w, h.z / h.w}, in);
            }
            if (!classified) return fail("could not classify a region by a ray test");
            ++im.stats.ray_tests;
            reg.inside = in ? kIn : kOut;
            queue.push_back(ri);
            if (!flood()) return fail(clash);
        }
        return true;
    }
};

} // namespace detail

/// Intersect the boundaries of two solids and split every face into regions
/// classified against the other solid. A and B must outlive the result.
// --- flat faces on rational carriers ----------------------------------------------
//
// A flat face can lie on a carrier that is not an affine parallelogram: the
// annular cap of a cylinder sector is rational in the angle. Against another
// flat face in the same plane the overlay needs an affine map between the two
// parameter planes, and there is none. Such a face is moved, before the
// Boolean, onto an exact parallelogram of its plane: the same face, the same
// edges and vertices, on another carrier. Its curved boundary edges become
// arc edges, the certified (transversal) intersection of the neighbouring
// surface with the new plane; its straight edges stay side edges, of a thin
// exact parallelogram built on the same side, so that they lie on a surface
// of the plane's class. Only operands where the conflict occurs are changed,
// and only copies of them.
namespace flat_detail {

/// Pieces of a re-certified curve: not a power of two, so that its joints
/// avoid the dyadic parameters where vertices usually sit (a piece starting
/// exactly at a joint would begin with an empty span).
constexpr int kRecarryPieces = 7;

/// The plane x_axis = h of a patch, exactly, for an axis-aligned flat patch.
inline bool axis_plane(const PolyPatchD& P, int& axis, double& h) {
    for (axis = 0; axis < 3; ++axis) {
        auto coord = [&](const Vec4d& c) { return axis == 0 ? c.x : (axis == 1 ? c.y : c.z); };
        const Vec4d& c0 = P.c.front();
        if (!(c0.w > 0.0)) return false;
        h = coord(c0) / c0.w;
        bool all = true;
        for (const Vec4d& c : P.c)
            all = all && c.w > 0.0 && (num::Expansion(coord(c)) - num::Expansion(h) * num::Expansion(c.w)).sign() == 0;
        if (all) return true;
    }
    return false;
}

inline Vec4d axis_point(int axis, double h, int a, double xa, int b, double xb) {
    double p[3] = {0.0, 0.0, 0.0};
    p[axis] = h;
    p[a] = xa;
    p[b] = xb;
    return {p[0], p[1], p[2], 1.0};
}

/// a + b, if exact.
inline bool exact_sum(double a, double b, double& s) {
    s = a + b;
    return (num::Expansion(s) - num::Expansion(a) - num::Expansion(b)).sign() == 0;
}

/// A thin parallelogram, uniform weight, whose side 0 is exactly the given
/// straight side (control points in order) and whose other side is that side
/// shifted by a dyadic step along an in-plane axis. False if the side's
/// weights differ or a shifted coordinate is not exact.
inline bool side_parallelogram(const std::vector<Vec4d>& side, int axis, PolyPatchD& Q) {
    const double w = side.front().w;
    for (const Vec4d& c : side)
        if (c.w != w) return false;
    const int n = static_cast<int>(side.size()) - 1;
    const Vec4d d = {side.back().x - side.front().x, side.back().y - side.front().y, side.back().z - side.front().z, 0.0};
    for (int dir = 0; dir < 3; ++dir) {
        if (dir == axis) continue;
        // Not along the side.
        const double along = dir == 0 ? d.x : (dir == 1 ? d.y : d.z);
        const double len = std::fabs(d.x) + std::fabs(d.y) + std::fabs(d.z);
        if (std::fabs(along) > 0.99 * len) continue;
        Q = PolyPatchD(n, 1);
        bool ok = true;
        for (int i = 0; i <= n && ok; ++i) {
            const Vec4d& c = side[static_cast<std::size_t>(i)];
            Vec4d s = c;
            double* coord = dir == 0 ? &s.x : (dir == 1 ? &s.y : &s.z);
            ok = exact_sum(*coord, 0.125 * w, *coord);
            Q.at(i, 0) = c;
            Q.at(i, 1) = s;
        }
        if (ok && planar_affine(Q)) return true;
    }
    return false;
}

/// Move the faces on flat carrier g of b onto an exact parallelogram of its
/// plane. Works on b in place; the caller keeps the original on failure.
inline bool recarry(Brep& b, int g, const solve::SsiOptions& sopt) {
    const PolyPatchD G = b.carriers[static_cast<std::size_t>(g)];
    int axis = -1;
    double h = 0.0;
    if (!axis_plane(G, axis, h)) return false;
    // The two in-plane axes, ordered so that the new normal has the
    // orientation of G's (Su x Sv), then a dyadic box around G's control points.
    const geom::Jet3 J = geom::BezierSurface(G).jet(0.5, 0.5);
    const Vec3d nrm = cross(J[static_cast<std::size_t>(geom::jet_index(1, 0))], J[static_cast<std::size_t>(geom::jet_index(0, 1))]);
    const double nax = axis == 0 ? nrm.x : (axis == 1 ? nrm.y : nrm.z);
    if (!(std::fabs(nax) > 0.0)) return false;
    int a = (axis + 1) % 3, c = (axis + 2) % 3; // e_a x e_c = +e_axis
    if (nax < 0.0) std::swap(a, c);
    double lo[3] = {1e300, 1e300, 1e300}, hi[3] = {-1e300, -1e300, -1e300};
    for (const Vec4d& p : G.c) {
        const double q[3] = {p.x / p.w, p.y / p.w, p.z / p.w};
        for (int i = 0; i < 3; ++i) { lo[i] = std::min(lo[i], q[i]); hi[i] = std::max(hi[i], q[i]); }
    }
    // Corners on a dyadic grid and a power-of-two extent, so that a vertex at
    // an exact dyadic point has exact parameters, (x - lo) 2^-k.
    const double grid = 16.0;
    double extent = 0.0;
    for (int i : {a, c}) {
        lo[i] = std::floor(lo[i] * grid) / grid - 1.0 / grid;
        extent = std::max(extent, std::ceil(hi[i] * grid) / grid + 1.0 / grid - lo[i]);
    }
    const double size = std::ldexp(1.0, static_cast<int>(std::ceil(std::log2(extent))));
    for (int i : {a, c}) hi[i] = lo[i] + size;
    PolyPatchD P(1, 1);
    P.at(0, 0) = axis_point(axis, h, a, lo[a], c, lo[c]);
    P.at(1, 0) = axis_point(axis, h, a, hi[a], c, lo[c]);
    P.at(0, 1) = axis_point(axis, h, a, lo[a], c, hi[c]);
    P.at(1, 1) = axis_point(axis, h, a, hi[a], c, hi[c]);
    if (!planar_affine(P)) return false;
    const int pg = static_cast<int>(b.carriers.size());
    b.carriers.push_back(P);
    const PolyPatchI Gi = geom::to_interval(G);

    // Faces and the edges they use.
    std::vector<int> edges;
    for (BFace& f : b.faces) {
        if (f.carrier != g) continue;
        f.carrier = pg;
        for (const auto& loop : f.loops)
            for (const BCoedge& ce : loop) edges.push_back(ce.edge);
    }
    std::sort(edges.begin(), edges.end());
    edges.erase(std::unique(edges.begin(), edges.end()), edges.end());

    // Vertices: parameters on the new plane, enclosing the exact ones.
    for (BVertex& v : b.vertices) {
        const Chart* cg = v.chart(g);
        if (!cg || v.chart(pg)) continue;
        const auto hp = Gi.eval(cg->u, cg->v);
        const Vec3<Interval> X{hp.x / hp.w, hp.y / hp.w, hp.z / hp.w};
        Chart m;
        m.carrier = pg;
        if (!plane_params(P, X, m.u, m.v)) return false;
        // A corner of G with a power-of-two weight is an exact point: its
        // parameters on the plane are exact where the arithmetic is.
        const bool corner = is_point(cg->u) && is_point(cg->v) && (cg->u.lo == 0.0 || cg->u.lo == 1.0) &&
                            (cg->v.lo == 0.0 || cg->v.lo == 1.0);
        if (corner) {
            const Vec4d& h4 = G.at(cg->u.lo == 1.0 ? G.p : 0, cg->v.lo == 1.0 ? G.q : 0);
            int ex = 0;
            if (h4.w > 0.0 && std::frexp(h4.w, &ex) == 0.5) {
                const double q[3] = {h4.x / h4.w, h4.y / h4.w, h4.z / h4.w}; // exact: division by a power of two
                double du = 0.0, dv = 0.0;
                if (exact_sum(q[a], -lo[a], du) && exact_sum(q[c], -lo[c], dv)) {
                    const double pu = du / size, pv = dv / size; // exact: division by a power of two
                    if (num::contains(m.u, pu) && num::contains(m.v, pv)) {
                        m.u = Interval(pu);
                        m.v = Interval(pv);
                    }
                }
            }
        }
        v.charts.push_back(m);
    }

    for (int ei : edges) {
        BEdge& e = b.edges[static_cast<std::size_t>(ei)];
        if (e.arc) {
            // The arc's carrier in the plane: g itself, or another flat
            // carrier of the same plane that holds no face (an earlier
            // operation left the arc there; the face is on g).
            auto in_plane = [&](int cc) {
                return cc == g || (cc != pg && !planar_affine(b.carriers[static_cast<std::size_t>(cc)]) &&
                                   same_plane_patches(b.carriers[static_cast<std::size_t>(cc)], G));
            };
            const int kk = in_plane(e.carrier[0]) ? 0 : (in_plane(e.carrier[1]) ? 1 : -1);
            if (kk < 0) continue;
            const int oc = e.carrier[kk];
            const PolyPatchI Oi = oc == g ? Gi : geom::to_interval(b.carriers[static_cast<std::size_t>(oc)]);
            // An intersection curve of that carrier with another surface m,
            // made by an earlier operation: certified again as the curve of
            // the plane and m, through points of the old arc.
            if (e.v0 == e.v1) return false;
            for (int end = 0; end < 2 && oc != g; ++end) {
                BVertex& v = b.vertices[static_cast<std::size_t>(end == 0 ? e.v0 : e.v1)];
                const Chart* co = v.chart(oc);
                if (v.chart(pg) || !co) continue;
                const auto hp = Oi.eval(co->u, co->v);
                Chart mc;
                mc.carrier = pg;
                if (!plane_params(P, {hp.x / hp.w, hp.y / hp.w, hp.z / hp.w}, mc.u, mc.v)) return false;
                v.charts.push_back(mc);
            }
            const ArcGeom& old = *e.geom;
            const SsiSystem& osys = *old.base->sys;
            const auto& osegs = old.base->curve->segs;
            const int m = e.carrier[1 - kk];
            const PolyPatchD& M = b.carriers[static_cast<std::size_t>(m)];
            auto sys = kk == 0 ? std::make_shared<SsiSystem>(P, M) : std::make_shared<SsiSystem>(M, P);
            const std::size_t gc = static_cast<std::size_t>(2 * kk), mc = static_cast<std::size_t>(2 * (1 - kk));
            const int ns = static_cast<int>(old.spans.size());
            SsiCurve curve;
            const bool ok = detail::known_curve(*sys, [&](double tau) {
                const double f = std::min(tau * ns, static_cast<double>(ns) - 1e-12);
                const int j = std::min(static_cast<int>(f), ns - 1);
                const Span& sp = old.spans[static_cast<std::size_t>(j)];
                const TubeSegment& seg = osegs[static_cast<std::size_t>(sp.seg)];
                const double c = tau >= 1.0 ? sp.cb : sp.ca + (f - j) * (sp.cb - sp.ca);
                Pt4 x;
                for (int i = 0; i < 4; ++i) x[i] = seg.a[i] + (c / seg.len) * (seg.b[i] - seg.a[i]);
                solve::newton_slice(osys, seg.m, seg.a, c, x);
                const auto hp = Oi.eval(Interval(x[gc]), Interval(x[gc + 1]));
                const Vec3<Interval> X{hp.x / hp.w, hp.y / hp.w, hp.z / hp.w};
                Interval pu(0.5), pv(0.5);
                plane_params(P, X, pu, pv);
                Pt4 y;
                y[gc] = num::mid(pu);
                y[gc + 1] = num::mid(pv);
                y[mc] = x[mc];
                y[mc + 1] = x[mc + 1];
                return y;
            }, sopt, curve, kRecarryPieces);
            if (!ok) return false;
            curve.known = false;
            auto base = make_curve_base(sys, std::make_shared<const SsiCurve>(std::move(curve)));
            CurvePos pos[2];
            for (int end = 0; end < 2; ++end) {
                const BVertex& v = b.vertices[static_cast<std::size_t>(end == 0 ? e.v0 : e.v1)];
                const Chart* cp = v.chart(pg);
                const Chart* cm = v.chart(m);
                if (!cp || !cm) return false;
                Box4 X;
                X[gc] = cp->u; X[gc + 1] = cp->v;
                X[mc] = cm->u; X[mc + 1] = cm->v;
                bool on = false;
                if (!locate_on_curve(*base, X, on, pos[end]) || !on) return false;
            }
            if (before(pos[0], pos[1]) != Tri::Yes) return false;
            e.carrier[kk] = pg;
            e.geom = make_arc(base, pos[0], pos[1], false);
            continue;
        }
        const int k = e.on[0].carrier == g ? 0 : (e.on[1].carrier == g ? 1 : -1);
        if (k < 0) continue; // a straight edge of another carrier lying in the plane
        const std::vector<Vec4d> S = side_control(G, e.on[k].side);
        bool straight = true;
        for (std::size_t i = 1; i + 1 < S.size(); ++i) straight = straight && collinear(S.front(), S[i], S.back());
        if (straight) {
            PolyPatchD Q;
            if (!side_parallelogram(S, axis, Q)) return false;
            const int q = static_cast<int>(b.carriers.size());
            b.carriers.push_back(Q);
            for (int end = 0; end < 2; ++end) {
                BVertex& v = b.vertices[static_cast<std::size_t>(end == 0 ? e.v0 : e.v1)];
                const Chart* cg = v.chart(g);
                if (!cg) return false;
                Interval tau;
                if (detail::on_side(*cg, e.on[k].side, tau) != Tri::Yes) return false;
                v.charts.push_back({q, tau, Interval(0.0)});
            }
            e.on[k].carrier = q;
            e.on[k].side = 0;
            continue;
        }
        // A curved side: the transversal intersection of the neighbouring
        // surface n with the plane, certified along n's side.
        const int o = 1 - k, n = e.on[o].carrier, ns = e.on[o].side;
        const PolyPatchD& N = b.carriers[static_cast<std::size_t>(n)];
        auto sys = std::make_shared<SsiSystem>(N, P);
        const PolyPatchI Ni = geom::to_interval(N);
        SsiCurve curve;
        const bool ok = detail::known_curve(*sys, [&](double tau) {
            double un, vn;
            side_uv(ns, tau, un, vn);
            const auto hp = Ni.eval(Interval(un), Interval(vn));
            const Vec3<Interval> X{hp.x / hp.w, hp.y / hp.w, hp.z / hp.w};
            Interval pu(0.5), pv(0.5);
            plane_params(P, X, pu, pv);
            return Pt4{un, vn, num::mid(pu), num::mid(pv)};
        }, sopt, curve, kRecarryPieces);
        if (!ok) return false;
        curve.known = false;
        auto base = make_curve_base(sys, std::make_shared<const SsiCurve>(std::move(curve)));
        CurvePos pos[2];
        for (int end = 0; end < 2; ++end) {
            const BVertex& v = b.vertices[static_cast<std::size_t>(end == 0 ? e.v0 : e.v1)];
            const Chart* cn = v.chart(n);
            const Chart* cp = v.chart(pg);
            if (!cn || !cp) return false;
            bool on = false;
            if (!locate_on_curve(*base, {cn->u, cn->v, cp->u, cp->v}, on, pos[end]) || !on) return false;
        }
        const Tri order = before(pos[0], pos[1]);
        if (order == Tri::Unknown) return false;
        if (order == Tri::No) {
            // The curve runs from v1 to v0: turn the edge round.
            std::swap(e.v0, e.v1);
            std::swap(pos[0], pos[1]);
            for (BFace& f : b.faces)
                for (auto& loop : f.loops)
                    for (BCoedge& ce : loop)
                        if (ce.edge == ei) ce.reversed = !ce.reversed;
        }
        e.arc = true;
        e.carrier[0] = n;
        e.carrier[1] = pg;
        e.geom = make_arc(base, pos[0], pos[1], false);
    }
    return true;
}

/// Flat carriers with faces that are not affine parallelograms.
inline bool flat_rational(const PolyPatchD& P) { return !planar_affine(P) && same_plane_patches(P, P); }

/// For both operands: the flat rational carriers to move, those in a plane
/// where a face of either solid lies on a carrier they are not affinely
/// related to (the overlay could not map one onto the other).
inline void conflicts(const Brep* s[2], std::vector<int> out[2]) {
    struct Item { int k, g; };
    std::vector<Item> flats, all;
    for (int k = 0; k < 2; ++k) {
        std::vector<char> used(s[k]->carriers.size(), 0);
        for (const BFace& f : s[k]->faces) used[static_cast<std::size_t>(f.carrier)] = 1;
        for (int g = 0; g < static_cast<int>(s[k]->carriers.size()); ++g) {
            if (!used[static_cast<std::size_t>(g)]) continue;
            all.push_back({k, g});
            if (flat_rational(s[k]->carriers[static_cast<std::size_t>(g)])) flats.push_back({k, g});
        }
    }
    for (const Item& f : flats) {
        const PolyPatchD& F = s[f.k]->carriers[static_cast<std::size_t>(f.g)];
        bool conflict = false, other_in_plane = false;
        for (const Item& o : all) {
            if (o.k == f.k) continue;
            const PolyPatchD& O = s[o.k]->carriers[static_cast<std::size_t>(o.g)];
            if (!same_plane_patches(O, F)) continue;
            other_in_plane = true;
            Affine m;
            if (!same_surface_map(F, O, m)) { conflict = true; break; }
        }
        // No face of the other solid in the plane, but a straight edge of it
        // lying there (a hill touching the plane along its base): the edge
        // would enter the cap's faces, where the rational parameters bend it.
        // (With a related face in the plane, the overlay maps that face's
        // edges, as before.) An edge along a side that this solid has too, by
        // construction, does not count: a result against its own operand
        // keeps that operand's surfaces, and their sides in the plane are the
        // cap's own boundary.
        Vec4d span[3];
        if (!conflict && !other_in_plane && flat_triple(F, span)) {
            const Brep& other = *s[1 - f.k];
            for (const BEdge& e : other.edges) {
                if (e.arc) continue;
                const PolyPatchD& X = other.carriers[static_cast<std::size_t>(e.on[0].carrier)];
                bool in_plane = true;
                for (const Vec4d& c : side_control(X, e.on[0].side))
                    in_plane = in_plane && c.w > 0.0 && planar_detail::coplanar4(span[0], span[1], span[2], c);
                if (!in_plane) continue;
                bool own = false;
                for (const PolyPatchD& C : s[f.k]->carriers)
                    for (int t = 0; t < 4 && !own; ++t) {
                        bool reversed = false;
                        own = detail::shared_side(X, e.on[0].side, C, t, reversed);
                    }
                if (!own) { conflict = true; break; }
            }
        }
        if (!conflict) continue;
        // Every flat rational carrier of either solid in that plane moves.
        for (const Item& x : flats)
            if (same_plane_patches(s[x.k]->carriers[static_cast<std::size_t>(x.g)], F)) out[x.k].push_back(x.g);
    }
    for (int k = 0; k < 2; ++k) {
        std::sort(out[k].begin(), out[k].end());
        out[k].erase(std::unique(out[k].begin(), out[k].end()), out[k].end());
    }
}

} // namespace flat_detail

inline Imprint imprint(const Brep& A, const Brep& B, const BooleanOptions& opt = {}) {
    Imprint im;
    im.solid[0] = &A;
    im.solid[1] = &B;
    if (opt.recarry_flat) {
        std::vector<int> move[2];
        flat_detail::conflicts(im.solid, move);
        for (int k = 0; k < 2; ++k) {
            if (move[k].empty()) continue;
            auto copy = std::make_shared<Brep>(*im.solid[k]);
            bool ok = true;
            for (int g : move[k]) ok = ok && flat_detail::recarry(*copy, g, opt.ssi);
            if (opt.trace) {
                std::string s = "solid " + std::to_string(k) + ": flat carriers to move onto planes:";
                for (int g : move[k]) s += " " + std::to_string(g);
                im.trace.push_back(s + (ok ? " (moved)" : " (could not be moved; left as they are)"));
            }
            if (!ok) continue; // left as it is; the Boolean reports the conflict
            im.owned[k] = copy;
            im.solid[k] = copy.get();
            im.stats.recarried_faces += static_cast<int>(move[k].size());
        }
    }
    detail::Builder builder(im, opt);
    im.ok = builder.run();
    return im;
}

struct BooleanResult {
    bool ok = false;
    bool empty = false; ///< the operation yields no solid
    std::string error;
    Op op = Op::Intersection;
    Brep solid;               ///< self-contained: does not refer to the imprint or the operands
    std::vector<int> regions; ///< parallel to solid.faces: the imprint region
    topo::AuditReport audit;
    /// False when parts of the result touch along an edge (four or more faces
    /// there), at a point (a pinched vertex), or along a line inside a face
    /// (`line_contacts`). Such a result is valid and can be an operand; its
    /// geometry around each edge with more than two faces was checked.
    bool manifold = true;
    /// Lines along which an edge of the result lies inside one of its faces,
    /// its two faces touching that face's surface from one side: the parts
    /// touch there without sharing an edge (the face is not split).
    int line_contacts = 0;
    /// Points at which a curved face of the result touches a flat face of it
    /// inside both (not at a vertex): the parts touch there.
    int point_contacts = 0;
};

/// Assemble the result of one Boolean operation from an imprint.
///
/// A region on the other solid's boundary is kept once, from solid 0, where
/// the operation keeps that boundary: with the same outward normal for union
/// and intersection, with the opposite one for difference.
inline BooleanResult select(const Imprint& im, Op op) {
    BooleanResult r;
    r.op = op;
    if (!im.ok) { r.error = im.error; return r; }
    Brep& out = r.solid;

    std::map<int, int> cmap, vmap;
    std::map<std::pair<bool, int>, int> emap;
    auto carrier = [&](int g) {
        const auto it = cmap.find(g);
        if (it != cmap.end()) return it->second;
        const int id = static_cast<int>(out.carriers.size());
        out.carriers.push_back(im.carriers[static_cast<std::size_t>(g)]);
        cmap[g] = id;
        return id;
    };
    auto vertex = [&](int v) {
        const auto it = vmap.find(v);
        if (it != vmap.end()) return it->second;
        const int id = static_cast<int>(out.vertices.size());
        out.vertices.push_back({im.vertices[static_cast<std::size_t>(v)].point, {}});
        out.vertices.back().exact = im.vertices[static_cast<std::size_t>(v)].exact;
        out.vertices.back().xp = im.vertices[static_cast<std::size_t>(v)].xp;
        vmap[v] = id;
        return id;
    };
    auto edge = [&](const Piece& p) {
        const std::pair<bool, int> key{p.fresh, p.id};
        const auto it = emap.find(key);
        if (it != emap.end()) return it->second;
        BEdge e;
        if (p.fresh) {
            const Imprint::NewArc& a = im.arcs[static_cast<std::size_t>(p.id)];
            e.v0 = vertex(a.v0);
            e.v1 = vertex(a.v1);
            e.arc = true;
            e.carrier[0] = carrier(a.carrier[0]);
            e.carrier[1] = carrier(a.carrier[1]);
            e.geom = a.geom;
        } else {
            const Imprint::OldPiece& pc = im.pieces[static_cast<std::size_t>(p.id)];
            const BEdge& src = im.solid[pc.solid]->edges[static_cast<std::size_t>(pc.edge)];
            const std::vector<int>& cof = im.carrier_of[pc.solid];
            e.v0 = vertex(pc.v0);
            e.v1 = vertex(pc.v1);
            e.arc = src.arc;
            if (src.arc) {
                e.carrier[0] = carrier(cof[static_cast<std::size_t>(src.carrier[0])]);
                e.carrier[1] = carrier(cof[static_cast<std::size_t>(src.carrier[1])]);
                e.geom = pc.geom;
            } else {
                for (int i = 0; i < 2; ++i)
                    e.on[i] = {carrier(cof[static_cast<std::size_t>(src.on[i].carrier)]), src.on[i].side, src.on[i].reversed};
                e.t0 = pc.t0;
                e.t1 = pc.t1;
            }
        }
        const int id = static_cast<int>(out.edges.size());
        out.edges.push_back(std::move(e));
        emap[key] = id;
        return id;
    };

    for (int ri = 0; ri < static_cast<int>(im.regions.size()); ++ri) {
        const Region& reg = im.regions[static_cast<std::size_t>(ri)];
        const int s = reg.solid, c = reg.inside;
        bool keep = false, flip = false;
        switch (op) {
        case Op::Union: keep = c == kOut || (c == kOnSame && s == 0); break;
        case Op::Intersection: keep = c == kIn || (c == kOnSame && s == 0); break;
        case Op::Difference:
            if (s == 0) keep = c == kOut || c == kOnOpposite;
            else { keep = c == kIn; flip = true; }
            break;
        }
        if (!keep) continue;
        const BFace& src = im.solid[s]->faces[static_cast<std::size_t>(reg.face)];
        BFace face;
        face.carrier = carrier(im.carrier_of[s][static_cast<std::size_t>(src.carrier)]);
        face.flipped = src.flipped != flip;
        for (const auto& loop : reg.loops) {
            std::vector<BCoedge> L;
            for (const Piece& p : loop) L.push_back({edge(p), !p.forward});
            face.loops.push_back(std::move(L));
        }
        out.faces.push_back(std::move(face));
        r.regions.push_back(ri);
    }
    // Vertex parameters on the surfaces the result uses.
    for (const auto& kv : vmap)
        for (const Chart& c : im.vertices[static_cast<std::size_t>(kv.first)].charts) {
            const auto it = cmap.find(c.carrier);
            if (it != cmap.end()) out.vertices[static_cast<std::size_t>(kv.second)].charts.push_back({it->second, c.u, c.v});
        }

    r.audit = audit(out);
    r.empty = out.faces.empty();
    r.manifold = r.audit.manifold();
    // A closed, oriented boundary; where it is not manifold (parts touching
    // along an edge or at a point), the faces around every edge with more
    // than two faces must alternate between material and empty space.
    std::string why;
    r.ok = r.empty || (r.audit.closed() && (r.manifold || check_edge_fans(out, why)));
    out.outward_known = r.ok;
    if (!r.ok)
        r.error = !why.empty() ? "the result was refused: " + why
                               : "the selected faces do not form a closed, oriented boundary";
    // Line contacts that the result keeps: the face's region and the edge, or
    // two regions touching along a line inside both.
    if (r.ok && !r.empty && (!im.line_contacts.empty() || !im.touching_regions.empty())) {
        std::vector<char> kept(im.regions.size(), 0);
        for (int ri : r.regions) kept[static_cast<std::size_t>(ri)] = 1;
        for (const auto& [r1, r2] : im.touching_regions)
            if (kept[static_cast<std::size_t>(r1)] && kept[static_cast<std::size_t>(r2)]) ++r.line_contacts;
        for (const Imprint::LineContact& c : im.line_contacts) {
            if (!kept[static_cast<std::size_t>(c.region)]) continue;
            bool edge_kept = false;
            for (int ri : r.regions)
                for (const auto& loop : im.regions[static_cast<std::size_t>(ri)].loops)
                    for (const Piece& p : loop) edge_kept = edge_kept || (!p.fresh && p.id == c.piece);
            if (edge_kept) ++r.line_contacts;
        }
        if (r.line_contacts > 0) r.manifold = false;
    }
    if (r.ok && !r.empty && !im.point_touching_regions.empty()) {
        std::vector<char> kept(im.regions.size(), 0);
        for (int ri : r.regions) kept[static_cast<std::size_t>(ri)] = 1;
        for (const auto& [r1, r2] : im.point_touching_regions)
            if (kept[static_cast<std::size_t>(r1)] && kept[static_cast<std::size_t>(r2)]) ++r.point_contacts;
        if (r.point_contacts > 0) r.manifold = false;
    }
    return r;
}

/// One Boolean operation. The result does not refer to the operands.
/// An empty operand (no faces) is the empty set: A or {} = A, A and {} = {},
/// A - {} = A, {} - B = {}.
inline BooleanResult boolean_op(const Brep& A, const Brep& B, Op op, const BooleanOptions& opt = {}) {
    if (A.faces.empty() || B.faces.empty()) {
        BooleanResult r;
        r.op = op;
        const Brep* keep = nullptr;
        if (op == Op::Union) keep = A.faces.empty() ? &B : &A;
        else if (op == Op::Difference && B.faces.empty()) keep = &A;
        if (keep && !keep->faces.empty()) {
            r.solid = *keep;
            r.audit = audit(r.solid);
            r.manifold = r.audit.manifold();
            std::string why;
            r.ok = r.audit.closed() && (r.manifold || check_edge_fans(r.solid, why));
            if (!r.ok) r.error = "the non-empty operand does not pass the audit";
        } else {
            r.ok = true;
            r.empty = true;
        }
        return r;
    }
    const Imprint im = imprint(A, B, opt);
    return select(im, op);
}

/// Certified classification of a parameter point against a region of an
/// imprint. Refines the region's loop geometry as a side effect.
inline Where locate(Imprint& im, int region, double u, double v, int max_depth = 40) {
    return locate(im.regions[static_cast<std::size_t>(region)].geom, Interval(u), Interval(v), max_depth);
}

/// Largest certified distance between the two surfaces, taken at one
/// certified point per segment of every intersection curve of the imprint.
inline double edge_gap_bound(const Imprint& im) {
    double worst = 0.0;
    for (const auto& base : im.curves) {
        const SsiSystem& sys = *base->sys;
        for (const TubeSegment& seg : base->curve->segs) {
            Box4 enc;
            if (!solve::curve_point(sys, seg, 0.5, enc)) continue;
            const auto pa = sys.point_a(enc), pb = sys.point_b(enc);
            const Interval d2 = num::sqr(pa.x - pb.x) + num::sqr(pa.y - pb.y) + num::sqr(pa.z - pb.z);
            worst = std::max(worst, num::sqrt(d2).hi);
        }
    }
    return worst;
}

} // namespace ak::op
