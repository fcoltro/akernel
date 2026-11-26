// akernel — certified self-intersection check of a solid (issue #6).
//
// The Booleans assume that each operand bounds a solid: its faces meet only
// along the edges and at the vertices the Brep records. The topological
// audit (op::audit) proves the incidence; this proves the geometry: no two
// faces meet anywhere else. The answer is one of
//   Valid      proven: every pair of faces meets at most along an edge they
//              share or at a corner they share;
//   Invalid    proven: two faces cross, with a witness (the faces, and a
//              point on a certified intersection curve, or an overlap of two
//              faces in one plane);
//   Undecided  neither could be proven for some pairs (listed, with the
//              reason for the first).
// Never Valid for an invalid solid; Undecided is the price of rigour.
//
// Each pair of faces is settled by the first test that applies:
//  1. Boxes. Interval enclosures of the two carriers that do not meet.
//  2. Exact tests with a flat face. A flat carrier (every control point in one
//     plane, any weights) has a plane; on which side of it each control point
//     of the other carrier lies is the sign of an exact 4x4 determinant of
//     homogeneous coordinates. By the convex hull property:
//       - the other face strictly on one side: they cannot meet;
//       - on one closed side, touching the plane only with the control points
//         of one straight side, and with the opposite row strictly off it:
//         its distance to the plane, a Bernstein polynomial with those signs,
//         is positive away from that side, so it meets the plane only there
//         (this is how a fillet is tangent to its plate). What is left is
//         collinear: the overlap of two segments on one line, decided by
//         exact comparisons of homogeneous coordinates;
//       - both faces in one plane: a line through a straight side of either,
//         with every control point of one on one closed side and of the other
//         on the other side, separates them; two parallelograms that no such
//         line separates overlap (Invalid).
//     A contact on a line must be a shared edge, a point must be a corner of
//     both faces. Two faces along a shared edge (curved, even tangent, as the
//     pieces of one smooth surface are) are also settled by a plane through
//     the edge with one face on each side, touching it only along the edge.
//  3. Certified surface/surface intersection (solve/ssi.hpp) over both whole
//     carriers, the edges the faces share given as known curves, each
//     certified as a tube. Complete, with no other curve: they meet only
//     along those edges. Another curve: Invalid. Incomplete (a tangency, or a
//     shared corner that subdivision cannot separate): undecided.
//
// Covered: Breps whose faces are whole carriers bounded by their four sides
// (the solids the kernel's constructors build, and fillets on them).
// Trimmed faces, procedural faces and edges along intersection curves are
// reported Undecided.
#pragma once

#include <algorithm>
#include <array>
#include <string>
#include <vector>

#include <akernel/num/expansion.hpp>
#include <akernel/op/boolean.hpp>
#include <akernel/op/planar.hpp>
#include <akernel/solve/ssi.hpp>

namespace ak::op {

struct Validity {
    enum State { Valid, Invalid, Undecided } state = Undecided;
    std::string why;
    int face_a = -1, face_b = -1; ///< the pair that decided Invalid, or the first undecided pair
    Vec3d where{};                ///< Invalid: approximately where they meet
    // How the pairs were settled.
    int pairs = 0, by_boxes = 0, by_planes = 0, by_ssi = 0, undecided = 0;
    std::vector<std::array<int, 2>> open_pairs; ///< the undecided pairs of faces
};

namespace validity_detail {

using num::Expansion;

/// Sign of the 4x4 determinant of homogeneous points (a, b, c, d): which side
/// of the plane through a, b, c the point d lies on (the weights positive).
inline int side(const Vec4d& a, const Vec4d& b, const Vec4d& c, const Vec4d& d) {
    using planar_detail::det3;
    const Expansion det = Expansion(d.w) * det3(a, b, c) - Expansion(c.w) * det3(a, b, d) +
                          Expansion(b.w) * det3(a, c, d) - Expansion(a.w) * det3(b, c, d);
    return det.sign();
}

/// 2D orientation of homogeneous points projected along `axis` (exact; the
/// weights are positive).
inline int orient2(const Vec4d& a, const Vec4d& b, const Vec4d& c, int axis) {
    auto x = [&](const Vec4d& p) { return axis == 0 ? p.y : p.x; };
    auto y = [&](const Vec4d& p) { return axis == 2 ? p.y : p.z; };
    using E = Expansion;
    const E d = E(x(a)) * (E(y(b)) * E(c.w) - E(b.w) * E(y(c))) - E(y(a)) * (E(x(b)) * E(c.w) - E(b.w) * E(x(c))) +
                E(a.w) * (E(x(b)) * E(y(c)) - E(y(b)) * E(x(c)));
    return d.sign();
}

/// Compare coordinate k of two homogeneous points (exact; weights positive).
inline int cmp(const Vec4d& a, const Vec4d& b, int k) {
    auto c = [&](const Vec4d& p) { return k == 0 ? p.x : k == 1 ? p.y : p.z; };
    return (Expansion(c(a)) * Expansion(b.w) - Expansion(c(b)) * Expansion(a.w)).sign();
}

inline bool same(const Vec4d& a, const Vec4d& b) { return planar_detail::same_point_h(a, b); }

/// The corners of a patch.
inline std::array<Vec4d, 4> corners(const PolyPatchD& P) {
    using planar_detail::corner;
    return {corner(P, 0, 0), corner(P, 1, 0), corner(P, 1, 1), corner(P, 0, 1)};
}
inline bool is_corner(const PolyPatchD& P, const Vec4d& x) {
    for (const Vec4d& c : corners(P))
        if (same(c, x)) return true;
    return false;
}

/// The controls of side s (patch_solid numbering: 0 v = 0, 1 u = 1, 2 v = 1,
/// 3 u = 0), or of the side opposite it.
inline std::vector<Vec4d> row(const PolyPatchD& P, int s, bool opposite) { return side_control(P, opposite ? (s + 2) % 4 : s); }

/// Overlap of two collinear point sets' hulls: empty, one point, or a segment.
struct Overlap1 {
    enum Kind { Empty, Point, Segment } kind = Empty;
    Vec4d at{};
};
inline Overlap1 overlap_on_line(const std::vector<Vec4d>& X, const std::vector<Vec4d>& Y) {
    Overlap1 r;
    if (X.empty() || Y.empty()) return r;
    // A coordinate along which the points are not all equal.
    int k = -1;
    for (int c = 0; c < 3 && k < 0; ++c)
        for (const auto* S : {&X, &Y})
            for (const Vec4d& p : *S)
                if (cmp(p, X[0], c) != 0) k = c;
    if (k < 0) {
        r.kind = Overlap1::Point;
        r.at = X[0];
        return r;
    }
    auto lo = [&](const std::vector<Vec4d>& S) {
        Vec4d m = S[0];
        for (const Vec4d& p : S)
            if (cmp(p, m, k) < 0) m = p;
        return m;
    };
    auto hi = [&](const std::vector<Vec4d>& S) {
        Vec4d m = S[0];
        for (const Vec4d& p : S)
            if (cmp(p, m, k) > 0) m = p;
        return m;
    };
    const Vec4d a = cmp(lo(X), lo(Y), k) > 0 ? lo(X) : lo(Y), b = cmp(hi(X), hi(Y), k) < 0 ? hi(X) : hi(Y);
    const int c = cmp(a, b, k);
    if (c > 0) return r;
    r.kind = c == 0 ? Overlap1::Point : Overlap1::Segment;
    r.at = a;
    return r;
}

enum class Pair { Disjoint, Touching, Overlap, Unknown };

/// The control points of B on the plane through t, and on which side the
/// rest lie (+1, -1, 0 if none is off the plane, 2 if both sides).
struct Sides {
    std::vector<Vec4d> on;
    int rest = 0;
};
inline Sides sides(const Vec4d t[3], const PolyPatchD& B) {
    Sides s;
    for (const Vec4d& c : B.c) {
        const int o = side(t[0], t[1], t[2], c);
        if (o == 0) s.on.push_back(c);
        else if (s.rest == 0) s.rest = o;
        else if (s.rest != o) s.rest = 2;
    }
    return s;
}

/// B on one closed side of a plane, touching it: the side of B whose controls
/// hold every control point on the plane, with the opposite row strictly off
/// it. Then B meets the plane exactly along that side. -1 if none.
inline int contact_side(const PolyPatchD& B, const Sides& s) {
    if (s.rest == 2 || s.rest == 0 || s.on.empty()) return -1;
    for (int sd = 0; sd < 4; ++sd) {
        const std::vector<Vec4d> r = row(B, sd, false);
        bool holds = true;
        for (const Vec4d& p : s.on) {
            bool in = false;
            for (const Vec4d& q : r) in = in || same(p, q);
            holds = holds && in;
        }
        if (!holds) continue;
        bool off = true;
        for (const Vec4d& q : row(B, sd, true))
            for (const Vec4d& p : s.on) off = off && !same(p, q);
        if (off) return sd;
    }
    return -1;
}

/// A projection axis in which the plane through t is not degenerate.
inline int axis_of(const Vec4d t[3]) {
    for (int k = 0; k < 3; ++k)
        if (orient2(t[0], t[1], t[2], k) != 0) return k;
    return -1;
}

/// Are the points collinear in the projection, and not all one point?
inline bool straight(const std::vector<Vec4d>& P, int axis) {
    if (cmp(P.front(), P.back(), 0) == 0 && cmp(P.front(), P.back(), 1) == 0 && cmp(P.front(), P.back(), 2) == 0) return false;
    for (const Vec4d& c : P)
        if (orient2(P.front(), P.back(), c, axis) != 0) return false;
    return true;
}

/// Sides of A and of B on edges the two faces share.
struct Shared {
    std::vector<int> sides_a, sides_b;
};
inline bool has(const std::vector<int>& v, int s) { return std::find(v.begin(), v.end(), s) != v.end(); }

/// Split the points by a line in the projection: the side the off-line ones
/// lie on (0 none, 2 both), and those on the line.
inline int split(const std::vector<Vec4d>& pts, const Vec4d& p, const Vec4d& q, int axis, std::vector<Vec4d>& on) {
    int s = 0;
    for (const Vec4d& c : pts) {
        const int o = orient2(p, q, c, axis);
        if (o == 0) on.push_back(c);
        else if (s == 0) s = o;
        else if (s != o) s = 2;
    }
    return s;
}

/// Two faces in one plane: a separating line through a straight side of
/// either, the contact on it at most a shared corner or a shared edge.
inline Pair coplanar(const PolyPatchD& A, const PolyPatchD& B, const Vec4d t[3], const Shared& sh) {
    const int axis = axis_of(t);
    if (axis < 0) return Pair::Unknown;
    for (int which = 0; which < 2; ++which) {
        const PolyPatchD& X = which == 0 ? A : B;
        const PolyPatchD& Y = which == 0 ? B : A;
        const std::vector<int>& xs = which == 0 ? sh.sides_a : sh.sides_b;
        for (int sd = 0; sd < 4; ++sd) {
            const std::vector<Vec4d> line = row(X, sd, false);
            if (!straight(line, axis)) continue;
            std::vector<Vec4d> zx, zy;
            const int sx = split(X.c, line.front(), line.back(), axis, zx);
            const int sy = split(Y.c, line.front(), line.back(), axis, zy);
            if (sx == 0 || sx == 2 || sy == 2 || sx == sy) continue;
            const Overlap1 o = overlap_on_line(zx, zy);
            if (o.kind == Overlap1::Empty) return Pair::Disjoint;
            if (o.kind == Overlap1::Point && is_corner(A, o.at) && is_corner(B, o.at)) return Pair::Touching;
            if (o.kind == Overlap1::Segment && has(xs, sd)) return Pair::Touching;
        }
    }
    // Two parallelograms that no edge line separates overlap.
    return planar_affine(A) && planar_affine(B) ? Pair::Overlap : Pair::Unknown;
}

/// A flat face A (its plane through tA) and any face B, decided exactly
/// where the signs allow.
inline Pair against_flat(const PolyPatchD& A, const Vec4d tA[3], const PolyPatchD& B, const Shared& sh) {
    for (const Vec4d& c : B.c)
        if (!(c.w > 0.0)) return Pair::Unknown;
    const Sides sB = sides(tA, B);
    if (sB.on.size() == B.c.size()) return coplanar(A, B, tA, sh);
    if (sB.on.empty() && sB.rest != 2) return Pair::Disjoint;
    const int cs = contact_side(B, sB);
    if (cs < 0) return Pair::Unknown;
    // B meets A's plane only along its side cs.
    if (has(sh.sides_b, cs)) return Pair::Touching; // the shared edge
    const int axis = axis_of(tA);
    const std::vector<Vec4d> seg = row(B, cs, false);
    if (axis < 0) return Pair::Unknown;
    // That side is a curve in A's plane, inside the hull of its controls:
    // separate the hull from A by a straight side line of A.
    for (int sd = 0; sd < 4; ++sd) {
        const std::vector<Vec4d> line = row(A, sd, false);
        if (!straight(line, axis)) continue;
        std::vector<Vec4d> za, zs;
        const int sa = split(A.c, line.front(), line.back(), axis, za);
        const int ss = split(seg, line.front(), line.back(), axis, zs);
        if (sa == 0 || sa == 2 || ss == 2 || ss == sa) continue;
        const Overlap1 o = overlap_on_line(za, zs);
        if (o.kind == Overlap1::Empty) return Pair::Disjoint;
        if (o.kind == Overlap1::Point && is_corner(A, o.at) && is_corner(B, o.at)) return Pair::Touching;
    }
    return Pair::Unknown;
}

/// Two faces along a shared side (A's side sa, B's side sb), curved or not,
/// tangent or not: a plane through the edge's end points that has A on one
/// closed side and B on the other, each meeting it only along the shared
/// side (contact_side). Then A and B meet only along the edge. Candidates:
/// the plane of a curved edge (an arc between a torus zone and a cylinder),
/// and a plane through a straight edge and A's normal at its middle, guessed
/// in doubles (it separates two pieces of one smooth surface cut along a
/// ruling). The test itself is exact, so any plane that passes is a proof.
inline bool separating_plane(const PolyPatchD& A, int sa, const PolyPatchD& B, int sb) {
    const std::vector<Vec4d> e = row(A, sa, false);
    const Vec4d &p = e.front(), &q = e.back();
    if (cmp(p, q, 0) == 0 && cmp(p, q, 1) == 0 && cmp(p, q, 2) == 0) return false;
    auto test = [&](const Vec4d t[3]) {
        const Sides x = sides(t, A), y = sides(t, B);
        if (x.rest == 2 || y.rest == 2 || x.rest == 0 || y.rest == 0 || x.rest == y.rest) return false;
        return contact_side(A, x) == sa && contact_side(B, y) == sb;
    };
    // A curved edge in a plane (an arc): that plane.
    for (const Vec4d& c : e) {
        const Vec4d t[3] = {p, c, q};
        if (axis_of(t) >= 0 && test(t)) return true;
    }
    double u = 0, v = 0;
    side_uv(sa, 0.5, u, v);
    const geom::Jet3 J = geom::BezierSurface(A).jet(u, v);
    const Vec3d n = cross(J[geom::jet_index(1, 0)], J[geom::jet_index(0, 1)]);
    const Vec3d m = J[0];
    for (const Vec3d& d : {n, n * -1.0}) {
        const Vec4d t[3] = {p, q, {m.x + d.x, m.y + d.y, m.z + d.z, 1.0}};
        if (test(t)) return true;
    }
    return false;
}

/// Do the interval boxes of two carriers miss each other?
inline bool boxes_apart(const PolyPatchD& A, const PolyPatchD& B) {
    const auto a = geom::projected_bounds(geom::to_interval(A)), b = geom::projected_bounds(geom::to_interval(B));
    return a.x.hi < b.x.lo || b.x.hi < a.x.lo || a.y.hi < b.y.lo || b.y.hi < a.y.lo || a.z.hi < b.z.lo || b.z.hi < a.z.lo;
}

} // namespace validity_detail

/// Is the solid free of self-intersection? See the file comment.
inline Validity check_self_intersection(const Brep& b, const solve::SsiOptions& opt = {}) {
    using namespace validity_detail;
    Validity res;
    const int nf = static_cast<int>(b.faces.size());
    auto undecided = [&](const std::string& why, int f, int g) {
        res.state = Validity::Undecided;
        res.why = why;
        res.face_a = f;
        res.face_b = g;
        return res;
    };
    auto pair_undecided = [&](const char* why, int f, int g) {
        res.open_pairs.push_back({f, g});
        if (res.undecided++ == 0) {
            res.why = why;
            res.face_a = f;
            res.face_b = g;
        }
    };
    auto invalid = [&](const char* why, int f, int g, const Vec3d& where) {
        res.state = Validity::Invalid;
        res.why = why;
        res.face_a = f;
        res.face_b = g;
        res.where = where;
        return res;
    };
    // Every face a whole carrier bounded by its four sides.
    for (int f = 0; f < nf; ++f) {
        const BFace& face = b.faces[static_cast<std::size_t>(f)];
        if (face.proc >= 0) return undecided("a procedural face", f, -1);
        if (face.loops.size() != 1 || face.loops[0].size() != 4) return undecided("a trimmed face", f, -1);
        for (const BCoedge& ce : face.loops[0]) {
            const BEdge& e = b.edges[static_cast<std::size_t>(ce.edge)];
            if (e.arc || e.proc || !side_on(e, face.carrier)) return undecided("a trimmed face", f, -1);
        }
        for (int g = 0; g < f; ++g)
            if (b.faces[static_cast<std::size_t>(g)].carrier == face.carrier) return undecided("two faces on one carrier", g, f);
    }
    for (int f = 0; f < nf; ++f)
        for (int g = f + 1; g < nf; ++g) {
            ++res.pairs;
            const int cf = b.faces[static_cast<std::size_t>(f)].carrier, cg = b.faces[static_cast<std::size_t>(g)].carrier;
            const PolyPatchD& A = b.carriers[static_cast<std::size_t>(cf)];
            const PolyPatchD& B = b.carriers[static_cast<std::size_t>(cg)];
            if (boxes_apart(A, B)) {
                ++res.by_boxes;
                continue;
            }
            // The edges the two faces share.
            Shared sh;
            std::vector<const BEdge*> shared_edges;
            for (const BEdge& e : b.edges) {
                if (e.arc || e.proc) continue;
                const BEdge::SideRef *ra = side_on(e, cf), *rb = side_on(e, cg);
                if (!ra || !rb || ra == rb) continue;
                sh.sides_a.push_back(ra->side);
                sh.sides_b.push_back(rb->side);
                shared_edges.push_back(&e);
            }
            // Exact tests against either face's plane.
            Vec4d tA[3], tB[3];
            Pair p = Pair::Unknown;
            if (flat_triple(A, tA)) p = against_flat(A, tA, B, sh);
            if (p == Pair::Unknown && flat_triple(B, tB)) p = against_flat(B, tB, A, Shared{sh.sides_b, sh.sides_a});
            if (p == Pair::Disjoint || p == Pair::Touching) {
                ++res.by_planes;
                continue;
            }
            if (p == Pair::Overlap) {
                const Vec4d c = A.eval(0.5, 0.5);
                return invalid("two faces in one plane overlap", f, g, {c.x / c.w, c.y / c.w, c.z / c.w});
            }
            // Two faces along one shared edge, separated by a plane through it.
            if (shared_edges.size() == 1 && separating_plane(A, sh.sides_a[0], B, sh.sides_b[0])) {
                ++res.by_planes;
                continue;
            }
            // Certified intersection, the shared edges as known curves.
            const solve::SsiSystem sys(A, B);
            std::vector<solve::SsiCurve> known;
            bool tangent_edge = false;
            for (const BEdge* e : shared_edges) {
                const BEdge::SideRef *ra = side_on(*e, cf), *rb = side_on(*e, cg);
                solve::SsiCurve c;
                if (!detail::side_curve(sys, ra->side, rb->side, ra->reversed != rb->reversed, opt, c)) {
                    tangent_edge = true;
                    break;
                }
                known.push_back(std::move(c));
            }
            if (tangent_edge) {
                pair_undecided("a shared edge not certified as a transversal intersection (the faces may be tangent there)", f, g);
                continue;
            }
            const solve::SsiResult r = solve::intersect(sys, opt, known);
            for (const solve::SsiCurve& c : r.curves) {
                if (c.known || c.segs.empty()) continue;
                const solve::TubeSegment& s = c.segs[c.segs.size() / 2];
                const Vec4d q = A.eval(s.a[0], s.a[1]);
                return invalid("two faces cross along a certified intersection curve", f, g, {q.x / q.w, q.y / q.w, q.z / q.w});
            }
            if (!r.complete) {
                pair_undecided("the intersection of two faces not proven complete (a tangency, or a shared corner)", f, g);
                continue;
            }
            ++res.by_ssi;
        }
    // Every pair is checked, so an Invalid pair is found even after an
    // undecided one; the first undecided pair is the one reported.
    res.state = res.undecided ? Validity::Undecided : Validity::Valid;
    return res;
}

} // namespace ak::op
