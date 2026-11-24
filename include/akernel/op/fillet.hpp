// akernel — fillets and chamfers applied to solids: trim the faces, replace
// the edge (docs/FEATURES.md, step F5).
//
// The solids here are swept: a planar cross-section extruded (a box) or
// revolved about an axis (a ring, by quarter turns up to a closed full turn).
// The edge being rounded or chamfered is a corner of the section, so the
// blend is the section's corner profile swept the same way: for a fillet the
// exact rolling-ball arc of geom/blend.hpp (a cylinder sector when extruded,
// a torus zone when revolved), for a chamfer a straight cut (a plane or a
// cone). A revolved edge is a closed chain of four edges with a tangent-
// continuous blend all round; the result is a closed PatchSolid either way:
// every face a whole rational patch, every edge one shared definition, so it
// is watertight by construction and can be an operand of the Booleans.
//
// The section. Section coordinates (x, y), the region x in [0, a], y in
// [0, b], the corner (a, b) cut. The supports are trimmed at the contacts
// A = (a, b - dA) on the side x = a and B = (a - dB, b) on the side y = b
// (dA = dB = r for a fillet). The section is cut into five ruled cells
// (s = dA / 2):
//
//      Y ----------------- B .
//      |        Qd         |    .  profile (arc or straight cut)
//     Yl ---------------- C1      .
//      |                   | Qcor   A
//      |        Qc         C2 ---/  |
//      |                   |   \    |
//      |                   | Qb |Qa |
//      O ---------------- Xr - Xs - X
//
//   C1 = (a - dB, b - dA - s), C2 = (a - dB/2, b - dA - s), Xr = (a - dB, 0),
//   Xs = (a - dB/2, 0), Yl = (0, b - dA - s).
//
// Why five cells. A fillet's profile is tangent to the trimmed sides at A and
// B, so a cell with a corner at A bounded by the side and the arc would have a
// straight angle there: a singular parametrisation. The corner cell Qcor is
// ruled between the profile A -> B and the segment C2 -> C1 instead, which
// gives interior angles below 180 degrees everywhere. The cells' edges cross
// the section boundary at Xr, Xs and Yl, so the side faces are split there to
// match.
//
// Weights. The exact right-angle arc has weights (1, 1, 2): contact A weight
// 1, contact B weight 2, no square root (geom/blend.hpp). Shared sides must
// have the same homogeneous controls, so every section point carries one
// weight wherever it is used: X and A weight 1, all others 2 (a uniform
// factor does not change a patch). The straight cut uses the same weights,
// its middle control the degree elevation of A -> B
// (section points are stored homogeneous so that this stays exact).
//
// Revolution. A section point (rho, h) with weight w at the circle control
// (c, s, k) is the homogeneous point (rho w c, rho w s, h w k, w k). The
// quarters alternate between the arc 0 -> 90 degrees with end weights (1, 2)
// and its mirror image with (2, 1): consecutive quarters share their end rows
// exactly, also across the seam of a full turn (geom/nurbs_curve.hpp's
// exact_circle scales by 2 per quarter instead and does not close).
//
// With dyadic dimensions and axis-aligned frames every control is exact.
#pragma once

#include <array>
#include <cmath>
#include <optional>
#include <string>
#include <vector>

#include <akernel/geom/bernstein.hpp>
#include <akernel/op/patch_solid.hpp>

namespace ak::op {

/// How a section corner is cut: a rolling-ball fillet of radius r (setbacks
/// dA = dB = r) or a straight chamfer with setbacks dA and dB.
struct EdgeCut {
    enum Kind { Round, Flat } kind = Round;
    double dA = 0.0, dB = 0.0;
    static EdgeCut round(double r) { return {Round, r, r}; }
    static EdgeCut flat(double dA, double dB) { return {Flat, dA, dB}; }
};

struct BlendedSolid {
    std::optional<PatchSolid> solid;
    std::string error;
};
using FilletedSolid = BlendedSolid;

namespace fillet_detail {

using geom::Vec4d;

struct SPt {
    double X = 0, Y = 0, w = 1; ///< homogeneous section point (x w, y w, w)
};
inline SPt spt(double x, double y, double w) { return {x * w, y * w, w}; }
/// Average of homogeneous points: the degree elevation of the segment p -> q.
inline SPt mid(const SPt& p, const SPt& q) { return {0.5 * (p.X + q.X), 0.5 * (p.Y + q.Y), 0.5 * (p.w + q.w)}; }
using SCurve = std::vector<SPt>;

struct Section {
    struct Cell {
        SCurve row0, row1;
    };
    struct Side {
        SCurve curve;
        double nx = 0, ny = 0; ///< outward normal in section coordinates
    };
    std::vector<Cell> cells;
    std::vector<Side> sides;
};

inline Section layout(double a, double b, const EdgeCut& cut) {
    const double dA = cut.dA, dB = cut.dB, s = 0.5 * dA;
    auto p2 = [](double x, double y) { return spt(x, y, 2.0); };
    const SPt O = p2(0, 0), Xr = p2(a - dB, 0), Xs = p2(a - 0.5 * dB, 0), X = spt(a, 0, 1);
    const SPt A = spt(a, b - dA, 1), B = p2(a - dB, b), Y = p2(0, b);
    const SPt C1 = p2(a - dB, b - dA - s), C2 = p2(a - 0.5 * dB, b - dA - s), Yl = p2(0, b - dA - s);
    // Round: where the tangents at A and B meet (geom/blend.hpp, fillet_arc).
    // Flat: the degree elevation of the homogeneous segment A -> B.
    const SPt M = cut.kind == EdgeCut::Round ? spt(a, b, 1) : mid(A, B);
    const SCurve profile{A, M, B}, line{C2, mid(C2, C1), C1};
    Section S;
    S.cells = {{profile, line},
               {{X, A}, {Xs, C2}},
               {{Xs, mid(Xr, Xs), Xr}, line},
               {{O, Xr}, {Yl, C1}},
               {{Yl, C1}, {Y, B}}};
    S.sides = {{{O, Xr}, 0, -1},       {{Xr, mid(Xr, Xs), Xs}, 0, -1}, {{Xs, X}, 0, -1},
               {{X, A}, 1, 0},         {profile, dA, dB},             {{B, Y}, 0, 1},
               {{Y, Yl}, -1, 0},       {{Yl, O}, -1, 0}};
    return S;
}

/// Patch of degree (p, 1) between two rows of p + 1 homogeneous controls.
inline PolyPatchD ruled(const std::vector<Vec4d>& r0, const std::vector<Vec4d>& r1) {
    const int p = static_cast<int>(r0.size()) - 1;
    PolyPatchD P(p, 1);
    for (int i = 0; i <= p; ++i) {
        P.at(i, 0) = r0[static_cast<std::size_t>(i)];
        P.at(i, 1) = r1[static_cast<std::size_t>(i)];
    }
    return P;
}

/// The same patch with u reversed, so that Su x Sv flips.
inline PolyPatchD reverse_u(const PolyPatchD& P) {
    PolyPatchD R(P.p, P.q);
    for (int i = 0; i <= P.p; ++i)
        for (int j = 0; j <= P.q; ++j) R.at(i, j) = P.at(P.p - i, j);
    return R;
}

/// Orient P so that Su x Sv points along `out` at the middle.
inline PolyPatchD outward(const PolyPatchD& P, const Vec3d& out) {
    const geom::Jet3 J = geom::BezierSurface(P).jet(0.5, 0.5);
    return dot(cross(J[1], J[2]), out) > 0.0 ? P : reverse_u(P);
}

inline bool check_cut(double a, double b, const EdgeCut& cut, std::string& error) {
    if (!(cut.dA > 0.0) || !(cut.dB > 0.0) || !(cut.dB < a) || !(1.5 * cut.dA < b)) {
        error = "setback out of range for this section (needs 0 < dB < a and 0 < 1.5 dA < b)";
        return false;
    }
    if (cut.kind == EdgeCut::Round && cut.dA != cut.dB) {
        error = "a fillet has equal setbacks";
        return false;
    }
    return true;
}

inline BlendedSolid finish(PatchSolid S) {
    BlendedSolid res;
    std::string err;
    if (!S.build(err)) {
        res.error = "the blended shell does not close: " + err;
        return res;
    }
    res.solid = std::move(S);
    return res;
}

/// Circle controls (c, s, k) of quarter q (see Revolution above).
inline std::array<std::array<double, 3>, 3> quarter(int q) {
    static constexpr double Q[4][3][3] = {{{1, 0, 1}, {1, 1, 1}, {0, 2, 2}},
                                          {{0, 2, 2}, {-1, 1, 1}, {-1, 0, 1}},
                                          {{-1, 0, 1}, {-1, -1, 1}, {0, -2, 2}},
                                          {{0, -2, 2}, {1, -1, 1}, {1, 0, 1}}};
    std::array<std::array<double, 3>, 3> r{};
    for (int i = 0; i < 3; ++i)
        for (int k = 0; k < 3; ++k) r[static_cast<std::size_t>(i)][static_cast<std::size_t>(k)] = Q[q][i][k];
    return r;
}

/// Concave edge: the L-shaped section [0, a] x [0, b1] united with
/// [0, a1] x [b1, b] (a step), its inside corner (a1, b1) filled by the
/// profile from A = (a1, b1 + dA) on the wall x = a1 to B = (a1 + dB, b1)
/// on the floor y = b1. Seven ruled cells:
//
//      Y -------- T -- Yt
//      |    K5a   | K5b |
//     Lc -------- A' -- A .
//      |          |      .  profile
//      |    K4    | Qcor  .
//     Lb -------- D ------- B -------- Xt
//      |          |   \__   ...  K3   |
//      |    K1    |      \_____       Ra      (D -> B -> Xt -> Ra -> D is K3)
//      |          |    K2              |
//      O -------- Xd ----------------- X
//
//   A' = (a1 - dB, b1 + dA), D = (a1 - dB, b1 - dA), T = (a1 - dB, b), Lc =
//   (0, b1 + dA), Lb = (0, b1 - dA), Xd = (a1 - dB, 0), Ra = (a, b1 - dA).
//   Qcor is ruled between the profile and the segment A' -> D (quadratic,
//   so K4 is quadratic along y). Only B carries weight 2.
inline Section layout_step(double a, double b, double a1, double b1, const EdgeCut& cut) {
    const double dA = cut.dA, dB = cut.dB;
    auto p1 = [](double x, double y) { return spt(x, y, 1.0); };
    const SPt O = p1(0, 0), Xd = p1(a1 - dB, 0), X = p1(a, 0), Ra = p1(a, b1 - dA), Xt = p1(a, b1);
    const SPt B = spt(a1 + dB, b1, 2.0), A = p1(a1, b1 + dA), Yt = p1(a1, b), T = p1(a1 - dB, b), Y = p1(0, b);
    const SPt Ap = p1(a1 - dB, b1 + dA), D = p1(a1 - dB, b1 - dA), Lc = p1(0, b1 + dA), Lb = p1(0, b1 - dA);
    const SPt M = cut.kind == EdgeCut::Round ? p1(a1, b1) : mid(A, B);
    const SCurve profile{A, M, B}, inner{Ap, mid(Ap, D), D};
    Section S;
    S.cells = {{profile, inner},
               {{O, Xd}, {Lb, D}},
               {{Xd, X}, {D, Ra}},
               {{D, Ra}, {B, Xt}},
               {{Lb, mid(Lb, Lc), Lc}, SCurve{D, mid(D, Ap), Ap}},
               {{Lc, Ap}, {Y, T}},
               {{Ap, A}, {T, Yt}}};
    S.sides = {{{O, Xd}, 0, -1},       {{Xd, X}, 0, -1},     {{X, Ra}, 1, 0},  {{Ra, Xt}, 1, 0},
               {{Xt, B}, 0, 1},        {profile, dA, dB}, {{A, Yt}, 1, 0},  {{Yt, T}, 0, 1},
               {{T, Y}, 0, 1},         {{Y, Lc}, -1, 0},     {{Lc, mid(Lc, Lb), Lb}, -1, 0}, {{Lb, O}, -1, 0}};
    return S;
}

inline bool check_step(double a, double b, double a1, double b1, const EdgeCut& cut, std::string& error) {
    if (!(cut.dA > 0.0) || !(cut.dB > 0.0) || !(cut.dB < a1) || !(cut.dA < b1) || !(a1 + cut.dB < a) ||
        !(b1 + cut.dA < b)) {
        error = "setback out of range for this step";
        return false;
    }
    if (cut.kind == EdgeCut::Round && cut.dA != cut.dB) {
        error = "a fillet has equal setbacks";
        return false;
    }
    return true;
}

/// The section extruded: o + x X + y Y + z Z, z in [0, c].
inline BlendedSolid extrude_section(const Section& sec, const Vec3d& o, const Vec3d& X, const Vec3d& Y,
                                    const Vec3d& Z, double c) {
    auto row = [&](const SCurve& cv, double z) {
        std::vector<Vec4d> r;
        for (const SPt& p : cv) {
            const Vec3d q = o * p.w + X * p.X + Y * p.Y + Z * (z * p.w);
            r.push_back({q.x, q.y, q.z, p.w});
        }
        return r;
    };
    PatchSolid S;
    for (int k = 0; k < 2; ++k)
        for (const auto& cell : sec.cells)
            S.faces.push_back(outward(ruled(row(cell.row0, k * c), row(cell.row1, k * c)), k ? Z : Z * -1.0));
    for (const auto& side : sec.sides)
        S.faces.push_back(outward(ruled(row(side.curve, 0.0), row(side.curve, c)), X * side.nx + Y * side.ny));
    return finish(std::move(S));
}

/// The section revolved about the z axis by `quarters` quarter turns from
/// +x, placed by (rho, h) = (p0r + ur x, p0h + vh y), ur, vh = +-1.
inline BlendedSolid revolve_section(const Section& sec, double p0r, double p0h, double ur, double vh,
                                    int quarters) {
    BlendedSolid res;
    for (const auto& cell : sec.cells)
        for (const SCurve* cv : {&cell.row0, &cell.row1})
            for (const SPt& p : *cv)
                if (!(p0r * p.w + ur * p.X > 0.0)) {
                    res.error = "the section reaches the axis";
                    return res;
                }
    auto hom = [&](const SPt& p, const std::array<double, 3>& cs) {
        const double rho = p0r * p.w + ur * p.X, h = p0h * p.w + vh * p.Y; // times w
        return Vec4d{rho * cs[0], rho * cs[1], h * cs[2], p.w * cs[2]};
    };
    auto row = [&](const SCurve& cv, const std::array<double, 3>& cs) {
        std::vector<Vec4d> r;
        for (const SPt& p : cv) r.push_back(hom(p, cs));
        return r;
    };
    PatchSolid S;
    for (int q = 0; q < quarters; ++q) {
        const auto Q = quarter(q);
        const double phi = (q + 0.5) * M_PI / 2;
        for (const auto& side : sec.sides) {
            PolyPatchD P(2, static_cast<int>(side.curve.size()) - 1);
            for (int i = 0; i < 3; ++i)
                for (std::size_t j = 0; j < side.curve.size(); ++j)
                    P.at(i, static_cast<int>(j)) = hom(side.curve[j], Q[static_cast<std::size_t>(i)]);
            const double nr = ur * side.nx, nh = vh * side.ny;
            S.faces.push_back(outward(P, {nr * std::cos(phi), nr * std::sin(phi), nh}));
        }
    }
    if (quarters < 4) {
        // End faces: radial half-planes at angle 0 and quarters * 90 degrees.
        const auto start = quarter(0)[0], end = quarter(quarters - 1)[2];
        const double th = quarters * M_PI / 2;
        const Vec3d out0{0, -1, 0}, out1{-std::sin(th), std::cos(th), 0};
        for (const auto& cell : sec.cells) {
            S.faces.push_back(outward(ruled(row(cell.row0, start), row(cell.row1, start)), out0));
            S.faces.push_back(outward(ruled(row(cell.row0, end), row(cell.row1, end)), out1));
        }
    }
    return finish(std::move(S));
}

} // namespace fillet_detail

/// Round (EdgeCut::round) or chamfer (EdgeCut::flat) the edge x = a, y = b
/// (along Z) of the box o + x X + y Y + z Z with sides a = |ex|, b = |ey|,
/// c = |ez|. ex, ey, ez must be orthogonal and right-handed. dA is the
/// setback on the face x = a, dB on the face y = b.
inline BlendedSolid blend_box_edge(const Vec3d& o, const Vec3d& ex, const Vec3d& ey, const Vec3d& ez,
                                   const EdgeCut& cut) {
    using namespace fillet_detail;
    BlendedSolid res;
    const double a = norm(ex), b = norm(ey), c = norm(ez);
    if (!check_cut(a, b, cut, res.error)) return res;
    const Vec3d X = ex / a, Y = ey / b, Z = ez / c;
    if (!(dot(cross(X, Y), Z) > 0.0)) {
        res.error = "the box axes must be right-handed";
        return res;
    }
    return extrude_section(layout(a, b, cut), o, X, Y, Z, c);
}

inline BlendedSolid fillet_box_edge(const Vec3d& o, const Vec3d& ex, const Vec3d& ey, const Vec3d& ez, double r) {
    return blend_box_edge(o, ex, ey, ez, EdgeCut::round(r));
}
inline BlendedSolid chamfer_box_edge(const Vec3d& o, const Vec3d& ex, const Vec3d& ey, const Vec3d& ez, double dA,
                                     double dB) {
    return blend_box_edge(o, ex, ey, ez, EdgeCut::flat(dA, dB));
}

/// Concave edge: the step block o + x X + y Y + z Z, [0, a] x [0, b] x [0, c]
/// (a = |ex|, b = |ey|, c = |ez|) without the notch x > a1, y > b1. Its inside
/// edge x = a1, y = b1 is filled by a fillet or chamfer (material is added);
/// dA is the setback on the wall x = a1, dB on the floor y = b1.
inline BlendedSolid blend_step_edge(const Vec3d& o, const Vec3d& ex, const Vec3d& ey, const Vec3d& ez, double a1,
                                    double b1, const EdgeCut& cut) {
    using namespace fillet_detail;
    BlendedSolid res;
    const double a = norm(ex), b = norm(ey), c = norm(ez);
    if (!check_step(a, b, a1, b1, cut, res.error)) return res;
    const Vec3d X = ex / a, Y = ey / b, Z = ez / c;
    if (!(dot(cross(X, Y), Z) > 0.0)) {
        res.error = "the box axes must be right-handed";
        return res;
    }
    return extrude_section(layout_step(a, b, a1, b1, cut), o, X, Y, Z, c);
}

/// The four circular edges of a ring.
enum class RingEdge { OuterTop, OuterBottom, InnerTop, InnerBottom };

/// Round or chamfer one circular edge of the ring r0 <= rho <= r1,
/// z0 <= z <= z1 about the z axis, swept by `quarters` quarter turns from
/// the +x axis (4: the closed ring, the edge a closed chain of four arcs).
/// dA is the setback on the cylinder, dB on the flat face; needs r0 > 0.
inline BlendedSolid blend_ring_edge(double r0, double r1, double z0, double z1, RingEdge edge, const EdgeCut& cut,
                                    int quarters = 4) {
    using namespace fillet_detail;
    BlendedSolid res;
    if (!(r0 > 0.0) || !(r1 > r0) || !(z1 > z0) || quarters < 1 || quarters > 4) {
        res.error = "bad ring";
        return res;
    }
    // Section frame: (rho, h) = p0 + x u + y v, the cut corner at (a, b).
    const bool outer = edge == RingEdge::OuterTop || edge == RingEdge::OuterBottom;
    const bool top = edge == RingEdge::OuterTop || edge == RingEdge::InnerTop;
    const double a = r1 - r0, b = z1 - z0;
    if (!check_cut(a, b, cut, res.error)) return res;
    return revolve_section(layout(a, b, cut), outer ? r0 : r1, top ? z0 : z1, outer ? 1.0 : -1.0, top ? 1.0 : -1.0,
                           quarters);
}

/// A stepped ring about the z axis, r0 <= rho <= r2, z0 <= z: the tall part
/// reaches z2 and the low part z1 < z2, divided at rho = r1. With
/// `tall_inside` the tall part is r0..r1 (a boss on a plate), otherwise
/// r1..r2 (a counterbore: a floor inside a wall). The inside circular edge
/// at (r1, z1) is filled by a fillet or chamfer: a torus zone or a cone that
/// adds material. dA is the setback on the cylinder, dB on the floor.
inline BlendedSolid blend_stepped_ring_edge(double r0, double r1, double r2, double z0, double z1, double z2,
                                            bool tall_inside, const EdgeCut& cut, int quarters = 4) {
    using namespace fillet_detail;
    BlendedSolid res;
    if (!(r0 > 0.0) || !(r0 < r1) || !(r1 < r2) || !(z0 < z1) || !(z1 < z2) || quarters < 1 || quarters > 4) {
        res.error = "bad stepped ring";
        return res;
    }
    const double a = r2 - r0, b = z2 - z0, a1 = tall_inside ? r1 - r0 : r2 - r1, b1 = z1 - z0;
    if (!check_step(a, b, a1, b1, cut, res.error)) return res;
    return revolve_section(layout_step(a, b, a1, b1, cut), tall_inside ? r0 : r2, z0, tall_inside ? 1.0 : -1.0,
                           1.0, quarters);
}

} // namespace ak::op
