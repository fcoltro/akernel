// akernel — applying a fillet to a solid: trim the faces, replace the edge
// (docs/FEATURES.md, step F5, first case).
//
// fillet_box_edge() rounds one edge of a box with the exact rolling-ball
// fillet of geom/blend.hpp and returns a closed PatchSolid: every face a whole
// rational patch, every edge one shared definition, so the result is
// watertight by construction and can be an operand of the Booleans.
//
// Box o + x X + y Y + z Z, x in [0, a], y in [0, b], z in [0, c], with (X, Y,
// Z) orthonormal and right-handed. The filleted edge is x = a, y = b, along
// Z. The supports are trimmed along the contact lines: the face x = a to
// y <= b - r, the face y = b to x <= a - r; the fillet (a cylinder sector)
// replaces the edge between them. Each end face (z = 0, z = c) becomes the
// rectangle with its corner rounded.
//
// Why the end faces are split. The end face's boundary is tangent-continuous
// where the arc meets the trimmed lines, so a single patch would have a
// corner with a straight angle there and a singular parametrisation. The
// end face is instead split into five regular patches (one cell layout,
// section coordinates, s = r / 2):
//
//      Y ----------------- B .
//      |        Qd         |    .  arc (the fillet's end)
//     Yl ---------------- C1      .
//      |                   | Qcor   A
//      |        Qc         C2 ---/  |
//      |                   |   \    |
//      |                   | Qb |Qa |
//      O ---------------- Xr - Xs - X
//
//   C1 = (a - r, b - 3s), C2 = (a - s, b - 3s): the corner patch Qcor is
//   ruled between the arc A -> B and the segment C2 -> C1, with interior
//   angles below 180 degrees at A and B. Cell edges continue across the side
//   faces, which are split to match (y = 0 at x = a - r and a - s, x = 0 at
//   y = b - 3s).
//
// Weights. The exact arc has weights (1, 1, 2) for a right angle, contact A
// weight 1 and contact B weight 2 (geom/blend.hpp: no square roots, every
// control round for round data). Sides shared by two patches must have the
// same homogeneous controls, so the patches on B's side carry weight 2
// throughout (a uniform factor does not change a patch), and the two
// patches next to x = a carry weight 1 on that side and 2 on the other.
// With dyadic a, b, c, r and axis-aligned (X, Y, Z) every control is exact.
#pragma once

#include <optional>
#include <string>
#include <vector>

#include <akernel/geom/blend.hpp>
#include <akernel/op/patch_solid.hpp>

namespace ak::op {

namespace fillet_detail {

using geom::Vec4d;

inline Vec4d hom(const Vec3d& p, double w) { return {p.x * w, p.y * w, p.z * w, w}; }
inline Vec3d point(const Vec4d& h) { return {h.x / h.w, h.y / h.w, h.z / h.w}; }

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

} // namespace fillet_detail

struct FilletedSolid {
    std::optional<PatchSolid> solid;
    geom::Fillet fillet; ///< the blend surface and its section
    std::string error;
};

/// Round the edge x = a, y = b (along Z) of the box o + x X + y Y + z Z with
/// sides a = |ex|, b = |ey|, c = |ez| by the rolling-ball fillet of radius r.
/// ex, ey, ez must be orthogonal and right-handed. Needs r < a and r < 2b/3
/// (room for the end-face layout above).
inline FilletedSolid fillet_box_edge(const Vec3d& o, const Vec3d& ex, const Vec3d& ey, const Vec3d& ez,
                                     double r) {
    using namespace fillet_detail;
    FilletedSolid res;
    const double a = norm(ex), b = norm(ey), c = norm(ez);
    if (!(r > 0.0) || !(r < a) || !(1.5 * r < b)) {
        res.error = "radius out of range for this box";
        return res;
    }
    const Vec3d X = ex / a, Y = ey / b, Z = ez / c;
    if (!(dot(cross(X, Y), Z) > 0.0)) {
        res.error = "the box axes must be right-handed";
        return res;
    }
    auto pt = [&](double x, double y, double z) { return o + X * x + Y * y + Z * z; };

    // The fillet: material x <= a, y <= b, so the ball side of each face is
    // into the material, -X and -Y.
    res.fillet = geom::plane_plane_fillet(pt(a, b, 0.0), ez, X * -1.0, Y * -1.0, r);
    if (res.fillet.status != geom::BlendStatus::Ok) {
        res.error = "no fillet of this radius";
        return res;
    }
    const auto pieces = res.fillet.surface.to_bezier();
    if (pieces.size() != 1) {
        res.error = "unexpected fillet surface";
        return res;
    }
    const PolyPatchD& F = pieces[0].patch; // u across (A at u = 0), v along Z

    const double s = 0.5 * r;
    PatchSolid S;
    std::vector<Vec4d> arc[2], line[2];
    // Per end layer k (z = 0, z = c): the cell points of the layout above.
    struct Layer {
        Vec4d A1, B2;                    // contacts, from the fillet (weights 1 and 2)
        Vec3d O, X_, Y_, Xr, Xs, C1, C2, Yl;
    } L[2];
    for (int k = 0; k < 2; ++k) {
        const double z = k ? c : 0.0;
        Layer& l = L[k];
        l.A1 = F.at(0, k);
        l.B2 = F.at(2, k);
        l.O = pt(0, 0, z);
        l.X_ = pt(a, 0, z);
        l.Y_ = pt(0, b, z);
        l.Xr = pt(a - r, 0, z);
        l.Xs = pt(a - s, 0, z);
        l.C1 = pt(a - r, b - 3 * s, z);
        l.C2 = pt(a - s, b - 3 * s, z);
        l.Yl = pt(0, b - 3 * s, z);
        arc[k] = {F.at(0, k), F.at(1, k), F.at(2, k)};
        line[k] = {hom(l.C2, 2), hom((l.C2 + l.C1) * 0.5, 2), hom(l.C1, 2)};
    }
    auto w2 = [](const Vec3d& p) { return hom(p, 2.0); };
    auto w1 = [](const Vec3d& p) { return hom(p, 1.0); };

    // End faces.
    for (int k = 0; k < 2; ++k) {
        const Layer& l = L[k];
        const Vec3d out = k ? Z : Z * -1.0;
        // Corner patch: arc A -> B, then C1 -> C2 back.
        S.faces.push_back(outward(ruled(arc[k], line[k]), out));
        // Qa: X, A (weight 1) to Xs, C2 (weight 2).
        S.faces.push_back(outward(ruled({w1(l.X_), l.A1}, {w2(l.Xs), w2(l.C2)}), out));
        // Qb: between x = a - r and a - s, degree 2 along x to match the corner's segment.
        S.faces.push_back(outward(ruled({w2(l.Xs), w2((l.Xr + l.Xs) * 0.5), w2(l.Xr)}, line[k]), out));
        // Qc and Qd.
        S.faces.push_back(outward(ruled({w2(l.O), w2(l.Xr)}, {w2(l.Yl), w2(l.C1)}), out));
        S.faces.push_back(outward(ruled({w2(l.Yl), w2(l.C1)}, {w2(l.Y_), l.B2}), out));
    }
    // Side faces between the two layers.
    const Layer &l0 = L[0], &l1 = L[1];
    S.faces.push_back(outward(F, X + Y));                                                  // the fillet
    S.faces.push_back(outward(ruled({w1(l0.X_), l0.A1}, {w1(l1.X_), l1.A1}), X));         // x = a, trimmed
    S.faces.push_back(outward(ruled({w2(l0.Y_), l0.B2}, {w2(l1.Y_), l1.B2}), Y));         // y = b, trimmed
    S.faces.push_back(outward(ruled({w2(l0.O), w2(l0.Xr)}, {w2(l1.O), w2(l1.Xr)}), Y * -1.0));
    S.faces.push_back(outward(ruled({w2(l0.Xr), w2((l0.Xr + l0.Xs) * 0.5), w2(l0.Xs)},
                                    {w2(l1.Xr), w2((l1.Xr + l1.Xs) * 0.5), w2(l1.Xs)}),
                              Y * -1.0));
    S.faces.push_back(outward(ruled({w2(l0.Xs), w1(l0.X_)}, {w2(l1.Xs), w1(l1.X_)}), Y * -1.0));
    S.faces.push_back(outward(ruled({w2(l0.O), w2(l0.Yl)}, {w2(l1.O), w2(l1.Yl)}), X * -1.0));
    S.faces.push_back(outward(ruled({w2(l0.Yl), w2(l0.Y_)}, {w2(l1.Yl), w2(l1.Y_)}), X * -1.0));

    std::string err;
    if (!S.build(err)) {
        res.error = "the filleted shell does not close: " + err;
        return res;
    }
    res.solid = std::move(S);
    return res;
}

} // namespace ak::op
