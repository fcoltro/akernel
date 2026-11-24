// akernel — constant-radius rolling-ball fillets whose surfaces are exact
// rational patches: between two planes (a cylinder sector), and between
// surfaces of revolution about one axis, planes perpendicular to it,
// cylinders, cones and spheres (a torus zone). docs/FEATURES.md, step F3a.
//
// The rolling ball. A ball of radius r touching two surfaces has its centre
// on both offsets at distance r (on the ball's side), and touches each
// surface at the foot of the centre. Sweeping it along the edge, the
// fillet is the envelope of the ball: in a plane section across the edge,
// a circular arc of radius r from contact A to contact B. With unit
// vectors u = (A - c) / r and v = (B - c) / r and k = u . v, the arc is the
// rational quadratic with controls
//     A (weight 1),   c + r (u + v) / (1 + k) (weight 1),   B (weight 2 / (1 + k)),
// the middle control being where the tangents at A and B meet (for two
// planes, the edge itself). The weights are those of the standard circular
// arc (1, cos(phi/2), 1) moved to a common factor: w1^2 / (w0 w2) =
// (1 + k) / 2 = cos^2(phi/2), which is the shape invariant of a conic arc.
// There is no square root in the weights, so with round data (right-angle
// edges, dyadic radii) every control is exact and the fillet lies on its
// cylinder or torus exactly, like the quadric patches of the Booleans.
//
// Revolved fillets are this arc in the meridian half-plane (x >= 0, y = 0,
// z along the axis), revolved by geom/sweep.hpp's exact quarter turns.
//
// What is checked here: existence (the ball fits between the surfaces in
// the section), the arc not crossing the axis, and a degenerate angle.
// Whether the ball fits on each curved support along the whole edge is the
// offset-validity predicate (geom/offset_validity.hpp), applied by
// support_fits() below. Trimming the supports and building the blend's
// topology are later steps (FEATURES F5, F6).
#pragma once

#include <cmath>

#include <akernel/geom/nurbs_curve.hpp>
#include <akernel/geom/offset_validity.hpp>
#include <akernel/geom/sweep.hpp>

namespace ak::geom {

enum class BlendStatus {
    Ok,
    Folded,      ///< the two contact directions are opposite (k = -1): no arc
    Tangent,     ///< the surfaces are tangent along the edge (k = 1): nothing to fillet
    NoBall,      ///< no ball of this radius touches both surfaces on the requested sides
    CrossesAxis, ///< a revolved fillet whose profile reaches the axis: not a regular surface
};

struct FilletArc {
    BlendStatus status = BlendStatus::NoBall;
    NurbsCurve arc;      ///< from contact A (t = 0) to contact B (t = 1)
    Vec3d centre{};
    Vec3d contact_a{}, contact_b{};
    Vec3d corner{};      ///< middle control: where the tangents at the contacts meet
    double setback = 0;  ///< |corner - contact_a| = |corner - contact_b| = r tan(phi / 2)
};

/// Circular arc of radius r about `centre` from centre + r u to centre + r v
/// (u, v unit), the shorter way round.
inline FilletArc fillet_arc(const Vec3d& centre, double r, const Vec3d& u, const Vec3d& v) {
    FilletArc f;
    const double k = dot(u, v);
    if (!(1.0 + k > 0.0)) { f.status = BlendStatus::Folded; return f; }
    if (!(1.0 - k > 0.0)) { f.status = BlendStatus::Tangent; return f; }
    f.centre = centre;
    f.contact_a = centre + u * r;
    f.contact_b = centre + v * r;
    f.corner = centre + (u + v) * (r / (1.0 + k));
    f.setback = norm(f.corner - f.contact_a);
    const double w2 = 2.0 / (1.0 + k);
    const Vec3d& a = f.contact_a;
    const Vec3d& m = f.corner;
    const Vec3d& b = f.contact_b;
    f.arc = bezier_curve({{a.x, a.y, a.z, 1.0}, {m.x, m.y, m.z, 1.0},
                          {b.x * w2, b.y * w2, b.z * w2, w2}});
    f.status = BlendStatus::Ok;
    return f;
}

struct Fillet {
    BlendStatus status = BlendStatus::NoBall;
    NurbsSurface surface; ///< u across the fillet (contact A at u = 0), v along the edge
    FilletArc section;    ///< the arc in the first section (v = 0)
};

/// Fillet of radius r along the straight edge from E to E + t between two
/// planes through it. mA and mB are normals of the planes pointing to the
/// side the ball is on: into the material for a convex edge (the fillet
/// removes material), into the empty side for a concave one (it adds).
/// They are normalised here; with axis-aligned normals and round E, t, r
/// the surface is an exact cylinder sector.
inline Fillet plane_plane_fillet(const Vec3d& E, const Vec3d& t, Vec3d mA, Vec3d mB, double r) {
    mA = normalized(mA);
    mB = normalized(mB);
    Fillet f;
    const double k = dot(mA, mB);
    if (!(1.0 + k > 0.0)) { f.status = BlendStatus::Folded; return f; }
    if (!(1.0 - k > 0.0)) { f.status = BlendStatus::Tangent; return f; }
    // Centre at distance r from both planes on the ball side, in the
    // section through E: c = E + x with x . mA = x . mB = r, x in span(mA, mB).
    const Vec3d c = E + (mA + mB) * (r / (1.0 + k));
    f.section = fillet_arc(c, r, -mA, -mB);
    f.status = f.section.status;
    if (f.status == BlendStatus::Ok) f.surface = extrude(f.section.arc, t);
    return f;
}

/// Fillet about the z axis between two supports whose meridians are
/// straight lines through the meridian point E = (x, 0, z): a plane z =
/// const has the meridian direction (1, 0, 0), a cylinder x = R the
/// direction (0, 0, 1), a cone any other. mA and mB are the meridian
/// normals (y = 0) pointing to the ball's side. The result is a torus zone
/// over `quarters` quarter turns, exact when the data are round.
inline Fillet revolved_fillet(const Vec3d& E, Vec3d mA, Vec3d mB, double r, int quarters) {
    mA = normalized(mA);
    mB = normalized(mB);
    Fillet f;
    const double k = dot(mA, mB);
    if (!(1.0 + k > 0.0)) { f.status = BlendStatus::Folded; return f; }
    if (!(1.0 - k > 0.0)) { f.status = BlendStatus::Tangent; return f; }
    const Vec3d c = E + (mA + mB) * (r / (1.0 + k));
    f.section = fillet_arc(c, r, -mA, -mB);
    f.status = f.section.status;
    if (f.status != BlendStatus::Ok) return f;
    // The arc lies in the convex hull of its controls (positive weights):
    // if they are all off the axis, so is the arc.
    if (!(f.section.contact_a.x > 0.0 && f.section.corner.x > 0.0 && f.section.contact_b.x > 0.0)) {
        f.status = BlendStatus::CrossesAxis;
        return f;
    }
    f.surface = revolve(f.section.arc, quarters);
    return f;
}

/// Fillet about the z axis between the plane z = h and the sphere of radius
/// Rs centred at (0, 0, zs). plane_side = +1 puts the ball above the plane
/// (z > h), -1 below; sphere_side = +1 outside the sphere, -1 inside. The
/// ball centre is at height h + plane_side r and at distance
/// Rs + sphere_side r from the sphere's centre; of the two such centres the
/// one with x > 0 is taken. The square root makes the data inexact in
/// general: the rounded controls are the fillet's definition.
inline Fillet plane_sphere_fillet(double h, int plane_side, double zs, double Rs, int sphere_side,
                                  double r, int quarters) {
    Fillet f;
    const double cz = h + plane_side * r;
    const double dist = Rs + sphere_side * r;
    const double dz = cz - zs;
    if (!(dist > 0.0) || !(dist * dist - dz * dz > 0.0)) { f.status = BlendStatus::NoBall; return f; }
    const Vec3d c{std::sqrt(dist * dist - dz * dz), 0.0, cz};
    const Vec3d u{0.0, 0.0, -static_cast<double>(plane_side)}; // towards the plane
    const Vec3d S{0.0, 0.0, zs};
    // Towards the sphere: away from its centre if inside, towards it if outside.
    const Vec3d v = (S - c) * (static_cast<double>(sphere_side) / dist);
    f.section = fillet_arc(c, r, u, v);
    f.status = f.section.status;
    if (f.status != BlendStatus::Ok) return f;
    if (!(f.section.contact_a.x > 0.0 && f.section.corner.x > 0.0 && f.section.contact_b.x > 0.0)) {
        f.status = BlendStatus::CrossesAxis;
        return f;
    }
    f.surface = revolve(f.section.arc, quarters);
    return f;
}

/// Does a ball of radius r fit on a support patch everywhere? The ball's
/// centre runs on the patch's offset at distance r on the ball's side, which
/// must be regular: 1 - d kappa_i > 0 with d = +r if the patch normal
/// (S_u x S_v) points to the ball, -r otherwise. Valid, Invalid with a
/// witness, or Undecided (geom/offset_validity.hpp).
inline OffsetCheck support_fits(const PolyPatchD& support, bool normal_towards_ball, double r) {
    return check_offset(support, normal_towards_ball ? r : -r);
}

} // namespace ak::geom
