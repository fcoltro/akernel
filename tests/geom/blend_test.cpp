// Exact rolling-ball fillets (geom/blend.hpp): the fillet lies on the
// rolling ball's cylinder or torus, touches each support along its contact
// curve with the support's tangent plane (G1), has the right setback, is
// exact for round data, and refuses what cannot be built. The radius limit
// on a curved support agrees with the certified offset test.
#include <akernel/geom/blend.hpp>
#include <akernel/geom/primitives.hpp>

#include <cmath>

#include "check.hpp"

using namespace ak::geom;

namespace {

Vec3d unit_normal(const NurbsSurface& s, double u, double v) {
    const Jet3 J = s.jet(u, v);
    return normalized(cross(J[1], J[2]));
}

/// |n1 x n2|: zero when the normals are parallel (either orientation).
double normal_gap(const Vec3d& a, const Vec3d& b) { return norm(cross(normalized(a), normalized(b))); }

/// Distance from p to the circle of radius cx at height cz about the z axis,
/// minus the tube radius: zero on the torus.
double torus_gap(const Vec3d& p, double cx, double cz, double r) {
    return std::hypot(std::hypot(p.x, p.y) - cx, p.z - cz) - r;
}

} // namespace

int main() {
    akt::Rng rng(31);

    // --- a box edge: right angle, round data, exact --------------------------
    {
        // Faces x = 0 and y = 0 of the material x, y > 0; edge along z.
        const Fillet f = plane_plane_fillet({0, 0, 0}, {0, 0, 5}, {1, 0, 0}, {0, 1, 0}, 0.5);
        AK_CHECK(f.status == BlendStatus::Ok);
        // Exact controls: contacts, corner and the weight 2.
        AK_CHECK(f.section.contact_a.x == 0.0 && f.section.contact_a.y == 0.5);
        AK_CHECK(f.section.contact_b.x == 0.5 && f.section.contact_b.y == 0.0);
        AK_CHECK(f.section.arc.cp[2].w == 2.0);
        AK_CHECK_NEAR(f.section.setback, 0.5, 0.0);
        for (int n = 0; n < 200; ++n) {
            const Vec3d p = f.surface.point(rng.unit(), rng.unit());
            AK_CHECK_NEAR(std::hypot(p.x - 0.5, p.y - 0.5), 0.5, 1e-15);
            AK_CHECK(p.x <= 0.5 && p.y <= 0.5); // on the corner side of the ball
        }
        // G1 with the supports along the contact lines.
        for (int n = 0; n <= 10; ++n) {
            const double v = n / 10.0;
            AK_CHECK(normal_gap(unit_normal(f.surface, 0.0, v), {1, 0, 0}) < 1e-15);
            AK_CHECK(normal_gap(unit_normal(f.surface, 1.0, v), {0, 1, 0}) < 1e-15);
        }
    }

    // --- random dihedral angles, convex and concave ---------------------------
    for (int trial = 0; trial < 200; ++trial) {
        const Vec3d t = normalized(Vec3d{rng.range(-1, 1), rng.range(-1, 1), rng.range(-1, 1)});
        // Two unit normals perpendicular to t, at a random angle in (5, 175) degrees.
        Vec3d e1 = normalized(cross(t, Vec3d{0.3, -0.7, 0.2}));
        const Vec3d e2 = cross(t, e1);
        const double a = rng.range(0.1, 3.0);
        const Vec3d mA = e1, mB = e1 * std::cos(a) + e2 * std::sin(a);
        const Vec3d E{rng.range(-2, 2), rng.range(-2, 2), rng.range(-2, 2)};
        const double r = rng.range(0.05, 2.0);
        const Fillet f = plane_plane_fillet(E, t * 3.0, mA, mB, r);
        AK_CHECK(f.status == BlendStatus::Ok);
        if (f.status != BlendStatus::Ok) continue;
        const Vec3d c = f.section.centre;
        // Centre at distance r from both planes; contacts on them.
        AK_CHECK_NEAR(dot(c - E, mA), r, 1e-13);
        AK_CHECK_NEAR(dot(c - E, mB), r, 1e-13);
        AK_CHECK(std::fabs(dot(f.section.contact_a - E, mA)) < 1e-13);
        AK_CHECK(std::fabs(dot(f.section.contact_b - E, mB)) < 1e-13);
        // Setback r tan(phi / 2), phi = a the angle the arc turns through.
        AK_CHECK_NEAR(f.section.setback, r * std::tan(0.5 * a), 1e-12);
        for (int n = 0; n < 20; ++n) {
            const double u = rng.unit(), v = rng.unit();
            const Vec3d p = f.surface.point(u, v);
            const Vec3d q = p - c - t * dot(p - c, t); // to the axis, across it
            AK_CHECK_NEAR(norm(q), r, 1e-13 * (1 + r));
        }
        AK_CHECK(normal_gap(unit_normal(f.surface, 0.0, 0.5), mA) < 1e-12);
        AK_CHECK(normal_gap(unit_normal(f.surface, 1.0, 0.5), mB) < 1e-12);
    }

    // --- refusals -------------------------------------------------------------
    AK_CHECK(plane_plane_fillet({0, 0, 0}, {0, 0, 1}, {1, 0, 0}, {-1, 0, 0}, 1).status ==
             BlendStatus::Folded);
    AK_CHECK(plane_plane_fillet({0, 0, 0}, {0, 0, 1}, {1, 0, 0}, {1, 0, 0}, 1).status ==
             BlendStatus::Tangent);

    // --- plane and cylinder: a boss on a plate and a hole through it ----------
    {
        const double R = 2.0, h = 1.0, r = 0.25;
        // Boss x <= R standing on the plate z <= h: concave edge, ball above the
        // plate and outside the boss.
        const Fillet boss = revolved_fillet({R, 0, h}, {0, 0, 1}, {1, 0, 0}, r, 4);
        AK_CHECK(boss.status == BlendStatus::Ok);
        // Exact data: every control is round, the weight is 2.
        AK_CHECK(boss.section.centre.x == R + r && boss.section.centre.z == h + r);
        double worst = 0.0;
        for (const auto& pc : boss.surface.to_bezier()) {
            const BezierSurface bz(pc.patch);
            for (int n = 0; n < 50; ++n) {
                const Vec3d p = bz.point(rng.unit(), rng.unit());
                worst = std::max(worst, std::fabs(torus_gap(p, R + r, h + r, r)));
            }
        }
        AK_CHECK(worst < 1e-14);
        // G1: the plane's normal at the contact circle, the cylinder's radial normal at the other.
        for (int n = 0; n <= 8; ++n) {
            const double u = n / 8.0;
            AK_CHECK(normal_gap(unit_normal(boss.surface, u, 0.0), {0, 0, 1}) < 1e-14);
            const Vec3d p = boss.surface.point(u, 1.0);
            AK_CHECK(normal_gap(unit_normal(boss.surface, u, 1.0), {p.x, p.y, 0}) < 1e-14);
            AK_CHECK_NEAR(std::hypot(p.x, p.y), R, 1e-15);
        }

        // Rim of a hole x >= R through the plate z <= h: convex edge, ball in
        // the material (below the top, outside the hole's radius).
        const Fillet rim = revolved_fillet({R, 0, h}, {0, 0, -1}, {1, 0, 0}, r, 4);
        AK_CHECK(rim.status == BlendStatus::Ok);
        AK_CHECK(rim.section.centre.x == R + r && rim.section.centre.z == h - r);

        // Bottom of a blind hole of radius R: concave edge between the hole's
        // wall and its floor, ball inside the hole. It fits only if r < R.
        const Fillet ok = revolved_fillet({R, 0, 0}, {0, 0, 1}, {-1, 0, 0}, 1.5, 4);
        const Fillet big = revolved_fillet({R, 0, 0}, {0, 0, 1}, {-1, 0, 0}, 2.5, 4);
        AK_CHECK(ok.status == BlendStatus::Ok);
        AK_CHECK(big.status == BlendStatus::CrossesAxis);
        // The certified offset test on the wall says the same: the wall's
        // normal points out of the hole, away from the ball.
        const PolyPatchD wall = make_cylinder_quarter(R, 0.0, 1.0);
        AK_CHECK(support_fits(wall, false, 1.5).verdict == OffsetVerdict::Valid);
        AK_CHECK(support_fits(wall, false, 2.5).verdict == OffsetVerdict::Invalid);
        // Outside a boss the wall is convex towards the ball: any radius fits.
        AK_CHECK(support_fits(wall, true, 50.0).verdict == OffsetVerdict::Valid);
    }

    // --- plane and cone -------------------------------------------------------
    {
        // Cone through (2, 0, 0) and (1, 0, 1) standing on the plane z = 0,
        // material inside; fillet in the concave corner outside the cone.
        const Vec3d E{2, 0, 0};
        const Vec3d along = normalized(Vec3d{-1, 0, 1});
        const Vec3d out{along.z, 0, -along.x}; // meridian normal pointing away from the axis
        const Fillet f = revolved_fillet(E, {0, 0, 1}, out, 0.3, 3);
        AK_CHECK(f.status == BlendStatus::Ok);
        const Vec3d c = f.section.centre;
        AK_CHECK_NEAR(c.z, 0.3, 1e-15);
        AK_CHECK_NEAR(dot(c - E, out), 0.3, 1e-15);
        for (const auto& pc : f.surface.to_bezier()) {
            const BezierSurface bz(pc.patch);
            for (int n = 0; n < 50; ++n)
                AK_CHECK(std::fabs(torus_gap(bz.point(rng.unit(), rng.unit()), c.x, c.z, 0.3)) < 1e-14);
        }
        // The contact circle on the cone lies on the cone x + z = 2 (in the meridian).
        const Vec3d b = f.section.contact_b;
        AK_CHECK_NEAR(b.x + b.z, 2.0, 1e-15);
    }

    // --- plane and sphere -----------------------------------------------------
    {
        // Sphere of radius 2 centred at height 1, cut by the plane z = 0; the
        // ball sits above the plane and outside the sphere.
        const Fillet f = plane_sphere_fillet(0.0, +1, 1.0, 2.0, +1, 0.4, 4);
        AK_CHECK(f.status == BlendStatus::Ok);
        const Vec3d c = f.section.centre;
        AK_CHECK_NEAR(c.z, 0.4, 0.0);
        AK_CHECK_NEAR(norm(c - Vec3d{0, 0, 1}), 2.4, 1e-15);
        AK_CHECK_NEAR(norm(f.section.contact_b - Vec3d{0, 0, 1}), 2.0, 1e-15);
        double worst = 0.0;
        for (const auto& pc : f.surface.to_bezier()) {
            const BezierSurface bz(pc.patch);
            for (int n = 0; n < 50; ++n)
                worst = std::max(worst, std::fabs(torus_gap(bz.point(rng.unit(), rng.unit()), c.x, c.z, 0.4)));
        }
        AK_CHECK(worst < 1e-14);
        // G1 with the sphere at the contact circle.
        for (int n = 0; n <= 8; ++n) {
            const double u = n / 8.0;
            const Vec3d p = f.surface.point(u, 1.0);
            AK_CHECK(normal_gap(unit_normal(f.surface, u, 1.0), p - Vec3d{0, 0, 1}) < 1e-13);
            AK_CHECK(normal_gap(unit_normal(f.surface, u, 0.0), {0, 0, 1}) < 1e-13);
        }
        // A plane far above a small sphere: the ball cannot touch both.
        AK_CHECK(plane_sphere_fillet(5.0, +1, 0.0, 1.0, +1, 0.5, 4).status == BlendStatus::NoBall);
    }

    return akt::finish("blend_test");
}
