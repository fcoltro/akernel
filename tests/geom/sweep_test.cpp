// Generative surfaces (geom/sweep.hpp, geom/nurbs_curve.hpp): the exact
// circle, extrusion, revolution (cylinder, cone, sphere zone, torus), the
// translational sweep, the ruled loft, the skinned loft, and the curve
// utilities they rest on (Bezier extraction, degree elevation, knot merging).
#include <akernel/geom/nurbs_curve.hpp>
#include <akernel/geom/sweep.hpp>

#include <cmath>
#include <vector>

#include "check.hpp"

using namespace ak::geom;

namespace {

double torus_residual(const Vec3d& p, double R, double r) {
    const double q = p.x * p.x + p.y * p.y + p.z * p.z + R * R - r * r;
    return q * q - 4.0 * R * R * (p.x * p.x + p.y * p.y);
}

/// Points of every Bezier piece of `s` satisfy f(p) ~ 0, and the pieces
/// agree with the NURBS evaluation.
template <class F>
void check_on(const NurbsSurface& s, F f, double tol, akt::Rng& rng) {
    for (const auto& piece : s.to_bezier()) {
        const BezierSurface bz(piece.patch);
        for (int n = 0; n < 40; ++n) {
            const double a = rng.unit(), b = rng.unit();
            const Vec3d p = bz.point(a, b);
            AK_CHECK_NEAR(f(p), 0.0, tol);
            const Vec3d q = s.point(piece.u0 + a * (piece.u1 - piece.u0),
                                    piece.v0 + b * (piece.v1 - piece.v0));
            AK_CHECK(norm(p - q) < 1e-12);
        }
    }
}

} // namespace

int main() {
    akt::Rng rng(11);

    // --- exact circle --------------------------------------------------------
    for (int k = 1; k <= 4; ++k) {
        const NurbsCurve c = exact_circle(3.0, 0.5, k);
        AK_CHECK(c.cp.size() == static_cast<std::size_t>(2 * k + 1));
        for (int n = 0; n <= 200; ++n) {
            const Vec3d p = c.point(n / 200.0);
            AK_CHECK_NEAR(p.x * p.x + p.y * p.y, 9.0, 1e-13);
            AK_CHECK(p.z == 0.5);
        }
        // Ends at angle k * 90 degrees.
        const Vec3d e = c.point(1.0);
        const double ang = std::atan2(e.y, e.x);
        const double want = k == 4 ? 0.0 : (k == 3 ? -0.5 : k * 0.5) * 3.14159265358979323846;
        AK_CHECK_NEAR(ang, want, 1e-14);
    }

    // --- curve utilities -----------------------------------------------------
    {
        // Random rational cubic B-spline, 7 controls.
        std::vector<Vec4d> cp;
        for (int i = 0; i < 7; ++i) {
            const double w = rng.range(0.5, 2.0);
            cp.push_back({w * rng.range(-1, 1), w * rng.range(-1, 1), w * rng.range(-1, 1), w});
        }
        const NurbsCurve c(3, open_uniform_knots(7, 3), cp);
        const auto segs = bezier_segments(c);
        AK_CHECK(segs.size() == 4u);
        const NurbsCurve e = elevate_degree(c, 2);
        AK_CHECK(e.p == 5);
        for (int n = 0; n <= 100; ++n) {
            const double t = n / 100.0;
            AK_CHECK(norm(e.point(t) - c.point(t)) < 1e-13);
        }
        // Each Bezier segment reproduces its span.
        const auto spans = std::vector<double>{0, 0.25, 0.5, 0.75, 1};
        for (std::size_t s = 0; s < segs.size(); ++s) {
            const NurbsCurve b = bezier_curve(segs[s]);
            for (int n = 0; n <= 20; ++n) {
                const double a = n / 20.0;
                const double t = spans[s] + a * (spans[s + 1] - spans[s]);
                AK_CHECK(norm(b.point(a) - c.point(t)) < 1e-13);
            }
        }
        // Compatibility with a quadratic on a different domain.
        NurbsCurve d(2, {2, 2, 2, 2.5, 3, 3, 3},
                     {{0, 0, 0, 1}, {1, 2, 0, 1}, {2, -1, 1, 1}, {3, 0, 0, 1}});
        const auto cs = make_compatible({c, d});
        AK_CHECK(cs[0].p == cs[1].p && cs[0].U == cs[1].U);
        AK_CHECK(cs[0].cp.size() == cs[1].cp.size());
        for (int n = 0; n <= 100; ++n) {
            const double t = n / 100.0;
            AK_CHECK(norm(cs[0].point(t) - c.point(t)) < 1e-13);
            AK_CHECK(norm(cs[1].point(t) - d.point(2.0 + t)) < 1e-13);
        }
    }

    // --- extrusion: a quarter of an exact circle gives a cylinder sector --------
    {
        const NurbsSurface s = extrude(exact_circle(2.0, -1.0, 4), {0, 0, 3});
        check_on(s, [](const Vec3d& p) { return p.x * p.x + p.y * p.y - 4.0; }, 1e-13, rng);
        AK_CHECK(s.to_bezier().size() == 4u);
    }

    // --- revolution -----------------------------------------------------------
    {
        // Torus: circle of radius 1 centred at (3, 0, 0) in the xz plane.
        // Profile: a full exact circle in the xz plane, built from the xy one.
        NurbsCurve prof = exact_circle(1.0, 0.0, 4);
        for (auto& c : prof.cp) c = {c.x + 3.0 * c.w, 0.0, c.y, c.w}; // (x, y) -> (x + 3, 0, y)
        const NurbsSurface torus = revolve(prof, 4);
        check_on(torus, [](const Vec3d& p) { return torus_residual(p, 3.0, 1.0); }, 1e-11, rng);
        AK_CHECK(torus.to_bezier().size() == 16u);

        // Pieces of one turn share their boundary rows bitwise (homogeneous).
        for (int j = 0; j < torus.nv; ++j)
            for (int q = 1; q < 4; ++q) {
                const Vec4d a = torus.at(2 * q, j);
                AK_CHECK(a.w > 0.0);
            }
        const auto pieces = torus.to_bezier();
        int shared = 0;
        for (const auto& A : pieces)
            for (const auto& B : pieces)
                if (A.u1 == B.u0 && A.v0 == B.v0)
                    for (int j = 0; j <= A.patch.q; ++j) {
                        const Vec4d a = A.patch.at(A.patch.p, j), b = B.patch.at(0, j);
                        AK_CHECK(a.x == b.x && a.y == b.y && a.z == b.z && a.w == b.w);
                        ++shared;
                    }
        AK_CHECK(shared == 3 * 4 * 3); // 3 interior joints x 4 profile pieces x 3 controls

        // Cone: a line from (1, 0, 0) to (0, 0, 2); x^2 + y^2 = (1 - z/2)^2.
        const NurbsCurve line = polynomial_curve(1, {0, 0, 1, 1}, {{1, 0, 0}, {0, 0, 2}});
        const NurbsSurface cone = revolve(line, 3);
        check_on(cone,
                 [](const Vec3d& p) { return p.x * p.x + p.y * p.y - (1 - p.z / 2) * (1 - p.z / 2); },
                 1e-14, rng);

        // Sphere zone: a quarter meridian arc that stops short of the pole.
        NurbsCurve arc = exact_circle(2.0, 0.0, 1);
        for (auto& c : arc.cp) c = {c.x, 0.0, c.y, c.w}; // (x, y) -> (x, 0, z)
        const NurbsSurface sphere = revolve(arc, 4);
        check_on(sphere, [](const Vec3d& p) { return dot(p, p) - 4.0; }, 1e-13, rng);
    }

    // --- translational sweep ---------------------------------------------------
    {
        // A rational profile (half circle in xz) moved along a cubic path.
        NurbsCurve prof = exact_circle(1.0, 0.0, 2);
        for (auto& c : prof.cp) c = {c.x, 0.0, c.y, c.w};
        const NurbsCurve path = polynomial_curve(
            3, {0, 0, 0, 0, 1, 1, 1, 1}, {{0, 0, 0}, {0, 1, 0.5}, {1, 2, 0}, {1, 3, 1}});
        const NurbsSurface s = translational_sweep(prof, path);
        for (int n = 0; n < 200; ++n) {
            const double u = rng.unit(), v = rng.unit();
            const Vec3d want = prof.point(u) + path.point(v) - path.point(0.0);
            AK_CHECK(norm(s.point(u, v) - want) < 1e-13);
        }
    }

    // --- ruled loft ---------------------------------------------------------------
    {
        const NurbsCurve a = exact_circle(1.0, 0.0, 4);
        const NurbsCurve b = polynomial_curve(
            3, {0, 0, 0, 0, 1, 1, 1, 1}, {{2, 0, 1}, {2, 2, 1}, {-2, 2, 1}, {-2, 0, 1}});
        const NurbsSurface s = ruled(a, b);
        for (int n = 0; n <= 50; ++n) {
            const double t = n / 50.0;
            AK_CHECK(norm(s.point(t, 0.0) - a.point(t)) < 1e-13);
            AK_CHECK(norm(s.point(t, 1.0) - b.point(t)) < 1e-13);
            // Rulings are straight: the midpoint lies on the chord (homogeneous
            // interpolation keeps v-lines straight even for rational rails).
            const Vec3d p0 = s.point(t, 0.0), p1 = s.point(t, 1.0), pm = s.point(t, 0.5);
            AK_CHECK(norm(cross(pm - p0, p1 - p0)) < 1e-12);
        }
    }

    // --- skinned loft through three circles ----------------------------------------
    {
        const std::vector<NurbsCurve> secs = {exact_circle(1.0, 0.0, 4), exact_circle(2.0, 1.0, 4),
                                              exact_circle(1.5, 3.0, 4)};
        const auto r = skin(secs, 2);
        AK_CHECK(r.has_value());
        if (r) {
            AK_CHECK(r->section_params.size() == 3u);
            for (std::size_t k = 0; k < secs.size(); ++k)
                for (int n = 0; n <= 80; ++n) {
                    const double t = n / 80.0;
                    AK_CHECK(norm(r->surface.point(t, r->section_params[k]) - secs[k].point(t)) <
                             1e-12);
                }
            // Same weights in every section: each v-line is the circle scaled
            // and lifted, so every point stays on a circle about z.
            for (int n = 0; n < 100; ++n) {
                const double u = rng.unit(), v = rng.unit();
                const Vec3d p = r->surface.point(u, v), q = r->surface.point(0.0, v);
                AK_CHECK_NEAR(p.x * p.x + p.y * p.y, q.x * q.x + q.y * q.y, 1e-12);
            }
        }
        // Coincident sections are refused, not lofted.
        AK_CHECK(!skin({secs[0], secs[0]}, 1).has_value());
    }

    return akt::finish("sweep_test");
}
