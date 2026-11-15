// Certified curve/surface intersection.
#include <akernel/geom/primitives.hpp>
#include <akernel/solve/csi.hpp>

#include "check.hpp"

using namespace ak;
using geom::Vec3d;
using solve::CsiResult;
using solve::CsiSystem;

namespace {

double hill(double h, double x, double y) { return 16.0 * h * x * (1.0 - x) * y * (1.0 - y); }

} // namespace

int main() {
    const auto plane = geom::make_plane({0, 0, 0}, {1, 0, 0}, {0, 1, 0});

    // 1. A line through a plane: one regular root at a known place.
    {
        const CsiSystem sys(solve::make_line({0.3, 0.4, -1.0}, {0.3, 0.4, 1.0}), plane);
        const CsiResult r = solve::intersect_curve_surface(sys);
        AK_CHECK(r.complete);
        AK_CHECK(r.roots.size() == 1);
        if (r.roots.size() == 1) {
            const auto& root = r.roots[0];
            AK_CHECK(root.interior);
            AK_CHECK(num::contains(root.enc[0], 0.5));
            AK_CHECK(num::contains(root.enc[1], 0.3) || std::fabs(num::mid(root.enc[1]) - 0.3) < 1e-15);
            AK_CHECK(std::fabs(num::mid(root.enc[2]) - 0.4) < 1e-15);
            for (int i = 0; i < 3; ++i) AK_CHECK(num::width(root.enc[i]) < 1e-12);
        }
    }

    // 2. A line that misses: proven empty.
    {
        const CsiSystem sys(solve::make_line({0.1, 0.2, 0.5}, {0.9, 0.7, 0.5}), plane);
        const CsiResult r = solve::intersect_curve_surface(sys);
        AK_CHECK(r.complete);
        AK_CHECK(r.roots.empty());
    }

    // 3. A line through the side of a hill: two roots, both on both objects.
    const double h = 0.5;
    const auto bump = geom::make_bump(h);
    {
        const CsiSystem sys(solve::make_line({0.0, 0.5, 0.25}, {1.0, 0.5, 0.25}), bump);
        const CsiResult r = solve::intersect_curve_surface(sys);
        AK_CHECK(r.complete);
        AK_CHECK(r.roots.size() == 2);
        for (const auto& root : r.roots) {
            AK_CHECK(root.interior);
            const double x = num::mid(root.enc[0]);
            AK_CHECK_NEAR(hill(h, x, 0.5), 0.25, 1e-13);
            AK_CHECK_NEAR(num::mid(root.enc[1]), x, 1e-13);
            for (int i = 0; i < 3; ++i) AK_CHECK(num::width(root.enc[i]) < 1e-11);
        }
    }

    // 4. Random lines against the hill, cross-checked by counting sign
    //    changes of (line height - hill height) on a dense sample.
    {
        akt::Rng rng(20261111);
        int complete = 0, agree = 0, total_roots = 0;
        const int trials = 300;
        for (int trial = 0; trial < trials; ++trial) {
            const Vec3d a{rng.range(0.02, 0.98), rng.range(0.02, 0.98), rng.range(-0.2, 0.7)};
            const Vec3d b{rng.range(0.02, 0.98), rng.range(0.02, 0.98), rng.range(-0.2, 0.7)};
            const CsiSystem sys(solve::make_line(a, b), bump);
            const CsiResult r = solve::intersect_curve_surface(sys);
            if (!r.complete) continue;
            ++complete;
            bool regular = true;
            for (const auto& root : r.roots) regular = regular && root.interior;
            int changes = 0;
            double prev = 0.0;
            const int samples = 20000;
            for (int i = 0; i <= samples; ++i) {
                const double t = static_cast<double>(i) / samples;
                const Vec3d p = a + (b - a) * t;
                const double g = p.z - hill(h, p.x, p.y);
                if (i > 0 && ((prev < 0.0) != (g < 0.0))) ++changes;
                prev = g;
            }
            if (regular && changes == static_cast<int>(r.roots.size())) ++agree;
            total_roots += static_cast<int>(r.roots.size());
        }
        std::printf("random lines vs hill: %d of %d proven complete, %d agree with sampling, %d roots\n",
                    complete, trials, agree, total_roots);
        AK_CHECK(complete >= trials - 3); // a near-tangent line may stay undecided
        AK_CHECK(agree == complete);
        AK_CHECK(total_roots > 100);
    }

    // 5. A line touching the top of the hill: not decided, and the
    //    unresolved region contains the contact point.
    {
        const CsiSystem sys(solve::make_line({0.0, 0.5, h}, {1.0, 0.5, h}), bump);
        solve::CsiOptions opt;
        opt.node_limit = 20000;
        const CsiResult r = solve::intersect_curve_surface(sys, opt);
        AK_CHECK(!r.complete);
        AK_CHECK(r.roots.empty());
        AK_CHECK(r.unverified_boxes > 0);
        for (int i = 0; i < 3; ++i) AK_CHECK(num::contains(r.unverified_hull[i], 0.5));
    }

    // 6. Rational surface: a radial line meets the unit cylinder at radius 1.
    {
        const auto cyl = geom::make_cylinder_quarter(1.0, 0.0, 1.0);
        const CsiSystem sys(solve::make_line({0.0, 0.0, 0.5}, {2.0, 2.0, 0.5}), cyl);
        const CsiResult r = solve::intersect_curve_surface(sys);
        AK_CHECK(r.complete);
        AK_CHECK(r.roots.size() == 1);
        if (r.roots.size() == 1) {
            AK_CHECK(r.roots[0].interior);
            AK_CHECK_NEAR(num::mid(r.roots[0].enc[0]), std::sqrt(0.5) / 2.0, 1e-14);
            AK_CHECK_NEAR(num::mid(r.roots[0].enc[2]), 0.5, 1e-14);
            const auto p = sys.curve_point(r.roots[0].enc[0]);
            AK_CHECK_NEAR(num::mid(p.x) * num::mid(p.x) + num::mid(p.y) * num::mid(p.y), 1.0, 1e-13);
        }
    }

    // 7. A cubic curve against a bicubic patch: every root is a real
    //    intersection, to rounding.
    {
        const auto patch = geom::fit_patch(
            [](double u, double v) { return Vec3d{u, v, 0.3 * u * u * u - 0.2 * u * v + 0.25 * v * v}; }, 3, 3);
        const auto curve = solve::make_curve({{0.1, 0.1, -0.2, 1.0},
                                              {0.4, 0.9, 0.9, 1.0},
                                              {0.6, 0.1, -0.7, 1.0},
                                              {0.9, 0.8, 0.5, 1.0}});
        const CsiSystem sys(curve, patch);
        const CsiResult r = solve::intersect_curve_surface(sys);
        AK_CHECK(r.complete);
        AK_CHECK(!r.roots.empty());
        for (const auto& root : r.roots) {
            const auto g = sys.G_at(root.x);
            for (int i = 0; i < 3; ++i) AK_CHECK(num::mag(g[static_cast<std::size_t>(i)]) < 1e-13);
        }
        std::printf("cubic vs bicubic: %d roots\n", static_cast<int>(r.roots.size()));
    }

    return akt::finish("csi_test");
}
