// Certified surface/surface intersection: curves, certificates, completeness.
#include <akernel/geom/primitives.hpp>
#include <akernel/solve/ssi.hpp>

#include <chrono>

#include "check.hpp"

using namespace ak::geom;
using namespace ak::solve;
using ak::num::Interval;

namespace {

const char* kind(EndKind k) {
    return k == EndKind::Boundary ? "boundary" : (k == EndKind::Closed ? "closed" : "singular");
}

struct Report {
    SsiResult res;
    double max_gap = 0.0;   // certified upper bound on |S1 - S2| at sampled curve points
    double max_width = 0.0; // widest certified parameter enclosure
    bool all_points_ok = true;
};

template <class OnPoint>
Report run(const char* name, const PolyPatchD& A, const PolyPatchD& B, OnPoint&& on_point,
           SsiOptions opt = {}) {
    const SsiSystem sys(A, B);
    const auto t0 = std::chrono::steady_clock::now();
    Report rep;
    rep.res = intersect(sys, opt);
    const double ms =
        std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();

    std::size_t segs = 0;
    for (const auto& c : rep.res.curves) {
        segs += c.segs.size();
        for (const auto& s : c.segs)
            for (double f : {0.0, 0.5, 1.0}) {
                Box4 enc;
                if (!curve_point(sys, s, f, enc)) { rep.all_points_ok = false; continue; }
                const auto pa = sys.point_a(enc), pb = sys.point_b(enc);
                const Interval d2 = sqr(pa.x - pb.x) + sqr(pa.y - pb.y) + sqr(pa.z - pb.z);
                rep.max_gap = std::max(rep.max_gap, std::sqrt(d2.hi));
                for (const auto& e : enc) rep.max_width = std::max(rep.max_width, width(e));
                on_point(pa, pb);
            }
    }
    std::printf("  %-28s curves %zu, segments %zu, complete %s, unverified %zu, "
                "gap <= %.1e, steps %zu/%zu rejected, verify nodes %zu, %.0f ms\n",
                name, rep.res.curves.size(), segs, rep.res.complete ? "yes" : "NO",
                rep.res.unverified_boxes, rep.max_gap, rep.res.stats.steps_rejected,
                rep.res.stats.steps_accepted + rep.res.stats.steps_rejected,
                rep.res.stats.verify_nodes, ms);
    for (const auto& c : rep.res.curves)
        std::printf("      curve: %zu segments, ends %s/%s, min sin(angle) %.3f\n", c.segs.size(),
                    kind(c.end_front), kind(c.end_back), c.min_sin_angle);
    return rep;
}

const auto no_check = [](const Vec3<Interval>&, const Vec3<Interval>&) {};

} // namespace

int main() {
    const PolyPatchD plane_z_half =
        make_plane({-0.5, -0.5, 0.5}, {2, 0, 0}, {0, 2, 0});

    // 1. Plane through a quarter cylinder: a quarter ellipse, boundary to boundary.
    {
        const PolyPatchD cyl = make_cylinder_quarter(1.0, 0.0, 2.0);
        const PolyPatchD plane = make_plane({-0.5, -0.5, 0.35}, {2, 0, 0.6}, {0, 2, 0}); // z = 0.5 + 0.3 x
        const Report r = run("plane x cylinder", cyl, plane,
            [](const Vec3<Interval>& pa, const Vec3<Interval>&) {
                // The certified point lies on the cylinder and on the plane.
                const Interval cylinder = sqr(pa.x) + sqr(pa.y) - Interval(1.0);
                const Interval plane_eq = pa.z - Interval(0.5) - Interval(0.3) * pa.x;
                AK_CHECK(contains_zero(cylinder) && width(cylinder) < 1e-11);
                AK_CHECK(std::fabs(mid(plane_eq)) < 1e-12);
            });
        AK_CHECK(r.res.curves.size() == 1u);
        AK_CHECK(r.res.complete);
        AK_CHECK(r.all_points_ok);
        AK_CHECK(r.max_gap < 1e-11);
        if (!r.res.curves.empty()) {
            AK_CHECK(!r.res.curves[0].closed);
            AK_CHECK(r.res.curves[0].end_front == EndKind::Boundary);
            AK_CHECK(r.res.curves[0].end_back == EndKind::Boundary);
        }
    }

    // 2. Plane through a hill: one closed loop strictly inside both domains.
    {
        const Report r = run("plane x hill (loop)", make_bump(1.0), plane_z_half,
            [](const Vec3<Interval>& pa, const Vec3<Interval>&) {
                AK_CHECK(std::fabs(mid(pa.z) - 0.5) < 1e-12);
            });
        AK_CHECK(r.res.curves.size() == 1u);
        AK_CHECK(r.res.complete);
        AK_CHECK(r.all_points_ok);
        if (!r.res.curves.empty()) AK_CHECK(r.res.curves[0].closed);
    }

    // 3. No intersection: proven empty.
    {
        const PolyPatchD high = make_plane({-0.5, -0.5, 2.0}, {2, 0, 0}, {0, 2, 0});
        const Report r = run("plane above hill (empty)", make_bump(1.0), high, no_check);
        AK_CHECK(r.res.curves.empty());
        AK_CHECK(r.res.complete);
    }

    // 4. Two planes: a straight line.
    {
        const PolyPatchD p1 = make_plane({0, 0, 0}, {1, 0, 0}, {0, 1, 0});
        const PolyPatchD p2 = make_plane({0.3, -0.2, -0.5}, {0.1, 1.3, 0.2}, {0.2, 0.1, 1.0});
        const Report r = run("plane x plane", p1, p2, no_check);
        AK_CHECK(r.res.curves.size() == 1u);
        AK_CHECK(r.res.complete);
        AK_CHECK(r.max_gap < 1e-12);
    }

    // 5. Two cylinders of different radius with perpendicular axes.
    {
        const PolyPatchD c1 = make_cylinder_quarter(1.0, 0.0, 2.0);
        const double swap_xz[9] = {0, 0, 1, 0, 1, 0, 1, 0, 0};
        const PolyPatchD c2 = transformed(make_cylinder_quarter(0.8, -1.0, 2.0), swap_xz, {0, 0, 0});
        const Report r = run("cylinder x cylinder", c1, c2,
            [](const Vec3<Interval>& pa, const Vec3<Interval>& pb) {
                const Interval e1 = sqr(pa.x) + sqr(pa.y) - Interval(1.0);
                const Interval e2 = sqr(pb.y) + sqr(pb.z) - Interval(0.64);
                AK_CHECK(contains_zero(e1) && width(e1) < 1e-11);
                AK_CHECK(std::fabs(mid(e2)) < 1e-11);
            });
        AK_CHECK(r.res.curves.size() == 1u);
        AK_CHECK(r.res.complete);
        AK_CHECK(r.all_points_ok);
    }

    // 6. Two freeform bicubic patches.
    {
        const PolyPatchD a = fit_patch([](double u, double v) {
            return Vec3d{u, v, 0.3 * u * u * u - 0.4 * u * v + 0.2 * v * v + 0.1};
        }, 3, 3);
        const PolyPatchD b = fit_patch([](double s, double t) {
            return Vec3d{0.1 + 0.8 * s + 0.1 * t * t, 0.9 * t + 0.05, 0.6 * s - 0.2 * t * s * s + 0.05 * t};
        }, 3, 3);
        const Report r = run("bicubic x bicubic", a, b, no_check);
        AK_CHECK(!r.res.curves.empty());
        AK_CHECK(r.res.complete);
        AK_CHECK(r.all_points_ok);
        AK_CHECK(r.max_gap < 1e-11);
    }

    // 7. Tangential contact: the plane touches the top of the hill in one
    // point. No certificate exists, and the kernel must say so instead of
    // returning "no intersection".
    {
        const PolyPatchD top = make_plane({-0.5, -0.5, 1.0}, {2, 0, 0}, {0, 2, 0});
        SsiOptions opt;
        opt.verify_depth = 24;
        opt.verify_node_limit = 200000;
        const Report r = run("plane tangent to hill", make_bump(1.0), top, no_check, opt);
        AK_CHECK(!r.res.complete);
        AK_CHECK(r.res.unverified_boxes > 0);
        // The unresolved region contains the true contact point (0.5, 0.5).
        AK_CHECK(contains(r.res.unverified_hull[0], 0.5) && contains(r.res.unverified_hull[1], 0.5));
        std::printf("      unresolved region: u in [%.3f, %.3f], v in [%.3f, %.3f]\n",
                    r.res.unverified_hull[0].lo, r.res.unverified_hull[0].hi,
                    r.res.unverified_hull[1].lo, r.res.unverified_hull[1].hi);
    }

    // 8 and 9. Regressions from the benchmark's hard family: two random
    // graphs over the same square (pairs number 1 and 139 of seed 4242).
    //   - Pair 1: a curve reaches the domain boundary with a very short last
    //     step. Its certified region used to be too small for the
    //     completeness pass, leaving one unresolved box at the end point.
    //   - Pair 139: the surfaces come within 1e-3 of touching (sine of the
    //     angle between normals 0.0007 there) without intersecting at that
    //     spot. Axis-aligned boxes could not separate them; the normal-
    //     direction test can.
    {
        akt::Rng rng(4242);
        for (int trial = 0; trial < 140; ++trial) {
            PolyPatchD a(3, 3), b(3, 3);
            for (int i = 0; i <= 3; ++i)
                for (int j = 0; j <= 3; ++j) {
                    a.at(i, j) = {i / 3.0, j / 3.0, rng.range(-0.4, 0.4), 1.0};
                    b.at(i, j) = {i / 3.0, j / 3.0, rng.range(-0.4, 0.4), 1.0};
                }
            if (trial != 1 && trial != 139) continue;
            const Report r = run(trial == 1 ? "graphs, short step to edge" : "graphs, near-tangent", a, b, no_check);
            AK_CHECK(r.res.complete);
            AK_CHECK(r.res.curves.size() == 2u);
            AK_CHECK(r.all_points_ok);
            AK_CHECK(r.max_gap < 1e-10);
        }
    }

    return akt::finish("ssi_test");
}
