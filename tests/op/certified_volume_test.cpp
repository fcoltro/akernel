// Certified volume enclosures (op/certified_volume.hpp).
//
// Each enclosure must contain the closed-form volume, be narrow (rounding
// only, under 1e-12: no quadrature error), agree with the double-precision measurement
// op::volume, and be refused, with a reason, outside the covered class.
#include <akernel/num/detmath.hpp>
#include <akernel/op/boolean.hpp>
#include <akernel/op/certified_volume.hpp>
#include <akernel/op/surfaces.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>

#include "check.hpp"

using namespace ak;
using geom::Vec3d;
using op::Brep;

namespace {

Brep box(double x0, double y0, double z0, double x1, double y1, double z1) {
    return op::to_brep(op::make_box({x0, y0, z0}, {x1, y1, z1}));
}

/// Encloses x, allowing for the rounding of a closed form computed in doubles.
bool encloses_value(const num::Interval& e, double x) {
    const double slack = 4e-16 * std::max(1.0, std::fabs(x));
    return e.lo <= x + slack && x - slack <= e.hi;
}

double g_widest = 0.0;
int g_certified = 0;

void check_closed_form(const char* name, const Brep& b, double exact) {
    std::string why;
    const auto v = op::certified_volume(b, &why);
    AK_CHECK(v.has_value());
    if (!v) {
        std::printf("FAIL %s: not certified (%s)\n", name, why.c_str());
        return;
    }
    ++g_certified;
    g_widest = std::max(g_widest, num::width(*v));
    std::printf("  %s: [%.17g, %.17g], width %.1e, closed form %.15f\n", name, v->lo, v->hi, num::width(*v), exact);
    AK_CHECK(encloses_value(*v, exact));
    AK_CHECK(num::width(*v) < 1e-12);
}

} // namespace

int main() {
    // Solids built directly.
    check_closed_form("unit cube", box(0, 0, 0, 1, 1, 1), 1.0);
    check_closed_form("offset box", box(-0.3, 0.25, 1.5, 0.7, 2.0, 2.125), 1.0 * 1.75 * 0.625);
    {
        const Vec3d o{0.1, -0.2, 0.3}, ex{1.0, 0.5, 0.0}, ey{-0.25, 0.75, 0.125}, ez{0.0, 0.25, 1.5};
        const double det = ex.x * (ey.y * ez.z - ey.z * ez.y) - ex.y * (ey.x * ez.z - ey.z * ez.x) + ex.z * (ey.x * ez.y - ey.y * ez.x);
        check_closed_form("parallelepiped", op::to_brep(op::make_box(o, ex, ey, ez)), det);
    }
    check_closed_form("box with a hill on top (degree 2 top)", op::to_brep(op::make_bump_box({0, 0, 0}, 1, 1, 1, 0.3)),
                      1.0 + 0.3 * 4.0 / 9.0);
    check_closed_form("larger hill box", op::to_brep(op::make_bump_box({-1, 0.5, 0}, 2, 0.5, 0.75, 0.25)),
                      2 * 0.5 * (0.75 + 0.25 * 4.0 / 9.0));

    // Booleans of flat-faced solids, including non-manifold results.
    {
        const Brep A = box(0, 0, 0, 1, 1, 1);
        struct Case { const char* name; Brep b; op::Op o; double v; };
        const Case cases[] = {
            {"cube - corner box", box(0.5, 0.5, 0.5, 1.5, 1.5, 1.5), op::Op::Difference, 0.875},
            {"cube + corner box", box(0.5, 0.5, 0.5, 1.5, 1.5, 1.5), op::Op::Union, 1.875},
            {"cube & corner box", box(0.5, 0.5, 0.5, 1.5, 1.5, 1.5), op::Op::Intersection, 0.125},
            {"cube - slot", box(0.25, -0.5, 0.75, 0.5, 1.5, 1.5), op::Op::Difference, 1.0 - 0.25 * 0.25},
            {"cube + box along an edge (non-manifold)", box(1, 1, 0, 2, 2, 1), op::Op::Union, 2.0},
            {"cube + box at a corner (non-manifold)", box(1, 1, 1, 2, 2, 2), op::Op::Union, 2.0},
        };
        for (const Case& c : cases) {
            const op::BooleanResult r = op::boolean_op(A, c.b, c.o);
            AK_CHECK(r.ok);
            if (!r.ok) {
                std::printf("FAIL %s: %s\n", c.name, r.error.c_str());
                continue;
            }
            check_closed_form(c.name, r.solid, c.v);
        }
    }

    // Random rotated boxes: Booleans with the unit cube. Every result has
    // flat faces only, so every one must be certified; the enclosure must
    // contain the double-precision measurement up to its quadrature error.
    {
        akt::Rng rng(25112026);
        const Brep A = box(0, 0, 0, 1, 1, 1);
        int results = 0, outside = 0;
        double worst_gap = 0.0;
        for (int trial = 0; trial < 60; ++trial) {
            const double ang = rng.range(0.2, 2.9);
            const double c = num::det::cos(ang), s = num::det::sin(ang);
            const double sx = rng.range(0.4, 1.1), sy = rng.range(0.4, 1.1), sz = rng.range(0.4, 1.1);
            const double ox = rng.range(-0.2, 0.8), oy = rng.range(-0.2, 0.8), oz = rng.range(-0.2, 0.8);
            // Rotation about z by ang, then about x by ang / 2.
            const double c2 = num::det::cos(0.5 * ang), s2 = num::det::sin(0.5 * ang);
            auto rot = [&](const Vec3d& p) {
                const Vec3d q{c * p.x - s * p.y, s * p.x + c * p.y, p.z};
                return Vec3d{q.x, c2 * q.y - s2 * q.z, s2 * q.y + c2 * q.z};
            };
            // Coordinates on a 2^-20 grid: the corner sums are then exact, so
            // each face is exactly flat and the intersection edges between
            // faces are exactly straight. (With unrounded coordinates the
            // faces are slightly twisted bilinear patches: see the refusals.)
            auto q = [](const Vec3d& p) {
                auto r = [](double x) { return std::ldexp(std::nearbyint(std::ldexp(x, 20)), -20); };
                return Vec3d{r(p.x), r(p.y), r(p.z)};
            };
            const Brep B =
                op::to_brep(op::make_box(q({ox, oy, oz}), q(rot({sx, 0, 0})), q(rot({0, sy, 0})), q(rot({0, 0, sz}))));
            const op::Op o = trial % 3 == 0 ? op::Op::Union : trial % 3 == 1 ? op::Op::Intersection : op::Op::Difference;
            const op::BooleanResult r = op::boolean_op(A, B, o);
            if (!r.ok || r.empty) continue;
            ++results;
            std::string why;
            const auto v = op::certified_volume(r.solid, &why);
            AK_CHECK(v.has_value());
            if (!v) {
                std::printf("FAIL random trial %d: %s\n", trial, why.c_str());
                continue;
            }
            ++g_certified;
            g_widest = std::max(g_widest, num::width(*v));
            const double m = op::volume(r.solid);
            const double gap = m < v->lo ? v->lo - m : m > v->hi ? m - v->hi : 0.0;
            worst_gap = std::max(worst_gap, gap);
            if (gap > 1e-12) ++outside;
        }
        std::printf("  random rotated boxes: %d results certified; double measurement outside the enclosure by at most %.1e\n",
                    results, worst_gap);
        AK_CHECK(results >= 50);
        AK_CHECK(outside == 0);
    }

    // Outside the covered class: refused with a reason, never a wrong bound.
    {
        std::string why;
        const Brep sector = op::to_brep(op::make_cylinder_sector(1.0, 2.0, 0.0, 1.0, 0.0, 1.0));
        AK_CHECK(!op::certified_volume(sector, &why).has_value());
        std::printf("  cylinder sector: not certified (%s)\n", why.c_str());
        AK_CHECK(why == "a rational carrier");
        // A rotated box whose corner sums round: its faces are twisted
        // bilinear patches, so the cube cuts them along curves.
        {
            const double c = num::det::cos(0.7), sn = num::det::sin(0.7);
            const Brep B = op::to_brep(op::make_box({0.3, 0.2, 0.1}, {0.6 * c, 0.6 * sn, 0.0}, {-0.6 * sn, 0.6 * c, 0.1},
                                                    {0.0, -0.1, 0.9}));
            const op::BooleanResult r = op::boolean_op(box(0, 0, 0, 1, 1, 1), B, op::Op::Difference);
            AK_CHECK(r.ok);
            if (r.ok) {
                AK_CHECK(!op::certified_volume(r.solid, &why).has_value());
                std::printf("  cube - twisted rotated box: not certified (%s)\n", why.c_str());
            }
        }
        const Brep hill = op::to_brep(op::make_bump_box({0, 0, 0}, 1, 1, 1, 0.3));
        const op::BooleanResult r = op::boolean_op(hill, box(-0.5, -0.5, 1.1, 1.5, 1.5, 2.0), op::Op::Difference);
        AK_CHECK(r.ok);
        if (r.ok) {
            AK_CHECK(!op::certified_volume(r.solid, &why).has_value());
            std::printf("  hill box cut through its top: not certified (%s)\n", why.c_str());
        }
    }

    std::printf("  %d enclosures, widest %.1e\n", g_certified, g_widest);
    return akt::finish("certified_volume_test");
}
