// Imprinting and Boolean operations on boundary representations.
//
// Every volume is checked against something the kernel did not compute:
// a closed form, a plain 2D quadrature, or a Monte Carlo count that uses
// only the analytic definition of the solids.
#include <akernel/geom/primitives.hpp>
#include <akernel/op/boolean.hpp>

#include <functional>
#include <memory>

#include "check.hpp"

using namespace ak;
using geom::Vec3d;
using op::BooleanResult;
using op::Brep;
using op::Imprint;
using op::Op;
using op::PatchSolid;

namespace {

using Member = std::function<bool(const Vec3d&)>;

struct Volumes {
    bool ok = false;
    double a = 0, b = 0, inter = 0, uni = 0, diff = 0;
    BooleanResult r_inter, r_uni, r_diff;
};

/// Checks every valid result must pass. `genus` is the expected total genus
/// of its shells.
void check_result(const BooleanResult& r, int genus = 0) {
    AK_CHECK(r.ok);
    if (!r.ok) std::printf("    result rejected: %s\n", r.error.c_str());
    if (r.ok && !r.empty) {
        AK_CHECK(r.audit.watertight());
        AK_CHECK(r.audit.ledger_entries == 0); // no tolerance was recorded anywhere
        int total = 0;
        for (const auto& sh : r.audit.shells) total += sh.genus;
        AK_CHECK(total == genus);
    }
}

/// Run all three operations on one imprint.
Volumes run_case(const char* name, const Brep& A, const Brep& B, Imprint& im) {
    Volumes v;
    op::BooleanOptions opt;
    opt.trace = true;
    im = op::imprint(A, B, opt);
    AK_CHECK(im.ok);
    if (!im.ok) {
        std::printf("%s: imprint failed: %s\n", name, im.error.c_str());
        for (const std::string& line : im.trace)
            if (line.rfind("  the curve", 0) != 0 && line.rfind("vertex ", 0) != 0) std::printf("    trace: %s\n", line.c_str());
        return v;
    }
    v.a = op::volume(A);
    v.b = op::volume(B);
    v.r_inter = op::select(im, Op::Intersection);
    v.r_uni = op::select(im, Op::Union);
    v.r_diff = op::select(im, Op::Difference);
    for (const BooleanResult* r : {&v.r_inter, &v.r_uni, &v.r_diff}) check_result(*r);
    v.inter = op::volume(v.r_inter.solid);
    v.uni = op::volume(v.r_uni.solid);
    v.diff = op::volume(v.r_diff.solid);
    // Inclusion-exclusion between three independently assembled results.
    AK_CHECK_NEAR(v.inter + v.uni, v.a + v.b, 1e-10);
    AK_CHECK_NEAR(v.diff, v.a - v.inter, 1e-10);
    AK_CHECK(op::edge_gap_bound(im) < 1e-10);
    std::printf("%s: curves %d (closed %d), new edges %d, new vertices %d (triple %d), regions %d, "
                "vertex width %.1e, gap bound %.1e\n"
                "    volumes: A %.12f  B %.12f  A&B %.12f  A|B %.12f  A-B %.12f\n",
                name, im.stats.curves, im.stats.closed_curves, im.stats.new_arcs, im.stats.new_vertices,
                im.stats.triple_vertices, im.stats.regions, im.stats.max_vertex_width, op::edge_gap_bound(im), v.a,
                v.b, v.inter, v.uni, v.diff);
    v.ok = true;
    return v;
}

/// Volume of a set given by a membership test, by Monte Carlo over a box.
double monte_carlo(const Vec3d& lo, const Vec3d& hi, const Member& in, int samples, double& sigma) {
    akt::Rng rng(777);
    int hits = 0;
    for (int i = 0; i < samples; ++i)
        if (in({rng.range(lo.x, hi.x), rng.range(lo.y, hi.y), rng.range(lo.z, hi.z)})) ++hits;
    const double box = (hi.x - lo.x) * (hi.y - lo.y) * (hi.z - lo.z);
    const double p = static_cast<double>(hits) / samples;
    sigma = box * std::sqrt(std::max(p * (1.0 - p), 1.0 / samples) / samples);
    return box * p;
}

Member in_box(const Vec3d& lo, const Vec3d& hi) {
    return [=](const Vec3d& p) {
        return p.x >= lo.x && p.x <= hi.x && p.y >= lo.y && p.y <= hi.y && p.z >= lo.z && p.z <= hi.z;
    };
}

} // namespace

int main() {
    const PatchSolid unit_patch = op::make_box({0, 0, 0}, {1, 1, 1});
    const Brep unit = op::to_brep(unit_patch);

    // --- the input model ---------------------------------------------------------
    {
        AK_CHECK(unit_patch.faces.size() == 6);
        AK_CHECK(unit_patch.edges.size() == 12);
        AK_CHECK(unit_patch.vertex_points.size() == 8);
        AK_CHECK_NEAR(op::volume(unit_patch), 1.0, 1e-14);
        AK_CHECK_NEAR(op::volume(unit), 1.0, 1e-14);
        AK_CHECK(op::audit(unit).watertight());

        const PatchSolid bump = op::make_bump_box({0, 0, 0}, 1, 1, 1, 0.3);
        AK_CHECK(bump.edges.size() == 12 && bump.vertex_points.size() == 8);
        AK_CHECK_NEAR(op::volume(op::to_brep(bump)), 1.0 + 0.3 * 4.0 / 9.0, 1e-13);

        const PatchSolid sector = op::make_annular_sector(1.0, 2.0, 0.0, 1.0);
        AK_CHECK(sector.edges.size() == 12 && sector.vertex_points.size() == 8);
        const double pi = 3.14159265358979323846;
        AK_CHECK_NEAR(op::volume(op::to_brep(sector)), pi / 4.0 * 3.0, 1e-9);

        // A face whose boundary does not match its neighbour bit for bit is
        // not accepted as closed.
        PatchSolid broken = unit_patch;
        broken.faces[0].at(0, 0).x += 1e-12;
        std::string err;
        AK_CHECK(!broken.build(err));
        AK_CHECK(!err.empty());
    }

    // --- placing points on a certified curve ------------------------------------
    {
        // The planes z = 0 (u = x, v = y) and x = 0.5 (s = y, t = z + 0.5)
        // meet in the segment x = 0.5, z = 0, y in [0, 1]: the points
        // (u, v, s, t) = (0.5, y, y, 0.5). Beyond y = 0 and y = 1 the same
        // line continues outside both parameter squares, and the certified
        // tubes reach a little way into that continuation.
        const auto pa = geom::make_plane({0, 0, 0}, {1, 0, 0}, {0, 1, 0});
        const auto pb = geom::make_plane({0.5, 0, -0.5}, {0, 1, 0}, {0, 0, 1});
        auto sys = std::make_shared<solve::SsiSystem>(pa, pb);
        const solve::SsiResult res = solve::intersect(*sys);
        AK_CHECK(res.complete && res.curves.size() == 1);
        if (res.curves.size() == 1) {
            const auto base = op::make_curve_base(sys, std::make_shared<const solve::SsiCurve>(res.curves[0]));
            auto at = [&](double y, bool& on) {
                const solve::Box4 P = {num::Interval(0.5), num::Interval(y), num::Interval(y), num::Interval(0.5)};
                op::CurvePos pos;
                return op::locate_on_curve(*base, P, on, pos);
            };
            // Beyond each end, the points of the continuation that lie inside
            // some certified region: the case that matters. How far the
            // margin reaches depends on the last step, so search for them.
            int tested = 0;
            for (double sign : {-1.0, 1.0})
                for (double d = 1e-6; d < 0.05; d *= 1.5) {
                    const double y = sign < 0.0 ? -d : 1.0 + d;
                    const solve::Box4 P = {num::Interval(0.5), num::Interval(y), num::Interval(y), num::Interval(0.5)};
                    bool in_some = false;
                    for (const auto& seg : res.curves[0].segs) in_some = in_some || solve::in_region(P, seg);
                    if (!in_some) continue;
                    bool on = true;
                    AK_CHECK(at(y, on));
                    AK_CHECK(!on); // on the continuation, not on the curve
                    ++tested;
                }
            std::printf("points on the continuation inside a certified margin: %d tested\n", tested);
            AK_CHECK(tested >= 10);
            for (double y : {0.0, 0.3, 0.75, 1.0}) {
                bool on = false;
                AK_CHECK(at(y, on));
                AK_CHECK(on);
            }
        }
    }

    // --- 1. two offset boxes: everything has a closed form ----------------------
    const Brep boxB = op::to_brep(op::make_box({0.5, 0.4, 0.3}, {1.5, 1.4, 1.3}));
    {
        Imprint im;
        const Volumes v = run_case("offset boxes", unit, boxB, im);
        if (v.ok) {
            AK_CHECK_NEAR(v.inter, 0.5 * 0.6 * 0.7, 1e-12);
            AK_CHECK_NEAR(v.uni, 2.0 - 0.21, 1e-12);
            AK_CHECK_NEAR(v.diff, 1.0 - 0.21, 1e-12);
            // The intersection is a box.
            AK_CHECK(v.r_inter.solid.faces.size() == 6);
            AK_CHECK(v.r_inter.solid.edges.size() == 12);
            AK_CHECK(v.r_inter.solid.vertices.size() == 8);
            AK_CHECK(v.r_inter.audit.shells.size() == 1);
            AK_CHECK(v.r_uni.audit.shells.size() == 1);
            AK_CHECK(v.r_diff.audit.shells.size() == 1);
            AK_CHECK(im.stats.ray_tests == 0);
            AK_CHECK(im.stats.new_vertices == 6 && im.stats.triple_vertices == 0);
        }
    }

    // --- 2. a rotated box: checked by Monte Carlo ---------------------------------
    {
        // Rotation by 0.6 rad about the axis (1, 2, 3) / sqrt(14).
        const double ang = 0.6, c = std::cos(ang), s = std::sin(ang), n = std::sqrt(14.0);
        const double ax = 1.0 / n, ay = 2.0 / n, az = 3.0 / n;
        const double R[3][3] = {{c + ax * ax * (1 - c), ax * ay * (1 - c) - az * s, ax * az * (1 - c) + ay * s},
                                {ay * ax * (1 - c) + az * s, c + ay * ay * (1 - c), ay * az * (1 - c) - ax * s},
                                {az * ax * (1 - c) - ay * s, az * ay * (1 - c) + ax * s, c + az * az * (1 - c)}};
        const double size = 0.9;
        const Vec3d o{0.35, 0.25, 0.3};
        const Vec3d ex{R[0][0] * size, R[1][0] * size, R[2][0] * size};
        const Vec3d ey{R[0][1] * size, R[1][1] * size, R[2][1] * size};
        const Vec3d ez{R[0][2] * size, R[1][2] * size, R[2][2] * size};
        const Brep B = op::to_brep(op::make_box(o, ex, ey, ez));
        AK_CHECK_NEAR(op::volume(B), size * size * size, 1e-13);
        Imprint im;
        const Volumes v = run_case("rotated box", unit, B, im);
        if (v.ok) {
            const Member in_a = in_box({0, 0, 0}, {1, 1, 1});
            auto in_b = [&](const Vec3d& p) {
                const Vec3d d = p - o;
                const double x = dot(d, ex) / (size * size), y = dot(d, ey) / (size * size),
                             z = dot(d, ez) / (size * size);
                return x >= 0 && x <= 1 && y >= 0 && y <= 1 && z >= 0 && z <= 1;
            };
            double sigma = 0.0;
            const double mc = monte_carlo({0, 0, 0}, {1, 1, 1}, [&](const Vec3d& p) { return in_a(p) && in_b(p); },
                                          2000000, sigma);
            std::printf("    Monte Carlo A&B %.6f +- %.6f, kernel %.6f\n", mc, sigma, v.inter);
            AK_CHECK(std::fabs(mc - v.inter) < 4.0 * sigma);
            AK_CHECK(v.inter > 0.05 && v.inter < 0.7);
        }
    }

    // --- 3. a closed intersection loop: the cap of a hill -------------------------
    const double bump = 0.3, cut = 1.15;
    const Brep hill = op::to_brep(op::make_bump_box({0, 0, 0}, 1, 1, 1, bump));
    const Brep slab = op::to_brep(op::make_box({-0.5, -0.5, cut}, {1.5, 1.5, 2.0}));
    auto top = [&](double x, double y) { return 1.0 + 16.0 * bump * x * (1 - x) * y * (1 - y); };
    Brep hill_cut; // the hill with its cap removed, used again below
    {
        Imprint im;
        const Volumes v = run_case("hill cap", hill, slab, im);
        if (v.ok) {
            hill_cut = v.r_diff.solid;
            // Cap volume by plain midpoint quadrature of max(0, top - cut).
            const int N = 2000;
            double cap = 0.0;
            for (int i = 0; i < N; ++i)
                for (int j = 0; j < N; ++j) cap += std::max(0.0, top((i + 0.5) / N, (j + 0.5) / N) - cut);
            cap /= static_cast<double>(N) * N;
            std::printf("    cap volume: quadrature %.9f, kernel %.9f\n", cap, v.inter);
            AK_CHECK_NEAR(v.inter, cap, 2e-7);
            AK_CHECK(im.stats.closed_curves == 1 && im.stats.curves == 1);

            // The cap: two faces, one closed edge, one seam vertex.
            AK_CHECK(v.r_inter.solid.faces.size() == 2);
            AK_CHECK(v.r_inter.solid.edges.size() == 1);
            AK_CHECK(v.r_inter.solid.vertices.size() == 1);
            // The remainder: the top face has a ring, and the cut face is flipped.
            AK_CHECK(v.r_diff.solid.faces.size() == 7);
            AK_CHECK(v.r_diff.audit.shells.size() == 1 && v.r_diff.audit.shells[0].rings == 1);
            int flipped = 0;
            for (const auto& f : v.r_diff.solid.faces) flipped += f.flipped ? 1 : 0;
            AK_CHECK(flipped == 1);

            // Certified point classification against the trimmed top face,
            // compared with the analytic loop.
            int region = -1;
            for (int ri = 0; ri < static_cast<int>(im.regions.size()); ++ri)
                if (im.regions[ri].solid == 0 && im.regions[ri].face == 5 && im.regions[ri].inside == 1) region = ri;
            AK_CHECK(region >= 0);
            if (region >= 0) {
                akt::Rng rng(99);
                int wrong = 0, undecided = 0, inside = 0;
                const int samples = 4000;
                for (int i = 0; i < samples; ++i) {
                    const double u = rng.unit(), w = rng.unit();
                    const op::Where where = op::locate(im, region, u, w);
                    if (where == op::Where::Undecided || where == op::Where::Boundary) { ++undecided; continue; }
                    if ((where == op::Where::Inside) != (top(u, w) > cut)) ++wrong;
                    if (where == op::Where::Inside) ++inside;
                }
                // Points very close to the loop.
                int near_wrong = 0, near_undecided = 0;
                for (int i = 0; i < 400; ++i) {
                    // Walk outward from the centre to the loop along a random ray by bisection.
                    const double th = rng.range(0.0, 6.283185307179586);
                    double lo = 0.0, hi = 0.49;
                    for (int it = 0; it < 60; ++it) {
                        const double m = 0.5 * (lo + hi);
                        if (top(0.5 + m * std::cos(th), 0.5 + m * std::sin(th)) > cut) lo = m; else hi = m;
                    }
                    for (double off : {-1e-6, 1e-6}) {
                        const double m = lo + off;
                        const op::Where where = op::locate(im, region, 0.5 + m * std::cos(th), 0.5 + m * std::sin(th));
                        if (where == op::Where::Undecided || where == op::Where::Boundary) { ++near_undecided; continue; }
                        if ((where == op::Where::Inside) != (off < 0.0)) ++near_wrong;
                    }
                }
                std::printf("    classification: %d random points, %d inside, %d wrong, %d undecided; "
                            "800 points 1e-6 from the loop: %d wrong, %d undecided\n",
                            samples, inside, wrong, undecided, near_wrong, near_undecided);
                AK_CHECK(wrong == 0);
                AK_CHECK(undecided == 0);
                AK_CHECK(near_wrong == 0);
                AK_CHECK(near_undecided == 0);
                AK_CHECK(inside > 100);
            }
        }
    }

    // --- 4. rational faces: a cylindrical shell cut by a box ---------------------
    const Brep shell = op::to_brep(op::make_annular_sector(1.0, 2.0, 0.0, 1.0));
    const Member in_shell = [](const Vec3d& p) {
        const double r2 = p.x * p.x + p.y * p.y;
        return p.x >= 0 && p.y >= 0 && r2 >= 1.0 && r2 <= 4.0 && p.z >= 0 && p.z <= 1;
    };
    {
        const Brep B = op::to_brep(op::make_box({0.5, 0.6, 0.25}, {1.3, 1.4, 0.75}));
        Imprint im;
        const Volumes v = run_case("cylinder shell and box", shell, B, im);
        if (v.ok) {
            // Part of the box footprint inside the unit circle.
            auto prim = [](double x) { return 0.5 * (x * std::sqrt(1.0 - x * x) + std::asin(x)); };
            const double inside_circle = prim(0.8) - prim(0.5) - 0.6 * 0.3;
            const double expect = (0.8 * 0.8 - inside_circle) * 0.5;
            std::printf("    closed form A&B %.12f, kernel %.12f\n", expect, v.inter);
            AK_CHECK_NEAR(v.inter, expect, 1e-9);
            AK_CHECK(v.r_inter.audit.shells.size() == 1);
            AK_CHECK(im.stats.new_vertices == 4);
            AK_CHECK(im.stats.new_arcs == 4);
        }
    }

    // --- 4b. two curved solids: a hill against the cylindrical shell ---------------
    {
        const Brep A = op::to_brep(op::make_bump_box({0.6, 0.6, -0.35}, 1, 1, 1, 0.3));
        Imprint im;
        const Volumes v = run_case("hill and cylinder shell", A, shell, im);
        if (v.ok) {
            auto in_a = [](const Vec3d& p) {
                const double x = p.x - 0.6, y = p.y - 0.6;
                return x >= 0 && x <= 1 && y >= 0 && y <= 1 && p.z >= -0.35 &&
                       p.z <= 0.65 + 4.8 * x * (1 - x) * y * (1 - y);
            };
            double sigma = 0.0;
            const double mc = monte_carlo({0.6, 0.6, -0.35}, {1.6, 1.6, 0.95},
                                          [&](const Vec3d& p) { return in_a(p) && in_shell(p); }, 2000000, sigma);
            std::printf("    Monte Carlo A&B %.6f +- %.6f, kernel %.6f\n", mc, sigma, v.inter);
            AK_CHECK(std::fabs(mc - v.inter) < 4.0 * sigma);
            AK_CHECK(v.inter > 0.1);
        }
    }

    // --- 5. solids that do not touch ------------------------------------------------
    {
        const Brep B = op::to_brep(op::make_box({2, 2, 2}, {3, 3, 3}));
        Imprint im;
        const Volumes v = run_case("disjoint boxes", unit, B, im);
        if (v.ok) {
            AK_CHECK(v.r_inter.empty);
            AK_CHECK(v.r_uni.audit.shells.size() == 2);
            AK_CHECK_NEAR(v.uni, 2.0, 1e-12);
            AK_CHECK_NEAR(v.diff, 1.0, 1e-12);
            AK_CHECK(im.stats.ray_tests == 2);
        }
    }
    {
        const Brep B = op::to_brep(op::make_box({0.3, 0.3, 0.3}, {0.6, 0.6, 0.6}));
        Imprint im;
        const Volumes v = run_case("nested boxes", unit, B, im);
        if (v.ok) {
            AK_CHECK_NEAR(v.inter, 0.027, 1e-12);
            AK_CHECK_NEAR(v.uni, 1.0, 1e-12);
            AK_CHECK_NEAR(v.diff, 1.0 - 0.027, 1e-12);
            AK_CHECK(v.r_diff.audit.shells.size() == 2); // a solid with a cavity
            AK_CHECK(v.r_inter.solid.faces.size() == 6);
        }
    }

    // --- 6. chained operations: a result as an operand ------------------------------
    {
        // Three boxes; every volume has a closed form.
        const Brep C = op::to_brep(op::make_box({0.2, 0.2, 0.6}, {0.8, 0.7, 1.5}));
        const BooleanResult a_minus_b = op::boolean_op(unit, boxB, Op::Difference);
        const BooleanResult a_or_b = op::boolean_op(unit, boxB, Op::Union);
        const BooleanResult a_and_b = op::boolean_op(unit, boxB, Op::Intersection);
        check_result(a_minus_b);
        check_result(a_or_b);
        check_result(a_and_b);
        if (a_minus_b.ok && a_or_b.ok && a_and_b.ok) {
            Imprint im;
            const Volumes v1 = run_case("(A - B) with C", a_minus_b.solid, C, im);
            if (v1.ok) {
                // A - B - C = A minus (B or C).
                AK_CHECK_NEAR(v1.diff, 1.0 - (0.21 + 0.12 - 0.036), 1e-11);
                AK_CHECK_NEAR(v1.inter, 0.12 - 0.036, 1e-11);
                AK_CHECK(im.stats.triple_vertices > 0); // C cuts edges made by the first operation
            }
            const Volumes v2 = run_case("(A | B) with C", a_or_b.solid, C, im);
            if (v2.ok) AK_CHECK_NEAR(v2.inter, 0.12 + 0.063 - 0.036, 1e-11);
            const Volumes v3 = run_case("(A & B) with C", a_and_b.solid, C, im);
            if (v3.ok) {
                AK_CHECK_NEAR(v3.uni, 0.21 + 0.27 - 0.036, 1e-11);
                AK_CHECK_NEAR(v3.inter, 0.036, 1e-11);
            }
            // Both operands may be results. E cuts a notch out of the part of
            // C that lies above A, so (A - B) & (C - E) = (A & C) - B.
            const Brep E = op::to_brep(op::make_box({0.1, 0.5, 1.05}, {0.45, 0.9, 1.7}));
            const BooleanResult c_minus_e = op::boolean_op(C, E, Op::Difference);
            check_result(c_minus_e);
            if (c_minus_e.ok) {
                AK_CHECK_NEAR(op::volume(c_minus_e.solid), 0.27 - 0.25 * 0.2 * 0.45, 1e-12);
                const Volumes v4 = run_case("(A - B) with (C - E)", a_minus_b.solid, c_minus_e.solid, im);
                if (v4.ok) AK_CHECK_NEAR(v4.inter, 0.12 - 0.036, 1e-11);
            }
        }
    }
    {
        // The hill with its cap removed, then cut by a box that crosses the
        // closed edge left by the first cut.
        const Vec3d lo{0.4, -0.5, 0.9}, hi{0.6, 1.5, 2.0};
        const Brep C = op::to_brep(op::make_box(lo, hi));
        if (!hill_cut.faces.empty()) {
            Imprint im;
            const Volumes v = run_case("(hill - slab) with a box across the loop", hill_cut, C, im);
            if (v.ok) {
                const Member in_c = in_box(lo, hi);
                const Member in_hill_cut = [&](const Vec3d& p) {
                    return p.x >= 0 && p.x <= 1 && p.y >= 0 && p.y <= 1 && p.z >= 0 && p.z <= top(p.x, p.y) && p.z <= cut;
                };
                double sigma = 0.0;
                const double mc = monte_carlo({0, 0, 0}, {1, 1, 1.3},
                                              [&](const Vec3d& p) { return in_hill_cut(p) && !in_c(p); }, 2000000, sigma);
                std::printf("    Monte Carlo (hill - slab) - box %.6f +- %.6f, kernel %.6f\n", mc, sigma, v.diff);
                AK_CHECK(std::fabs(mc - v.diff) < 4.0 * sigma);
                // The two faces x = 0.4 and x = 0.6 each cross the closed edge twice.
                AK_CHECK(im.stats.triple_vertices == 4);
            }
        }
    }
    {
        // A longer chain with fixed operands, each step checked by Monte Carlo.
        struct Step {
            Vec3d lo, hi;
            Op op;
            int genus;
        };
        // Step 1 cuts a notch out of a top corner. Step 3 cuts a pocket up
        // from the bottom that runs into the notch, which makes a tunnel:
        // from then on the solid has genus 1.
        const Step steps[] = {{{0.55, -0.2, 0.35}, {1.3, 0.45, 1.2}, Op::Difference, 0},
                              {{-0.3, 0.6, 0.15}, {0.4, 1.25, 0.7}, Op::Union, 0},
                              {{0.1, 0.1, -0.4}, {0.9, 0.85, 0.5}, Op::Difference, 1},
                              {{-0.6, -0.5, 0.05}, {1.6, 1.7, 0.95}, Op::Intersection, 1},
                              {{0.3, 0.3, 0.6}, {0.75, 1.5, 1.4}, Op::Union, 1}};
        Brep current = unit;
        Member member = in_box({0, 0, 0}, {1, 1, 1});
        int depth = 0;
        for (const Step& st : steps) {
            const Brep B = op::to_brep(op::make_box(st.lo, st.hi));
            const BooleanResult r = op::boolean_op(current, B, st.op);
            check_result(r, st.genus);
            if (!r.ok || r.empty) break;
            const Member prev = member, in_b = in_box(st.lo, st.hi);
            if (st.op == Op::Difference) member = [=](const Vec3d& p) { return prev(p) && !in_b(p); };
            else if (st.op == Op::Union) member = [=](const Vec3d& p) { return prev(p) || in_b(p); };
            else member = [=](const Vec3d& p) { return prev(p) && in_b(p); };
            double sigma = 0.0;
            const double mc = monte_carlo({-0.7, -0.6, -0.5}, {1.7, 1.8, 1.5}, member, 1000000, sigma);
            const double vol = op::volume(r.solid);
            ++depth;
            std::printf("chain step %d: %d faces, %d edges, %d vertices, volume %.6f, Monte Carlo %.6f +- %.6f, "
                        "edge gap bound %.1e, vertex width %.1e\n",
                        depth, static_cast<int>(r.solid.faces.size()), static_cast<int>(r.solid.edges.size()),
                        static_cast<int>(r.solid.vertices.size()), vol, mc, sigma, op::edge_gap_bound(r.solid),
                        op::vertex_width(r.solid));
            AK_CHECK(std::fabs(mc - vol) < 4.0 * sigma);
            AK_CHECK(op::edge_gap_bound(r.solid) < 1e-10);
            AK_CHECK(op::vertex_width(r.solid) < 1e-10);
            current = r.solid;
        }
        AK_CHECK(depth == 5);
    }

    // --- 7. not in general position: reported, never answered wrongly ------------
    {
        // An imprint that is not solved gives no result.
        const Brep B = op::to_brep(op::make_box({-0.5, -0.5, 1.3}, {1.5, 1.5, 2.0}));
        const Imprint im = op::imprint(op::to_brep(op::make_bump_box({0, 0, 0}, 1, 1, 1, 0.3)), B);
        AK_CHECK(!im.ok);
        const BooleanResult r = op::select(im, Op::Union);
        AK_CHECK(!r.ok && r.solid.faces.empty());
    }
    {
        // A box whose bottom touches the top of a hill in one point.
        const Brep B = op::to_brep(op::make_box({-0.5, -0.5, 1.3}, {1.5, 1.5, 2.0}));
        const Imprint im = op::imprint(hill, B);
        AK_CHECK(!im.ok);
        std::printf("box touching a hill: %s\n", im.error.c_str());
    }
    {
        // A tilted box whose edge from (0.5, 0.5, 0.5) to (0.5, 1.5, 1.5)
        // leaves the unit box exactly through its edge y = 1, z = 1. That edge
        // lies in the plane y = z of one of B's faces, and the unit box's two
        // faces there lie on the two sides of it (z = 1 inside B, y = 1
        // outside). Solved since the classification no longer crosses an edge
        // lying in a face plane of the other solid. B is 0.5 <= x <= 1.3,
        // 1 <= y + z <= 3, 0 <= z - y <= 2: the intersection is x in [0.5, 1]
        // times the triangle z >= y, y + z >= 1 of the unit square, 1/8.
        const Brep B = op::to_brep(op::make_box({0.5, 0.5, 0.5}, {0.8, 0, 0}, {0, 1, 1}, {0, -1, 1}));
        AK_CHECK_NEAR(op::volume(B), 1.6, 1e-12);
        const Imprint im = op::imprint(unit, B);
        AK_CHECK(im.ok);
        double v[3] = {0.0, 0.0, 0.0};
        const Op ops[3] = {Op::Intersection, Op::Union, Op::Difference};
        for (int i = 0; i < 3 && im.ok; ++i) {
            const BooleanResult r = op::select(im, ops[i]);
            AK_CHECK(r.ok && !r.empty && r.audit.closed());
            if (r.ok && !r.empty) v[i] = op::volume(r.solid);
        }
        AK_CHECK_NEAR(v[0], 0.125, 1e-12);
        AK_CHECK_NEAR(v[1], 2.475, 1e-12);
        AK_CHECK_NEAR(v[2], 0.875, 1e-12);
        std::printf("edge leaving through an edge: %s; volumes %.12f %.12f %.12f\n", im.ok ? "solved" : im.error.c_str(), v[0],
                    v[1], v[2]);
    }
    {
        // An edge of B crossing the vertical edge x = 1, y = 1 of the unit
        // box at (1, 1, 0.5), with every face pair otherwise transversal.
        // Solved since the crossing point is exact (op/boolean.hpp, "edges
        // through edges"). In B's coordinates o + a e1 + b e2 + c e3 the
        // intersection is 0.5 + 0.4(b - c) <= a <= 0.5 - 0.4(b - c) with
        // 0.4 b + 0.8 c <= 0.5: the integral of 0.8 (c - b) over that region,
        // 25/1152, times |det| = 0.96, so 1/48.
        const Brep B = op::to_brep(op::make_box({0.5, 1.5, 0.5}, {1, -1, 0}, {0.4, 0.4, 0.4}, {-0.4, -0.4, 0.8}));
        AK_CHECK_NEAR(op::volume(B), 0.96, 1e-12);
        const Imprint im = op::imprint(unit, B);
        AK_CHECK(im.ok);
        AK_CHECK(im.stats.edge_crossings == 1);
        double v[3] = {0.0, 0.0, 0.0};
        const Op ops[3] = {Op::Intersection, Op::Union, Op::Difference};
        for (int i = 0; i < 3 && im.ok; ++i) {
            const BooleanResult r = op::select(im, ops[i]);
            AK_CHECK(r.ok && !r.empty && r.audit.closed());
            if (r.ok && !r.empty) v[i] = op::volume(r.solid);
        }
        AK_CHECK_NEAR(v[0], 1.0 / 48.0, 1e-12);
        AK_CHECK_NEAR(v[1], 1.96 - 1.0 / 48.0, 1e-12);
        AK_CHECK_NEAR(v[2], 1.0 - 1.0 / 48.0, 1e-12);
        std::printf("edge crossing edge: %s; volumes %.12f %.12f %.12f\n", im.ok ? "solved" : im.error.c_str(), v[0], v[1], v[2]);
    }
    // --- 7b. coincidence decided exactly: boxes in contact ----------------------
    // Different patches in one plane, each parameterised its own way; the kernel
    // decides exactly that they are one plane.
    struct Contact {
        const char* name;
        Brep b;
        double inter, uni, diff;
        int faces_uni; // -1: not checked
    };
    const Contact contacts[] = {
        {"boxes face to face", op::to_brep(op::make_box({1, 0, 0}, {2, 1, 1})), 0.0, 2.0, 1.0, 10},
        {"face to face, offset", op::to_brep(op::make_box({1, 0.2, 0.3}, {2, 1.3, 0.9})), 0.0, 1.0 + 0.66, 1.0, -1},
        {"stacked", op::to_brep(op::make_box({0.3, 0.2, 1}, {0.8, 0.7, 1.5})), 0.0, 1.125, 1.0, 11},
        {"overlapping, two planes shared", op::to_brep(op::make_box({0.5, 0, 0.3}, {1.5, 1, 0.8})), 0.25, 1.25, 0.75, -1},
        {"inside, touching one face", op::to_brep(op::make_box({0.25, 0.25, 0}, {0.75, 0.75, 0.5})), 0.125, 1.0, 0.875, -1},
        // The unit cube again, built mirrored: every face is the same square,
        // parameterised the other way round.
        {"the same cube, built another way", op::to_brep(op::make_box({1, 0, 1}, {-1, 0, 0}, {0, 1, 0}, {0, 0, -1})), 1.0, 1.0, 0.0, 6},
    };
    for (const Contact& c : contacts) {
        Imprint im;
        const Volumes v = run_case(c.name, unit, c.b, im);
        if (!v.ok) continue;
        std::printf("    exactly coplanar surface pairs: %d; shared planes: %d; regions on the other's boundary: %d\n",
                    im.stats.coplanar_pairs, im.stats.shared_carriers, im.stats.on_regions);
        AK_CHECK(im.stats.coplanar_pairs > 0);
        AK_CHECK_NEAR(v.inter, c.inter, 1e-12);
        AK_CHECK_NEAR(v.uni, c.uni, 1e-12);
        AK_CHECK_NEAR(v.diff, c.diff, 1e-12);
        AK_CHECK(c.inter > 0.0 || v.r_inter.empty);
        AK_CHECK(c.diff > 0.0 || v.r_diff.empty);
        if (c.faces_uni >= 0) AK_CHECK(static_cast<int>(v.r_uni.solid.faces.size()) == c.faces_uni);
        AK_CHECK(v.r_uni.audit.shells.size() == 1);
    }
    {
        // Results of contact operations chain: a row of three boxes, then a
        // notch cut through the seams between them.
        const BooleanResult ab = op::boolean_op(unit, op::to_brep(op::make_box({1, 0, 0}, {2, 1, 1})), Op::Union);
        check_result(ab);
        if (ab.ok) {
            const BooleanResult abc = op::boolean_op(ab.solid, op::to_brep(op::make_box({2, 0, 0}, {3, 1, 1})), Op::Union);
            check_result(abc);
            if (abc.ok) {
                AK_CHECK_NEAR(op::volume(abc.solid), 3.0, 1e-12);
                const BooleanResult notch = op::boolean_op(abc.solid, op::to_brep(op::make_box({0.5, -0.5, 0.5}, {2.5, 0.5, 1.5})), Op::Difference);
                check_result(notch);
                if (notch.ok) {
                    std::printf("row of three boxes with a notch: %d faces, volume %.12f\n",
                                static_cast<int>(notch.solid.faces.size()), op::volume(notch.solid));
                    AK_CHECK_NEAR(op::volume(notch.solid), 3.0 - 2.0 * 0.5 * 0.5, 1e-12);
                }
            }
        }
    }

    // --- 8. shared by construction: a result with its own operand ---------------
    {
        // (A - B) shares three faces, their edges and a corner with B.
        const BooleanResult a_minus_b = op::boolean_op(unit, boxB, Op::Difference);
        check_result(a_minus_b);
        if (a_minus_b.ok) {
            Imprint im;
            const Volumes v = run_case("(A - B) with B", a_minus_b.solid, boxB, im);
            if (v.ok) {
                AK_CHECK_NEAR(v.uni, 1.79, 1e-12);  // = A | B
                AK_CHECK(v.r_inter.empty);           // they only touch
                AK_CHECK_NEAR(v.diff, 0.79, 1e-12);  // = A - B
                AK_CHECK(v.r_diff.solid.faces.size() == a_minus_b.solid.faces.size());
                AK_CHECK(v.r_uni.audit.shells.size() == 1);
                AK_CHECK(im.stats.shared_carriers == 3);
                AK_CHECK(im.stats.on_regions == 6); // three faces of each solid
                AK_CHECK(im.stats.new_vertices == 0); // every vertex existed already
                std::printf("    shared: %d surfaces, %d sides, %d vertices, %d edge pieces; %d regions on the "
                            "other's boundary\n",
                            im.stats.shared_carriers, im.stats.shared_sides, im.stats.merged_vertices,
                            im.stats.unified_pieces, im.stats.on_regions);
            }
            // The union, used again: (A - B) | B is A | B, so cutting B out
            // once more gives A - B back.
            const BooleanResult back = op::boolean_op(v.r_uni.solid, boxB, Op::Difference);
            check_result(back);
            if (back.ok) AK_CHECK_NEAR(op::volume(back.solid), 0.79, 1e-12);
        }
    }
    {
        // A solid with itself: every surface, edge and vertex is shared.
        Imprint im;
        const Volumes v = run_case("A with A", unit, unit, im);
        if (v.ok) {
            AK_CHECK_NEAR(v.uni, 1.0, 1e-12);
            AK_CHECK_NEAR(v.inter, 1.0, 1e-12);
            AK_CHECK(v.r_diff.empty);
            AK_CHECK(v.r_uni.solid.faces.size() == 6);
            AK_CHECK(v.r_inter.solid.faces.size() == 6);
            AK_CHECK(im.stats.shared_carriers == 6);
            AK_CHECK(im.stats.unified_pieces == 12);
        }
    }
    {
        // A rotated box cut out and put back: shared surfaces that are not
        // axis-aligned (their corners are rounded, so they are not exactly
        // planar).
        const double c = std::cos(0.5), s = std::sin(0.5);
        const Brep R = op::to_brep(op::make_box({0.4, 0.3, 0.35}, {0.8 * c, 0.8 * s, 0}, {-0.7 * s, 0.7 * c, 0.1}, {0.05, -0.05, 0.75}));
        const BooleanResult a_minus_r = op::boolean_op(unit, R, Op::Difference);
        check_result(a_minus_r);
        if (a_minus_r.ok) {
            Imprint im;
            const Volumes v = run_case("(A - rotated box) with the rotated box", a_minus_r.solid, R, im);
            if (v.ok) {
                const BooleanResult a_or_r = op::boolean_op(unit, R, Op::Union);
                check_result(a_or_r);
                AK_CHECK_NEAR(v.uni, op::volume(a_or_r.solid), 1e-11);
                AK_CHECK(v.r_inter.empty);
                AK_CHECK_NEAR(v.diff, op::volume(a_minus_r.solid), 1e-11);
            }
        }
    }
    {
        // Curved: the hill with its cap removed, combined with the slab that
        // removed it. They share the slab's bottom surface and the closed edge.
        if (!hill_cut.faces.empty()) {
            Imprint im;
            const Volumes v = run_case("(hill - slab) with slab", hill_cut, slab, im);
            if (v.ok) {
                const BooleanResult hill_or_slab = op::boolean_op(hill, slab, Op::Union);
                check_result(hill_or_slab);
                AK_CHECK_NEAR(v.uni, op::volume(hill_or_slab.solid), 1e-11);
                AK_CHECK(v.r_inter.empty);
                AK_CHECK_NEAR(v.diff, op::volume(hill_cut), 1e-11);
            }
        }
    }

    return akt::finish("boolean_test");
}
