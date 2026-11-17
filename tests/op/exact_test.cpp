// Exact points (op/boolean.hpp, "exact points"; op/surfaces.hpp,
// exact_point_at). A root whose parameters are dyadic is known exactly, and
// a vertex with an exact point is decided exactly against a surface it lies
// on: the coincidences of models with round dimensions.
#include <akernel/op/boolean.hpp>
#include <akernel/op/surfaces.hpp>
#include <akernel/op/tessellate.hpp>

#include <cmath>
#include <cstdio>

#include "check.hpp"

using namespace ak;
using geom::Vec3d;
using geom::Vec4d;
using op::Brep;
using op::Op;

namespace {

const double kPi = 3.14159265358979323846;

/// Integral of sqrt(R^2 - x^2) from 0 to x.
double seg(double R, double x) { return 0.5 * (x * std::sqrt(R * R - x * x) + R * R * std::asin(x / R)); }

/// All three operations from one imprint: each result passes the audit and
/// tessellates closed, and the three volumes, from three separately assembled
/// results, satisfy inclusion-exclusion against the operands' volumes. `v`
/// receives the volumes of intersection, union and difference. False if the
/// operation was not solved.
bool consistent(const char* name, const Brep& a, const Brep& b, op::Imprint& im, double v[3]) {
    im = op::imprint(a, b);
    AK_CHECK(im.ok);
    if (!im.ok) {
        std::printf("FAIL %s: %s\n", name, im.error.c_str());
        op::BooleanOptions traced;
        traced.trace = true;
        for (const std::string& line : op::imprint(a, b, traced).trace) std::printf("  %s\n", line.c_str());
        return false;
    }
    v[0] = v[1] = v[2] = 0.0;
    const Op ops[3] = {Op::Intersection, Op::Union, Op::Difference};
    for (int i = 0; i < 3; ++i) {
        const op::BooleanResult r = op::select(im, ops[i]);
        AK_CHECK(r.ok);
        if (!r.ok) { std::printf("FAIL %s: %s\n", name, r.error.c_str()); return false; }
        if (r.empty) continue;
        v[i] = op::volume(r.solid);
        AK_CHECK(r.audit.closed());
        const op::TessResult t = op::tessellate(r.solid);
        AK_CHECK(t.ok);
        if (t.ok) {
            const op::MeshAudit m = op::audit(t.mesh);
            AK_CHECK(m.watertight() && m.euler == r.audit.euler);
        }
    }
    const double va = op::volume(a), vb = op::volume(b);
    AK_CHECK_NEAR(v[0] + v[2], va, 1e-10);
    AK_CHECK_NEAR(v[0] + v[1], va + vb, 1e-10);
    std::printf("%s: and %.12f, or %.12f, minus %.12f; exact roots %d, roots at an exact end vertex %d, edges in a "
                "curved surface %d\n",
                name, v[0], v[1], v[2], im.stats.exact_roots, im.stats.exact_ends, im.stats.iso_edges);
    return true;
}

} // namespace

int main() {
    std::setvbuf(stdout, nullptr, _IONBF, 0);

    // --- candidates ----------------------------------------------------------------
    {
        double a = 0.0, b = 0.0;
        AK_CHECK(op::rational_near(1.0 / 3.0, a, b) && a == 1.0 && b == 3.0);
        AK_CHECK(op::rational_near(0.5, a, b) && a == 1.0 && b == 2.0);
        AK_CHECK(op::rational_near(5.0 / 7.0 + 3e-14, a, b) && a == 5.0 && b == 7.0);
        AK_CHECK(op::rational_near(-0.25, a, b) && a == -1.0 && b == 4.0);
        AK_CHECK(op::rational_near(2.0, a, b) && a == 2.0 && b == 1.0);
        AK_CHECK(!op::rational_near(kPi, a, b));
        AK_CHECK(op::snap_dyadic(0.375 + 1e-13) == 0.375);
    }

    // --- the exact evaluator against the ordinary one ------------------------------
    const Brep sector = op::to_brep(op::make_cylinder_sector(1.25, 2.0, 0.0, 1.0, 0.0, 1.0));
    {
        const Brep hill = op::to_brep(op::make_bump_box({0.25, 0.5, 0.0}, 1.0, 0.75, 0.5, 0.3));
        const double params[][4] = {{1, 2, 1, 4}, {0, 1, 1, 1}, {1, 3, 2, 3}, {3, 8, 5, 7}, {1, 1, 1, 2}};
        int tested = 0, representable = 0;
        for (const Brep* s : {&sector, &hill})
            for (const geom::PolyPatchD& P : s->carriers)
                for (const auto& q : params) {
                    ++tested;
                    Vec4d e;
                    if (!op::exact_point_at(P, q[0], q[1], q[2], q[3], e)) continue;
                    ++representable;
                    const Vec4d d = P.eval(q[0] / q[1], q[2] / q[3]);
                    AK_CHECK_NEAR(e.x / e.w, d.x / d.w, 1e-14);
                    AK_CHECK_NEAR(e.y / e.w, d.y / d.w, 1e-14);
                    AK_CHECK_NEAR(e.z / e.w, d.z / d.w, 1e-14);
                }
        // Not every exact point is a vector of doubles (sevenths on a degree-2
        // patch need more bits); such a point is simply not used.
        AK_CHECK(representable > tested * 3 / 4);
        std::printf("exact evaluation: %d points, %d representable\n", tested, representable);
        // A point of the circle r = 1.25 that is rational but not dyadic in
        // the parameter: (1, 0.75) at t = 1/3, exactly on the circle.
        bool found = false;
        for (const geom::PolyPatchD& P : sector.carriers) {
            Vec4d e;
            if (!op::exact_point_at(P, 1, 3, 0, 1, e) && !op::exact_point_at(P, 0, 1, 1, 3, e)) continue;
            for (int swap = 0; swap < 2 && !found; ++swap) {
                if (!(swap == 0 ? op::exact_point_at(P, 1, 3, 0, 1, e) : op::exact_point_at(P, 0, 1, 1, 3, e))) continue;
                found = op::planar_detail::same_point_h(e, {1.0, 0.75, 0.0, 1.0}) ||
                        op::planar_detail::same_point_h(e, {1.0, 0.75, 1.0, 1.0});
            }
            if (found) break;
        }
        AK_CHECK(found);
    }

    // --- a vertex with an exact point, in a face of the next operand ----------------
    {
        // A notch in the sector: its edge along y (x = 0.75, z = 0.5) pierces
        // the cylinder r = 1.25 at (0.75, 1, 0.5), a 3-4-5 point. All its
        // parameters are dyadic, so the new vertex is known exactly.
        const Brep notch = op::to_brep(op::make_box({0.75, -1.0, 0.5}, {1.75, 3.0, 2.0}));
        const op::Imprint im1 = op::imprint(sector, notch);
        AK_CHECK(im1.ok);
        const op::BooleanResult r1 = op::select(im1, Op::Difference);
        AK_CHECK(r1.ok && !r1.empty);
        if (r1.ok && !r1.empty) {
            AK_CHECK(im1.stats.exact_roots > 0);
            const double area = seg(2.0, 1.75) - seg(2.0, 0.75) - (seg(1.25, 1.25) - seg(1.25, 0.75));
            const double v1 = 0.25 * kPi * (4.0 - 1.5625) - 0.5 * area;
            AK_CHECK_NEAR(op::volume(r1.solid), v1, 1e-11);
            std::printf("notched sector: volume %.12f (closed form %.12f), exact roots %d of %d new vertices\n",
                        op::volume(r1.solid), v1, im1.stats.exact_roots, im1.stats.new_vertices);
            // A tilted box whose face y + z = 1.5 passes exactly through that
            // vertex, across the three edges that end there. The first of them
            // gives a root that cannot be separated from its end; the vertex's
            // exact point lies exactly in the plane, so it is that root, and
            // the vertex gets parameters on the plane, which settle the other
            // two edges. Without exact points the operation is refused ("the
            // surface passes through a vertex of the edge").
            // (Its other five faces lie outside the sector.)
            const Brep tilted = op::to_brep(op::make_box({-0.75, -2.0, 3.5}, {3.0, 0.0, 0.0}, {0.0, 4.0, -4.0}, {0.0, 1.0, 1.0}));
            op::Imprint im2;
            double v2[3];
            consistent("notched sector and a face through its exact vertex", r1.solid, tilted, im2, v2);
            AK_CHECK(im2.ok && im2.stats.exact_ends > 0);
        }
    }
    // --- an edge in a curved surface, along one of its straight lines ---------------
    {
        // Boxes with a vertical edge through (0.75, 1) or (1, 0.75), points
        // of the circle r = 1.25, so that the edge lies in the sector's inner
        // cylinder along a line of constant angle (parameter 1/2 or 1/3).
        const double sv = 0.25 * kPi * (4.0 - 1.5625);
        struct Case {
            const char* name;
            Vec3d lo, hi;
            double in; // volume of the box inside the sector
        };
        // The part of [0.75, 1.5] x [0.5, 1] with r >= 1.25, times the height.
        const double through = 0.5 * (0.375 - (seg(1.25, 1.0) - seg(1.25, 0.5) - 0.375));
        const Case cases[] = {
            // The box is inside the sector and touches its inner wall along the edge.
            {"box against the wall, from inside the sector", {0.75, 1.0, 0.25}, {1.25, 1.5, 0.75}, 0.125},
            {"the same at a third of the circle", {1.0, 0.75, 0.25}, {1.5, 1.25, 0.75}, 0.125},
            // The box is in the hole and touches the wall along the edge from the other side.
            {"box against the wall, from the hole", {0.25, 0.5, 0.25}, {0.75, 1.0, 0.75}, 0.0},
            // The wall passes through the box's corner region: the edge is where it enters.
            {"wall through the box, along its edge", {0.75, 0.5, 0.25}, {1.5, 1.0, 0.75}, through},
        };
        for (const Case& c : cases) {
            const Brep box = op::to_brep(op::make_box(c.lo, c.hi));
            const double bv = (c.hi.x - c.lo.x) * (c.hi.y - c.lo.y) * (c.hi.z - c.lo.z);
            op::Imprint im;
            double v[3] = {0.0, 0.0, 0.0};
            if (!consistent(c.name, sector, box, im, v)) continue;
            AK_CHECK(im.stats.iso_edges > 0);
            // Where the box only touches the wall (the first three cases), the
            // result that keeps both the wall and the box's edge has its parts
            // touching along that line: reported, and not manifold.
            int contacts = 0;
            for (Op o : {Op::Intersection, Op::Union, Op::Difference}) {
                const op::BooleanResult r = op::select(im, o);
                if (!r.ok || r.empty) continue;
                contacts += r.line_contacts;
                if (r.line_contacts > 0) AK_CHECK(!r.manifold);
            }
            AK_CHECK((contacts > 0) == (c.in == 0.0 || c.in == 0.125));
            std::printf("  line contacts over the three results: %d\n", contacts);
            AK_CHECK_NEAR(v[0], c.in, 1e-11);
            AK_CHECK_NEAR(v[1], sv + bv - c.in, 1e-11);
            AK_CHECK_NEAR(v[2], sv - c.in, 1e-11);
        }
    }
    // --- one vertex reached from two edges ---------------------------------------------
    {
        // Two chains of the benchmark (§11, chains 30 and 28), whose operands
        // lie on dyadic grids. In the last step shown, the vertex (0.75, 1, z)
        // where a hill box's edge meets the sector's inner cylinder r = 1.25
        // is found again from another edge, on other surfaces. Both roots are
        // the same exact point, so they are one vertex; before, the two could
        // not be ordered along their curve and the step was refused. Each
        // result is checked by inclusion-exclusion and against a Monte Carlo
        // count that uses only the analytic definitions of the operands.
        struct Sector { double r0, r1, t0, t1, z0, z1; };
        struct Hill { Vec3d o; double sx, sy, sz, bump; };
        auto in_sector = [](const Sector& s, const Vec3d& p) {
            const double rho = std::sqrt(p.x * p.x + p.y * p.y), th = std::atan2(p.y, p.x);
            return rho >= s.r0 && rho <= s.r1 && p.z >= s.z0 && p.z <= s.z1 && th >= 2.0 * std::atan(s.t0) &&
                   th <= 2.0 * std::atan(s.t1);
        };
        auto in_hill = [](const Hill& h, const Vec3d& p) {
            const double x = (p.x - h.o.x) / h.sx, y = (p.y - h.o.y) / h.sy;
            if (x < 0 || x > 1 || y < 0 || y > 1 || p.z < h.o.z) return false;
            return p.z <= h.o.z + h.sz + 16.0 * h.bump * x * (1 - x) * y * (1 - y);
        };
        struct Step { int kind; Sector s; Hill h; Op op; }; // kind 0 sector, 1 hill box
        const Sector shell{1.0, 2.0, 0.0, 1.0, 0.0, 1.0};
        const std::vector<std::vector<Step>> chains = {
            // chain 30, steps 1 to 3
            {{1, {}, {{0.0, 0.75, 0.25}, 1.0, 0.5, 0.5, 0.2}, Op::Difference},
             {1, {}, {{0.75, 0.25, -0.25}, 1.0, 1.0, 0.75, 0.1}, Op::Difference},
             {0, {1.0, 1.25, 0.5, 0.75, 0.25, 1.25}, {}, Op::Intersection}},
            // chain 28, steps 1 to 5 (step 3 is a union beyond the shell)
            {{0, {1.0, 1.5, 0.5, 1.0, -0.25, 0.25}, {}, Op::Union},
             {1, {}, {{0.0, 0.0, -0.25}, 0.75, 0.75, 0.5, 0.3}, Op::Difference},
             {0, {1.5, 2.0, -0.25, 0.25, 0.0, 1.0}, {}, Op::Union},
             {1, {}, {{0.75, 0.25, 0.25}, 1.0, 0.75, 0.5, 0.1}, Op::Intersection},
             {0, {1.25, 1.75, 0.5, 1.0, 0.5, 0.75}, {}, Op::Difference}},
            // chain 19, steps 1 to 6: in the last, a vertex with an exact point
            // decided on a side whose parameter across it was an enclosure
            {{1, {}, {{0.0, 0.25, 0.25}, 0.5, 0.5, 0.5, 0.3}, Op::Difference},
             {1, {}, {{0.0, 0.0, 0.25}, 0.75, 1.0, 0.5, 0.1}, Op::Difference},
             {0, {0.75, 1.75, 0.25, 0.75, -0.25, 0.25}, {}, Op::Difference},
             {1, {}, {{-0.25, 0.0, 0.0}, 0.75, 0.75, 0.5, 0.3}, Op::Difference},
             {1, {}, {{0.5, 0.0, 0.25}, 1.0, 1.0, 0.75, 0.3}, Op::Union},
             {1, {}, {{0.5, 0.0, -0.25}, 0.75, 0.5, 0.75, 0.1}, Op::Intersection}},
            // chain 38, steps 1 to 3: in the last, a cap is moved onto a plane
            // and an arc on another flat carrier of that plane moves with it
            {{0, {0.75, 1.25, 0.5, 0.75, 0.25, 0.75}, {}, Op::Union},
             {0, {0.75, 1.75, 0.25, 0.75, 0.0, 0.25}, {}, Op::Intersection},
             {1, {}, {{0.0, 0.75, -0.25}, 0.5, 0.75, 0.5, 0.3}, Op::Difference}},
            // chain 23, steps 1 to 4: in the last, the sector's inner bottom
            // circle (r = 0.75, z = 0.5) starts tangent to the wall x = 0.75
            // left by the hill box, and passes over that hill's top, which
            // meets the plane z = 0.5 only along its sides: both arcs are
            // decided by sides-only contacts (no singular root is searched)
            {{0, {1.25, 1.5, -0.25, 0.0, -0.25, 0.0}, {}, Op::Difference},
             {0, {1.5, 1.75, 0.25, 0.5, 0.0, 0.25}, {}, Op::Difference},
             {1, {}, {{0.25, -0.25, -0.25}, 0.5, 1.0, 0.75, 0.1}, Op::Difference},
             {0, {0.75, 1.0, 0.0, 0.5, 0.5, 0.75}, {}, Op::Difference}},
        };
        const char* names[] = {"chain 30, step 3 (sector)", "chain 28, step 5 (sector)", "chain 19, step 6 (hill box)",
                               "chain 38, step 3 (hill box)", "chain 23, step 4 (sector)"};
        for (std::size_t ci = 0; ci < chains.size(); ++ci) {
            Brep current = op::to_brep(op::make_cylinder_sector(shell.r0, shell.r1, shell.t0, shell.t1, shell.z0, shell.z1));
            std::vector<Step> done;
            bool ok = true;
            for (std::size_t si = 0; si < chains[ci].size() && ok; ++si) {
                const Step& st = chains[ci][si];
                const Brep B = op::to_brep(st.kind == 0 ? op::make_cylinder_sector(st.s.r0, st.s.r1, st.s.t0, st.s.t1, st.s.z0, st.s.z1)
                                                        : op::make_bump_box(st.h.o, st.h.sx, st.h.sy, st.h.sz, st.h.bump));
                const bool last = si + 1 == chains[ci].size();
                op::BooleanResult r;
                if (last) {
                    op::Imprint im;
                    double v[3];
                    ok = consistent(names[ci], current, B, im, v);
                    if (!ok) break;
                    if (ci < 2) AK_CHECK(im.stats.merged_vertices > 0);
                    if (ci == 3) AK_CHECK(im.stats.recarried_faces > 0);
                    r = op::select(im, st.op);
                } else {
                    r = op::boolean_op(current, B, st.op);
                }
                AK_CHECK(r.ok && !r.empty);
                if (!r.ok || r.empty) { ok = false; break; }
                current = r.solid;
                done.push_back(st);
            }
            if (!ok) continue;
            // Monte Carlo over a box holding every operand.
            auto member = [&](const Vec3d& p) {
                bool in = in_sector(shell, p);
                for (const Step& st : done) {
                    const bool b = st.kind == 0 ? in_sector(st.s, p) : in_hill(st.h, p);
                    in = st.op == Op::Union ? (in || b) : st.op == Op::Intersection ? (in && b) : (in && !b);
                }
                return in;
            };
            // Around every operand; for chain 30, whose result is a small part of
            // its last operand (an intersection), around that operand only.
            const Vec3d lo = ci == 0 ? Vec3d{0.25, 0.75, 0.25} : Vec3d{-0.5, -1.5, -0.5};
            const Vec3d hi = ci == 0 ? Vec3d{0.8, 1.25, 1.25} : Vec3d{2.75, 2.75, 1.75};
            akt::Rng rng(4711 + ci);
            const int n = 400000;
            int hits = 0;
            for (int i = 0; i < n; ++i)
                if (member({lo.x + rng.unit() * (hi.x - lo.x), lo.y + rng.unit() * (hi.y - lo.y), lo.z + rng.unit() * (hi.z - lo.z)})) ++hits;
            const double box = (hi.x - lo.x) * (hi.y - lo.y) * (hi.z - lo.z), p = static_cast<double>(hits) / n;
            const double mc = box * p, sigma = box * std::sqrt(p * (1.0 - p) / n), vol = op::volume(current);
            AK_CHECK(std::fabs(vol - mc) < 4.0 * sigma);
            std::printf("  result volume %.9f, Monte Carlo %.6f +- %.6f (%.2f sigma)\n", vol, mc, sigma, std::fabs(vol - mc) / sigma);
        }
    }
    // --- a plane tangent to a cylinder along a ruling inside the patch ---------------
    {
        // Sectors over t in [-1/4, 1/4] (angles within about 28 degrees of the
        // x axis) and a box whose face x = 1.25 is tangent to the cylinder
        // r = 1.25 along its ruling y = 0, at the middle of the patch. Not
        // along a side, so only the line contact of op/contact.hpp decides it.
        auto sv = [](double r0, double r1, double t0, double t1, double z0, double z1) {
            return 0.5 * (2.0 * std::atan(t1) - 2.0 * std::atan(t0)) * (r1 * r1 - r0 * r0) * (z1 - z0);
        };
        // (T1) The sector r in [0.75, 1.25]: its outer cylinder touches the box
        // from inside; the box is outside the sector.
        {
            const Brep S = op::to_brep(op::make_cylinder_sector(0.75, 1.25, -0.25, 0.25, 0.0, 1.0));
            const Brep box = op::to_brep(op::make_box({1.25, -0.25, 0.25}, {2.0, 0.25, 0.75}));
            op::Imprint im;
            double v[3];
            if (consistent("box tangent to a sector's outer cylinder along a ruling", S, box, im, v)) {
                AK_CHECK(im.stats.line_tangencies > 0);
                const double s = sv(0.75, 1.25, -0.25, 0.25, 0.0, 1.0);
                AK_CHECK_NEAR(v[0], 0.0, 1e-12);
                AK_CHECK_NEAR(v[1], s + 0.1875, 1e-11);
                AK_CHECK_NEAR(v[2], s, 1e-11);
                // The union keeps both faces along the ruling: a line contact.
                const op::BooleanResult u = op::select(im, Op::Union);
                AK_CHECK(u.ok && u.line_contacts == 1 && !u.manifold);
            }
        }
        // (T2) The sector r in [1.25, 2]: the tangent plane runs through the
        // sector's material (points of x = 1.25 off the ruling have r > 1.25).
        // The box x in [0.5, 1.25] then meets the sector where r >= 1.25:
        // for each y, x from sqrt(1.5625 - y^2) to 1.25, over |y| <= 1/4, times
        // the height 1/2.
        {
            const Brep S = op::to_brep(op::make_cylinder_sector(1.25, 2.0, -0.25, 0.25, 0.0, 1.0));
            const Brep box = op::to_brep(op::make_box({0.5, -0.25, 0.25}, {1.25, 0.25, 0.75}));
            op::Imprint im;
            double v[3];
            if (consistent("box tangent to a sector's inner cylinder along a ruling", S, box, im, v)) {
                AK_CHECK(im.stats.line_tangencies > 0);
                const double s = sv(1.25, 2.0, -0.25, 0.25, 0.0, 1.0);
                const double in = 0.5 * (1.25 * 0.5 - 2.0 * seg(1.25, 0.25));
                AK_CHECK_NEAR(v[0], in, 1e-11);
                AK_CHECK_NEAR(v[1], s + 0.1875 - in, 1e-11);
                AK_CHECK_NEAR(v[2], s - in, 1e-11);
            }
        }
    }

    // --- an edge through an edge, in space -----------------------------------------
    {
        // The unit box A, and a parallelepiped B with dyadic edge vectors
        // whose edge from (0.75, 0.5, 0.75) along (0.5, 0, 0.5) passes
        // through A's edge x = 1, z = 1 at (1, 0.5, 1). B's two faces there
        // are the planes y = 0.5 and x + y - z = 0.5; neither contains A's
        // edge, and A's faces there (x = 1, z = 1) do not contain B's. The
        // crossing lies on four planes. Volumes are checked by
        // inclusion-exclusion and against Monte Carlo with the analytic
        // definition of B.
        const Brep A = op::to_brep(op::make_box({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}));
        const Vec3d o{0.75, 0.5, 0.75}, e1{0.5, 0.0, 0.5}, e2{0.0, 0.5, 0.5}, e3{-0.5, 0.0, 0.5};
        const Brep B = op::to_brep(op::make_box(o, e1, e2, e3));
        op::Imprint im;
        double v[3];
        if (consistent("box and a parallelepiped whose edge crosses its edge", A, B, im, v)) {
            AK_CHECK(im.stats.edge_crossings > 0);
            // B's volume |det(e1, e2, e3)| = 0.25. In B's coordinates
            // (a, b, c), x <= 1 is a - c <= 1/2 and z <= 1 is a + b + c <= 1/2,
            // which implies the first: A and B is that simplex, of volume
            // (1/2)^3 / 6 = 1/48 in (a, b, c), so 1/192.
            AK_CHECK_NEAR(op::volume(B), 0.25, 1e-12);
            AK_CHECK_NEAR(v[0], 1.0 / 192.0, 1e-12);
            AK_CHECK_NEAR(v[1], 1.25 - 1.0 / 192.0, 1e-12);
            AK_CHECK_NEAR(v[2], 1.0 - 1.0 / 192.0, 1e-12);
            // Monte Carlo of A and B over B's bounding box.
            auto in_B = [&](const Vec3d& p) {
                // p = o + a e1 + b e2 + c e3, all of a, b, c in [0, 1].
                const Vec3d d = p - o;
                const double det = e1.x * (e2.y * e3.z - e2.z * e3.y) - e1.y * (e2.x * e3.z - e2.z * e3.x) +
                                   e1.z * (e2.x * e3.y - e2.y * e3.x);
                const double a = (d.x * (e2.y * e3.z - e2.z * e3.y) - d.y * (e2.x * e3.z - e2.z * e3.x) +
                                  d.z * (e2.x * e3.y - e2.y * e3.x)) / det;
                const double b = (e1.x * (d.y * e3.z - d.z * e3.y) - e1.y * (d.x * e3.z - d.z * e3.x) +
                                  e1.z * (d.x * e3.y - d.y * e3.x)) / det;
                const double c = (e1.x * (e2.y * d.z - e2.z * d.y) - e1.y * (e2.x * d.z - e2.z * d.x) +
                                  e1.z * (e2.x * d.y - e2.y * d.x)) / det;
                return a >= 0 && a <= 1 && b >= 0 && b <= 1 && c >= 0 && c <= 1;
            };
            const Vec3d lo{0.25, 0.5, 0.75}, hi{1.25, 1.0, 1.75};
            akt::Rng rng(31337);
            const int n = 400000;
            int hits = 0;
            for (int i = 0; i < n; ++i) {
                const Vec3d p{lo.x + rng.unit() * (hi.x - lo.x), lo.y + rng.unit() * (hi.y - lo.y), lo.z + rng.unit() * (hi.z - lo.z)};
                if (p.x <= 1.0 && p.z <= 1.0 && in_B(p)) ++hits;
            }
            const double box = (hi.x - lo.x) * (hi.y - lo.y) * (hi.z - lo.z), p = static_cast<double>(hits) / n;
            const double mc = box * p, sigma = box * std::sqrt(p * (1.0 - p) / n);
            AK_CHECK(std::fabs(v[0] - mc) < 4.0 * sigma);
            std::printf("  A and B: %.9f, Monte Carlo %.6f +- %.6f (%.2f sigma); crossings %d\n", v[0], mc, sigma,
                        std::fabs(v[0] - mc) / sigma, im.stats.edge_crossings);
        }
    }

    // --- two hill tops crossing along a line they share ------------------------------
    {
        // Two hill boxes with the same top height 0.5 whose tops share the
        // line y = 0.5, x in [0, 1]: both bumps vanish there. A's top is
        // 0.5 + 0.2 * 16 a(1-a) b(1-b) over x = a, y = 0.5 + 0.5 b; B's is
        // 0.5 + 0.1 * 16 c(1-c) d(1-d) over x = c, y = 0.25 + 0.5 d. Near the
        // line both rise with y; the tops cross along a curve that leaves the
        // shared line, found by dividing the line out of their height
        // difference exactly (graph_crossing).
        const Brep A = op::to_brep(op::make_bump_box({0.5, 0.5, -0.25}, 0.5, 0.5, 0.75, 0.2));
        const Brep B = op::to_brep(op::make_bump_box({0.5, 0.25, 0.0}, 0.5, 0.5, 0.5, 0.1));
        op::Imprint im;
        double v[3];
        if (consistent("two hill tops crossing along a shared line", A, B, im, v)) {
            AK_CHECK(im.stats.graph_crossings > 0);
            // Closed form. Below z = 0.5 the two meet in the slab
            // x in [0.5, 1], y in [0.5, 0.75], z in [0, 0.5]: 0.0625. Above,
            // over that rectangle, the lower of the two bumps. Both have the
            // x factor 16 a(1-a) (integral 16/6 times 0.5 in x); in y the
            // lower of 0.2 b(1-b) and 0.1 d(1-d), which cross at
            // b = c := 1 - sqrt(0.75); I is that y integral, by parts.
            const double c = 1.0 - std::sqrt(0.75);
            const double I = 3.2 * (c * c / 2.0 - c * c * c / 3.0) + 1.6 * (0.25 * (0.5 - c) - (0.125 - c * c * c) / 3.0);
            const double inter = 0.0625 + 0.25 * (1.0 / 6.0) * I;
            const double va = 0.5 * 0.5 * 0.75 + 0.2 * 16.0 / 36.0 * 0.25, vb = 0.5 * 0.5 * 0.5 + 0.1 * 16.0 / 36.0 * 0.25;
            AK_CHECK_NEAR(op::volume(A), va, 1e-12);
            AK_CHECK_NEAR(op::volume(B), vb, 1e-12);
            AK_CHECK_NEAR(v[0], inter, 1e-11);
            AK_CHECK_NEAR(v[1], va + vb - inter, 1e-11);
            AK_CHECK_NEAR(v[2], va - inter, 1e-11);
            std::printf("  closed form: and %.12f, or %.12f, minus %.12f; graph crossings %d\n", inter, va + vb - inter, va - inter,
                        im.stats.graph_crossings);
        }
    }
    // --- two hill tops whose crossing curve touches the shared line --------------------
    {
        // A over x in [0.5, 1], y in [0.5, 1], bump 0.25; B over x in
        // [0.25, 1.25], y in [0.5, 1.25], bump 0.375. Along y = 0.5 both rise
        // with slope 6.4 x'(1-x') / ... equal at x = 0.75 (the middle of both)
        // and nowhere else: their difference has a double root there, and the
        // crossing curve touches the shared line instead of leaving through
        // it. Dyadic data, so decided exactly and refused with that reason.
        const Brep A = op::to_brep(op::make_bump_box({0.5, 0.5, -0.25}, 0.5, 0.5, 0.75, 0.25));
        const Brep B = op::to_brep(op::make_bump_box({0.25, 0.5, -0.25}, 1.0, 0.75, 0.75, 0.375));
        const op::Imprint im = op::imprint(A, B);
        AK_CHECK(!im.ok);
        AK_CHECK(im.error.find("touches the boundary line they share") != std::string::npos);
        std::printf("crossing curve tangent to the shared line: %s\n", im.ok ? "solved (unexpected)" : im.error.c_str());
    }
    return akt::finish("exact_test");
}
