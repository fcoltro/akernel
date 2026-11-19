// Multi-patch faces: a box whose top is a quadratic B-spline graph, split
// into its Bezier pieces (one face each, C1 across the seams), and the
// bottom and sides below them likewise. Booleans with boxes, a cylinder
// sector and a second such box, checked by inclusion-exclusion, against a
// Monte Carlo count from the analytic definition, and by closed meshes.
#include <akernel/op/boolean.hpp>
#include <akernel/op/tessellate.hpp>

#include <cmath>
#include <cstdio>
#include <functional>

#include "check.hpp"

using namespace ak;
using geom::Vec3d;
using op::Brep;
using op::Op;

namespace {

struct Spline {
    Vec3d o;
    double sx, sy, sz;
    int m, n;
    std::vector<double> h;
    geom::NurbsSurface top;
    Spline(Vec3d o_, double sx_, double sy_, double sz_, int m_, int n_, std::vector<double> h_)
        : o(o_), sx(sx_), sy(sy_), sz(sz_), m(m_), n(n_), h(std::move(h_)),
          top(2, 2, m_ + 2, n_ + 2, geom::open_uniform_knots(m_ + 2, 2), geom::open_uniform_knots(n_ + 2, 2)) {
        const auto gu = geom::greville(m + 2, 2, top.U), gv = geom::greville(n + 2, 2, top.V);
        for (int i = 0; i < m + 2; ++i)
            for (int j = 0; j < n + 2; ++j)
                top.at(i, j) = {gu[static_cast<std::size_t>(i)], gv[static_cast<std::size_t>(j)], h[static_cast<std::size_t>(i * (n + 2) + j)], 1.0};
    }
    op::PatchSolid solid() const { return op::make_spline_box(o, sx, sy, sz, m, n, h); }
    bool inside(const Vec3d& p) const {
        const double x = (p.x - o.x) / sx, y = (p.y - o.y) / sy;
        if (x < 0 || x > 1 || y < 0 || y > 1 || p.z < o.z) return false;
        return p.z <= o.z + sz + top.point(x, y).z;
    }
    /// Closed form: the box, plus the integral of the spline, sum h_ij I_i
    /// I_j sx sy with I_i the integral of the i-th basis function,
    /// (U[i+3] - U[i]) / 3.
    double volume() const {
        double s = 0.0;
        for (int i = 0; i < m + 2; ++i)
            for (int j = 0; j < n + 2; ++j)
                s += h[static_cast<std::size_t>(i * (n + 2) + j)] * (top.U[static_cast<std::size_t>(i + 3)] - top.U[static_cast<std::size_t>(i)]) / 3.0 *
                     (top.V[static_cast<std::size_t>(j + 3)] - top.V[static_cast<std::size_t>(j)]) / 3.0;
        return sx * sy * (sz + s);
    }
};

void run(const char* name, const Brep& A, const Brep& B, const std::function<bool(const Vec3d&)>& inA,
         const std::function<bool(const Vec3d&)>& inB, const Vec3d& lo, const Vec3d& hi, std::uint64_t seed) {
    const op::Imprint im = op::imprint(A, B);
    AK_CHECK(im.ok);
    if (!im.ok) {
        std::printf("FAIL %s: %s\n", name, im.error.c_str());
        return;
    }
    const double va = op::volume(A), vb = op::volume(B);
    double v[3];
    int k = 0;
    for (Op o : {Op::Intersection, Op::Union, Op::Difference}) {
        const op::BooleanResult r = op::select(im, o);
        AK_CHECK(r.ok);
        v[k] = r.ok && !r.empty ? op::volume(r.solid) : 0.0;
        if (r.ok && !r.empty) {
            const op::TessResult t = op::tessellate(r.solid);
            AK_CHECK(t.ok && op::audit(t.mesh).watertight());
        }
        akt::Rng rng(seed + static_cast<std::uint64_t>(k));
        const int N = 100000;
        int hits = 0;
        for (int i = 0; i < N; ++i) {
            const Vec3d q{rng.range(lo.x, hi.x), rng.range(lo.y, hi.y), rng.range(lo.z, hi.z)};
            const bool a = inA(q), b = inB(q);
            if (o == Op::Intersection ? (a && b) : o == Op::Union ? (a || b) : (a && !b)) ++hits;
        }
        const double box = (hi.x - lo.x) * (hi.y - lo.y) * (hi.z - lo.z), p = static_cast<double>(hits) / N;
        const double mc = box * p, sg = box * std::sqrt(std::max(p * (1 - p), 1.0 / N) / N);
        AK_CHECK(std::fabs(v[k] - mc) < 4.0 * sg);
        std::printf("  %s, op %d: %.9f, Monte Carlo %.6f (%.2f sigma), %d faces\n", name, static_cast<int>(o), v[k], mc,
                    std::fabs(v[k] - mc) / sg, r.ok ? static_cast<int>(r.solid.faces.size()) : -1);
        ++k;
    }
    AK_CHECK_NEAR(v[0] + v[2], va, 1e-12);
    AK_CHECK_NEAR(v[0] + v[1], va + vb, 1e-12);
}

} // namespace

int main() {
    const Spline s({0, 0, 0}, 1, 1, 0.5, 2, 2, {0, 0, 0, 0, 0, 0.25, 0.25, 0, 0, 0.25, 0.25, 0, 0, 0, 0, 0});
    const op::PatchSolid ps = s.solid();
    // 2 x 2 pieces on top and bottom, two on each side: 16 faces.
    AK_CHECK(ps.faces.size() == 16u);
    AK_CHECK_NEAR(op::volume(ps), s.volume(), 1e-12);
    AK_CHECK_NEAR(s.volume(), 0.5 + 0.25 * 4.0 / 9.0, 1e-15);
    {
        // 4 x 2 spans with heights on the border too (curved top edges of the
        // sides).
        std::vector<double> h;
        for (int i = 0; i < 6 * 4; ++i) h.push_back(0.0625 * static_cast<double>((i * 7) % 5));
        const Spline t({0.25, -0.5, 0.0}, 0.75, 1.0, 0.25, 4, 2, h);
        const op::PatchSolid pt = t.solid();
        AK_CHECK(pt.faces.size() == static_cast<std::size_t>(2 * 4 * 2 + 2 * (4 + 2)));
        AK_CHECK_NEAR(op::volume(pt), t.volume(), 1e-12);
    }
    const Brep A = op::to_brep(ps);
    auto inA = [&](const Vec3d& p) { return s.inside(p); };
    // A box across both seams of the top (x = 0.5 and y = 0.5).
    run("box across both seams", A, op::to_brep(op::make_box({0.25, -0.5, 0.25}, {0.75, 1.5, 1.0})), inA,
        [](const Vec3d& p) { return p.x >= 0.25 && p.x <= 0.75 && p.y >= -0.5 && p.y <= 1.5 && p.z >= 0.25 && p.z <= 1.0; },
        {-0.1, -0.6, -0.1}, {1.1, 1.6, 1.1}, 101);
    // A box through the top: closed curves over four pieces.
    run("box through the top", A, op::to_brep(op::make_box({0.3, 0.3, 0.55}, {0.8, 0.9, 0.7})), inA,
        [](const Vec3d& p) { return p.x >= 0.3 && p.x <= 0.8 && p.y >= 0.3 && p.y <= 0.9 && p.z >= 0.55 && p.z <= 0.7; },
        {-0.1, -0.1, -0.1}, {1.1, 1.1, 1.1}, 202);
    // A sector whose inner cylinder r = 0.5 is tangent to the plane x = 0.5
    // of the top's seam: the seam edge lies in that plane, which the cylinder
    // meets only along its first ruling (in_plane_meeting_sides).
    run("sector through the top", A, op::to_brep(op::make_cylinder_sector(0.5, 1.0, 0.0, 1.0, 0.25, 1.0)), inA,
        [](const Vec3d& p) {
            const double r = std::sqrt(p.x * p.x + p.y * p.y);
            return r >= 0.5 && r <= 1.0 && p.x >= 0 && p.y >= 0 && p.z >= 0.25 && p.z <= 1.0;
        },
        {-0.1, -0.1, -0.1}, {1.1, 1.1, 1.1}, 303);
    // Two spline boxes with coplanar bottoms whose seams differ: a region of
    // one bottom whose first interior point lies on a seam of the other.
    {
        const Spline t({0.375, 0.25, 0.0}, 1, 1, 0.5, 2, 2, {0, 0, 0, 0, 0, 0.25, 0.125, 0, 0, 0.125, 0.25, 0, 0, 0, 0, 0});
        run("two spline boxes", A, op::to_brep(t.solid()), inA, [&](const Vec3d& p) { return t.inside(p); }, {-0.1, -0.1, -0.1},
            {1.5, 1.4, 1.1}, 404);
    }
    // The quarter shell r in [1, 2] against a spline box whose bottom seam
    // y = 1, z = 0.25 is tangent to the cylinder r = 1 at (0, 1, 0.25), inside
    // the seam: a straight edge, in the plane y = 1 through it, which the
    // cylinder meets only along its last ruling (in_plane_meeting_sides, a
    // plane through the line and an axis).
    {
        const Spline t({-0.5, 0.5, 0.25}, 1, 1, 0.5, 2, 2, {0, 0, 0, 0, 0, 0.25, 0.125, 0, 0, 0.125, 0.25, 0, 0, 0, 0, 0});
        run("shell and a spline box tangent inside a seam", op::to_brep(op::make_cylinder_sector(1.0, 2.0, 0.0, 1.0, 0.0, 1.0)),
            op::to_brep(t.solid()),
            [](const Vec3d& p) {
                const double r = std::sqrt(p.x * p.x + p.y * p.y);
                return r >= 1.0 && r <= 2.0 && p.x >= 0 && p.y >= 0 && p.z >= 0.0 && p.z <= 1.0;
            },
            [&](const Vec3d& p) { return t.inside(p); }, {-0.6, -0.1, -0.1}, {2.1, 2.1, 1.1}, 505);
    }
    // A sector over angles 0 to 2 atan(1/4), r in [1.25, 2.25], z in [0.25,
    // 0.75], against a spline box whose top corner (1.25, 0.5, 0.25) lies
    // exactly on the sector's bottom cap (a flat face on the rational annulus
    // carrier), inside it; the box's vertical edge there ends on the cap. The
    // cap is moved onto an exact parallelogram of its plane first (a vertex of
    // the other solid in its plane), where the end is decided exactly
    // (benchmark section 13, chain 8, step 5).
    {
        const Spline t({0.75, 0.0, -0.25}, 0.5, 0.5, 0.5, 2, 2,
                       {0.25, 0, 0.125, 0, 0.0625, 0.125, 0.125, 0.25, 0.0625, 0.25, 0.25, 0.0625, 0.0625, 0.125, 0.125, 0});
        const double t1 = 0.25;
        run("spline box with a corner on a sector's cap", op::to_brep(op::make_cylinder_sector(1.25, 2.25, 0.0, t1, 0.25, 0.75)),
            op::to_brep(t.solid()),
            [&](const Vec3d& p) {
                const double r = std::sqrt(p.x * p.x + p.y * p.y), th = std::atan2(p.y, p.x);
                return r >= 1.25 && r <= 2.25 && th >= 0.0 && th <= 2.0 * std::atan(t1) && p.z >= 0.25 && p.z <= 0.75;
            },
            [&](const Vec3d& p) { return t.inside(p); }, {0.6, -0.1, -0.35}, {2.35, 1.2, 0.85}, 606);
    }
    return akt::finish("multipatch_test");
}
