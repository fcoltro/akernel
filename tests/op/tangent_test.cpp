// Two curved surfaces touching along a line, decided exactly
// (op/quadric_contact.hpp, and its use in op/boolean.hpp): cylinder sectors
// side by side touching along a ruling, inside one another, in one another's
// hole, and overlapping with two tangencies at once. Before, every one of
// these was refused as "not in general position".
#include <akernel/op/boolean.hpp>
#include <akernel/op/brep_transform.hpp>
#include <akernel/op/quadric_contact.hpp>
#include <akernel/op/surfaces.hpp>
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

/// A cylinder sector as make_cylinder_sector makes it, placed by a signed
/// permutation: for point membership in Monte Carlo checks.
struct Sector {
    double r0, r1, t0, t1, z0, z1;
    op::SignedPermutation m;
    Brep brep() const { return op::transformed(op::to_brep(op::make_cylinder_sector(r0, r1, t0, t1, z0, z1)), m); }
    bool contains(const Vec3d& p) const {
        // Undo m: p = s * x_perm + shift with perm the identity here.
        const double x = (p.x - m.shift.x) * m.sign[0], y = (p.y - m.shift.y) * m.sign[1], z = (p.z - m.shift.z) * m.sign[2];
        const double r = std::sqrt(x * x + y * y);
        if (r < r0 || r > r1 || z < z0 || z > z1) return false;
        // The angle's half tangent, t = y / (r + x), within [t0, t1].
        if (r + x <= 0.0) return false;
        const double t = y / (r + x);
        return t >= t0 && t <= t1;
    }
};

/// All three operations from one imprint, each result audited closed and
/// tessellated watertight; the volumes satisfy inclusion-exclusion and the
/// intersection matches Monte Carlo over the box [lo, hi]. The imprint
/// refers to the operands, so the caller keeps a and b (made by sa and sb).
bool consistent(const char* name, const Sector& sa, const Sector& sb, const Brep& a, const Brep& b, const Vec3d& lo,
                const Vec3d& hi, op::Imprint& im, double v[3], int seed) {
    im = op::imprint(a, b);
    AK_CHECK(im.ok);
    if (!im.ok) {
        std::printf("FAIL %s: %s\n", name, im.error.c_str());
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
            const op::MeshAudit ma = op::audit(t.mesh);
            AK_CHECK(ma.watertight() && ma.euler == r.audit.euler);
        }
    }
    const double va = op::volume(a), vb = op::volume(b);
    AK_CHECK_NEAR(v[0] + v[2], va, 1e-10);
    AK_CHECK_NEAR(v[0] + v[1], va + vb, 1e-10);
    akt::Rng rng(static_cast<std::uint64_t>(seed));
    const int n = 400000;
    int hits = 0;
    for (int i = 0; i < n; ++i) {
        const Vec3d p{lo.x + rng.unit() * (hi.x - lo.x), lo.y + rng.unit() * (hi.y - lo.y), lo.z + rng.unit() * (hi.z - lo.z)};
        if (sa.contains(p) && sb.contains(p)) ++hits;
    }
    const double box = (hi.x - lo.x) * (hi.y - lo.y) * (hi.z - lo.z), p = static_cast<double>(hits) / n;
    const double mc = box * p, sigma = std::max(box * std::sqrt(p * (1.0 - p) / n), 1e-12);
    AK_CHECK(std::fabs(v[0] - mc) < 4.0 * sigma + 1e-12);
    std::printf("%s: and %.12f (Monte Carlo %.6f +- %.6f), or %.12f, minus %.12f; touching lines %d\n", name, v[0], mc, sigma,
                v[1], v[2], im.stats.curved_line_tangencies);
    return true;
}

op::SignedPermutation place(double sx, double sy, const Vec3d& shift) {
    op::SignedPermutation m;
    m.sign = {sx, sy, 1.0};
    m.shift = shift;
    return m;
}

} // namespace

int main() {
    std::setvbuf(stdout, nullptr, _IONBF, 0);

    // --- the decision on its own -----------------------------------------------------
    {
        const Brep A = op::to_brep(op::make_cylinder_sector(0.75, 1.25, -0.25, 0.25, 0.0, 1.0));
        // The outer cylinder of a sector about (cx, 0): its patch of degree
        // (2, 1) at distance 5/4 from the axis.
        auto outer = [](const Brep& B, double cx) {
            for (const auto& c : B.carriers) {
                const auto h = c.at(0, 0);
                const double x = h.x / h.w - cx, y = h.y / h.w;
                if (c.p == 2 && c.q == 1 && std::fabs(x * x + y * y - 1.5625) < 1e-12) return c;
            }
            AK_CHECK(false);
            return B.carriers[0];
        };
        const geom::PolyPatchD X = outer(A, 0.0);
        {
            op::Quadric Q;
            AK_CHECK(op::implicit_quadric(X, Q));
            // x^2 + y^2 - r^2 w^2, up to a constant: no z, no cross terms.
            AK_CHECK(Q.A[2][2] == 0.0 && Q.A[0][1] == 0.0 && Q.A[0][0] == Q.A[1][1] && Q.A[0][0] != 0.0);
            AK_CHECK(Q.A[3][3] == -1.5625 * Q.A[0][0]);
        }
        // The same sector turned half round about (1.25, 0): touching along
        // the ruling at angle 0, which is u = 1/2 on both.
        for (double gap : {0.0, std::ldexp(1.0, -30), -std::ldexp(1.0, -30), -0.125}) {
            const Brep B = op::transformed(op::to_brep(op::make_cylinder_sector(0.75, 1.25, -0.25, 0.25, 0.0, 1.0)),
                                           place(-1.0, -1.0, {2.5 + gap, 0.0, 0.0}));
            const geom::PolyPatchD Y = outer(B, 2.5 + gap);
            const op::QuadricLineContact c = op::quadric_line_contact(X, Y);
            if (gap == 0.0) {
                AK_CHECK(c.found && c.fixed_u && c.r.size() == 1 && c.r[0].lo == 0.5 && c.r[0].hi == 0.5);
            } else {
                // Apart by 2^-30, or crossing: no contact is claimed.
                AK_CHECK(!c.found);
            }
        }
        // A contact at a parameter that is not dyadic: the sector over
        // t in [1/4, 1] meets t = 1/2 (angle 2 atan(1/2), direction (3/5, 4/5))
        // at u = 1/3, and the sector over [-5/2, -7/4] about (3/2, 2) at u = 2/3.
        {
            const Brep A2 = op::to_brep(op::make_cylinder_sector(0.75, 1.25, 0.25, 1.0, 0.0, 1.0));
            const Brep B2 = op::transformed(op::to_brep(op::make_cylinder_sector(0.75, 1.25, -2.5, -1.75, 0.0, 1.0)),
                                            place(1.0, 1.0, {1.5, 2.0, 0.0}));
            int found = 0;
            for (const auto& x : A2.carriers)
                for (const auto& y : B2.carriers) {
                    const op::QuadricLineContact c = op::quadric_line_contact(x, y);
                    if (!c.found) continue;
                    ++found;
                    AK_CHECK(c.r.size() == 1);
                    // One of 1/3 (on A2) or 2/3 (on B2), enclosed tightly.
                    const double r = num::mid(c.r[0]);
                    AK_CHECK(std::fabs(r - 1.0 / 3.0) < 1e-14 || std::fabs(r - 2.0 / 3.0) < 1e-14);
                    AK_CHECK(c.r[0].hi - c.r[0].lo < 1e-14);
                }
            AK_CHECK(found == 1);
            const op::Imprint im = op::imprint(A2, B2);
            AK_CHECK(im.ok && im.stats.curved_line_tangencies == 1);
            if (!im.ok) std::printf("FAIL contact at u = 1/3: %s\n", im.error.c_str());
        }
    }

    // --- Booleans ----------------------------------------------------------------------
    const Sector A{0.75, 1.25, -0.25, 0.25, 0.0, 1.0, op::SignedPermutation{}};
    const Vec3d lo{-0.5, -1.0, -0.5}, hi{3.75, 1.0, 1.5};
    op::Imprint im;
    double v[3];
    // (C1) Side by side: the same sector turned half round about (1.25, 0, 0)
    // and raised: the two outer cylinders touch along x = 1.25, y = 0.
    {
        const Sector B{0.75, 1.25, -0.25, 0.25, 0.0, 0.5, place(-1.0, -1.0, {2.5, 0.0, 0.25})};
        const Brep oa = A.brep(), ob = B.brep();
        if (consistent("side by side", A, B, oa, ob, lo, hi, im, v, 11)) {
            AK_CHECK(im.stats.curved_line_tangencies == 1);
            AK_CHECK_NEAR(v[0], 0.0, 1e-12);
            const op::BooleanResult u = op::select(im, Op::Union);
            AK_CHECK(u.ok && u.line_contacts == 1 && !u.manifold);
        }
    }
    // (C2) The same, taller than A: A's flat ends cut B's cylinder through
    // the touching line.
    {
        const Sector B{0.75, 1.25, -0.25, 0.25, 0.0, 2.0, place(-1.0, -1.0, {2.5, 0.0, -0.5})};
        const Brep oa = A.brep(), ob = B.brep();
        if (consistent("side by side, taller", A, B, oa, ob, lo, hi, im, v, 12)) {
            AK_CHECK(im.stats.curved_line_tangencies == 1);
            AK_CHECK_NEAR(v[0], 0.0, 1e-12);
        }
    }
    // (C3) Inside one another: B's inner cylinder (r = 3/2 about (-1/4, 0))
    // holds A's outer one (r = 5/4) and touches it along x = 5/4, y = 0.
    {
        const Sector B{1.5, 2.0, -0.25, 0.25, 0.0, 0.5, place(1.0, 1.0, {-0.25, 0.0, 0.25})};
        const Brep oa = A.brep(), ob = B.brep();
        if (consistent("inside one another", A, B, oa, ob, lo, hi, im, v, 13)) {
            AK_CHECK(im.stats.curved_line_tangencies == 1);
            AK_CHECK_NEAR(v[0], 0.0, 1e-12);
        }
    }
    // (C4) In A's hole: B's outer cylinder (r = 1/2 about (1/4, 0)) touches
    // A's inner one (r = 3/4) from inside along x = 3/4, y = 0.
    {
        const Sector B{0.25, 0.5, -0.25, 0.25, 0.0, 0.5, place(1.0, 1.0, {0.25, 0.0, 0.25})};
        const Brep oa = A.brep(), ob = B.brep();
        if (consistent("in the hole", A, B, oa, ob, lo, hi, im, v, 14)) {
            AK_CHECK(im.stats.curved_line_tangencies == 1);
            AK_CHECK_NEAR(v[0], 0.0, 1e-12);
        }
    }
    // (C5) Overlapping, with two tangencies: B about (-1/4, 0) with radii 1
    // and 3/2 touches A's inner cylinder along x = 3/4 and its outer one
    // along x = 5/4, and holds the material of A between them near y = 0.
    {
        const Sector B{1.0, 1.5, -0.25, 0.25, 0.0, 0.5, place(1.0, 1.0, {-0.25, 0.0, 0.25})};
        const Brep oa = A.brep(), ob = B.brep();
        if (consistent("overlapping, touching twice", A, B, oa, ob, lo, hi, im, v, 15)) {
            AK_CHECK(im.stats.curved_line_tangencies == 2);
            AK_CHECK(v[0] > 0.2);
        }
    }
    // (C6) A contact at a parameter that is not dyadic (u = 1/3 and 2/3).
    {
        const Sector A2{0.75, 1.25, 0.25, 1.0, 0.0, 1.0, op::SignedPermutation{}};
        const Sector B2{0.75, 1.25, -2.5, -1.75, 0.25, 0.75, place(1.0, 1.0, {1.5, 2.0, 0.0})};
        const Brep oa = A2.brep(), ob = B2.brep();
        if (consistent("touching at a third", A2, B2, oa, ob, {-0.5, -0.5, -0.5}, {3.0, 3.5, 1.5}, im, v, 16)) {
            AK_CHECK(im.stats.curved_line_tangencies == 1);
            AK_CHECK_NEAR(v[0], 0.0, 1e-12);
        }
    }
    return akt::finish("tangent_test");
}
