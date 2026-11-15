// Exact quadrics: patches that lie exactly on cylinders, spheres and cones,
// and Boolean operations between solids whose curved faces lie on the same
// cylinder without being the same patch (decided exactly, op/surfaces.hpp).
// Every volume is checked against a closed form.
#include <akernel/geom/primitives.hpp>
#include <akernel/op/boolean.hpp>
#include <akernel/op/compare.hpp>
#include <akernel/op/surfaces.hpp>
#include <akernel/op/tessellate.hpp>

#include <cmath>
#include <cstdio>

#include "check.hpp"

using namespace ak;
using geom::PolyPatchD;
using geom::Vec4d;
using op::Brep;
using op::Op;
using op::Affine;
using num::Interval;

namespace {

/// Homogeneous Bezier form of the unit circle (1 - t^2, 2t, 1 + t^2) on [0, 1].
const Vec4d kCircle[3] = {{1, 0, 0, 1}, {1, 1, 0, 1}, {0, 2, 0, 2}};

/// Volume of make_cylinder_sector(r0, r1, t0, t1, z0, z1).
double sector_volume(double r0, double r1, double t0, double t1, double z0, double z1) {
    return 0.5 * (2.0 * std::atan(t1) - 2.0 * std::atan(t0)) * (r1 * r1 - r0 * r0) * (z1 - z0);
}

Brep sector(double r0, double r1, double t0, double t1, double z0, double z1) {
    return op::to_brep(op::make_cylinder_sector(r0, r1, t0, t1, z0, z1));
}

/// Run an operation, check volume, audit, tessellation; return the result.
Brep run(const char* name, const Brep& a, const Brep& b, Op o, double volume) {
    const op::BooleanResult r = op::boolean_op(a, b, o);
    AK_CHECK(r.ok);
    if (!r.ok) {
        std::printf("FAIL %s: %s\n", name, r.error.c_str());
        return Brep{};
    }
    const double v = r.empty ? 0.0 : op::volume(r.solid);
    AK_CHECK_NEAR(v, volume, 1e-12);
    if (!r.empty) {
        AK_CHECK(r.audit.closed());
        const op::TessResult t = op::tessellate(r.solid);
        AK_CHECK(t.ok);
        if (t.ok) {
            const op::MeshAudit m = op::audit(t.mesh);
            AK_CHECK(m.watertight() && m.euler == r.audit.euler);
            AK_CHECK(op::mesh_deviation(t.mesh, r.solid) < 2e-3);
        }
    }
    std::printf("%s: volume %.12f (closed form %.12f), %d faces%s\n", name, v, volume,
                static_cast<int>(r.solid.faces.size()), r.empty ? ", empty" : "");
    return r.solid;
}

} // namespace

int main() {
    // --- exact incidence on quadrics -------------------------------------------------
    {
        // The factories return the intended matrices (see the note on Quadric).
        const op::Quadric c3 = op::Quadric::cylinder_z(0.5, -0.25, 2.0);
        const double want[4][4] = {{1, 0, 0, -0.5}, {0, 1, 0, 0.25}, {0, 0, 0, 0}, {-0.5, 0.25, 0, 0.25 + 0.0625 - 4.0}};
        bool exact = true;
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j) exact = exact && c3.A[i][j] == want[i][j];
        AK_CHECK(exact);
        const op::Quadric unit_cyl = op::Quadric::cylinder_z(0, 0, 1), cyl2 = op::Quadric::cylinder_z(0, 0, 2);
        AK_CHECK(op::on_quadric(geom::make_cylinder_quarter(1.0, 0.0, 1.0), unit_cyl));
        AK_CHECK(op::on_quadric(geom::make_cylinder_quarter(2.0, -1.0, 3.0), cyl2));
        AK_CHECK(!op::on_quadric(geom::make_cylinder_quarter(2.0, -1.0, 3.0), unit_cyl));
        // The usual quarter circle with middle weight sqrt(2)/2, rounded: a
        // nearby conic, not the circle.
        PolyPatchD usual(2, 1);
        const double w = std::sqrt(0.5);
        for (int j = 0; j < 2; ++j) {
            usual.at(0, j) = {1, 0, double(j), 1};
            usual.at(1, j) = {w, w, w * j, w};
            usual.at(2, j) = {0, 1, double(j), 1};
        }
        AK_CHECK(!op::on_quadric(usual, unit_cyl));
        // A sector of every circle range lies on its cylinders exactly.
        const op::PatchSolid s = op::make_cylinder_sector(1.0, 2.0, -0.75, 0.5, 0.0, 1.0);
        int on = 0;
        for (const PolyPatchD& f : s.faces) on += op::on_quadric(f, unit_cyl) || op::on_quadric(f, cyl2) ? 1 : 0;
        AK_CHECK(on == 2);
        // The old quarter-shell constructor is the sector over [0, 1].
        const op::PatchSolid old = op::make_annular_sector(1.0, 2.0, 0.0, 1.0);
        const op::PatchSolid neu = op::make_cylinder_sector(1.0, 2.0, 0.0, 1.0, 0.0, 1.0);
        bool same = old.faces.size() == neu.faces.size();
        for (std::size_t f = 0; same && f < old.faces.size(); ++f) {
            same = old.faces[f].c.size() == neu.faces[f].c.size();
            for (std::size_t i = 0; same && i < old.faces[f].c.size(); ++i)
                same = old.faces[f].c[i].x == neu.faces[f].c[i].x && old.faces[f].c[i].y == neu.faces[f].c[i].y &&
                       old.faces[f].c[i].z == neu.faces[f].c[i].z && old.faces[f].c[i].w == neu.faces[f].c[i].w;
        }
        AK_CHECK(same);

        // Sphere: the product of two circles, (x_u X_v, y_u X_v, w_u Y_v, w_u W_v).
        PolyPatchD sph(2, 2);
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) {
                const Vec4d& a = kCircle[i];
                const Vec4d& b = kCircle[j];
                sph.at(i, j) = {a.x * b.x, a.y * b.x, a.w * b.y, a.w * b.w};
            }
        AK_CHECK(op::on_quadric(sph, op::Quadric::sphere(0, 0, 0, 1)));
        AK_CHECK(!op::on_quadric(sph, op::Quadric::sphere(0, 0, 0, 2)));
        // Cone x^2 + y^2 = z^2: s * (circle, 1) for s in [1, 2].
        PolyPatchD cone(2, 1);
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 2; ++j) {
                const double sc = 1.0 + j;
                const Vec4d& a = kCircle[i];
                cone.at(i, j) = {sc * a.x, sc * a.y, sc * a.w, a.w};
            }
        AK_CHECK(op::on_quadric(cone, op::Quadric::cone_z(1.0)));
        AK_CHECK(!op::on_quadric(cone, op::Quadric::cone_z(2.0)));
    }

    // --- one cylinder, two parameterisations ----------------------------------------------
    {
        // faces[1] is the outer cylinder (make_hex order: x0, x1, y0, y1, z0, z1).
        const PolyPatchD a = op::make_cylinder_sector(1, 2, 0, 1, 0, 1).faces[1];
        const PolyPatchD b = op::make_cylinder_sector(1, 2, 0.5, 1, 0.25, 0.75).faces[1];
        AK_CHECK(op::on_quadric(a, op::Quadric::cylinder_z(0, 0, 2)) && op::on_quadric(b, op::Quadric::cylinder_z(0, 0, 2)));
        Affine m;
        AK_CHECK(op::affine_reparam(a, b, m));
        // b's (u, v) = (0, 0) is a's (0.5, 0.25); (1, 1) is a's (1, 0.75).
        Interval u, v;
        m.apply(Interval(0.0), Interval(0.0), u, v);
        AK_CHECK(u.lo == 0.5 && u.hi == 0.5 && v.lo == 0.25 && v.hi == 0.25);
        m.apply(Interval(1.0), Interval(1.0), u, v);
        AK_CHECK(u.lo == 1.0 && u.hi == 1.0 && v.lo == 0.75 && v.hi == 0.75);
        // Another radius is another cylinder.
        AK_CHECK(!op::affine_reparam(a, op::make_cylinder_sector(1, 2.5, 0.5, 1, 0, 1).faces[1], m));
        // The radial faces at t = 1 (weight 2 on every corner) are flat and in
        // one plane, decided exactly.
        const PolyPatchD ra = op::make_cylinder_sector(1, 2, 0, 1, 0, 1).faces[3];
        const PolyPatchD rb = op::make_cylinder_sector(1.5, 2, 0.5, 1, 0.25, 0.75).faces[3];
        AK_CHECK(op::planar_affine(ra) && op::planar_affine(rb) && op::coplanar(ra, rb));
    }

    // --- Booleans with faces on the same cylinder ----------------------------------------
    const double VA = sector_volume(1, 2, 0, 1, 0, 1);
    const Brep A = sector(1, 2, 0, 1, 0, 1);

    // B inside A, its outer face on A's outer cylinder (r = 2), over part of it.
    const Brep B = sector(1.5, 2, 0.5, 1, 0.25, 0.75);
    const double VB = sector_volume(1.5, 2, 0.5, 1, 0.25, 0.75);
    {
        const Brep u = run("A or B, B against A's outer wall", A, B, Op::Union, VA);
        const Brep i = run("A and B", A, B, Op::Intersection, VB);
        const Brep d = run("A - B: a pocket in the outer wall", A, B, Op::Difference, VA - VB);
        if (!u.faces.empty()) AK_CHECK(op::compare_solids(u, A).same);
        if (!i.faces.empty()) AK_CHECK(op::compare_solids(i, B).same);
        if (!d.faces.empty()) {
            // Putting B back gives A again.
            const Brep back = run("(A - B) or B", d, B, Op::Union, VA);
            if (!back.faces.empty()) {
                const op::SolidComparison c = op::compare_solids(back, A);
                AK_CHECK(c.same);
                if (!c.same) std::printf("(A - B) or B vs A: %s\n", c.reason.c_str());
            }
        }
    }
    // Two sectors side by side around the axis: they share a radial plane and
    // both their cylinders.
    {
        const Brep L = sector(1, 2, 0, 0.5, 0, 1), R = sector(1, 2, 0.5, 1, 0, 1);
        const Brep u = run("left or right half of the quarter", L, R, Op::Union, VA);
        if (!u.faces.empty()) {
            const op::SolidComparison c = op::compare_solids(u, A);
            AK_CHECK(c.same);
            if (!c.same) std::printf("glued quarter vs A: %s\n", c.reason.c_str());
        }
        run("left and right (touching)", L, R, Op::Intersection, 0.0);
        run("left - right", L, R, Op::Difference, sector_volume(1, 2, 0, 0.5, 0, 1));
    }
    // Overlapping sectors of the same cylinders. The parameter relation must be
    // a dyadic shift with a power-of-two scale (here Q is P shifted by half its
    // range) for it to be proven.
    {
        const Brep P = sector(1, 2, 0, 0.5, 0, 1), Q = sector(1, 2, 0.25, 0.75, 0, 0.5);
        const double VP = sector_volume(1, 2, 0, 0.5, 0, 1), VQ = sector_volume(1, 2, 0.25, 0.75, 0, 0.5);
        const double VPQ = sector_volume(1, 2, 0.25, 0.5, 0, 0.5);
        const Brep i = run("overlapping sectors, and", P, Q, Op::Intersection, VPQ);
        run("overlapping sectors, or", P, Q, Op::Union, VP + VQ - VPQ);
        run("overlapping sectors, minus", P, Q, Op::Difference, VP - VPQ);
        if (!i.faces.empty()) {
            const op::SolidComparison c = op::compare_solids(i, sector(1, 2, 0.25, 0.5, 0, 0.5));
            AK_CHECK(c.same);
            if (!c.same) std::printf("intersection vs closed-form sector: %s\n", c.reason.c_str());
        }
    }
    // Found by the curved-chain benchmark (ak_bench §11): a sector over the
    // upper part of the right half, on the same inner and outer cylinders and
    // the same radial plane. The cut where its radial plane meets A's top runs
    // between two vertices exactly on the ends of the intersection curve.
    run("A - sector over the upper right part", A, sector(1, 2, 0.5, 1, 0.5, 1.5), Op::Difference,
        VA - sector_volume(1, 2, 0.5, 1, 0.5, 1));
    // A notch sharing four surfaces with A, each on a different patch: the
    // maps between the coincident planes must keep exact parameters exact.
    run("A - notch on four shared surfaces", A, sector(1.5, 2, 0, 0.25, 0, 1), Op::Difference,
        VA - sector_volume(1.5, 2, 0, 0.25, 0, 1));
    // The flat faces of a box with a hill on top are degree-elevated
    // parallelograms, and are recognised as flat; the hill is not.
    {
        const op::PatchSolid hb = op::make_bump_box({0, 0, 0}, 1, 1, 1, 0.3);
        int flat = 0;
        for (const PolyPatchD& f : hb.faces) flat += op::planar_affine(f) ? 1 : 0;
        AK_CHECK(flat == 5);
    }

    // A relation with scale 2/3 is not accepted: the operation is reported as
    // not solved, never answered from intervals that cannot decide.
    {
        const op::BooleanResult r = op::boolean_op(sector(1, 2, 0, 0.75, 0, 1), sector(1, 2, 0.25, 1, 0, 0.5), Op::Union);
        AK_CHECK(!r.ok);
        std::printf("scale 2/3 between two sectors of one cylinder: %s\n", r.ok ? "solved" : r.error.c_str());
    }
    // Sectors over negative circle parameters (the fourth quadrant), joined
    // with the first.
    {
        const Brep N = sector(1, 2, -0.5, 0, 0, 1);
        run("fourth-quadrant sector or first-quadrant sector", N, A, Op::Union,
            sector_volume(1, 2, -0.5, 1, 0, 1));
    }

    return akt::finish("quadric_test");
}
