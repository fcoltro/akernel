// Tangential contact with planes, decided exactly (op/contact.hpp), and
// Boolean operations in which a box face touches a cylinder along one of its
// rulings. Every volume is checked against a closed form.
#include <akernel/op/boolean.hpp>
#include <akernel/op/contact.hpp>
#include <akernel/op/surfaces.hpp>
#include <akernel/op/tessellate.hpp>

#include <cmath>
#include <cstdio>

#include "check.hpp"

using namespace ak;
using geom::PolyPatchD;
using op::Brep;
using op::Op;

namespace {

const double kPi = 3.14159265358979323846;

/// Area between the unit circle and its tangent at the top, over x in [0, a]:
/// the integral of 1 - sqrt(1 - x^2).
double cap_area(double a) { return a - 0.5 * (a * std::sqrt(1.0 - a * a) + std::asin(a)); }

/// Run an operation, check volume, audit and tessellation.
void run(const char* name, const Brep& a, const Brep& b, Op o, double volume, bool expect_manifold = true) {
    const op::BooleanResult r = op::boolean_op(a, b, o);
    AK_CHECK(r.ok);
    if (!r.ok) {
        std::printf("FAIL %s: %s\n", name, r.error.c_str());
        return;
    }
    const double v = r.empty ? 0.0 : op::volume(r.solid);
    AK_CHECK_NEAR(v, volume, 1e-11);
    if (!r.empty) {
        AK_CHECK(r.audit.closed());
        AK_CHECK(r.manifold == expect_manifold);
        const op::TessResult t = op::tessellate(r.solid);
        AK_CHECK(t.ok);
        if (t.ok) {
            const op::MeshAudit m = op::audit(t.mesh);
            AK_CHECK(m.watertight() && m.euler == r.audit.euler);
            AK_CHECK(op::mesh_deviation(t.mesh, r.solid) < 2e-3);
        }
    }
    std::printf("%s: volume %.12f (closed form %.12f), %d faces%s%s\n", name, v, volume,
                static_cast<int>(r.solid.faces.size()), r.empty ? ", empty" : "", r.manifold ? "" : ", not manifold");
}

} // namespace

int main() {
    // The quarter shell 1 <= r <= 2, 0 <= angle <= 90 degrees, 0 <= z <= 1, on
    // the exact circle parameterisation. Its faces: 0 inner cylinder, 1 outer
    // cylinder, 2 and 3 the radial planes y = 0 and x = 0, 4 and 5 the caps.
    const op::PatchSolid shell_ps = op::make_cylinder_sector(1.0, 2.0, 0.0, 1.0, 0.0, 1.0);
    const PolyPatchD& inner = shell_ps.faces[0];
    const PolyPatchD& outer = shell_ps.faces[1];

    // --- the exact contact tests ------------------------------------------------------
    {
        const op::PatchSolid top = op::make_box({-1.0, 0.0, -1.0}, {1.0, 1.0, 2.0}); // its face y = 1
        const op::PatchSolid right = op::make_box({0.0, -1.0, -1.0}, {1.0, 1.0, 2.0}); // its face x = 1
        const PolyPatchD* plane_y1 = nullptr;
        const PolyPatchD* plane_x1 = nullptr;
        for (const PolyPatchD& f : top.faces)
            if (f.c[0].y == 1.0 && f.c.back().y == 1.0 && f.c[1].y == 1.0) plane_y1 = &f;
        for (const PolyPatchD& f : right.faces)
            if (f.c[0].x == 1.0 && f.c.back().x == 1.0 && f.c[1].x == 1.0) plane_x1 = &f;
        AK_CHECK(plane_y1 && plane_x1);
        if (plane_y1 && plane_x1) {
            // The inner cylinder touches y = 1 along its side at 90 degrees and
            // x = 1 along its side at 0 degrees: order 2, and nowhere else.
            int side_y = -1, side_x = -1;
            for (int s = 0; s < 4; ++s) {
                if (op::side_in_plane(inner, s, *plane_y1)) side_y = s;
                if (op::side_in_plane(inner, s, *plane_x1)) side_x = s;
            }
            AK_CHECK(side_y >= 0 && side_x >= 0 && side_x != side_y);
            const op::SideContact cy = op::side_contact(inner, side_y, *plane_y1);
            const op::SideContact cx = op::side_contact(inner, side_x, *plane_x1);
            AK_CHECK(cy.only_side && cy.order == 2);
            AK_CHECK(cx.only_side && cx.order == 2);
            // On the enlarged square too (the circle stays on one side of its tangent).
            AK_CHECK(op::side_contact(inner, side_y, *plane_y1, 12, 1.0 / 64.0).only_side);
            // The radial plane x = 0 meets y = 1 transversally along its side: order 1.
            const PolyPatchD& radial = shell_ps.faces[3];
            int side_r = -1;
            for (int s = 0; s < 4; ++s)
                if (op::side_in_plane(radial, s, *plane_y1)) side_r = s;
            AK_CHECK(side_r >= 0);
            if (side_r >= 0) {
                const op::SideContact cr = op::side_contact(radial, side_r, *plane_y1);
                AK_CHECK(cr.only_side && cr.order == 1);
            }
            // The circle at the bottom of the inner cylinder touches x = 1 at
            // one end, to second order, and meets it nowhere else.
            for (int s = 0; s < 4; ++s) {
                const auto cp = op::side_control(inner, s);
                if (cp.size() != 3 || cp[0].z != 0.0 || cp[2].z != 0.0) continue;
                const op::EndContact ec = op::end_contact(solve::make_curve(cp), *plane_x1);
                AK_CHECK(ec.decided);
                AK_CHECK((ec.order0 == 2 && ec.order1 == 0) || (ec.order0 == 0 && ec.order1 == 2));
            }
        }
        // A plane through the cylinder's interior is not a contact along a side.
        const op::PatchSolid mid = op::make_box({-1.0, 0.0, -1.0}, {1.0, 0.5, 2.0});
        for (const PolyPatchD& f : mid.faces) {
            if (!(f.c[0].y == 0.5 && f.c.back().y == 0.5 && f.c[1].y == 0.5)) continue;
            for (int s = 0; s < 4; ++s) AK_CHECK(!op::side_contact(inner, s, f).only_side);
        }
        // A hill touches the plane of its base along all four sides: no single
        // side is the whole contact.
        const op::PatchSolid hill = op::make_bump_box({0, 0, 0}, 1, 1, 1, 0.3);
        const op::PatchSolid lid = op::make_box({-1.0, -1.0, 0.0}, {2.0, 2.0, 1.0});
        for (const PolyPatchD& h : hill.faces) {
            bool raised = false; // the hill top: control points above z = 1
            for (const geom::Vec4d& c : h.c) raised = raised || c.z > 1.0 * c.w;
            if (!raised) continue;
            for (const PolyPatchD& f : lid.faces) {
                if (!(f.c[0].z == 1.0 && f.c.back().z == 1.0 && f.c[1].z == 1.0)) continue;
                int sides = 0;
                for (int s = 0; s < 4; ++s) {
                    sides += op::side_in_plane(h, s, f) ? 1 : 0;
                    AK_CHECK(!op::side_contact(h, s, f).only_side);
                }
                AK_CHECK(sides == 4);
            }
        }
        (void)outer;
    }

    // --- Boolean operations with a contact along a ruling -----------------------------
    const Brep shell = op::to_brep(shell_ps);
    const double shell_volume = 0.75 * kPi;
    {
        // A box whose face y = 1 touches the inner cylinder along its line at
        // 90 degrees, inside that face (x = 0 between -0.25 and 0.25). (Its
        // bottom is kept out of the plane of the shell's annular cap: flat faces
        // in one plane that are not affinely related are not supported.)
        const Brep box = op::to_brep(op::make_box({-0.25, 0.0, 0.25}, {0.25, 1.0, 0.75}));
        const double in = 0.5 * cap_area(0.25), box_volume = 0.25;
        run("shell and box touching y = 1", shell, box, Op::Intersection, in);
        run("shell or box touching y = 1", shell, box, Op::Union, shell_volume + box_volume - in);
        run("shell minus box touching y = 1", shell, box, Op::Difference, shell_volume - in);
    }
    {
        // The same with the face x = 1, at 0 degrees, and the box partly below
        // the shell.
        const Brep box = op::to_brep(op::make_box({0.0, -0.25, 0.25}, {1.0, 0.75, 0.75}));
        const double in = 0.5 * cap_area(0.75), box_volume = 0.5;
        run("shell and box touching x = 1", shell, box, Op::Intersection, in);
        run("shell or box touching x = 1", shell, box, Op::Union, shell_volume + box_volume - in);
        run("shell minus box touching x = 1", shell, box, Op::Difference, shell_volume - in);
    }
    {
        // A box outside the shell whose face x = 2 touches the outer cylinder
        // along its line at 0 degrees: the solids meet only along that line.
        const Brep box = op::to_brep(op::make_box({2.0, -0.25, 0.25}, {2.5, 0.25, 0.75}));
        run("shell and box touching outside", shell, box, Op::Intersection, 0.0);
        run("shell or box touching outside", shell, box, Op::Union, shell_volume + 0.125, false);
        run("shell minus box touching outside", shell, box, Op::Difference, shell_volume);
    }

    return akt::finish("contact_test");
}
