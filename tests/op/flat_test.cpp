// Flat faces in one plane on carriers that are not affinely related: the
// annular cap of a cylinder sector (rational in the angle) against a box face
// or against edges of the other solid that lie in its plane. The imprint
// moves such faces onto exact parallelograms of their plane (op/boolean.hpp,
// flat_detail). Two caps whose parameterisations are swapped are related by
// an affine map that exchanges the parameters and are not moved. Also the
// contacts that moving a cap leaves behind: a vertex of the other solid in
// the cap's plane at the end of an edge, and a circle that is a side of a
// cylinder patch in one solid and an intersection curve in the other. Every
// volume is checked against a closed form.
#include <akernel/op/boolean.hpp>
#include <akernel/op/surfaces.hpp>
#include <akernel/op/tessellate.hpp>

#include <cmath>
#include <cstdio>

#include "check.hpp"

using namespace ak;
using op::Brep;
using op::Op;

namespace {

const double kPi = 3.14159265358979323846;

double sector_volume(double r0, double r1, double t0, double t1, double z0, double z1) {
    return 0.5 * (2.0 * std::atan(t1) - 2.0 * std::atan(t0)) * (r1 * r1 - r0 * r0) * (z1 - z0);
}

/// Area between the unit circle and its tangent at the top, over x in [0, a].
double cap_area(double a) { return a - 0.5 * (a * std::sqrt(1.0 - a * a) + std::asin(a)); }

enum Recarry { kMoved, kNotMoved, kEither };

/// One operation with a closed-form volume. The result must pass the audit,
/// tessellate into a closed mesh and work as an operand. Returns the result
/// (empty if the operation was not solved).
Brep run(const char* name, const Brep& a, const Brep& b, Op o, double volume, Recarry expect = kMoved) {
    std::printf("%s ...\n", name);
    const op::Imprint im = op::imprint(a, b);
    AK_CHECK(im.ok);
    if (!im.ok) {
        std::printf("FAIL %s: %s\n", name, im.error.c_str());
        op::BooleanOptions traced;
        traced.trace = true;
        for (const std::string& line : op::imprint(a, b, traced).trace) std::printf("  %s\n", line.c_str());
        return {};
    }
    if (expect != kEither) AK_CHECK((im.stats.recarried_faces > 0) == (expect == kMoved));
    const op::BooleanResult r = op::select(im, o);
    AK_CHECK(r.ok);
    if (!r.ok) {
        std::printf("FAIL %s: %s\n", name, r.error.c_str());
        return {};
    }
    const double v = r.empty ? 0.0 : op::volume(r.solid);
    AK_CHECK_NEAR(v, volume, 1e-11);
    if (!r.empty) {
        AK_CHECK(r.audit.closed());
        const op::TessResult t = op::tessellate(r.solid);
        AK_CHECK(t.ok);
        if (t.ok) {
            const op::MeshAudit m = op::audit(t.mesh);
            AK_CHECK(m.watertight() && m.euler == r.audit.euler);
            AK_CHECK(op::mesh_deviation(t.mesh, r.solid) < 2e-3);
        }
        // The result is an operand: combined with a box far away it comes back whole.
        const Brep far = op::to_brep(op::make_box({5.0, 5.0, 5.0}, {6.0, 6.0, 6.0}));
        const op::BooleanResult again = op::boolean_op(r.solid, far, Op::Union);
        AK_CHECK(again.ok);
        if (again.ok) AK_CHECK_NEAR(op::volume(again.solid), volume + 1.0, 1e-11);
    }
    std::printf("%s: volume %.12f (closed form %.12f), %d faces, %d carriers moved%s\n", name, v, volume,
                static_cast<int>(r.solid.faces.size()), im.stats.recarried_faces, r.empty ? ", empty" : "");
    return r.solid;
}

/// Volume of a box with a hill on top (op::make_bump_box): the hill is
/// 16 b x (1 - x) y (1 - y) over the unit square of the top.
double bump_box_volume(double sx, double sy, double sz, double bump) { return sx * sy * sz + 4.0 / 9.0 * bump * sx * sy; }

} // namespace

int main() {
    std::setvbuf(stdout, nullptr, _IONBF, 0); // progress is visible even when redirected
    const Brep shell = op::to_brep(op::make_cylinder_sector(1.0, 2.0, 0.0, 1.0, 0.0, 1.0));
    const double shell_volume = 0.75 * kPi;
    {
        // A box standing on the plane of the shell's bottom cap, inside the shell.
        const Brep box = op::to_brep(op::make_box({1.25, 0.0, 0.0}, {1.75, 0.5, 0.5}));
        run("shell and box on its cap", shell, box, Op::Intersection, 0.125);
        run("shell or box on its cap", shell, box, Op::Union, shell_volume);
        run("shell minus box on its cap", shell, box, Op::Difference, shell_volume - 0.125);
    }
    {
        // A box on the cap's plane, partly outside the shell, whose face y = 1
        // also touches the inner cylinder along its line (both contacts at once).
        const Brep box = op::to_brep(op::make_box({-0.25, 0.0, 0.0}, {0.25, 1.0, 0.5}));
        const double in = 0.5 * cap_area(0.25);
        run("shell and touching box on its cap", shell, box, Op::Intersection, in);
        run("shell or touching box on its cap", shell, box, Op::Union, shell_volume + 0.25 - in);
        run("shell minus touching box on its cap", shell, box, Op::Difference, shell_volume - in);
    }
    {
        // A sector below the shell: its top cap and the shell's bottom cap lie in
        // z = 0 with swapped parameters (radius and angle exchanged). An affine
        // map that exchanges the parameters relates them exactly, so the
        // overlay applies and nothing is moved.
        const Brep below = op::to_brep(op::make_cylinder_sector(1.5, 2.0, 0.0, 0.25, -0.5, 0.0));
        const double sv = sector_volume(1.5, 2.0, 0.0, 0.25, -0.5, 0.0);
        run("shell and sector below", shell, below, Op::Intersection, 0.0, kNotMoved);
        run("shell or sector below", shell, below, Op::Union, shell_volume + sv, kNotMoved);
        run("shell minus sector below", shell, below, Op::Difference, shell_volume, kNotMoved);
        // The same sector together with a box beside it, between the same two
        // planes: the sector's caps share their planes with the box's faces
        // and are moved, although the two solids do not meet. The union (two
        // shells, caps on parallelograms) is then put under the shell: z = 0
        // holds the shell's cap, which is moved as well, the sector's moved
        // cap and the box's face.
        const Brep box = op::to_brep(op::make_box({0.0, 1.25, -0.5}, {0.5, 1.75, 0.0}));
        const Brep both = run("sector below or box below", below, box, Op::Union, sv + 0.125);
        if (!both.faces.empty()) {
            run("shell and (sector, box) below", shell, both, Op::Intersection, 0.0);
            run("shell or (sector, box) below", shell, both, Op::Union, shell_volume + sv + 0.125);
            run("shell minus (sector, box) below", shell, both, Op::Difference, shell_volume);
        }
    }
    {
        // A cap that does not conflict with anything is left where it is.
        const Brep box = op::to_brep(op::make_box({1.25, 0.25, 0.25}, {1.75, 0.5, 0.5}));
        run("shell and box inside, no shared plane", shell, box, Op::Intersection, 0.03125, kNotMoved);
    }
    {
        // A box with a hill on top whose base (z = 1) is the plane of the
        // shell's top cap: no face of the box lies in that plane, but the four
        // edges around the hill do, and they cross the cap. The cap is moved;
        // the hill touches the plane along all four of its sides (exact), and
        // the box's face y = 1 touches the inner cylinder along its line.
        const Brep hill = op::to_brep(op::make_bump_box({-0.25, 0.0, 0.25}, 0.75, 1.0, 0.75, 0.3));
        const double hv = bump_box_volume(0.75, 1.0, 0.75, 0.3), in = 0.75 * cap_area(0.5);
        run("shell and hill box up to its cap", shell, hill, Op::Intersection, in);
        run("shell or hill box up to its cap", shell, hill, Op::Union, shell_volume + hv - in);
        run("shell minus hill box up to its cap", shell, hill, Op::Difference, shell_volume - in);
    }
    {
        // A sector in the shell's hole, against its inner cylinder: the two
        // faces on that cylinder are parameterised the other way round (an
        // outer and an inner face), which an affine map exchanging the
        // parameters relates. The difference is the shell again, with the
        // outline of the contact left on its inner face.
        const Brep plug = op::to_brep(op::make_cylinder_sector(0.75, 1.0, 0.25, 0.5, 0.5, 0.75));
        const double pv = sector_volume(0.75, 1.0, 0.25, 0.5, 0.5, 0.75);
        run("shell and plug in its hole", shell, plug, Op::Intersection, 0.0, kEither);
        run("shell or plug in its hole", shell, plug, Op::Union, shell_volume + pv, kEither);
        const Brep marked = run("shell minus plug in its hole", shell, plug, Op::Difference, shell_volume, kEither);
        // A slot through the marked shell whose top (z = 0.5) holds the lower
        // edge of that outline and whose inner face runs along it. The slot's
        // cap is moved; a vertex of the outline then lies in the cap's plane at
        // the end of an edge (decided exactly), and the outline's edge, a side
        // of a cylinder patch, is one curve with the cap's circle, which is
        // now the intersection curve of the cylinder and the plane.
        if (!marked.faces.empty()) {
            const Brep slot = op::to_brep(op::make_cylinder_sector(1.0, 2.0, 0.25, 0.75, 0.25, 0.5));
            const double sv = sector_volume(1.0, 2.0, 0.25, 0.75, 0.25, 0.5);
            run("marked shell and slot", marked, slot, Op::Intersection, sv, kEither);
            run("marked shell or slot", marked, slot, Op::Union, shell_volume, kEither);
            run("marked shell minus slot", marked, slot, Op::Difference, shell_volume - sv, kEither);
        }
    }
    {
        // A pocket with a hill-shaped roof whose base edges lie in z = 0.75,
        // then a sector pocket up to that plane, then a sector put back into
        // it. The last sector's top cap lies in z = 0.75 and is crossed there
        // by the base edges of the roof, which meet at a corner inside the
        // cap; the cap's straight side runs diagonally past that corner, with
        // end points known only as enclosures. Every volume differs from the
        // previous one by a closed form.
        const Brep hill = op::to_brep(op::make_bump_box({0.5, 0.75, 0.0}, 1.0, 0.75, 0.75, 0.3));
        const op::BooleanResult s1 = op::boolean_op(shell, hill, Op::Difference);
        AK_CHECK(s1.ok && !s1.empty);
        if (s1.ok && !s1.empty) {
            const Brep pocket = op::to_brep(op::make_cylinder_sector(1.5, 2.0, 0.0, 0.5, -0.25, 0.75));
            const op::BooleanResult s2 = op::boolean_op(s1.solid, pocket, Op::Difference);
            AK_CHECK(s2.ok && !s2.empty);
            if (s2.ok && !s2.empty) {
                const double v2 = op::volume(s2.solid);
                const Brep back = op::to_brep(op::make_cylinder_sector(1.5, 2.0, 0.0, 0.25, 0.25, 0.75));
                const double bv = sector_volume(1.5, 2.0, 0.0, 0.25, 0.25, 0.75);
                run("pocketed shell and sector put back", s2.solid, back, Op::Intersection, 0.0, kEither);
                run("pocketed shell or sector put back", s2.solid, back, Op::Union, v2 + bv, kEither);
                run("pocketed shell minus sector put back", s2.solid, back, Op::Difference, v2, kEither);
            }
        }
    }
    {
        // A box with a hill on top, face to face with a plain box. The hill
        // box's flat faces have degree 2, so the four edges around the common
        // square are degree-2 sides in one solid and degree-1 sides in the
        // other. They are straight all the same (control points exactly on one
        // line) and are one edge each, fixed by their end vertices.
        const Brep hill = op::to_brep(op::make_bump_box({0.0, 0.0, 0.0}, 1.0, 1.0, 1.0, 0.3));
        const Brep box = op::to_brep(op::make_box({1.0, 0.0, 0.0}, {2.0, 1.0, 1.0}));
        const double hv = bump_box_volume(1.0, 1.0, 1.0, 0.3);
        run("hill box and box face to face", hill, box, Op::Intersection, 0.0, kNotMoved);
        run("hill box or box face to face", hill, box, Op::Union, hv + 1.0, kNotMoved);
        run("hill box minus box face to face", hill, box, Op::Difference, hv, kNotMoved);
    }
    return akt::finish("flat_test");
}
