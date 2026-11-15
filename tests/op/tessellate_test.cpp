// Watertight tessellation of Breps: input solids and Boolean results.
//
// What must hold exactly: the mesh is closed (every edge used once in each
// direction), has no degenerate triangle, and its Euler characteristic and
// number of components match the Brep's integer audit. What is measured:
// the deviation from the surfaces and the volume, both of which must shrink
// with the chord tolerance.
#include <akernel/op/boolean.hpp>
#include <akernel/op/tessellate.hpp>

#include <cmath>
#include <cstdio>

#include "check.hpp"

using namespace ak;
using geom::Vec3d;
using op::Brep;
using op::Op;

namespace {

struct MeshCheck {
    bool ok = false;
    op::MeshAudit audit;
    double deviation = 0.0, volume_error = 0.0;
    op::TessResult tess;
};

MeshCheck check_mesh(const char* name, const Brep& b, double chord, bool quiet = false) {
    MeshCheck mc;
    op::TessOptions opt;
    opt.chord = chord;
    mc.tess = op::tessellate(b, opt);
    AK_CHECK(mc.tess.ok);
    if (!mc.tess.ok) {
        std::printf("%s: tessellation failed: %s\n", name, mc.tess.error.c_str());
        return mc;
    }
    const op::Mesh& m = mc.tess.mesh;
    const op::MeshAudit a = op::audit(m);
    const topo::AuditReport br = op::audit(b);
    const int euler = br.euler;
    AK_CHECK(a.watertight());
    AK_CHECK(a.degenerate == 0);
    AK_CHECK(a.unused_points == 0);
    AK_CHECK(a.components == static_cast<int>(br.shells.size()));
    AK_CHECK(a.euler == euler);
    // Every edge of the Brep appears in the mesh as its sampled polyline.
    AK_CHECK(m.edges.size() == b.edges.size());
    mc.deviation = op::mesh_deviation(m, b);
    mc.volume_error = std::fabs(a.volume - op::volume(b));
    AK_CHECK(mc.deviation <= 2.0 * chord);
    AK_CHECK(mc.volume_error <= 2.0 * chord * a.area);
    if (!quiet)
        std::printf("%s (chord %.0e): %d points, %d triangles, %d shells, euler %d, watertight %s, "
                    "deviation %.2e, volume error %.2e, %.1f ms\n",
                    name, chord, a.points, a.triangles, a.components, a.euler, a.watertight() ? "yes" : "no",
                    mc.deviation, mc.volume_error, mc.tess.stats.ms);
    mc.audit = a;
    mc.ok = a.watertight();
    return mc;
}

} // namespace

int main() {
    const Brep unit = op::to_brep(op::make_box({0, 0, 0}, {1, 1, 1}));

    // --- a box: two triangles per face, nothing else ------------------------------
    {
        const MeshCheck c = check_mesh("unit box", unit, 1e-3);
        AK_CHECK(c.audit.points == 8);
        AK_CHECK(c.audit.triangles == 12);
        AK_CHECK_NEAR(c.audit.volume, 1.0, 1e-15);
        AK_CHECK_NEAR(c.audit.area, 6.0, 1e-15);
    }

    // --- curved solids, and convergence with the chord tolerance ------------------
    const Brep hill = op::to_brep(op::make_bump_box({0, 0, 0}, 1, 1, 1, 0.3));
    const Brep shell = op::to_brep(op::make_annular_sector(1.0, 2.0, 0.0, 1.0));
    for (const auto& [name, b] : {std::pair<const char*, const Brep*>{"hill", &hill}, {"cylinder shell", &shell}}) {
        double prev_dev = 1e9, prev_vol = 1e9;
        int prev_tris = 0;
        for (double chord : {1e-2, 1e-3, 1e-4}) {
            const MeshCheck c = check_mesh(name, *b, chord);
            AK_CHECK(c.deviation < prev_dev);
            AK_CHECK(c.volume_error < prev_vol || c.volume_error < 1e-9);
            AK_CHECK(c.audit.triangles > prev_tris);
            prev_dev = c.deviation;
            prev_vol = c.volume_error;
            prev_tris = c.audit.triangles;
        }
    }

    // --- Boolean results ----------------------------------------------------------
    auto check_op = [](const char* name, const Brep& a, const Brep& b, Op o, double chord = 1e-3) {
        const op::BooleanResult r = op::boolean_op(a, b, o);
        AK_CHECK(r.ok);
        if (!r.ok) { std::printf("%s: Boolean failed: %s\n", name, r.error.c_str()); return Brep{}; }
        if (!r.empty) check_mesh(name, r.solid, chord);
        return r.solid;
    };
    const Brep boxB = op::to_brep(op::make_box({0.5, 0.4, 0.3}, {1.5, 1.4, 1.3}));
    check_op("offset boxes, union", unit, boxB, Op::Union);
    check_op("offset boxes, difference", unit, boxB, Op::Difference);
    check_op("box with a cavity", unit, op::to_brep(op::make_box({0.3, 0.3, 0.3}, {0.6, 0.6, 0.6})), Op::Difference);

    const Brep slab = op::to_brep(op::make_box({-0.5, -0.5, 1.15}, {1.5, 1.5, 2.0}));
    const Brep hill_cut = check_op("hill minus slab (a ring face)", hill, slab, Op::Difference);
    check_op("hill cap", hill, slab, Op::Intersection);
    check_op("cylinder shell and box, intersection", shell,
             op::to_brep(op::make_box({0.5, 0.6, 0.25}, {1.3, 1.4, 0.75})), Op::Intersection);
    const Brep hill2 = op::to_brep(op::make_bump_box({0.6, 0.6, -0.35}, 1, 1, 1, 0.3));
    check_op("hill and cylinder shell, union", hill2, shell, Op::Union);
    check_op("hill and cylinder shell, difference", hill2, shell, Op::Difference);
    if (!hill_cut.faces.empty())
        check_op("(hill - slab) - box across the loop", hill_cut,
                 op::to_brep(op::make_box({0.2, 0.2, 0.6}, {0.8, 0.7, 1.5})), Op::Difference);

    // Faces in a common plane: edges drawn on carriers they are not defined on.
    check_op("boxes face to face, union", unit, op::to_brep(op::make_box({1, 0.2, 0.3}, {2, 1.3, 0.9})), Op::Union);
    check_op("stacked boxes, union", unit, op::to_brep(op::make_box({0.3, 0.2, 1}, {0.8, 0.7, 1.5})), Op::Union);

    // The demo: a hill cut by a tilted box, then drilled through (genus 1).
    {
        const double cx = std::cos(0.5), sx = std::sin(0.5), cz = std::cos(0.4), sz = std::sin(0.4);
        const Vec3d ex{0.9 * cz, 0.9 * sz, 0.0};
        const Vec3d ey{-0.8 * sz * cx, 0.8 * cz * cx, 0.8 * sx};
        const Vec3d ez{0.7 * sz * sx, -0.7 * cz * sx, 0.7 * cx};
        const Brep B = op::to_brep(op::make_box({0.45, 0.2, 0.75}, ex, ey, ez));
        const Brep diff = check_op("hill minus tilted box", hill, B, Op::Difference);
        if (!diff.faces.empty()) {
            const op::BooleanResult drilled =
                op::boolean_op(diff, op::to_brep(op::make_box({0.3, 0.55, -0.5}, {0.62, 0.8, 2.0})), Op::Difference);
            AK_CHECK(drilled.ok);
            if (drilled.ok) {
                const MeshCheck c = check_mesh("drilled", drilled.solid, 1e-3);
                AK_CHECK(c.audit.euler == 0); // genus 1

                // Determinism: the same input gives the same mesh, bit for bit.
                const op::TessResult again = op::tessellate(drilled.solid, op::TessOptions{1e-3});
                AK_CHECK(again.ok && again.mesh.tris == c.tess.mesh.tris);
                bool same = again.ok && again.mesh.points.size() == c.tess.mesh.points.size();
                for (std::size_t i = 0; same && i < again.mesh.points.size(); ++i) {
                    const Vec3d &p = again.mesh.points[i], &q = c.tess.mesh.points[i];
                    same = p.x == q.x && p.y == q.y && p.z == q.z;
                }
                AK_CHECK(same);

                // A length limit on top of the chord.
                op::TessOptions fine;
                fine.chord = 1e-3;
                fine.max_edge = 0.05;
                const op::TessResult r = op::tessellate(drilled.solid, fine);
                AK_CHECK(r.ok);
                if (r.ok) {
                    const op::MeshAudit a = op::audit(r.mesh);
                    AK_CHECK(a.watertight() && a.euler == 0);
                    AK_CHECK(a.triangles > c.audit.triangles);
                    std::printf("drilled with max edge 0.05: %d triangles, watertight %s\n", a.triangles,
                                a.watertight() ? "yes" : "no");
                }
            }
        }
    }

    return akt::finish("tessellate_test");
}
