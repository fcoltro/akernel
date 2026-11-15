// Non-manifold results: solids whose parts touch along an edge or at a point.
//
// Such a result is accepted when its boundary is closed and oriented and the
// faces around every edge with more than two faces alternate between material
// and empty space (op::check_edge_fans). It must then work like any other
// solid: volume, tessellation, comparison, and as an operand of further
// operations. Every volume here has a closed form.
#include <akernel/op/boolean.hpp>
#include <akernel/op/compare.hpp>
#include <akernel/op/tessellate.hpp>

#include <cstdio>
#include <string>

#include "check.hpp"

using namespace ak;
using op::Brep;
using op::Op;

namespace {

Brep box(double x0, double y0, double z0, double x1, double y1, double z1) {
    return op::to_brep(op::make_box({x0, y0, z0}, {x1, y1, z1}));
}

/// Run one operation and check the result against a closed-form volume and
/// the expected manifoldness. Returns the result (empty Brep on failure).
Brep run(const char* name, const Brep& a, const Brep& b, Op o, double volume, bool manifold, int shells = -1) {
    const op::BooleanResult r = op::boolean_op(a, b, o);
    AK_CHECK(r.ok);
    if (!r.ok) {
        std::printf("FAIL %s: %s\n", name, r.error.c_str());
        return Brep{};
    }
    AK_CHECK(!r.empty);
    AK_CHECK(r.audit.closed());
    AK_CHECK(r.manifold == manifold);
    const double v = op::volume(r.solid);
    AK_CHECK_NEAR(v, volume, 1e-12);
    if (shells >= 0) AK_CHECK(static_cast<int>(r.audit.shells.size()) == shells);
    // The mesh must be closed and have the Brep's Euler characteristic.
    const op::TessResult t = op::tessellate(r.solid);
    AK_CHECK(t.ok);
    if (t.ok) {
        const op::MeshAudit m = op::audit(t.mesh);
        AK_CHECK(m.watertight() && m.degenerate == 0);
        AK_CHECK(m.euler == r.audit.euler);
        AK_CHECK(m.components == static_cast<int>(r.audit.shells.size()));
        AK_CHECK_NEAR(m.volume, volume, 1e-12);
        std::printf("%s: volume %.12f, %s, %d faces, %d shell(s), non-manifold edges %d, pinched vertices %d, "
                    "Euler characteristic %d (mesh %d, %s)\n",
                    name, v, r.manifold ? "manifold" : "not manifold", static_cast<int>(r.solid.faces.size()),
                    static_cast<int>(r.audit.shells.size()), r.audit.nonmanifold_edges, r.audit.nonmanifold_vertices,
                    r.audit.euler, m.euler, m.manifold() ? "manifold mesh" : "non-manifold mesh");
    }
    return r.solid;
}

} // namespace

int main() {
    const Brep A = box(0, 0, 0, 1, 1, 1);

    // --- touching along a whole edge ------------------------------------------------
    const Brep edge_pair = run("edge contact, union", A, box(1, 1, 0, 2, 2, 1), Op::Union, 2.0, false, 1);
    {
        const op::BooleanResult i = op::boolean_op(A, box(1, 1, 0, 2, 2, 1), Op::Intersection);
        AK_CHECK(i.ok && i.empty);
        const op::BooleanResult d = op::boolean_op(A, box(1, 1, 0, 2, 2, 1), Op::Difference);
        AK_CHECK(d.ok && d.manifold && op::compare_solids(d.solid, A).same);
    }
    // Along part of an edge: the shared piece is the middle of A's edge.
    run("edge contact over part of the edge, union", A, box(1, 1, 0.25, 2, 2, 0.75), Op::Union, 1.5, false, 1);
    // At one point.
    const Brep point_pair = run("point contact, union", A, box(1, 1, 1, 2, 2, 2), Op::Union, 2.0, false, 2);

    // --- non-manifold operands ----------------------------------------------------------
    if (!edge_pair.faces.empty()) {
        // A slot across the shared edge cuts both parts; they still touch above
        // and below it.
        run("(edge pair) - slot across the edge", edge_pair, box(0.5, 0.5, 0.4, 1.5, 1.5, 0.6), Op::Difference,
            2.0 - 0.1, false, 1);
        // Filling the slot region joins the parts there; the edge is still
        // shared above and below.
        run("(edge pair) or block across the edge", edge_pair, box(0.5, 0.5, 0.4, 1.5, 1.5, 0.6), Op::Union,
            2.0 + 0.2 - 0.1, false, 1);
        // The part of the block inside the pair: two pieces touching along a
        // short piece of the edge.
        run("(edge pair) and block across the edge", edge_pair, box(0.5, 0.5, 0.4, 1.5, 1.5, 0.6), Op::Intersection,
            0.1, false, 1);
        // Removing a column around the edge separates the parts: manifold again.
        run("(edge pair) - column around the edge", edge_pair, box(0.8, 0.8, -1, 1.2, 1.2, 2), Op::Difference,
            2.0 - 2 * 0.04, true, 2);
        // The non-manifold solid as the second operand: a slab with two square
        // holes whose walls meet along a piece of the shared edge, so the
        // result is not manifold either.
        run("slab - (edge pair)", box(-0.5, -0.5, 0.25, 2.5, 2.5, 0.75), edge_pair, Op::Difference,
            3.0 * 3.0 * 0.5 - 2 * 0.5, false, 1);
        // An identity with a non-manifold operand.
        const op::BooleanResult x = op::boolean_op(edge_pair, box(0.5, 0.5, 0.4, 1.5, 1.5, 0.6), Op::Intersection);
        const op::BooleanResult y = op::boolean_op(box(0.5, 0.5, 0.4, 1.5, 1.5, 0.6), edge_pair, Op::Intersection);
        AK_CHECK(x.ok && y.ok);
        if (x.ok && y.ok) {
            const op::SolidComparison c = op::compare_solids(x.solid, y.solid);
            AK_CHECK(c.same);
            if (!c.same) std::printf("commutativity with a non-manifold operand: %s\n", c.reason.c_str());
        }
    }
    if (!point_pair.faces.empty()) {
        // A box over the touching corner joins the two parts into one manifold solid.
        run("(point pair) or box over the corner", point_pair, box(0.75, 0.75, 0.75, 1.25, 1.25, 1.25), Op::Union,
            2.0 + 0.125 - 2 * 0.015625, true, 1);
        // A box that misses the corner leaves it pinched.
        run("(point pair) - box away from the corner", point_pair, box(0.2, 0.2, 0.2, 0.4, 0.4, 0.4), Op::Difference,
            2.0 - 0.008, false, 3);
    }

    // --- an invalid non-manifold boundary is refused --------------------------------------
    // Turning the second part inside out keeps every edge balanced in the
    // integer audit, but around the shared edge the faces no longer alternate
    // between material and empty space.
    if (!edge_pair.faces.empty()) {
        Brep bad = edge_pair;
        for (op::BFace& f : bad.faces) {
            const geom::BezierSurface S(bad.carriers[static_cast<std::size_t>(f.carrier)]);
            const geom::Vec3d c = S.point(0.5, 0.5);
            if (c.x > 1.0 + 1e-9 || c.y > 1.0 + 1e-9) f.flipped = !f.flipped;
        }
        bad.outward_known = false;
        AK_CHECK(op::audit(bad).closed());
        std::string why;
        AK_CHECK(!op::check_edge_fans(bad, why));
        std::printf("inside-out part refused: %s\n", why.c_str());
        std::string ok_why;
        AK_CHECK(op::check_edge_fans(edge_pair, ok_why));
    }

    return akt::finish("nonmanifold_test");
}
