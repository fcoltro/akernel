// Body types (solid, sheet, wire, acorn, general) and attributes.
//
// The type of a body is decided on its incidence structure alone, by the
// integer audit; the wire edges' geometry is checked exactly. Attributes
// must survive copying, rigid maps (including mirror images), merging and
// Boolean operations.
#include <akernel/geom/nurbs_curve.hpp>
#include <akernel/op/boolean.hpp>
#include <akernel/op/brep_transform.hpp>
#include <akernel/topo/complex.hpp>

#include <string>
#include <vector>

#include "check.hpp"

using namespace ak;
using geom::Vec3d;
using geom::Vec4d;
using op::Attribute;
using op::Brep;
using topo::BodyKind;
using topo::Complex;

namespace {

// Unit cube faces, outward; vertex index bits are (x, y, z).
const std::vector<std::vector<int>> kCubeQuads = {
    {0, 2, 3, 1}, {4, 5, 7, 6}, {0, 1, 5, 4}, {2, 6, 7, 3}, {0, 4, 6, 2}, {1, 3, 7, 5},
};

Complex cube_faces(std::vector<int> which) {
    Complex c;
    for (int i = 0; i < 8; ++i) c.add_vertex();
    for (int f : which) c.add_polygon_face(f, {kCubeQuads[static_cast<std::size_t>(f)]});
    return c;
}

// An open tube: the four side faces of the cube (no top, no bottom).
Complex tube() { return cube_faces({2, 3, 4, 5}); }

void topology_classification() {
    // Solid.
    {
        const topo::AuditReport r = cube_faces({0, 1, 2, 3, 4, 5}).audit();
        AK_CHECK(r.kind() == BodyKind::Solid);
        AK_CHECK(!r.sheet());
        AK_CHECK(r.shells[0].boundary_loops == 0 && r.shells[0].genus == 0);
    }
    // Open box: a disc, one boundary cycle of four edges.
    {
        const topo::AuditReport r = cube_faces({0, 2, 3, 4, 5}).audit();
        AK_CHECK(r.kind() == BodyKind::Sheet);
        AK_CHECK(r.sheet());
        AK_CHECK(r.shells.size() == 1u);
        AK_CHECK(r.shells[0].boundary_edges == 4 && r.shells[0].boundary_loops == 1);
        AK_CHECK(r.shells[0].genus == 0);
    }
    // One quad (only its own four vertices: the cube's others would be acorns).
    {
        Complex c;
        for (int i = 0; i < 4; ++i) c.add_vertex();
        c.add_polygon_face(0, {{0, 1, 2, 3}});
        const topo::AuditReport r = c.audit();
        AK_CHECK(r.kind() == BodyKind::Sheet);
        AK_CHECK(r.shells[0].boundary_loops == 1 && r.shells[0].genus == 0);
    }
    // A tube: an annulus, two boundary cycles, genus 0.
    {
        const topo::AuditReport r = tube().audit();
        AK_CHECK(r.kind() == BodyKind::Sheet);
        AK_CHECK(r.shells[0].boundary_loops == 2 && r.shells[0].genus == 0);
    }
    // Two disjoint quads: a sheet with two shells.
    {
        const topo::AuditReport r = cube_faces({0, 1}).audit();
        AK_CHECK(r.kind() == BodyKind::Sheet);
        AK_CHECK(r.shells.size() == 2u);
    }
    // A torus with one face removed: genus 1, one boundary cycle.
    {
        Complex c;
        for (int i = 0; i < 9; ++i) c.add_vertex();
        auto v = [](int i, int j) { return ((i + 3) % 3) * 3 + ((j + 3) % 3); };
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                if (i || j) c.add_polygon_face(0, {{v(i, j), v(i + 1, j), v(i + 1, j + 1), v(i, j + 1)}});
        const topo::AuditReport r = c.audit();
        AK_CHECK(r.kind() == BodyKind::Sheet);
        AK_CHECK(r.shells[0].boundary_loops == 1 && r.shells[0].genus == 1);
    }
    // Two quads sharing only a corner: a pinch, so a general body.
    {
        Complex c;
        for (int i = 0; i < 7; ++i) c.add_vertex();
        c.add_polygon_face(0, {{0, 1, 2, 3}});
        c.add_polygon_face(1, {{0, 4, 5, 6}});
        const topo::AuditReport r = c.audit();
        AK_CHECK(r.nonmanifold_vertices == 1);
        AK_CHECK(r.kind() == BodyKind::General);
    }
    // Three quads on one edge: a general body (the edge's uses cannot balance).
    {
        Complex c;
        for (int i = 0; i < 8; ++i) c.add_vertex();
        c.add_polygon_face(0, {{0, 1, 2, 3}});
        c.add_polygon_face(1, {{1, 0, 4, 5}});
        c.add_polygon_face(2, {{1, 0, 6, 7}});
        const topo::AuditReport r = c.audit();
        AK_CHECK(r.nonmanifold_edges == 1 && r.same_way_edges == 0);
        AK_CHECK(r.kind() == BodyKind::General);
    }
    // A face turned over inside a sheet: an edge used twice the same way.
    {
        Complex d;
        for (int i = 0; i < 8; ++i) d.add_vertex();
        for (int f : {0, 2, 3, 4, 5}) {
            std::vector<int> q = kCubeQuads[static_cast<std::size_t>(f)];
            if (f == 2) q = {q[3], q[2], q[1], q[0]};
            d.add_polygon_face(f, {q});
        }
        const topo::AuditReport r = d.audit();
        AK_CHECK(r.misoriented_edges > 0 && r.same_way_edges > 0);
        AK_CHECK(r.kind() == BodyKind::Invalid);
    }
    // Wires: an open path, a closed cycle, a ring edge, a star.
    {
        Complex c;
        for (int i = 0; i < 4; ++i) c.add_vertex();
        c.add_edge(0, 1);
        c.add_edge(1, 2);
        c.add_edge(2, 3);
        topo::AuditReport r = c.audit();
        AK_CHECK(r.kind() == BodyKind::Wire);
        AK_CHECK(r.wires.size() == 1u && !r.wires[0].closed && r.wires[0].vertices == 4 && r.wires[0].edges == 3);
        c.add_edge(3, 0);
        r = c.audit();
        AK_CHECK(r.kind() == BodyKind::Wire && r.wires[0].closed);
        c.add_edge(0, 2); // a chord: vertices of degree 3
        r = c.audit();
        AK_CHECK(r.kind() == BodyKind::General && r.wires[0].max_degree == 3);
    }
    {
        Complex c;
        c.add_vertex();
        c.add_edge(0, 0);
        const topo::AuditReport r = c.audit();
        AK_CHECK(r.kind() == BodyKind::Wire && r.wires[0].closed && r.wires[0].max_degree == 2);
    }
    {
        Complex c;
        for (int i = 0; i < 4; ++i) c.add_vertex();
        c.add_edge(0, 1);
        c.add_edge(2, 3);
        const topo::AuditReport r = c.audit();
        AK_CHECK(r.kind() == BodyKind::Wire && r.wires.size() == 2u);
    }
    // Acorn, empty, and a vertex next to a wire.
    {
        Complex c;
        c.add_vertex();
        AK_CHECK(c.audit().kind() == BodyKind::Acorn);
        AK_CHECK(Complex{}.audit().kind() == BodyKind::Empty);
        c.add_vertex();
        c.add_vertex();
        c.add_edge(1, 2);
        AK_CHECK(c.audit().kind() == BodyKind::General);
    }
    // A solid with a wire edge hanging off it: general.
    {
        Complex c = cube_faces({0, 1, 2, 3, 4, 5});
        c.add_vertex();
        c.add_edge(0, 8);
        const topo::AuditReport r = c.audit();
        AK_CHECK(!r.closed());
        AK_CHECK(r.kind() == BodyKind::General);
    }
}

Brep unit_box() { return op::to_brep(op::make_box({0, 0, 0}, {1, 1, 1})); }

// The polyline 0 -> a -> b ... as degree-1 curves.
std::vector<geom::NurbsCurve> segments(const std::vector<Vec3d>& pts) {
    std::vector<geom::NurbsCurve> cs;
    for (std::size_t i = 0; i + 1 < pts.size(); ++i)
        cs.push_back(geom::bezier_curve({{pts[i].x, pts[i].y, pts[i].z, 1.0}, {pts[i + 1].x, pts[i + 1].y, pts[i + 1].z, 1.0}}));
    return cs;
}

void sheets() {
    const op::PatchSolid box = op::make_box({0, 0, 0}, {1, 1, 1});
    // The box without one face.
    std::vector<geom::PolyPatchD> open(box.faces.begin(), box.faces.end() - 1);
    Brep s;
    std::string error;
    AK_CHECK(op::make_sheet(open, s, error));
    const topo::AuditReport r = op::audit(s);
    AK_CHECK(op::body_kind(s) == BodyKind::Sheet);
    AK_CHECK(s.faces.size() == 5u && s.edges.size() == 12u && s.vertices.size() == 8u);
    AK_CHECK(r.shells[0].boundary_edges == 4 && r.shells[0].boundary_loops == 1 && r.shells[0].genus == 0);
    // Every boundary edge lies on one carrier only.
    int boundary = 0;
    for (const op::BEdge& e : s.edges) boundary += e.on[1].carrier < 0;
    AK_CHECK(boundary == 4);
    // Face geometry for point location works on sheet faces as on solid faces.
    op::FaceGeom g;
    AK_CHECK(op::face_geometry(s, 0, g, error));
    AK_CHECK(op::locate(g, num::Interval(0.5), num::Interval(0.5)) == op::Where::Inside);

    // One patch is a sheet.
    Brep one;
    AK_CHECK(op::make_sheet({box.faces[0]}, one, error));
    AK_CHECK(op::body_kind(one) == BodyKind::Sheet);
    // All six patches: a solid, through the same constructor.
    Brep all;
    AK_CHECK(op::make_sheet(box.faces, all, error));
    AK_CHECK(op::body_kind(all) == BodyKind::Solid);
    // A sheet is refused as a Boolean operand, with a reason.
    const op::BooleanResult br = op::boolean_op(s, unit_box(), op::Op::Union);
    AK_CHECK(!br.ok && !br.error.empty());

    // A mirror image of a sheet is still a sheet with the same edge count.
    op::SignedPermutation mirror;
    mirror.sign = {-1.0, 1.0, 1.0};
    const Brep m = op::transformed(s, mirror);
    AK_CHECK(op::body_kind(m) == BodyKind::Sheet);
}

void wires() {
    std::string error;
    // An L-shaped open polyline and a closed square.
    Brep L;
    AK_CHECK(op::make_wire(segments({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}}), L, error));
    AK_CHECK(op::body_kind(L) == BodyKind::Wire);
    AK_CHECK(L.vertices.size() == 3u && L.edges.size() == 2u);
    AK_CHECK(op::check_wire_geometry(L, error));

    Brep sq;
    AK_CHECK(op::make_wire(segments({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0, 0, 0}}), sq, error));
    const topo::AuditReport r = op::audit(sq);
    AK_CHECK(r.kind() == BodyKind::Wire && r.wires.size() == 1u && r.wires[0].closed);
    AK_CHECK(sq.vertices.size() == 4u);

    // The exact full circle ends at 16 * (r, 0, 0, 1): the same point as its
    // start in homogeneous coordinates, so it is one ring edge.
    Brep ring;
    AK_CHECK(op::make_wire({geom::exact_circle(2.0, 0.5, 4)}, ring, error));
    AK_CHECK(ring.vertices.size() == 1u && ring.edges.size() == 1u);
    AK_CHECK(ring.edges[0].v0 == ring.edges[0].v1);
    AK_CHECK(op::body_kind(ring) == BodyKind::Wire);
    AK_CHECK(op::check_wire_geometry(ring, error));

    // Three quarters of a circle and the two radii back to the centre: a
    // closed wire of a curved edge and two straight ones.
    {
        std::vector<geom::NurbsCurve> cs = {geom::exact_circle(1.0, 0.0, 3)};
        const Vec4d end = cs[0].cp.back();
        cs.push_back(geom::bezier_curve({end, {0, 0, 0, end.w}}));
        cs.push_back(geom::bezier_curve({{0, 0, 0, 1}, {1, 0, 0, 1}}));
        Brep sector;
        AK_CHECK(op::make_wire(cs, sector, error));
        AK_CHECK(op::body_kind(sector) == BodyKind::Wire);
        AK_CHECK(op::audit(sector).wires[0].closed && sector.vertices.size() == 3u);
        AK_CHECK(op::check_wire_geometry(sector, error));
    }

    // Moving a wire keeps it exact: curves and vertices map the same way.
    const op::SignedPermutation turn = op::SignedPermutation::quarter_turn_z({0.25, -3, 8});
    const Brep moved = op::transformed(ring, turn);
    AK_CHECK(op::check_wire_geometry(moved, error));
    AK_CHECK(op::body_kind(moved) == BodyKind::Wire);

    // The geometry check catches a vertex that is not the curve's end.
    Brep bad = L;
    bad.vertices[1].xp = {1.0, 0.0, 0.0, 2.0};
    AK_CHECK(!op::check_wire_geometry(bad, error));
    // And a knot vector that is not clamped is refused at construction.
    geom::NurbsCurve open_knots(1, {0.0, 0.5, 1.0, 1.0}, {{0, 0, 0, 1}, {1, 0, 0, 1}});
    Brep none;
    AK_CHECK(!op::make_wire({open_knots}, none, error));
    // A wire is refused as a Boolean operand, not taken for the empty solid.
    for (op::Op o : {op::Op::Union, op::Op::Intersection, op::Op::Difference}) {
        AK_CHECK(!op::boolean_op(L, unit_box(), o).ok);
        AK_CHECK(!op::boolean_op(unit_box(), L, o).ok);
    }
    AK_CHECK(op::boolean_op(Brep{}, unit_box(), op::Op::Union).ok); // the empty body still is
}

Attribute colour(std::int64_t rgb) {
    Attribute a;
    a.name = "colour";
    a.ints = {rgb};
    return a;
}

void attributes() {
    Brep b = unit_box();
    Attribute datum;
    datum.name = "datum";
    datum.text = "top";
    datum.reals = {0.125};
    datum.points = {{0.5, 0.5, 1.0}};
    datum.directions = {{0.0, 0.0, 1.0}};
    op::set_attribute(b.faces[0].attributes, datum);
    op::set_attribute(b.faces[0].attributes, colour(0xff0000));
    op::set_attribute(b.faces[0].attributes, colour(0x00ff00)); // replaces
    op::set_attribute(b.edges[3].attributes, colour(0x0000ff));
    op::set_attribute(b.vertices[2].attributes, colour(7));
    op::set_attribute(b.attributes, colour(1));
    AK_CHECK(b.faces[0].attributes.size() == 2u);
    AK_CHECK(op::find_attribute(b.faces[0].attributes, "colour")->ints[0] == 0x00ff00);

    // Copy.
    const Brep c = b;
    AK_CHECK(op::same_attribute(*op::find_attribute(c.faces[0].attributes, "datum"), datum));
    AK_CHECK(op::find_attribute(c.edges[3].attributes, "colour")->ints[0] == 0x0000ff);
    AK_CHECK(op::find_attribute(c.vertices[2].attributes, "colour")->ints[0] == 7);

    // A quarter turn about z and a shift: the point moves, the direction turns.
    {
        const op::SignedPermutation m = op::SignedPermutation::quarter_turn_z({10, 0, 0});
        const Brep t = op::transformed(b, m);
        const Attribute* a = op::find_attribute(t.faces[0].attributes, "datum");
        AK_CHECK(a != nullptr);
        if (a) {
            AK_CHECK(a->text == "top" && a->reals == std::vector<double>{0.125});
            AK_CHECK(op::same_vec(a->points[0], {9.5, 0.5, 1.0}));
            AK_CHECK(op::same_vec(a->directions[0], {0.0, 0.0, 1.0}));
        }
        AK_CHECK(op::find_attribute(t.attributes, "colour")->ints[0] == 1);
        AK_CHECK(op::find_attribute(t.edges[3].attributes, "colour") != nullptr);
    }
    // A mirror in z: the direction is mirrored too, and the faces turn over.
    {
        op::SignedPermutation m;
        m.sign = {1.0, 1.0, -1.0};
        const Brep t = op::transformed(b, m);
        const Attribute* a = op::find_attribute(t.faces[0].attributes, "datum");
        AK_CHECK(a && op::same_vec(a->points[0], {0.5, 0.5, -1.0}) && op::same_vec(a->directions[0], {0.0, 0.0, -1.0}));
        AK_CHECK(op::body_kind(t) == BodyKind::Solid);
    }
    // Remove.
    {
        Brep d = b;
        AK_CHECK(op::remove_attribute(d.faces[0].attributes, "colour"));
        AK_CHECK(!op::remove_attribute(d.faces[0].attributes, "colour"));
        AK_CHECK(d.faces[0].attributes.size() == 1u);
    }

    // Merging keeps every entity's attributes and re-indexes the geometry.
    {
        Brep w;
        std::string error;
        AK_CHECK(op::make_wire(segments({{3, 0, 0}, {4, 0, 0}}), w, error));
        op::set_attribute(w.edges[0].attributes, colour(42));
        op::set_attribute(w.attributes, colour(2));
        Attribute tag;
        tag.name = "tag";
        op::set_attribute(w.attributes, tag);
        const Brep m = op::merged(b, w);
        AK_CHECK(op::body_kind(m) == BodyKind::General);
        AK_CHECK(m.edges.back().curve == 0 && m.edges.back().v0 == static_cast<int>(b.vertices.size()));
        AK_CHECK(op::find_attribute(m.edges.back().attributes, "colour")->ints[0] == 42);
        AK_CHECK(op::find_attribute(m.attributes, "colour")->ints[0] == 1); // A's wins
        AK_CHECK(op::find_attribute(m.attributes, "tag") != nullptr);
        AK_CHECK(op::check_wire_geometry(m, error));
        // Two boxes apart: one solid body with two shells.
        const Brep two = op::merged(unit_box(), op::transformed(unit_box(), op::SignedPermutation::translation({3, 0, 0})));
        const topo::AuditReport r = op::audit(two);
        AK_CHECK(r.kind() == BodyKind::Solid && r.shells.size() == 2u);
    }

    // Booleans: regions of faces keep their face's attributes, pieces of
    // edges keep their edge's.
    {
        Brep A = unit_box(), B = op::transformed(unit_box(), op::SignedPermutation::translation({0.5, 0.25, 0.125}));
        for (op::BFace& f : A.faces) op::set_attribute(f.attributes, colour(0xA));
        for (op::BFace& f : B.faces) op::set_attribute(f.attributes, colour(0xB));
        for (op::BEdge& e : A.edges) op::set_attribute(e.attributes, colour(0xAE));
        for (op::Op o : {op::Op::Union, op::Op::Intersection, op::Op::Difference}) {
            const op::BooleanResult r = op::boolean_op(A, B, o);
            AK_CHECK(r.ok);
            if (!r.ok) continue;
            AK_CHECK(op::body_kind(r.solid) == BodyKind::Solid);
            int a = 0, bb = 0, untagged = 0;
            for (const op::BFace& f : r.solid.faces) {
                const Attribute* at = op::find_attribute(f.attributes, "colour");
                if (!at) ++untagged;
                else if (at->ints[0] == 0xA) ++a;
                else if (at->ints[0] == 0xB) ++bb;
            }
            AK_CHECK(untagged == 0);
            AK_CHECK(a > 0 && bb > 0);
            int tagged_edges = 0;
            for (const op::BEdge& e : r.solid.edges) tagged_edges += op::find_attribute(e.attributes, "colour") != nullptr;
            AK_CHECK(tagged_edges > 0);
        }
    }
}

} // namespace

int main() {
    topology_classification();
    sheets();
    wires();
    attributes();
    return akt::finish("body_test");
}
