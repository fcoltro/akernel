// Topology audit: valid solids pass, each kind of defect is named.
#include <algorithm>
#include <akernel/topo/complex.hpp>

#include "check.hpp"

using namespace ak::topo;

namespace {

// Unit cube; vertex index bits are (x, y, z). Faces are oriented outward.
Complex cube(bool skip_top = false, bool flip_bottom = false) {
    Complex c;
    for (int i = 0; i < 8; ++i) c.add_vertex();
    std::vector<std::vector<int>> quads = {
        {0, 2, 3, 1}, // z = 0, normal -z
        {4, 5, 7, 6}, // z = 1, normal +z
        {0, 1, 5, 4}, // y = 0
        {2, 6, 7, 3}, // y = 1
        {0, 4, 6, 2}, // x = 0
        {1, 3, 7, 5}, // x = 1
    };
    if (flip_bottom) std::reverse(quads[0].begin(), quads[0].end());
    for (std::size_t f = 0; f < quads.size(); ++f) {
        if (skip_top && f == 1) continue;
        c.add_polygon_face(static_cast<int>(f), {quads[f]});
    }
    return c;
}

// Torus as a 3 x 3 grid of quads with periodic identification.
Complex torus() {
    Complex c;
    for (int i = 0; i < 9; ++i) c.add_vertex();
    auto v = [](int i, int j) { return ((i + 3) % 3) * 3 + ((j + 3) % 3); };
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            c.add_polygon_face(0, {{v(i, j), v(i + 1, j), v(i + 1, j + 1), v(i, j + 1)}});
    return c;
}

// Square block with a square hole through it: top and bottom faces have a ring.
Complex block_with_hole() {
    Complex c;
    // outer bottom 0-3, outer top 4-7, inner bottom 8-11, inner top 12-15,
    // each counter-clockwise seen from above.
    for (int i = 0; i < 16; ++i) c.add_vertex();
    auto ob = [](int i) { return i % 4; };
    auto ot = [](int i) { return 4 + i % 4; };
    auto ib = [](int i) { return 8 + i % 4; };
    auto it = [](int i) { return 12 + i % 4; };
    c.add_polygon_face(0, {{ot(0), ot(1), ot(2), ot(3)}, {it(0), it(3), it(2), it(1)}}); // top
    c.add_polygon_face(1, {{ob(0), ob(3), ob(2), ob(1)}, {ib(0), ib(1), ib(2), ib(3)}}); // bottom
    for (int i = 0; i < 4; ++i) {
        c.add_polygon_face(2 + i, {{ob(i), ob(i + 1), ot(i + 1), ot(i)}});     // outer wall
        c.add_polygon_face(6 + i, {{ib(i + 1), ib(i), it(i), it(i + 1)}});     // hole wall
    }
    return c;
}

} // namespace

int main() {
    // Cube: one closed shell of genus 0.
    {
        const AuditReport r = cube().audit();
        AK_CHECK(r.watertight());
        AK_CHECK(r.shells.size() == 1u);
        AK_CHECK(r.shells[0].vertices == 8 && r.shells[0].edges == 12 && r.shells[0].faces == 6);
        AK_CHECK(r.shells[0].euler == 2 && r.shells[0].genus == 0);
        AK_CHECK(r.ledger_entries == 0);
    }

    // Missing face: four boundary edges, not watertight, but still a valid chain complex.
    {
        const AuditReport r = cube(true).audit();
        AK_CHECK(!r.watertight());
        AK_CHECK(r.chain_condition && r.loops_closed);
        AK_CHECK(r.boundary_edges == 4);
        AK_CHECK(!r.shells[0].closed);
    }

    // One face flipped: every one of its edges is used twice the same way.
    {
        const AuditReport r = cube(false, true).audit();
        AK_CHECK(!r.watertight());
        AK_CHECK(r.misoriented_edges == 4);
        AK_CHECK(r.boundary_edges == 0);
    }

    // Torus: genus 1.
    {
        const AuditReport r = torus().audit();
        AK_CHECK(r.watertight());
        AK_CHECK(r.shells[0].vertices == 9 && r.shells[0].edges == 18 && r.shells[0].faces == 9);
        AK_CHECK(r.shells[0].genus == 1);
    }

    // Block with a through hole: rings enter Euler-Poincare, genus 1.
    {
        const AuditReport r = block_with_hole().audit();
        AK_CHECK(r.watertight());
        AK_CHECK(r.shells.size() == 1u);
        AK_CHECK(r.shells[0].vertices == 16 && r.shells[0].edges == 24 && r.shells[0].faces == 10);
        AK_CHECK(r.shells[0].rings == 2);
        AK_CHECK(r.shells[0].genus == 1);
    }

    // Two disjoint cubes: two shells.
    {
        Complex c = cube();
        const Complex d = cube();
        const int v_off = static_cast<int>(c.vertices.size()), e_off = static_cast<int>(c.edges.size());
        for (const Vertex& v : d.vertices) c.vertices.push_back(v);
        for (Edge e : d.edges) { e.v0 += v_off; e.v1 += v_off; c.edges.push_back(e); }
        for (Face f : d.faces) {
            for (Loop& l : f.loops)
                for (Coedge& ce : l) ce.edge += e_off;
            c.faces.push_back(f);
        }
        const AuditReport r = c.audit();
        AK_CHECK(r.watertight());
        AK_CHECK(r.shells.size() == 2u);
    }

    // Valid solids have one fan of faces at every vertex.
    AK_CHECK(cube().audit().nonmanifold_vertices == 0);
    AK_CHECK(torus().audit().nonmanifold_vertices == 0);
    AK_CHECK(block_with_hole().audit().nonmanifold_vertices == 0);

    // Two cubes sharing one corner vertex: every edge is fine, but the
    // boundary is pinched at that vertex (two fans), so it is not a manifold
    // solid. The per-shell Euler counts alone would not show it.
    {
        Complex c = cube();
        const Complex d = cube();
        // Vertex 7 of the second cube (1,1,1) becomes vertex 0 of the first
        // when the second cube is placed at (1,1,1) - (1,1,1) = corner contact.
        std::vector<int> map(8);
        for (int i = 0; i < 8; ++i) map[static_cast<std::size_t>(i)] = i == 0 ? 7 : c.add_vertex();
        const int e_off = static_cast<int>(c.edges.size());
        for (Edge e : d.edges) {
            e.v0 = map[static_cast<std::size_t>(e.v0)];
            e.v1 = map[static_cast<std::size_t>(e.v1)];
            c.edges.push_back(e);
        }
        for (Face f : d.faces) {
            for (Loop& l : f.loops)
                for (Coedge& ce : l) ce.edge += e_off;
            c.faces.push_back(f);
        }
        const AuditReport r = c.audit();
        AK_CHECK(r.boundary_edges == 0 && r.nonmanifold_edges == 0 && r.misoriented_edges == 0);
        AK_CHECK(r.shells.size() == 2u && r.shells[0].genus == 0 && r.shells[1].genus == 0);
        AK_CHECK(r.nonmanifold_vertices == 1);
        AK_CHECK(!r.watertight());
    }

    // A third face on an edge: non-manifold.
    {
        Complex c = cube();
        const int extra = c.add_vertex();
        c.add_polygon_face(9, {{0, 1, extra}});
        const AuditReport r = c.audit();
        AK_CHECK(!r.watertight());
        AK_CHECK(r.nonmanifold_edges == 1);
    }

    // A loop that does not close breaks the chain condition.
    {
        Complex c;
        for (int i = 0; i < 4; ++i) c.add_vertex();
        const int e0 = c.add_edge(0, 1), e1 = c.add_edge(1, 2), e2 = c.add_edge(3, 0);
        c.add_face(0, {{{e0, false}, {e1, false}, {e2, false}}});
        const AuditReport r = c.audit();
        AK_CHECK(!r.loops_closed);
        AK_CHECK(!r.chain_condition);
        AK_CHECK(!r.watertight());
    }

    // Tolerance ledger: recorded merges are visible in the audit.
    {
        Complex c = cube();
        c.edges[3].tol = 2e-6;
        c.vertices[5].tol = 5e-7;
        const AuditReport r = c.audit();
        AK_CHECK(r.watertight());
        AK_CHECK(r.ledger_entries == 2);
        AK_CHECK(r.max_tolerance == 2e-6);
    }

    return akt::finish("complex_test");
}
