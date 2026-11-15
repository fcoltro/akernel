// Boolean demo: a box with a hill on top, cut by a tilted box, and then the
// result cut again by a second box (a chained operation).
//
// Writes Wavefront OBJ files into the directory given as the first argument
// (default: current directory). Each file holds the watertight tessellation
// of one result (op/tessellate.hpp), one group per face, and the model's
// edges as polylines. The two faces along an edge share its points, so the
// mesh has no cracks; the program checks that before writing.
#include <akernel/op/boolean.hpp>
#include <akernel/op/tessellate.hpp>

#include <cstdio>
#include <fstream>
#include <string>

using namespace ak;
using geom::Vec3d;

namespace {

bool save_mesh(const std::string& path, const op::Brep& b, double chord) {
    op::TessOptions opt;
    opt.chord = chord;
    const op::TessResult t = op::tessellate(b, opt);
    if (!t.ok) {
        std::printf("  %s: not tessellated: %s\n", path.c_str(), t.error.c_str());
        return false;
    }
    const op::MeshAudit a = op::audit(t.mesh);
    std::ofstream out(path);
    op::write_obj(out, t.mesh);
    std::printf("  %s: %d triangles, %d points, closed %s, Euler characteristic %d, deviation %.1e, %.0f ms\n",
                path.c_str(), a.triangles, a.points, a.watertight() ? "yes" : "NO", a.euler,
                op::mesh_deviation(t.mesh, b), t.stats.ms);
    return a.watertight();
}

void report(const char* name, const op::BooleanResult& r) {
    if (!r.ok) {
        std::printf("%s: not solved: %s\n", name, r.error.c_str());
        return;
    }
    std::printf("%s: %s, %d faces, %d edges, %d vertices, %d shell(s), volume %.9f\n", name,
                r.empty ? "empty" : "closed and oriented", static_cast<int>(r.solid.faces.size()),
                static_cast<int>(r.solid.edges.size()), static_cast<int>(r.solid.vertices.size()),
                static_cast<int>(r.audit.shells.size()), op::volume(r.solid));
}

} // namespace

int main(int argc, char** argv) {
    const std::string dir = argc > 1 ? std::string(argv[1]) + "/" : std::string();

    const op::Brep A = op::to_brep(op::make_bump_box({0, 0, 0}, 1, 1, 1, 0.3));
    // A box tilted by 0.5 rad about the x axis and 0.4 rad about the z axis.
    const double cx = std::cos(0.5), sx = std::sin(0.5), cz = std::cos(0.4), sz = std::sin(0.4);
    const Vec3d ex{0.9 * cz, 0.9 * sz, 0.0};
    const Vec3d ey{-0.8 * sz * cx, 0.8 * cz * cx, 0.8 * sx};
    const Vec3d ez{0.7 * sz * sx, -0.7 * cz * sx, 0.7 * cx};
    const op::Brep B = op::to_brep(op::make_box({0.45, 0.2, 0.75}, ex, ey, ez));

    const op::Imprint im = op::imprint(A, B);
    if (!im.ok) {
        std::printf("The kernel did not solve this configuration: %s\n", im.error.c_str());
        return 1;
    }
    std::printf("Imprint: %d intersection curves, %d new edges, %d new vertices, %d regions.\n", im.stats.curves,
                im.stats.new_arcs, im.stats.new_vertices, im.stats.regions);
    std::printf("Certified bound on the gap between the two surfaces along the new edges: %.1e\n",
                op::edge_gap_bound(im));

    const struct { op::Op o; const char* name; } ops[3] = {{op::Op::Intersection, "intersection"},
                                                             {op::Op::Union, "union"},
                                                             {op::Op::Difference, "difference"}};
    double vol[3] = {0, 0, 0};
    op::Brep difference;
    for (int k = 0; k < 3; ++k) {
        const op::BooleanResult r = op::select(im, ops[k].o);
        report(ops[k].name, r);
        if (!r.ok) return 1;
        vol[k] = op::volume(r.solid);
        if (!r.empty) save_mesh(dir + "boolean_" + ops[k].name + ".obj", r.solid, 1e-3);
        if (ops[k].o == op::Op::Difference) difference = r.solid;
    }
    const double va = op::volume(A), vb = op::volume(B);
    std::printf("Check: vol(A) + vol(B) - vol(A and B) - vol(A or B) = %.1e;  vol(A) - vol(A and B) - vol(A minus B) = %.1e\n",
                va + vb - vol[0] - vol[1], va - vol[0] - vol[2]);

    // Chained: drill a square hole through the difference. The second cut
    // crosses edges that the first cut created.
    const op::Brep C = op::to_brep(op::make_box({0.3, 0.55, -0.5}, {0.62, 0.8, 2.0}));
    const op::BooleanResult chained = op::boolean_op(difference, C, op::Op::Difference);
    report("(A minus B) minus C", chained);
    if (chained.ok) {
        std::printf("  certified edge gap bound %.1e, widest vertex enclosure %.1e\n", op::edge_gap_bound(chained.solid),
                    op::vertex_width(chained.solid));
        save_mesh(dir + "boolean_chained.obj", chained.solid, 1e-3);
    }
    return chained.ok ? 0 : 1;
}
