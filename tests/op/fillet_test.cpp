// A fillet applied to a solid (op/fillet.hpp): one edge of a box rounded by
// the exact rolling-ball fillet, with the supports trimmed along the contact
// lines. The result must be a closed, regular shell with the exact volume,
// tessellate watertight, and work as an operand of the Booleans (cuts
// through the fillet included), each volume against a closed form.
#include <akernel/op/boolean.hpp>
#include <akernel/op/fillet.hpp>
#include <akernel/op/tessellate.hpp>

#include <cmath>
#include <cstdio>

#include "check.hpp"

using namespace ak;
using geom::Vec3d;

namespace {

/// Area of the corner square [0, r]^2 outside the quarter disc of radius r
/// about its far corner, for the part with offset t in [0, w] from the
/// disc's centre line: integral of r - sqrt(r^2 - t^2).
double corner_area(double r, double w) {
    return r * w - (0.5 * w * std::sqrt(r * r - w * w) + 0.5 * r * r * std::asin(w / r));
}

void check_solid(const char* name, const op::Brep& b, double volume) {
    const auto rep = op::audit(b);
    AK_CHECK(rep.closed() && rep.euler == 2);
    const double v = op::volume(b);
    AK_CHECK_NEAR(v, volume, 1e-12);
    const op::TessResult t = op::tessellate(b);
    AK_CHECK(t.ok);
    double dev = -1;
    if (t.ok) {
        const op::MeshAudit m = op::audit(t.mesh);
        AK_CHECK(m.watertight());
        dev = op::mesh_deviation(t.mesh, b);
        AK_CHECK(dev < 2e-3);
    }
    std::printf("%s: %zu faces, volume %.15f (closed form %.15f), mesh deviation %.2e\n", name, b.faces.size(), v,
                volume, dev);
}

} // namespace

int main() {
    const double a = 2.0, b = 1.5, c = 1.0, r = 0.5;
    const double V = a * b * c - c * r * r * (1.0 - M_PI / 4.0);

    const op::FilletedSolid f = op::fillet_box_edge({0, 0, 0}, {a, 0, 0}, {0, b, 0}, {0, 0, c}, r);
    AK_CHECK(f.solid.has_value());
    if (!f.solid) {
        std::printf("fillet_box_edge: %s\n", f.error.c_str());
        return akt::finish("fillet_test");
    }
    const op::Brep B = op::to_brep(*f.solid);

    // Exact data: every control of every face is a dyadic, so the fillet
    // and the trimmed faces meet exactly (one shared definition per edge).
    for (const auto& P : f.solid->faces)
        for (const auto& h : P.c)
            for (double x : {h.x, h.y, h.z, h.w}) AK_CHECK(std::ldexp(x, 20) == std::round(std::ldexp(x, 20)));
    // Regular everywhere: no face has a singular corner.
    double jmin = 1e300;
    for (const auto& P : f.solid->faces) {
        const geom::BezierSurface S(P);
        for (int i = 0; i <= 20; ++i)
            for (int j = 0; j <= 20; ++j) {
                const geom::Jet3 J = S.jet(i / 20.0, j / 20.0);
                jmin = std::min(jmin, geom::norm(geom::cross(J[1], J[2])));
            }
    }
    std::printf("filleted box: min |Su x Sv| over all faces %.3g\n", jmin);
    AK_CHECK(jmin > 1e-3);
    std::string why;
    AK_CHECK(op::check_edge_fans(B, why));
    check_solid("filleted box", B, V);

    // As an operand.
    struct Case {
        const char* name;
        op::Brep tool;
        op::Op o;
        double volume;
    };
    const double w = 0.25; // a cut at x = a - r + w crosses the fillet
    const Case cases[] = {
        {"half (z <= 1/2): plane across the fillet", op::to_brep(op::make_box({-1, -1, -1}, {3, 3, 0.5})),
         op::Op::Intersection, V / 2},
        {"cut x >= 7/4: plane along the fillet's axis", op::to_brep(op::make_box({a - r + w, -1, -1}, {3, 3, 2})),
         op::Op::Difference, ((a - r + w) * b - corner_area(r, w)) * c},
        {"union with a block on top", op::to_brep(op::make_box({0, 0, 1}, {1, 1, 1.5})), op::Op::Union, V + 0.5},
    };
    for (const Case& k : cases) {
        const op::BooleanResult res = op::boolean_op(B, k.tool, k.o);
        AK_CHECK(res.ok);
        if (!res.ok) {
            std::printf("FAIL %s: %s\n", k.name, res.error.c_str());
            continue;
        }
        check_solid(k.name, res.solid, k.volume);
    }

    // Refusals.
    AK_CHECK(!op::fillet_box_edge({0, 0, 0}, {a, 0, 0}, {0, b, 0}, {0, 0, c}, 1.2).solid); // r > 2b/3
    AK_CHECK(!op::fillet_box_edge({0, 0, 0}, {a, 0, 0}, {0, b, 0}, {0, 0, c}, 0.0).solid);
    AK_CHECK(!op::fillet_box_edge({0, 0, 0}, {0, b, 0}, {a, 0, 0}, {0, 0, c}, r).solid); // left-handed

    return akt::finish("fillet_test");
}
