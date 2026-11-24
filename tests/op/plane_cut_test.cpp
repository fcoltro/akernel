// Plane cuts of solids, including solids with procedural faces
// (op/plane_cut.hpp): the first step of Booleans on procedural faces. Every
// result must be a closed shell that meshes watertight, with the volume of
// a closed form or of an independent slice integration, and must not move
// with the solid (closure).
#include <akernel/op/plane_cut.hpp>
#include <akernel/op/proc_fillet.hpp>
#include <akernel/op/tessellate.hpp>

#include <cmath>
#include <cstdio>

#include "check.hpp"

using namespace ak;
using geom::Vec3d;

namespace {

op::PolyPatchD top_patch(double a, double L, double h, double c, double x0 = 0, double y0 = 0) {
    op::PolyPatchD T(1, 3);
    static constexpr double s[4] = {0, 0, 1, 1};
    for (int i = 0; i <= 1; ++i)
        for (int j = 0; j <= 3; ++j) T.at(i, j) = {x0 + a * i, y0 + L * j / 3.0, h + c * a * i * s[j], 1.0};
    return T;
}

/// Check a cut: closed, watertight, and the volume returned.
double check_cut(const char* name, const op::PlaneCutResult& r, int faces) {
    if (!r.solid) std::printf("%s: refused: %s\n", name, r.error.c_str());
    AK_CHECK(r.solid.has_value());
    if (!r.solid) return 0.0;
    const auto rep = op::audit(*r.solid);
    AK_CHECK(rep.closed() && rep.euler == 2);
    AK_CHECK(static_cast<int>(r.solid->faces.size()) == faces);
    const double v = op::volume(*r.solid);
    const op::TessResult t = op::tessellate(*r.solid);
    if (!t.ok) std::printf("  tessellation failed: %s\n", t.error.c_str());
    AK_CHECK(t.ok);
    double dev = -1.0;
    if (t.ok) {
        AK_CHECK(op::audit(t.mesh).watertight());
        dev = op::mesh_deviation(t.mesh, *r.solid);
        AK_CHECK(dev < 2e-3);
    }
    std::printf("%s: %zu faces, %d section curves, closed %d, volume %.15f, mesh deviation %.2e, %.2f s\n", name,
                r.solid->faces.size(), r.section_curves, static_cast<int>(rep.closed()), v, dev, r.seconds);
    return v;
}

/// Removed volume of the rounded block over y in [0, y_end] by slices (as in
/// proc_fillet_test).
double removed_by_slices(const solve::BlendCarrier& blend, double a, double L, double h, double c, double y_end) {
    auto cut_area = [&](double y0, int n) {
        const double t = y0 / L, sy = t * t * (3 - 2 * t);
        std::vector<Vec3d> poly;
        for (int i = 0; i <= n; ++i) {
            const double th = static_cast<double>(i) / n;
            Vec3d P, Ps, Pt;
            auto f = [&](double sg) {
                blend.eval(sg, th, P, Ps, Pt);
                return P.y - y0;
            };
            double lo = 0.0, hi = 1.0, flo = f(lo), fhi = f(hi);
            int side = 0;
            for (int it = 0; it < 100 && hi - lo > 1e-15; ++it) {
                const double m = (lo * fhi - hi * flo) / (fhi - flo), fm = f(m);
                if (fm == 0.0) { lo = hi = m; break; }
                if ((fm > 0.0) == (flo > 0.0)) { lo = m; flo = fm; if (side == -1) fhi /= 2; side = -1; }
                else { hi = m; fhi = fm; if (side == 1) flo /= 2; side = 1; }
                if (std::fabs(fm) < 1e-14) { lo = hi = m; break; }
            }
            f(0.5 * (lo + hi));
            poly.push_back(P);
        }
        poly.push_back({a, y0, h + c * a * sy});
        double A = 0.0;
        for (std::size_t i = 0; i < poly.size(); ++i) {
            const Vec3d& p = poly[i];
            const Vec3d& q = poly[(i + 1) % poly.size()];
            A += p.x * q.z - q.x * p.z;
        }
        return std::fabs(A) / 2;
    };
    static constexpr double gx[4] = {0.1834346424956498, 0.5255324099163290, 0.7966664774136267, 0.9602898564975363};
    static constexpr double gw[4] = {0.3626837833783620, 0.3137066458778873, 0.2223810344533745, 0.1012285362903763};
    const int cells = 8;
    double total = 0.0;
    for (int k = 0; k < cells; ++k)
        for (int i = 0; i < 4; ++i)
            for (int sgn = -1; sgn <= 1; sgn += 2) {
                const double y = y_end * (k + 0.5 + 0.5 * sgn * gx[i]) / cells;
                const double a1 = cut_area(y, 300), a2 = cut_area(y, 600);
                total += 0.5 * gw[i] * (y_end / cells) * (a2 + (a2 - a1) / 3);
            }
    return total;
}

} // namespace

int main() {
    // A unit cube cut through its centre by x + y + z <= 3/2: half the cube,
    // a hexagonal cap.
    {
        const op::Brep cube = op::to_brep(op::make_box(Vec3d{0, 0, 0}, Vec3d{1, 1, 1}));
        const double v = check_cut("cube, x + y + z <= 3/2", op::cut_by_plane(cube, {1, 1, 1}, 1.5), 7);
        AK_CHECK_NEAR(v, 0.5, 1e-12);
    }

    // The rounded freeform block (a = 1, L = 2, h = 3/4, r = 1/4, c = 1/2) cut
    // across the edge at y = 1: the section runs through the blend face.
    const double a = 1.0, L = 2.0, h = 0.75, r = 0.25, c = 0.5;
    const op::ProcFilletResult rb = op::round_block_edge(top_patch(a, L, h, c), r);
    AK_CHECK(rb.solid.has_value());
    if (rb.solid) {
        const double y1 = 1.0, T = y1 / L;
        const double block = a * y1 * h + c * a * a / 2 * L * (T * T * T - T * T * T * T / 2);
        const double v = check_cut("rounded freeform block, y <= 1", op::cut_by_plane(*rb.solid, {0, 1, 0}, y1), 7);
        const double slices = removed_by_slices(*rb.blend, a, L, h, c, y1);
        std::printf("  removed %.12f, by slices %.12f, difference %.2e\n", block - v, slices, block - v - slices);
        AK_CHECK_NEAR(block - v, slices, 1e-9);

        // A slanted cut through the blend, the top and the side: the same
        // volume when solid and plane are moved together.
        const Vec3d nrm{0.3, 1.0, 0.4};
        const double d = 1.2;
        const double v1 = check_cut("rounded freeform block, slanted cut", op::cut_by_plane(*rb.solid, nrm, d), 7);
        const Vec3d shift{0.5, -0.5, 0.0};
        const op::ProcFilletResult moved = op::round_block_edge(top_patch(a, L, h, c, shift.x, shift.y), r);
        if (moved.solid) {
            const double v2 = check_cut("same, moved", op::cut_by_plane(*moved.solid, nrm, d + dot(nrm, shift)), 7);
            std::printf("  difference %.2e\n", v2 - v1);
            AK_CHECK_NEAR(v1, v2, 1e-10);
        }
    }
    return akt::finish("plane_cut_test");
}
