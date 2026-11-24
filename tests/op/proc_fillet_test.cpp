// A freeform fillet inside a solid (op/proc_fillet.hpp): a block whose top is
// a polynomial graph has its edge with a side rounded by a rolling ball; the
// blend is a procedural face (the certified spine plus the section arc), the
// supports are trimmed along the contact curves. The solid must be a closed
// shell that meshes watertight; its volume is checked against the closed form when the top is flat
// (the blend is then an exact quarter cylinder), and on a freeform top
// against translation (the divergence sum of a closed surface does not move
// with the origin) and against an independent slice-by-slice integration.
#include <akernel/op/boolean.hpp>
#include <akernel/op/proc_fillet.hpp>
#include <akernel/op/tessellate.hpp>

#include <chrono>
#include <cmath>
#include <cstdio>

#include "check.hpp"

using namespace ak;
using geom::Vec3d;

namespace {

/// Top over [x0, x0 + a] x [y0, y0 + L], degree (1, 3): z = h + c (x - x0) s(y)
/// with s the smoothstep 3t^2 - 2t^3 (Bernstein coefficients 0, 0, 1, 1),
/// so the slope across the end planes is zero.
op::PolyPatchD top_patch(double a, double L, double h, double c, double x0 = 0, double y0 = 0) {
    op::PolyPatchD T(1, 3);
    static constexpr double s[4] = {0, 0, 1, 1};
    for (int i = 0; i <= 1; ++i)
        for (int j = 0; j <= 3; ++j) T.at(i, j) = {x0 + a * i, y0 + L * j / 3.0, h + c * a * i * s[j], 1.0};
    return T;
}

double seconds_since(std::chrono::steady_clock::time_point t0) {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
}

op::ProcFilletResult build(const op::PolyPatchD& top, double r, const char* name, double cut = 1.0) {
    const auto t0 = std::chrono::steady_clock::now();
    op::ProcFilletResult res = op::round_block_edge(top, r, {}, cut);
    if (!res.solid) std::printf("%s: refused: %s\n", name, res.error.c_str());
    AK_CHECK(res.solid.has_value());
    if (res.solid) {
        const auto rep = op::audit(*res.solid);
        AK_CHECK(rep.closed() && rep.euler == 2);
        std::printf("%s: %zu faces, %zu edges, %zu vertices, spine %zu segments, %.2f s\n", name,
                    res.solid->faces.size(), res.solid->edges.size(), res.solid->vertices.size(),
                    res.blend->segments(), seconds_since(t0));
        // The mesh: watertight, consistently oriented, close to the surfaces.
        const auto t1 = std::chrono::steady_clock::now();
        const op::TessResult t = op::tessellate(*res.solid);
        if (!t.ok) std::printf("  tessellation failed: %s\n", t.error.c_str());
        AK_CHECK(t.ok);
        if (t.ok) {
            const op::MeshAudit m = op::audit(t.mesh);
            AK_CHECK(m.watertight());
            const double dev = op::mesh_deviation(t.mesh, *res.solid);
            AK_CHECK(dev < 2e-3);
            std::printf("  mesh: %d triangles, watertight %d, deviation %.2e, %.2f s\n", m.triangles,
                        static_cast<int>(m.watertight()), dev, seconds_since(t1));
        }
    }
    return res;
}

/// The removed volume measured independently of the divergence sum: the
/// area cut from each section y = y0 (between the top, which is straight in
/// x there, the side x = a and the blend's section curve, found by solving
/// P(sigma, theta).y = y0 per theta), by a shoelace on n chords with
/// Richardson extrapolation, integrated over y by Gauss 8 on 8 cells.
double removed_by_slices(const solve::BlendCarrier& blend, double a, double L, double h, double c, double y_end) {
    auto cut_area = [&](double y0, int n) {
        const double t = y0 / L, sy = t * t * (3 - 2 * t);
        std::vector<Vec3d> poly;
        for (int i = 0; i <= n; ++i) {
            const double th = static_cast<double>(i) / n;
            // y is monotone in sigma along each section: regula falsi
            // (Illinois) on [0, 1], which does not need the speed in sigma
            // (uniform per spine segment, so it jumps between them).
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
                if ((fm > 0.0) == (flo > 0.0)) {
                    lo = m; flo = fm;
                    if (side == -1) fhi /= 2;
                    side = -1;
                } else {
                    hi = m; fhi = fm;
                    if (side == 1) flo /= 2;
                    side = 1;
                }
                if (std::fabs(fm) < 1e-14) { lo = hi = m; break; }
            }
            f(0.5 * (lo + hi));
            poly.push_back(P);
        }
        poly.push_back({a, y0, h + c * a * sy}); // the corner the blend replaces
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
    const double a = 1.0, L = 2.0, h = 0.75, r = 0.25;

    // Flat top: the blend is a quarter cylinder along y.
    {
        auto res = build(top_patch(a, L, h, 0.0), r, "flat top");
        if (res.solid) {
            const double v = op::volume(*res.solid), exact = a * L * h - (1 - M_PI / 4) * r * r * L;
            std::printf("  volume %.15f, closed form %.15f, error %.2e\n", v, exact, v - exact);
            AK_CHECK_NEAR(v, exact, 1e-10);
        }
    }

    // Freeform top: the dihedral angle at the edge varies along it.
    const double c = 0.5;
    double v_free = 0.0;
    {
        auto res = build(top_patch(a, L, h, c), r, "freeform top");
        if (res.solid) {
            v_free = op::volume(*res.solid);
            const double block = a * L * h + c * a * a / 2 * L / 2; // integral of s over [0, L] is L / 2
            const op::Brep unrounded = op::to_brep(op::freeform_block(top_patch(a, L, h, c)));
            const auto rep = op::audit(unrounded);
            AK_CHECK(rep.closed() && rep.euler == 2);
            AK_CHECK_NEAR(op::volume(unrounded), block, 1e-12);
            std::printf("  volume %.15f, block %.15f, removed %.6f\n", v_free, block, block - v_free);
            AK_CHECK(v_free < block && block - v_free < r * r * L);
            const double slices = removed_by_slices(*res.blend, a, L, h, c, L);
            std::printf("  removed by slices %.12f, difference %.2e\n", slices, block - v_free - slices);
            AK_CHECK_NEAR(block - v_free, slices, 1e-9);
            // Booleans do not read procedural faces yet: they refuse.
            const op::Brep box = op::to_brep(op::freeform_block(top_patch(a, L, h, 0.0)));
            AK_CHECK(!op::boolean_op(*res.solid, box, op::Op::Union).error.empty());
        }
    }
    // The same solid moved in x and y: the same volume.
    {
        auto res = build(top_patch(a, L, h, c, 3.0, -5.0), r, "freeform top, moved");
        if (res.solid) {
            const double v = op::volume(*res.solid);
            std::printf("  volume %.15f, difference %.2e\n", v, v - v_free);
            AK_CHECK_NEAR(v, v_free, 1e-10);
        }
    }
    // The block cut at y = 3L/4, where the top is not flat across the end:
    // the blend leaves the end plane and is trimmed by a traced curve.
    const double cut = 0.75, y1 = cut * L;
    {
        auto res = build(top_patch(a, L, h, 0.0), r, "flat top, cut at 3L/4", cut);
        if (res.solid) {
            const double v = op::volume(*res.solid), exact = a * y1 * h - (1 - M_PI / 4) * r * r * y1;
            std::printf("  volume %.15f, closed form %.15f, error %.2e\n", v, exact, v - exact);
            AK_CHECK_NEAR(v, exact, 1e-10);
        }
    }
    double v_cut = 0.0;
    {
        auto res = build(top_patch(a, L, h, c), r, "freeform top, cut at 3L/4", cut);
        if (res.solid) {
            v_cut = op::volume(*res.solid);
            const double T = cut, block = a * y1 * h + c * a * a / 2 * L * (T * T * T - T * T * T * T / 2);
            const double slices = removed_by_slices(*res.blend, a, L, h, c, y1);
            std::printf("  volume %.15f, block %.15f, removed %.12f, by slices %.12f, difference %.2e\n", v_cut, block,
                        block - v_cut, slices, block - v_cut - slices);
            AK_CHECK_NEAR(block - v_cut, slices, 1e-9);
        }
    }
    {
        auto res = build(top_patch(a, L, h, c, 0.5, -0.5), r, "freeform top, cut at 3L/4, moved", cut);
        if (res.solid) {
            const double v = op::volume(*res.solid);
            std::printf("  volume %.15f, difference %.2e\n", v, v - v_cut);
            AK_CHECK_NEAR(v, v_cut, 1e-10);
        }
    }
    // Refusals.
    {
        op::PolyPatchD T = top_patch(a, L, h, c);
        T.at(1, 1).z += 0.1; // slope across y = 0
        AK_CHECK(!op::round_block_edge(T, r).solid);
    }
    return akt::finish("proc_fillet_test");
}
