// Fillets on general solids (op/proc_fillet.hpp, round_edge): any convex
// side edge of a Brep rounded by a rolling ball, with the blend trimmed at
// both ends by its certified curves with the end faces. Each result must be
// a closed shell that meshes watertight; volumes are checked against closed
// forms on a box, against round_block_edge on the same freeform edge,
// against a mirrored block, and against translation.
#include <akernel/op/boolean.hpp>
#include <akernel/op/proc_fillet.hpp>
#include <akernel/op/tessellate.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>

#include "check.hpp"

using namespace ak;

namespace {

constexpr double kSmooth[4] = {0, 0, 1, 1};
constexpr double kLinear[4] = {0, 1.0 / 3, 2.0 / 3, 1};

/// Top over [x0, x0 + a] x [y0, y0 + L], degree (1, 3): z = h + c X(x) s(y),
/// X = (x - x0) / a, or 1 - X when mirrored.
op::PolyPatchD top_patch(double a, double L, double h, double c, double x0 = 0, double y0 = 0,
                         const double (&s)[4] = kSmooth, bool mirrored = false) {
    op::PolyPatchD T(1, 3);
    for (int i = 0; i <= 1; ++i)
        for (int j = 0; j <= 3; ++j)
            T.at(i, j) = {x0 + a * i, y0 + L * j / 3.0, h + c * a * (mirrored ? 1 - i : i) * s[j], 1.0};
    return T;
}

/// The edge shared by carriers c0 and c1 (freeform_block numbering: 0 bottom,
/// 1 side x = 0, 2 side x = a, 3 top, 4 end y = 0, 5 end y = L).
int edge_between(const op::Brep& b, int c0, int c1) {
    for (std::size_t i = 0; i < b.edges.size(); ++i) {
        const op::BEdge& e = b.edges[i];
        if ((e.on[0].carrier == c0 && e.on[1].carrier == c1) || (e.on[0].carrier == c1 && e.on[1].carrier == c0))
            return static_cast<int>(i);
    }
    return -1;
}

double seconds_since(std::chrono::steady_clock::time_point t0) {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
}

/// The same solid with carrier c reparameterised by u -> 1 - u and its face
/// marked flipped: sides, edge references, charts and the loop follow.
op::Brep flip_carrier(op::Brep b, int c) {
    op::PolyPatchD& P = b.carriers[static_cast<std::size_t>(c)];
    const op::PolyPatchD Q = P;
    for (int i = 0; i <= P.p; ++i)
        for (int j = 0; j <= P.q; ++j) P.at(i, j) = Q.at(P.p - i, j);
    for (op::BEdge& e : b.edges)
        for (auto& ref : e.on)
            if (ref.carrier == c) {
                if (ref.side == 1 || ref.side == 3) ref.side = 4 - ref.side;
                else ref.reversed = !ref.reversed;
            }
    for (op::BVertex& v : b.vertices)
        for (op::Chart& ch : v.charts)
            if (ch.carrier == c) ch.u = num::Interval(1.0) - ch.u;
    for (op::BFace& f : b.faces)
        if (f.carrier == c && f.proc < 0) {
            f.flipped = !f.flipped;
            for (auto& L : f.loops) {
                std::reverse(L.begin(), L.end());
                for (op::BCoedge& ce : L) ce.reversed = !ce.reversed;
            }
        }
    return b;
}

/// Round the edge between carriers c0 and c1 of the block under `top`; check
/// it and return its volume (NaN if refused).
double rounded(const op::Brep& block, int c0, int c1, double r, const char* name) {
    const auto t0 = std::chrono::steady_clock::now();
    const int e = edge_between(block, c0, c1);
    AK_CHECK(e >= 0);
    const op::ProcFilletResult res = op::round_edge(block, e, r);
    if (!res.solid) std::printf("%s: refused: %s\n", name, res.error.c_str());
    AK_CHECK(res.solid.has_value());
    if (!res.solid) return std::nan("");
    const auto rep = op::audit(*res.solid);
    AK_CHECK(rep.closed() && rep.euler == 2);
    const double v = op::volume(*res.solid);
    std::printf("%s: %zu faces, %zu edges, %zu vertices, spine %zu segments, %.2f s\n", name, res.solid->faces.size(),
                res.solid->edges.size(), res.solid->vertices.size(), res.blend->segments(), seconds_since(t0));
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
    return v;
}

double rounded(const op::PolyPatchD& top, int c0, int c1, double r, const char* name) {
    return rounded(op::to_brep(op::freeform_block(top)), c0, c1, r, name);
}

} // namespace

int main() {
    const double a = 1.0, L = 2.0, h = 0.75, r = 0.25, c = 0.5;
    const double q = (1 - M_PI / 4) * r * r; // cross-section removed by a quarter-round

    // A box: a vertical edge, a top edge along x, a top edge along y.
    {
        const op::PolyPatchD box = top_patch(a, L, h, 0.0);
        const double V = a * L * h;
        double v = rounded(box, 2, 4, r, "box, vertical edge x = a, y = 0");
        std::printf("  volume %.15f, closed form %.15f, error %.2e\n", v, V - q * h, v - (V - q * h));
        AK_CHECK_NEAR(v, V - q * h, 1e-10);
        v = rounded(box, 3, 5, r, "box, top edge at y = L");
        std::printf("  volume %.15f, closed form %.15f, error %.2e\n", v, V - q * a, v - (V - q * a));
        AK_CHECK_NEAR(v, V - q * a, 1e-10);
        v = rounded(box, 1, 3, r, "box, top edge at x = 0");
        std::printf("  volume %.15f, closed form %.15f, error %.2e\n", v, V - q * L, v - (V - q * L));
        AK_CHECK_NEAR(v, V - q * L, 1e-10);
    }
    // A face whose outward normal is -(Su x Sv): the side x = a reversed in u.
    {
        const op::Brep flipped = flip_carrier(op::to_brep(op::freeform_block(top_patch(a, L, h, 0.0))), 2);
        const auto rep = op::audit(flipped);
        AK_CHECK(rep.closed() && rep.euler == 2);
        const double V = a * L * h, v = rounded(flipped, 2, 4, r, "box, side x = a flipped, vertical edge");
        std::printf("  volume %.15f, closed form %.15f, error %.2e\n", v, V - q * h, v - (V - q * h));
        AK_CHECK_NEAR(v, V - q * h, 1e-10);
    }
    // Two edges in turn: the second rounding reads the first one's result.
    {
        const op::Brep box = op::to_brep(op::freeform_block(top_patch(a, L, h, 0.0)));
        const op::ProcFilletResult first = op::round_edge(box, edge_between(box, 2, 4), r);
        AK_CHECK(first.solid.has_value());
        if (first.solid) {
            const double V = a * L * h, v = rounded(*first.solid, 1, 5, r, "box, two vertical edges");
            std::printf("  volume %.15f, closed form %.15f, error %.2e\n", v, V - 2 * q * h, v - (V - 2 * q * h));
            AK_CHECK_NEAR(v, V - 2 * q * h, 1e-10);
        }
    }
    // The freeform edge of proc_fillet_test, through round_edge: the same volume.
    {
        const op::ProcFilletResult ref = op::round_block_edge(top_patch(a, L, h, c), r);
        AK_CHECK(ref.solid.has_value());
        const double v = rounded(top_patch(a, L, h, c), 2, 3, r, "freeform top, edge x = a");
        if (ref.solid) {
            const double w = op::volume(*ref.solid);
            std::printf("  volume %.15f, round_block_edge %.15f, difference %.2e\n", v, w, v - w);
            AK_CHECK_NEAR(v, w, 1e-10);
        }
    }
    // The freeform edge at x = 0 against the mirrored block's edge at x = a.
    {
        const double v = rounded(top_patch(a, L, h, c), 1, 3, r, "freeform top, edge x = 0");
        const op::ProcFilletResult ref = op::round_block_edge(top_patch(a, L, h, c, 0, 0, kSmooth, true), r);
        AK_CHECK(ref.solid.has_value());
        if (ref.solid) {
            const double w = op::volume(*ref.solid);
            std::printf("  volume %.15f, mirrored %.15f, difference %.2e\n", v, w, v - w);
            AK_CHECK_NEAR(v, w, 1e-10);
        }
    }
    // A sloped top (z = h + c x y / L): its edge with the end y = L, where the
    // dihedral angle varies along x, and the same solid moved.
    {
        const double v = rounded(top_patch(a, L, h, c, 0, 0, kLinear), 3, 5, r, "sloped top, edge y = L");
        const double w = rounded(top_patch(a, L, h, c, 2.0, -1.5, kLinear), 3, 5, r, "  moved");
        const double block = a * L * h + c * a * a / 2 * L / 2;
        std::printf("  volume %.15f, moved %.15f, difference %.2e, removed %.6f\n", v, w, v - w, block - v);
        AK_CHECK_NEAR(v, w, 1e-10);
        AK_CHECK(v < block && block - v < r * r * a * 1.5);
    }
    // Refusals: a concave edge does not exist on a block, so a radius too
    // large for the face.
    {
        const op::Brep block = op::to_brep(op::freeform_block(top_patch(a, L, h, 0.0)));
        AK_CHECK(!op::round_edge(block, edge_between(block, 2, 4), 1.2).solid);
    }
    return akt::finish("round_edge_test");
}
