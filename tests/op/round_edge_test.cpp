// Fillets on general solids (op/proc_fillet.hpp, round_edge): any convex
// side edge of a Brep rounded by a rolling ball, with the blend trimmed at
// both ends by its certified curves with the end faces. Each result must be
// a closed shell that meshes watertight; volumes are checked against closed
// forms on a box, against round_block_edge on the same freeform edge,
// against a mirrored block, and against translation.
#include <akernel/num/detmath.hpp>
#include <akernel/op/boolean.hpp>
#include <akernel/op/proc_fillet.hpp>
#include <akernel/op/tessellate.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <numbers>
#include <vector>

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
        if (e.proc || e.arc) continue;
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
    AK_CHECK(rep.closed() && rep.euler == 2 * static_cast<int>(rep.shells.size()));
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

/// The solid b with the box c as a cavity: c's faces appended, reversed.
op::Brep with_cavity(op::Brep b, const op::Brep& c) {
    const int nc = static_cast<int>(b.carriers.size()), ne = static_cast<int>(b.edges.size()),
              nv = static_cast<int>(b.vertices.size());
    for (const op::PolyPatchD& x : c.carriers) b.carriers.push_back(x);
    for (op::BVertex v : c.vertices) {
        for (op::Chart& ch : v.charts) ch.carrier += nc;
        b.vertices.push_back(v);
    }
    for (op::BEdge e : c.edges) {
        e.v0 += nv;
        e.v1 += nv;
        for (auto& ref : e.on)
            if (ref.carrier >= 0) ref.carrier += nc;
        for (int& k : e.carrier)
            if (k >= 0) k += nc;
        b.edges.push_back(e);
    }
    for (op::BFace f : c.faces) {
        f.carrier += nc;
        f.flipped = !f.flipped;
        for (auto& L : f.loops)
            for (op::BCoedge& ce : L) ce.edge += ne;
        b.faces.push_back(f);
    }
    return b;
}

} // namespace

int main() {
    const double a = 1.0, L = 2.0, h = 0.75, r = 0.25, c = 0.5;
    const double q = (1 - std::numbers::pi / 4) * r * r; // cross-section removed by a quarter-round

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
    // Mitre corners: one edge rounded with r1, then an edge at its end with
    // r2. The two fillets meet in a mitre through the point where both touch
    // their shared face. Smaller r2: the mitre ends on the earlier fillet's
    // end curve; equal: at its end vertex; larger: on its far contact, and
    // the new fillet goes on to meet the face past it. The volume is the box
    // less both quarter-round regions plus their overlap near the corner,
    // integrated independently along the axis normal to the shared face
    // (Gauss-Legendre after x = a - rm sin^2 t, rm the smaller radius, which
    // takes out the square roots at both ends).
    {
        auto overlap = [&](double r1, double r2) {
            static constexpr double gx[4] = {0.1834346424956498, 0.5255324099163290, 0.7966664774136267, 0.9602898564975363};
            static constexpr double gw[4] = {0.3626837833783620, 0.3137066458778873, 0.2223810344533745, 0.1012285362903763};
            const double rm = std::min(r1, r2);
            auto f = [&](double x) {
                const double d1 = x - (a - r1), d2 = x - (a - r2);
                return (r1 - std::sqrt(std::max(0.0, r1 * r1 - d1 * d1))) * (r2 - std::sqrt(std::max(0.0, r2 * r2 - d2 * d2)));
            };
            double I = 0.0;
            const int n = 400;
            for (int k = 0; k < n; ++k) {
                const double t0 = std::numbers::pi / 2 * k / n, t1 = std::numbers::pi / 2 * (k + 1) / n;
                for (int i = 0; i < 4; ++i)
                    for (int sg = -1; sg <= 1; sg += 2) {
                        const double t = 0.5 * (t0 + t1) + 0.5 * sg * gx[i] * (t1 - t0);
                        I += 0.5 * (t1 - t0) * gw[i] * f(a - rm * num::det::sin(t) * num::det::sin(t)) * 2 * rm * num::det::sin(t) * num::det::cos(t);
                    }
            }
            return I;
        };
        const op::Brep box = op::to_brep(op::freeform_block(top_patch(a, L, h, 0.0)));
        // Edge (c0, c1) of length len1 with r1, then (d0, d1) of length len2.
        auto corner = [&](int c0, int c1, double len1, double r1, int d0, int d1, double len2, double r2, const char* what) {
            const op::ProcFilletResult first = op::round_edge(box, edge_between(box, c0, c1), r1);
            AK_CHECK(first.solid.has_value());
            if (!first.solid) return;
            char name[96];
            std::snprintf(name, sizeof name, "box corner, %s, r1 = %.2f, r2 = %.2f", what, r1, r2);
            const double V = a * L * h, q1 = (1 - std::numbers::pi / 4) * r1 * r1, q2 = (1 - std::numbers::pi / 4) * r2 * r2;
            const double exact = V - q1 * len1 - q2 * len2 + overlap(r1, r2);
            const double v = rounded(*first.solid, d0, d1, r2, name);
            std::printf("  volume %.15f, exact %.15f, error %.2e\n", v, exact, v - exact);
            AK_CHECK_NEAR(v, exact, 1e-10);
        };
        // The top edge at x = a, then the vertical edge at x = a, y = 0.
        for (const double r2 : {0.15, 0.1, 0.2, 0.25, 0.35}) corner(2, 3, L, r, 2, 4, h, r2, "top then vertical");
        // The other order, and the shared face on the new blend's other side.
        corner(2, 4, h, 0.15, 2, 3, L, r, "vertical then top");
        corner(3, 4, a, r, 2, 4, h, 0.15, "front top then vertical");
        corner(3, 4, a, r, 2, 4, h, r, "front top then vertical");
        corner(3, 4, a, 0.15, 2, 4, h, r, "front top then vertical");
    }
    // Three fillets at one corner: the top edge at x = a (r1), the vertical
    // edge at x = a, y = 0 (r2), then the top edge at y = 0 (r3), whose mitre
    // crosses both earlier fillets through the point where all three meet.
    // The volume by inclusion and exclusion of the three quarter-round
    // regions R1, R2, R3: pairwise overlaps as above, and R1 & R2 & R3 as the
    // integral over R2's section of min(g1(u), g3(v)), g the depth of a
    // region at distance t from its face, g(r, t) = r - sqrt(t (2r - t)),
    // after t = r (1 - cos p).
    {
        static constexpr double gx[4] = {0.1834346424956498, 0.5255324099163290, 0.7966664774136267, 0.9602898564975363};
        static constexpr double gw[4] = {0.3626837833783620, 0.3137066458778873, 0.2223810344533745, 0.1012285362903763};
        auto quad = [&](const auto& f, double lo, double hi, int n) {
            double I = 0.0;
            for (int k = 0; k < n; ++k) {
                const double t0 = lo + (hi - lo) * k / n, t1 = lo + (hi - lo) * (k + 1) / n;
                for (int i = 0; i < 4; ++i)
                    for (int sg = -1; sg <= 1; sg += 2) I += 0.5 * (t1 - t0) * gw[i] * f(0.5 * (t0 + t1) + 0.5 * sg * gx[i] * (t1 - t0));
            }
            return I;
        };
        auto g = [](double rr, double t) { return t >= rr ? 0.0 : rr - std::sqrt(std::max(0.0, t * (2 * rr - t))); };
        auto ginv = [](double rr, double y) { return rr - std::sqrt(std::max(0.0, rr * rr - (rr - y) * (rr - y))); };
        const double hp = std::numbers::pi / 2;
        auto ov = [&](double ri, double rj) {
            const double rm = std::min(ri, rj);
            return quad([&](double p) { return g(ri, rm * (1 - num::det::cos(p))) * g(rj, rm * (1 - num::det::cos(p))) * rm * num::det::sin(p); },
                        0.0, hp, 200);
        };
        auto triple = [&](double r1, double r2, double r3) {
            auto inner = [&](double u) {
                const double G1 = g(r1, u), vmax = std::min(g(r2, u), r3);
                if (vmax <= 0.0 || G1 <= 0.0) return 0.0;
                const double vs = G1 < r3 ? ginv(r3, G1) : 0.0, lo = std::min(vs, vmax);
                double tot = G1 * lo;
                if (vmax > lo) {
                    const double ca = std::acos(1 - lo / r3), cb = std::acos(std::min(1.0, 1 - vmax / r3));
                    tot += quad([&](double ch) { return g(r3, r3 * (1 - num::det::cos(ch))) * r3 * num::det::sin(ch); }, ca, cb, 20);
                }
                return tot;
            };
            // The outer integrand's kinks: g1 = r3, g2 = r3, and where the
            // switch point reaches the section's edge.
            std::vector<double> ks = {0.0, hp};
            for (double uk : {r3 < r1 ? ginv(r1, r3) : -1.0, r3 < r2 ? ginv(r2, r3) : -1.0})
                if (uk > 0.0 && uk < r2) ks.push_back(std::acos(1 - uk / r2));
            auto sw = [&](double p) {
                const double u = r2 * (1 - num::det::cos(p)), G1 = g(r1, u);
                return (G1 < r3 ? ginv(r3, G1) : 0.0) - std::min(g(r2, u), r3);
            };
            const int N = 2000;
            for (int i = 1; i <= N; ++i) {
                double lo = hp * (i - 1) / N, hi = hp * i / N;
                const double flo = sw(lo);
                if ((flo > 0.0) == (sw(hi) > 0.0)) continue;
                for (int it = 0; it < 80; ++it) {
                    const double m = 0.5 * (lo + hi);
                    if ((sw(m) > 0.0) == (flo > 0.0)) lo = m;
                    else hi = m;
                }
                ks.push_back(0.5 * (lo + hi));
            }
            std::sort(ks.begin(), ks.end());
            double tot = 0.0;
            for (std::size_t i = 0; i + 1 < ks.size(); ++i)
                tot += quad([&](double p) { return inner(r2 * (1 - num::det::cos(p))) * r2 * num::det::sin(p); }, ks[i], ks[i + 1], 60);
            return tot;
        };
        const op::Brep box = op::to_brep(op::freeform_block(top_patch(a, L, h, 0.0)));
        const double r1 = r, r2 = 0.15;
        const op::ProcFilletResult first = op::round_edge(box, edge_between(box, 2, 3), r1);
        AK_CHECK(first.solid.has_value());
        const op::ProcFilletResult second = first.solid ? op::round_edge(*first.solid, edge_between(*first.solid, 2, 4), r2)
                                                        : op::ProcFilletResult{};
        AK_CHECK(second.solid.has_value());
        for (const double r3 : {0.01, 0.05, 0.1, 0.2}) {
            if (!second.solid) break;
            char name[96];
            std::snprintf(name, sizeof name, "box corner, three fillets, r = %.2f, %.2f, %.2f", r1, r2, r3);
            auto q3 = [](double rr) { return (1 - std::numbers::pi / 4) * rr * rr; };
            const double exact = a * L * h - q3(r1) * L - q3(r2) * h - q3(r3) * a + ov(r1, r2) + ov(r1, r3) + ov(r2, r3) -
                                 triple(r1, r2, r3);
            const double v = rounded(*second.solid, 3, 4, r3, name);
            std::printf("  volume %.15f, exact %.15f, error %.2e\n", v, exact, v - exact);
            AK_CHECK_NEAR(v, exact, 1e-10);
        }
    }
    // A mitre on a freeform top: its edge at x = a rounded, then the vertical
    // edge at y = 0, and the same solid moved.
    {
        auto corner = [&](double x0, double y0, double r2, const char* name) {
            const op::Brep block = op::to_brep(op::freeform_block(top_patch(a, L, h, c, x0, y0)));
            const op::ProcFilletResult first = op::round_edge(block, edge_between(block, 2, 3), r);
            AK_CHECK(first.solid.has_value());
            return first.solid ? rounded(*first.solid, 2, 4, r2, name) : std::nan("");
        };
        for (const double r2 : {0.15, 0.35}) {
            char name[64];
            std::snprintf(name, sizeof name, "freeform top corner, mitre r2 = %.2f", r2);
            const double v = corner(0, 0, r2, name), w = corner(2.0, -1.5, r2, "  moved");
            std::printf("  volume %.15f, moved %.15f, difference %.2e\n", v, w, v - w);
            AK_CHECK_NEAR(v, w, 1e-10);
        }
    }
    // Clearance: a small cavity inside the corner region of the edge x = a
    // on top. With r = 0.25 the ball would cut it: refused. With r = 0.1 the
    // blend passes it (their carriers proven apart) and the volume is the
    // closed form less the cavity.
    {
        const op::Brep block = with_cavity(op::to_brep(op::freeform_block(top_patch(a, L, h, 0.0))),
                                           op::to_brep(op::make_box({0.9, 0.9, 0.6}, {0.97, 1.1, 0.7})));
        const auto rep = op::audit(block);
        AK_CHECK(rep.closed() && rep.euler == 4);
        const op::ProcFilletResult big = op::round_edge(block, edge_between(block, 2, 3), r);
        std::printf("cavity in the corner region, r = 0.25: %s\n", big.solid ? "built" : big.error.c_str());
        AK_CHECK(!big.solid);
        const double r0 = 0.1, cav = 0.07 * 0.2 * 0.1;
        const double v = rounded(block, 2, 3, r0, "cavity beside the blend, r = 0.10");
        const double exact = a * L * h - cav - (1 - std::numbers::pi / 4) * r0 * r0 * L;
        std::printf("  volume %.15f, closed form %.15f, error %.2e\n", v, exact, v - exact);
        AK_CHECK_NEAR(v, exact, 1e-10);
    }
    // Refusals: a concave edge does not exist on a block, so a radius too
    // large for the face.
    {
        const op::Brep block = op::to_brep(op::freeform_block(top_patch(a, L, h, 0.0)));
        AK_CHECK(!op::round_edge(block, edge_between(block, 2, 4), 1.2).solid);
    }
    return akt::finish("round_edge_test");
}
