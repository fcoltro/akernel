// Fillets and chamfers applied to solids (op/fillet.hpp): an edge of a box,
// circular edges of rings (plane meeting cylinder, the blend a torus zone or
// a cone; a full turn is a closed chain of four edges), and concave edges
// where the blend adds material (a step, a boss on a plate, a counterbore),
// and a box with every edge rounded (sphere-octant vertex blends),
// with the supports
// trimmed along the contact curves. Every result must be a closed, regular
// shell with exact (dyadic) controls and the closed-form volume, tessellate
// watertight, and work as an operand of the Booleans.
#include <akernel/op/boolean.hpp>
#include <akernel/op/fillet.hpp>
#include <akernel/op/tessellate.hpp>

#include <cmath>
#include <numbers>
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

void check_solid(const char* name, const op::Brep& b, double volume, int euler = 2) {
    const auto rep = op::audit(b);
    AK_CHECK(rep.closed() && rep.euler == euler);
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

/// Area and first moments of the section [0, a] x [0, b] with its corner
/// (a, b) cut by `cut`.
struct SectionMoments {
    double area, mx, my; ///< integrals of 1, x and y
};
SectionMoments section(double a, double b, const op::EdgeCut& cut) {
    SectionMoments m{a * b, a * b * a / 2, a * b * b / 2};
    double ar, cx, cy; // the removed corner piece and its centroid
    if (cut.kind == op::EdgeCut::Round) {
        const double r = cut.dA, q = std::numbers::pi * r * r / 4, e = 4 * r / (3 * std::numbers::pi);
        ar = r * r - q;
        cx = (r * r * (a - r / 2) - q * (a - r + e)) / ar;
        cy = (r * r * (b - r / 2) - q * (b - r + e)) / ar;
    } else {
        ar = cut.dA * cut.dB / 2;
        cx = a - cut.dB / 3;
        cy = b - cut.dA / 3;
    }
    m.area -= ar;
    m.mx -= ar * cx;
    m.my -= ar * cy;
    return m;
}

/// The same for the step section [0, a] x [0, b1] + [0, a1] x [b1, b] with
/// its inside corner (a1, b1) filled by `cut`.
SectionMoments step_section(double a, double b, double a1, double b1, const op::EdgeCut& cut) {
    SectionMoments m{a * b1 + a1 * (b - b1), a * b1 * a / 2 + a1 * (b - b1) * a1 / 2,
                     a * b1 * b1 / 2 + a1 * (b - b1) * (b1 + b) / 2};
    double ar, cx, cy; // the added corner piece and its centroid
    if (cut.kind == op::EdgeCut::Round) {
        const double r = cut.dA, q = std::numbers::pi * r * r / 4, e = 4 * r / (3 * std::numbers::pi);
        ar = r * r - q;
        cx = (r * r * (a1 + r / 2) - q * (a1 + r - e)) / ar;
        cy = (r * r * (b1 + r / 2) - q * (b1 + r - e)) / ar;
    } else {
        ar = cut.dA * cut.dB / 2;
        cx = a1 + cut.dB / 3;
        cy = b1 + cut.dA / 3;
    }
    m.area += ar;
    m.mx += ar * cx;
    m.my += ar * cy;
    return m;
}

/// Exact controls: every control of every face a dyadic rational; and no
/// face with a singular point on a 21 x 21 grid.
void check_faces(const char* name, const op::PatchSolid& s) {
    bool dyadic = true;
    double jmin = 1e300;
    for (const auto& P : s.faces) {
        for (const auto& h : P.c)
            for (double x : {h.x, h.y, h.z, h.w}) dyadic = dyadic && std::ldexp(x, 20) == std::round(std::ldexp(x, 20));
        const geom::BezierSurface S(P);
        for (int i = 0; i <= 20; ++i)
            for (int j = 0; j <= 20; ++j) {
                const geom::Jet3 J = S.jet(i / 20.0, j / 20.0);
                jmin = std::min(jmin, geom::norm(geom::cross(J[1], J[2])));
            }
    }
    AK_CHECK(dyadic);
    AK_CHECK(jmin > 1e-3);
    std::printf("%s: controls dyadic %d, min |Su x Sv| %.3g\n", name, static_cast<int>(dyadic), jmin);
}

} // namespace

int main() {
    const double a = 2.0, b = 1.5, c = 1.0, r = 0.5;
    const double V = a * b * c - c * r * r * (1.0 - std::numbers::pi / 4.0);

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

    // --- chamfers on the box ---------------------------------------------------
    for (const op::EdgeCut cut : {op::EdgeCut::flat(0.5, 0.5), op::EdgeCut::flat(0.5, 0.75)}) {
        const auto ch = op::blend_box_edge({0, 0, 0}, {a, 0, 0}, {0, b, 0}, {0, 0, c}, cut);
        AK_CHECK(ch.solid.has_value());
        if (!ch.solid) continue;
        check_faces("chamfered box", *ch.solid);
        check_solid("chamfered box", op::to_brep(*ch.solid), c * section(a, b, cut).area);
    }

    // --- rings: a circular edge between a flat face and a cylinder ------------
    // The blend is a torus zone (fillet) or a cone (chamfer); with four
    // quarters the edge is a closed chain and the ring has genus 1.
    {
        const double r0 = 1.0, r1 = 2.0, z0 = 0.0, z1 = 1.0;
        struct Ring {
            const char* name;
            op::RingEdge edge;
            op::EdgeCut cut;
            int quarters;
        };
        const Ring rings[] = {
            {"ring, outer top edge rounded, full turn", op::RingEdge::OuterTop, op::EdgeCut::round(0.25), 4},
            {"ring, outer top edge rounded, one quarter", op::RingEdge::OuterTop, op::EdgeCut::round(0.25), 1},
            {"ring, inner top edge rounded, three quarters", op::RingEdge::InnerTop, op::EdgeCut::round(0.5), 3},
            {"ring, outer bottom edge rounded, half", op::RingEdge::OuterBottom, op::EdgeCut::round(0.25), 2},
            {"ring, inner bottom edge chamfered, full turn", op::RingEdge::InnerBottom, op::EdgeCut::flat(0.25, 0.5), 4},
        };
        for (const Ring& rg : rings) {
            const auto res = op::blend_ring_edge(r0, r1, z0, z1, rg.edge, rg.cut, rg.quarters);
            AK_CHECK(res.solid.has_value());
            if (!res.solid) {
                std::printf("FAIL %s: %s\n", rg.name, res.error.c_str());
                continue;
            }
            const bool outer = rg.edge == op::RingEdge::OuterTop || rg.edge == op::RingEdge::OuterBottom;
            const SectionMoments m = section(r1 - r0, z1 - z0, rg.cut);
            // Pappus: volume = angle * integral of rho over the section.
            const double rho_moment = outer ? r0 * m.area + m.mx : r1 * m.area - m.mx;
            const double volume = rg.quarters * (std::numbers::pi / 2) * rho_moment;
            check_faces(rg.name, *res.solid);
            check_solid(rg.name, op::to_brep(*res.solid), volume, rg.quarters == 4 ? 0 : 2);
        }
        // The rounded ring as an operand: keep the half y >= 0.
        const auto full = op::blend_ring_edge(r0, r1, z0, z1, op::RingEdge::OuterTop, op::EdgeCut::round(0.25), 4);
        if (full.solid) {
            const auto res = op::boolean_op(op::to_brep(*full.solid),
                                            op::to_brep(op::make_box({-3, 0.5, -1}, {3, 3, 2})),
                                            op::Op::Intersection);
            AK_CHECK(res.ok);
            if (res.ok) {
                const double v = op::volume(res.solid);
                std::printf("rounded ring cut by y >= 1/2: %zu faces, volume %.12f, closed %d\n",
                            res.solid.faces.size(), v, static_cast<int>(res.audit.closed()));
                AK_CHECK(res.audit.closed());
                // Reference by quadrature over the section: the angle range with
                // rho sin(theta) >= 1/2 is pi - 2 asin(1 / (2 rho)).
                const double rr = 0.25, A = r1 - r0, H = z1 - z0;
                auto f = [&](double x) {
                    const double rho = r0 + x;
                    return rho * (std::numbers::pi - 2 * std::asin(0.5 / rho));
                };
                auto simpson = [](auto g, double lo, double hi, int n) {
                    double sum = g(lo) + g(hi);
                    for (int i = 1; i < n; ++i) sum += g(lo + (hi - lo) * i / n) * (i % 2 ? 4 : 2);
                    return sum * (hi - lo) / (3 * n);
                };
                const double ref = simpson([&](double x) { return f(x) * H; }, 0, A - rr, 2000) +
                                   simpson([&](double t) {
                                       const double x = A - rr + rr * std::sin(t);
                                       return f(x) * (H - rr + rr * std::cos(t)) * rr * std::cos(t);
                                   }, 0, std::numbers::pi / 2, 2000);
                std::printf("   reference by quadrature %.12f\n", ref);
                AK_CHECK_NEAR(v, ref, 1e-9);
            } else {
                std::printf("FAIL rounded ring cut: %s\n", res.error.c_str());
            }
        }
    }

    // --- concave edges: the blend adds material ---------------------------------
    {
        // Step block 2 x 1.5 x 1 without the notch x > 1, y > 3/4.
        const double A = 2.0, Bh = 1.5, a1 = 1.0, b1 = 0.75;
        for (const op::EdgeCut cut : {op::EdgeCut::round(0.25), op::EdgeCut::flat(0.25, 0.5)}) {
            const auto st = op::blend_step_edge({0, 0, 0}, {A, 0, 0}, {0, Bh, 0}, {0, 0, c}, a1, b1, cut);
            AK_CHECK(st.solid.has_value());
            if (!st.solid) continue;
            const char* name = cut.kind == op::EdgeCut::Round ? "step, inside edge rounded" : "step, inside edge chamfered";
            check_faces(name, *st.solid);
            const double vol = c * step_section(A, Bh, a1, b1, cut).area;
            const op::Brep sb = op::to_brep(*st.solid);
            check_solid(name, sb, vol);
            // Halved across the blend by z <= 1/2.
            const auto res = op::boolean_op(sb, op::to_brep(op::make_box({-1, -1, -1}, {3, 3, 0.5})),
                                            op::Op::Intersection);
            AK_CHECK(res.ok);
            if (res.ok) check_solid("   halved by z <= 1/2", res.solid, vol / 2);
            else std::printf("FAIL halving: %s\n", res.error.c_str());
        }
        // A boss on a plate and a counterbore, revolved: torus and cone blends
        // that add material; full turns are closed chains (genus 1).
        struct Stepped {
            const char* name;
            bool tall_inside;
            op::EdgeCut cut;
            int quarters;
        };
        const double r0 = 1.0, r1 = 1.75, r2 = 2.5, z0 = 0.0, z1 = 0.5, z2 = 1.25;
        const Stepped cases[] = {
            {"boss on a plate, base rounded, full turn", true, op::EdgeCut::round(0.25), 4},
            {"boss on a plate, base chamfered, one quarter", true, op::EdgeCut::flat(0.25, 0.5), 1},
            {"counterbore, floor edge rounded, full turn", false, op::EdgeCut::round(0.25), 4},
            {"counterbore, floor edge chamfered, half", false, op::EdgeCut::flat(0.375, 0.5), 2},
        };
        for (const Stepped& k : cases) {
            const auto res = op::blend_stepped_ring_edge(r0, r1, r2, z0, z1, z2, k.tall_inside, k.cut, k.quarters);
            AK_CHECK(res.solid.has_value());
            if (!res.solid) {
                std::printf("FAIL %s: %s\n", k.name, res.error.c_str());
                continue;
            }
            const double a = r2 - r0, b = z2 - z0, a1 = k.tall_inside ? r1 - r0 : r2 - r1;
            const SectionMoments m = step_section(a, b, a1, z1 - z0, k.cut);
            const double rho_moment = k.tall_inside ? r0 * m.area + m.mx : r2 * m.area - m.mx;
            check_faces(k.name, *res.solid);
            check_solid(k.name, op::to_brep(*res.solid), k.quarters * (std::numbers::pi / 2) * rho_moment,
                        k.quarters == 4 ? 0 : 2);
        }
    }

    // --- all edges rounded: three fillets meet at every corner ----------------
    {
        const double A = 2.0, Bh = 1.5, Ch = 1.0, rr = 0.25;
        const auto rb = op::rounded_box({0, 0, 0}, A, Bh, Ch, rr);
        AK_CHECK(rb.solid.has_value());
        if (rb.solid) {
            // The corner patches lie on the unit sphere exactly in exact
            // arithmetic; in doubles to rounding.
            double worst = 0.0;
            for (const auto& P : op::fillet_detail::octant_patches()) {
                const geom::BezierSurface S(P);
                for (int i = 0; i <= 20; ++i)
                    for (int j = 0; j <= 20; ++j) {
                        const Vec3d p = S.point(i / 20.0, j / 20.0);
                        worst = std::max(worst, std::fabs(geom::norm(p) - 1.0));
                        AK_CHECK(p.x > -1e-15 && p.y > -1e-15 && p.z > -1e-15);
                    }
            }
            std::printf("sphere octant patches: off the unit sphere by %.2e\n", worst);
            AK_CHECK(worst < 1e-14);
            check_faces("rounded box", *rb.solid);
            // Core box plus face slabs, quarter cylinders and the ball.
            const double a0 = A - 2 * rr, b0 = Bh - 2 * rr, c0 = Ch - 2 * rr;
            const double vol = a0 * b0 * c0 + 2 * rr * (a0 * b0 + b0 * c0 + a0 * c0) +
                               std::numbers::pi * rr * rr * (a0 + b0 + c0) + 4.0 / 3.0 * std::numbers::pi * rr * rr * rr;
            const op::Brep rbb = op::to_brep(*rb.solid);
            check_solid("rounded box (all 12 edges)", rbb, vol);
            // Cut through four corners: keep x <= a - r/2. The removed slab
            // is the half-ball sweep of the core's x face from depth r/2 to r.
            const double cut_vol = b0 * c0 * rr / 2 +
                                   2 * (b0 + c0) * (std::numbers::pi * rr * rr / 6 - std::sqrt(3.0) * rr * rr / 8) +
                                   std::numbers::pi * 5 * rr * rr * rr / 24;
            const auto res = op::boolean_op(rbb, op::to_brep(op::make_box({-1, -1, -1}, {A - rr / 2, 3, 3})),
                                            op::Op::Intersection);
            AK_CHECK(res.ok);
            if (res.ok) check_solid("   cut x <= a - r/2 through the corners", res.solid, vol - cut_vol);
            else std::printf("FAIL rounded box cut: %s\n", res.error.c_str());
        }
        AK_CHECK(!op::rounded_box({0, 0, 0}, 1, 1, 1, 0.5).solid);
    }

    // Refusals.
    AK_CHECK(!op::fillet_box_edge({0, 0, 0}, {a, 0, 0}, {0, b, 0}, {0, 0, c}, 1.2).solid); // r > 2b/3
    AK_CHECK(!op::fillet_box_edge({0, 0, 0}, {a, 0, 0}, {0, b, 0}, {0, 0, c}, 0.0).solid);
    AK_CHECK(!op::fillet_box_edge({0, 0, 0}, {0, b, 0}, {a, 0, 0}, {0, 0, c}, r).solid); // left-handed
    AK_CHECK(!op::blend_ring_edge(0.0, 1.0, 0.0, 1.0, op::RingEdge::OuterTop, op::EdgeCut::round(0.25)).solid);
    AK_CHECK(!op::blend_box_edge({0, 0, 0}, {a, 0, 0}, {0, b, 0}, {0, 0, c}, op::EdgeCut{op::EdgeCut::Round, 0.25, 0.5}).solid);

    return akt::finish("fillet_test");
}
