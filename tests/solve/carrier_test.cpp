// Procedural carriers (geom/carrier.hpp) and the certified spine of a
// rolling-ball blend (solve/carrier_ssi.hpp). Enclosures must contain every
// sampled value (soundness) and shrink with the box (convergence); the
// offset of a cylinder is a cylinder; spines traced between a cylinder and
// a plane, and between a freeform hill and a wall, are certified and match
// the rolling-ball conditions.
#include <akernel/geom/carrier.hpp>
#include <akernel/geom/primitives.hpp>
#include <akernel/solve/carrier_ssi.hpp>

#include <chrono>
#include <cmath>
#include <cstdio>

#include "check.hpp"

using namespace ak::geom;
using ak::num::Interval;
using ak::solve::Pt4;

namespace {

bool inside(const IVec3& box, const Vec3d& p, double slack) {
    auto in = [&](const Interval& i, double x) { return i.lo - slack <= x && x <= i.hi + slack; };
    return in(box.x, p.x) && in(box.y, p.y) && in(box.z, p.z);
}
double width(const IVec3& b) { return std::max({b.x.hi - b.x.lo, b.y.hi - b.y.lo, b.z.hi - b.z.lo}); }

/// Soundness and convergence of a carrier's enclosures on random boxes.
void check_enclosures(const Carrier& c, akt::Rng& rng, double& w_coarse, double& w_fine) {
    w_coarse = w_fine = 0.0;
    for (int trial = 0; trial < 30; ++trial) {
        for (double h : {1.0 / 8, 1.0 / 64}) {
            const double u0 = rng.range(0, 1 - h), v0 = rng.range(0, 1 - h);
            IVec3 P, Pu, Pv;
            AK_CHECK(c.jet1_over(u0, u0 + h, v0, v0 + h, P, Pu, Pv));
            for (int n = 0; n < 10; ++n) {
                const double u = u0 + h * rng.unit(), v = v0 + h * rng.unit();
                Vec3d p, pu, pv;
                c.eval(u, v, p, pu, pv);
                AK_CHECK(inside(P, p, 1e-12) && inside(Pu, pu, 1e-11) && inside(Pv, pv, 1e-11));
                AK_CHECK(inside(c.point_at(u, v), p, 1e-12));
            }
            (h > 0.1 ? w_coarse : w_fine) = std::max(h > 0.1 ? w_coarse : w_fine, width(Pu));
        }
    }
}

} // namespace

int main() {
    akt::Rng rng(41);

    // --- the offset of a cylinder is a cylinder -------------------------------
    {
        const PolyPatchD cyl = make_cylinder_quarter(2.0, 0.0, 1.0);
        const OffsetCarrier off(cyl, 0.5);
        for (int n = 0; n < 200; ++n) {
            const IVec3 p = off.point_at(rng.unit(), rng.unit());
            const Interval r2 = ak::num::sqr(p.x) + ak::num::sqr(p.y);
            AK_CHECK(ak::num::contains(r2, 6.25)); // proven to hold the exact value
            AK_CHECK(r2.hi - r2.lo < 1e-12);
        }
        double wc, wf;
        check_enclosures(off, rng, wc, wf);
        std::printf("offset cylinder: max width of the S_u enclosure %.3g (h = 1/8), %.3g (h = 1/64)\n", wc, wf);
        AK_CHECK(wf < wc / 4); // converges
        const PatchCarrier pc(cyl);
        check_enclosures(pc, rng, wc, wf);
    }

    // --- a rational freeform patch and its offset ------------------------------
    {
        PolyPatchD p(3, 2);
        for (int i = 0; i <= 3; ++i)
            for (int j = 0; j <= 2; ++j) {
                const double w = rng.range(0.7, 1.4);
                p.at(i, j) = {w * i / 3.0, w * j / 2.0, w * rng.range(-0.2, 0.2), w};
            }
        double wc, wf;
        check_enclosures(OffsetCarrier(p, 0.15), rng, wc, wf);
        AK_CHECK(wf < wc / 4);
    }

    // --- spine between a cylinder and a plane: two parallel lines --------------
    {
        const PolyPatchD cyl = make_cylinder_quarter(2.0, 0.0, 1.0); // normal outward
        const PolyPatchD wall = make_plane({1, 0, -0.5}, {0, 3, 0}, {0, 0, 2}); // x = 1, normal +x
        const double r = 0.25;
        const double yc = std::sqrt(2.25 * 2.25 - 1.25 * 1.25);
        const double th = std::atan2(yc, 1.25);
        const Pt4 x0{std::tan(0.5 * th) + 0.01, 0.5, yc / 3.0 - 0.01, 0.5};
        const auto sp = ak::solve::trace_spine(cyl, wall, r, x0);
        AK_CHECK(sp.has_value());
        if (sp) {
            AK_CHECK(sp->end_front == ak::solve::EndKind::Boundary);
            AK_CHECK(sp->end_back == ak::solve::EndKind::Boundary);
            const OffsetCarrier oc(cyl, r);
            double vmin = 1, vmax = 0, worst = 0;
            for (const auto& s : sp->segs)
                for (const Pt4& q : {s.a, s.b}) {
                    Vec3d c, cu, cv;
                    oc.eval(q[0], q[1], c, cu, cv);
                    worst = std::max({worst, std::fabs(c.x - 1.25), std::fabs(c.y - yc)});
                    vmin = std::min(vmin, q[1]);
                    vmax = std::max(vmax, q[1]);
                }
            std::printf("cylinder/plane spine: %zu certified segments, centre off the line by %.2e\n",
                        sp->segs.size(), worst);
            AK_CHECK(worst < 1e-12);
            AK_CHECK(vmin < 1e-12 && vmax > 1 - 1e-12); // runs the cylinder's full height
        }
    }

    // --- spine between a freeform hill and a wall ------------------------------
    {
        const PolyPatchD hill = make_bump(0.25); // normal up
        const PolyPatchD wall = make_plane({0.6, -0.5, -0.5}, {0, 0, 2}, {0, 2, 0}); // x = 0.6, normal -x
        const double r = 0.1;
        AK_CHECK(support_fits(hill, true, r).verdict == OffsetVerdict::Valid);
        const Pt4 x0{0.5, 0.5, 0.35, 0.5};
        const auto sp = ak::solve::trace_spine(hill, wall, r, x0);
        AK_CHECK(sp.has_value());
        if (sp) {
            AK_CHECK(sp->end_front == ak::solve::EndKind::Boundary);
            AK_CHECK(sp->end_back == ak::solve::EndKind::Boundary);
            const BezierSurface hb(hill);
            double worst = 0, vmin = 1, vmax = 0;
            int arcs = 0;
            for (const auto& s : sp->segs)
                for (const Pt4& q : {s.a, s.b}) {
                const auto sec = ak::solve::spine_section(hill, wall, r, q);
                AK_CHECK(sec.status == BlendStatus::Ok);
                arcs += sec.status == BlendStatus::Ok;
                // Ball centre at distance r from the wall and from the hill along its normal.
                worst = std::max({worst, std::fabs(sec.centre.x - 0.5),
                                  std::fabs(norm(sec.centre - hb.point(q[0], q[1])) - r),
                                  std::fabs(sec.contact_b.x - 0.6)});
                vmin = std::min(vmin, q[1]);
                vmax = std::max(vmax, q[1]);
                }
            std::printf("hill/wall spine: %zu certified segments, worst residual %.2e, worst sin angle %.3f\n",
                        sp->segs.size(), worst, sp->min_sin_angle);
            AK_CHECK(worst < 1e-12);
            AK_CHECK(arcs == 2 * static_cast<int>(sp->segs.size()));
            AK_CHECK(vmin < 1e-12 && vmax > 1 - 1e-12); // across the whole hill
        }
    }

    // --- spines found without a guess, with the completeness proof ------------
    {
        struct Case { const char* name; PolyPatchD a, b; double r; std::size_t curves; };
        const Case cases[] = {
            {"cylinder/plane", make_cylinder_quarter(2.0, 0.0, 1.0),
             make_plane({1, 0, -0.5}, {0, 3, 0}, {0, 0, 2}), 0.25, 1},
            {"hill/wall", make_bump(0.25), make_plane({0.6, -0.5, -0.5}, {0, 0, 2}, {0, 2, 0}), 0.1, 1},
            {"hill/far wall", make_bump(0.25), make_plane({3.0, -0.5, -0.5}, {0, 0, 2}, {0, 2, 0}), 0.1, 0},
        };
        for (const Case& c : cases) {
            const auto t0 = std::chrono::steady_clock::now();
            const auto res = ak::solve::find_spines(c.a, c.b, c.r);
            const double secs =
                std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
            std::printf("find_spines %s: %zu curve(s), complete %d, %zu verifier nodes, %.2f s\n", c.name,
                        res.curves.size(), static_cast<int>(res.complete), res.stats.verify_nodes, secs);
            AK_CHECK(res.complete);
            AK_CHECK(res.curves.size() == c.curves);
            for (const auto& cv : res.curves) {
                AK_CHECK(cv.end_front == ak::solve::EndKind::Boundary);
                AK_CHECK(cv.end_back == ak::solve::EndKind::Boundary);
            }
        }
    }

    // --- the blend carrier: a surface other code can intersect and mesh -----
    {
        // Cylinder of radius 2 and the plane x = 1, r = 0.25: the ball centre
        // runs along the line x = 1.25, y = yc, so the blend is a piece of the
        // cylinder of radius r about that line.
        const PolyPatchD cyl = make_cylinder_quarter(2.0, 0.0, 1.0);
        const PolyPatchD wall = make_plane({1, 0, -0.5}, {0, 3, 0}, {0, 0, 2});
        const double r = 0.25;
        const double yc = std::sqrt(2.25 * 2.25 - 1.25 * 1.25);
        const auto sp = ak::solve::find_spines(cyl, wall, r);
        AK_CHECK(sp.complete && sp.curves.size() == 1);
        if (sp.curves.size() == 1) {
            const ak::solve::BlendCarrier bl(cyl, wall, r, sp.curves[0].segs);
            // On the round tube about the spine line, in doubles and certified.
            double worst = 0.0;
            for (int n = 0; n < 300; ++n) {
                const double s = rng.unit(), t = rng.unit();
                Vec3d p, ps, pt;
                bl.eval(s, t, p, ps, pt);
                worst = std::max(worst, std::fabs(std::hypot(p.x - 1.25, p.y - yc) - r));
                const IVec3 P = bl.point_at(s, t);
                AK_CHECK(inside(P, p, 1e-12));
                const Interval d2 = ak::num::sqr(P.x - Interval(1.25)) + ak::num::sqr(P.y - Interval(yc));
                AK_CHECK(d2.lo <= r * r + 1e-9 && d2.hi >= r * r - 1e-9);
                AK_CHECK(width(P) < 1e-9);
            }
            std::printf("blend carrier (cylinder/plane): %zu spine segments, off the tube by %.2e\n",
                        bl.segments(), worst);
            AK_CHECK(worst < 1e-12);
            // Box enclosures hold the sampled points and derivatives (the
            // doubles' derivatives are central differences: slack 1e-6), and
            // shrink with the box.
            double wc = 0.0, wf = 0.0;
            for (int trial = 0; trial < 20; ++trial)
                for (double h : {1.0 / 8, 1.0 / 64}) {
                    const double s0 = rng.range(0, 1 - h), t0 = rng.range(0, 1 - h);
                    IVec3 P, Ps, Pt;
                    AK_CHECK(bl.jet1_over(s0, s0 + h, t0, t0 + h, P, Ps, Pt));
                    for (int n = 0; n < 10; ++n) {
                        Vec3d p, ps, pt;
                        bl.eval(s0 + h * rng.unit(), t0 + h * rng.unit(), p, ps, pt);
                        AK_CHECK(inside(P, p, 1e-12) && inside(Ps, ps, 1e-6) && inside(Pt, pt, 1e-6));
                    }
                    (h > 0.1 ? wc : wf) = std::max(h > 0.1 ? wc : wf, width(P));
                }
            std::printf("blend carrier: max width of the point enclosure %.3g (h = 1/8), %.3g (h = 1/64)\n",
                        wc, wf);
            AK_CHECK(wf < wc / 4);
            // G1 with the supports along the contact curves.
            for (int n = 0; n <= 10; ++n) {
                const double s = n / 10.0;
                Vec3d p, ps, pt;
                bl.eval(s, 0.0, p, ps, pt);
                const Vec3d nb = cross(ps, pt);
                AK_CHECK(std::fabs(std::hypot(p.x, p.y) - 2.0) < 1e-12);
                AK_CHECK(norm(cross(normalized(nb), normalized(Vec3d{p.x, p.y, 0}))) < 1e-6);
                bl.eval(s, 1.0, p, ps, pt);
                AK_CHECK(std::fabs(p.x - 1.0) < 1e-12);
                AK_CHECK(norm(cross(normalized(cross(ps, pt)), Vec3d{1, 0, 0})) < 1e-6);
            }
            // Intersected with the plane z = 0.5 like any other carrier: one
            // section arc, found without a guess, with the completeness proof.
            const PatchCarrier cut(make_plane({0, 0, 0.5}, {3, 0, 0}, {0, 3, 0}));
            const auto t0 = std::chrono::steady_clock::now();
            const auto res = ak::solve::intersect_carriers(bl, cut);
            const double secs =
                std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
            std::printf("blend x plane: %zu curve(s), complete %d, %zu verifier nodes, %.2f s\n",
                        res.curves.size(), static_cast<int>(res.complete), res.stats.verify_nodes, secs);
            AK_CHECK(res.complete);
            AK_CHECK(res.curves.size() == 1);
            double off = 0.0;
            for (const auto& cv : res.curves)
                for (const auto& sg : cv.segs)
                    for (const Pt4& q : {sg.a, sg.b}) {
                        Vec3d p, ps, pt;
                        bl.eval(q[0], q[1], p, ps, pt);
                        off = std::max({off, std::fabs(p.z - 0.5),
                                        std::fabs(std::hypot(p.x - 1.25, p.y - yc) - r)});
                    }
            AK_CHECK(off < 1e-10);
        }
    }

    // --- a freeform blend: hill and wall ---------------------------------------
    {
        const PolyPatchD hill = make_bump(0.25);
        const PolyPatchD wall = make_plane({0.6, -0.5, -0.5}, {0, 0, 2}, {0, 2, 0});
        const double r = 0.1;
        const auto sp = ak::solve::find_spines(hill, wall, r);
        AK_CHECK(sp.complete && sp.curves.size() == 1);
        if (sp.curves.size() == 1) {
            const ak::solve::BlendCarrier bl(hill, wall, r, sp.curves[0].segs);
            const BezierSurface hb(hill);
            for (int n = 0; n < 100; ++n) {
                const double s = rng.unit(), t = rng.unit();
                Vec3d p, ps, pt;
                bl.eval(s, t, p, ps, pt);
                AK_CHECK(inside(bl.point_at(s, t), p, 1e-12));
                // Every blend point is at distance r from the ball centre.
                Pt4 x;
                AK_CHECK(bl.spine_at(s, x));
                const auto sec = ak::solve::spine_section(hill, wall, r, x);
                AK_CHECK_NEAR(norm(p - sec.centre), r, 1e-12);
            }
            for (int trial = 0; trial < 20; ++trial) {
                const double h = 1.0 / 32, s0 = rng.range(0, 1 - h), t0 = rng.range(0, 1 - h);
                IVec3 P, Ps, Pt;
                AK_CHECK(bl.jet1_over(s0, s0 + h, t0, t0 + h, P, Ps, Pt));
                for (int n = 0; n < 10; ++n) {
                    Vec3d p, ps, pt;
                    bl.eval(s0 + h * rng.unit(), t0 + h * rng.unit(), p, ps, pt);
                    AK_CHECK(inside(P, p, 1e-12) && inside(Ps, ps, 1e-6) && inside(Pt, pt, 1e-6));
                }
            }
            // G1: tangent to the wall at theta = 1, to the hill at theta = 0.
            for (int n = 0; n <= 10; ++n) {
                Vec3d p, ps, pt;
                Pt4 x;
                bl.spine_at(n / 10.0, x);
                bl.eval(n / 10.0, 0.0, p, ps, pt);
                const Jet3 J = BezierSurface(hill).jet(x[0], x[1]);
                AK_CHECK(norm(cross(normalized(cross(ps, pt)), normalized(cross(J[1], J[2])))) < 1e-6);
                bl.eval(n / 10.0, 1.0, p, ps, pt);
                AK_CHECK(norm(cross(normalized(cross(ps, pt)), Vec3d{1, 0, 0})) < 1e-6);
            }
            std::printf("blend carrier (hill/wall): %zu spine segments, sampled and enclosed\n",
                        bl.segments());
        }
    }

    return akt::finish("carrier_test");
}
