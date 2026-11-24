// Procedural carriers (geom/carrier.hpp) and the certified spine of a
// rolling-ball blend (solve/carrier_ssi.hpp). Enclosures must contain every
// sampled value (soundness) and shrink with the box (convergence); the
// offset of a cylinder is a cylinder; spines traced between a cylinder and
// a plane, and between a freeform hill and a wall, are certified and match
// the rolling-ball conditions.
#include <akernel/geom/carrier.hpp>
#include <akernel/geom/primitives.hpp>
#include <akernel/solve/carrier_ssi.hpp>

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

    return akt::finish("carrier_test");
}
