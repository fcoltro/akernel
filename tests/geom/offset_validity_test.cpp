// Certified offset validity and curvature bounds (geom/offset_validity.hpp).
// Soundness: a Valid verdict is never contradicted by sampling, an Invalid
// witness really violates 1 - d kappa_i > 0, and nothing on the focal
// boundary is called Valid. Tightness: verdicts flip at the analytic radius
// on cylinders, tori and saddles, and the curvature bound converges to the
// sampled maximum.
#include <akernel/geom/continuity.hpp>
#include <akernel/geom/offset_validity.hpp>
#include <akernel/geom/primitives.hpp>
#include <akernel/geom/sweep.hpp>

#include <cmath>
#include <vector>

#include "check.hpp"

using namespace ak::geom;

namespace {

/// Principal curvatures at (u, v) in doubles, sampled.
void kappas(const BezierSurface& s, double u, double v, double& k1, double& k2) {
    const Jet3 J = s.jet(u, v);
    const SurfaceFrame f = surface_frame(J);
    const double H = mean_curvature(f), K = gauss_curvature(f);
    const double r = std::sqrt(std::max(0.0, H * H - K));
    k1 = H - r;
    k2 = H + r;
}

double sampled_max_kappa(const PolyPatchD& p, int n) {
    const BezierSurface s(p);
    double m = 0.0;
    for (int i = 0; i <= n; ++i)
        for (int j = 0; j <= n; ++j) {
            double k1, k2;
            kappas(s, static_cast<double>(i) / n, static_cast<double>(j) / n, k1, k2);
            m = std::max({m, std::fabs(k1), std::fabs(k2)});
        }
    return m;
}

/// Smallest sampled value of min(1 - d k1, 1 - d k2): positive where the
/// offset is regular at every sample.
double sampled_margin(const PolyPatchD& p, double d, int n) {
    const BezierSurface s(p);
    double m = 1e300;
    for (int i = 0; i <= n; ++i)
        for (int j = 0; j <= n; ++j) {
            double k1, k2;
            kappas(s, static_cast<double>(i) / n, static_cast<double>(j) / n, k1, k2);
            m = std::min({m, 1 - d * k1, 1 - d * k2});
        }
    return m;
}

} // namespace

int main() {
    akt::Rng rng(23);

    // --- cylinder of radius 2: outward normal, kappa = -1/2 ------------------
    {
        const PolyPatchD cyl = make_cylinder_quarter(2.0, 0.0, 1.0);
        AK_CHECK(check_offset(cyl, -1.9).verdict == OffsetVerdict::Valid);
        AK_CHECK(check_offset(cyl, -1.99).verdict == OffsetVerdict::Valid);
        AK_CHECK(check_offset(cyl, 100.0).verdict == OffsetVerdict::Valid);
        AK_CHECK(check_offset(cyl, -2.1).verdict == OffsetVerdict::Invalid);
        // Exactly on the axis: the offset collapses. Never called valid.
        AK_CHECK(check_offset(cyl, -2.0, 8).verdict != OffsetVerdict::Valid);
        // The bound converges to 1/2 from above as the grid is refined.
        const double kb4 = curvature_bound(cyl, 4), kb6 = curvature_bound(cyl, 6);
        AK_CHECK(kb4 >= 0.5 && kb6 >= 0.5 && kb6 <= kb4);
        AK_CHECK(kb6 < 0.5 * 1.05);
    }

    // --- torus (R = 3, r = 1) from revolution -----------------------------------
    // One side degenerates at d = 1 (the tube's centre circle), the other at
    // d = 2 (the axis, reached first at the inner equator).
    {
        NurbsCurve prof = exact_circle(1.0, 0.0, 4);
        for (auto& c : prof.cp) c = {c.x + 3.0 * c.w, 0.0, c.y, c.w};
        // The four quarter turns are congruent by exact rotations, so the
        // pieces of the first quarter (the whole profile) cover the torus.
        std::vector<NurbsSurface::BezierPiece> pieces;
        for (const auto& pc : revolve(prof, 4).to_bezier())
            if (pc.u0 == 0.0) pieces.push_back(pc);
        AK_CHECK(pieces.size() == 4u);
        // Orientation: is the normal of piece 0 pointing away from the tube centre?
        const BezierSurface s0(pieces[0].patch);
        const Jet3 J = s0.jet(0.5, 0.5);
        const Vec3d N = cross(J[1], J[2]), P = J[0];
        const Vec3d centre = Vec3d{P.x, P.y, 0.0} * (3.0 / std::hypot(P.x, P.y));
        const double out = dot(N, P - centre) > 0 ? 1.0 : -1.0;
        auto all = [&](double d, OffsetVerdict want) {
            int hits = 0;
            for (const auto& pc : pieces) hits += check_offset(pc.patch, d).verdict == want;
            return hits;
        };
        auto any = [&](double d, OffsetVerdict want) {
            for (const auto& pc : pieces)
                if (check_offset(pc.patch, d).verdict == want) return true;
            return false;
        };
        // Towards the tube centre: valid below 1, invalid beyond.
        AK_CHECK(all(-0.9 * out, OffsetVerdict::Valid) == 4);
        AK_CHECK(any(-1.1 * out, OffsetVerdict::Invalid));
        // Away from the tube centre: valid below 2, invalid beyond.
        AK_CHECK(all(1.9 * out, OffsetVerdict::Valid) == 4);
        AK_CHECK(all(1.1 * out, OffsetVerdict::Valid) == 4);
        AK_CHECK(any(2.1 * out, OffsetVerdict::Invalid));
    }

    // --- the hill z = 16 h u(1-u) v(1-v): saddles at the corners ---------------
    {
        const double h = 0.25;
        const PolyPatchD hill = make_bump(h);
        const double sampled = sampled_max_kappa(hill, 256);
        const double kb4 = curvature_bound(hill, 4), kb6 = curvature_bound(hill, 6);
        std::printf("hill: sampled max |kappa| %.6f, bound %.6f (16x16), %.6f (64x64)\n", sampled,
                    kb4, kb6);
        AK_CHECK(kb4 >= sampled && kb6 >= sampled);
        AK_CHECK(kb6 <= kb4);
        AK_CHECK(kb6 < sampled * 1.10);
        // Offsets below 1/bound are certified valid on both sides.
        AK_CHECK(check_offset(hill, 0.95 / kb6).verdict == OffsetVerdict::Valid);
        AK_CHECK(check_offset(hill, -0.95 / kb6).verdict == OffsetVerdict::Valid);
        // Beyond the sampled maximum, invalid on both sides (saddle corners).
        AK_CHECK(check_offset(hill, 1.1 / sampled).verdict == OffsetVerdict::Invalid);
        AK_CHECK(check_offset(hill, -1.1 / sampled).verdict == OffsetVerdict::Invalid);
    }

    // --- random rational patches: verdicts never contradicted ------------------
    int valid = 0, invalid = 0, undecided = 0;
    for (int trial = 0; trial < 40; ++trial) {
        PolyPatchD p(2, 2);
        for (int i = 0; i <= 2; ++i)
            for (int j = 0; j <= 2; ++j) {
                const double w = trial % 2 ? rng.range(0.6, 1.6) : 1.0;
                const double z = rng.range(-0.3, 0.3);
                p.at(i, j) = {w * 0.5 * i, w * 0.5 * j, w * z, w};
            }
        const double d = rng.range(-1.5, 1.5);
        const OffsetCheck c = check_offset(p, d, 10, 20000);
        const double margin = sampled_margin(p, d, 64);
        if (c.verdict == OffsetVerdict::Valid) {
            ++valid;
            AK_CHECK(margin > 0.0);
        } else if (c.verdict == OffsetVerdict::Invalid) {
            ++invalid;
            double k1, k2;
            kappas(BezierSurface(p), c.u, c.v, k1, k2);
            AK_CHECK(!(1 - d * k1 > 1e-9 && 1 - d * k2 > 1e-9));
        } else {
            // Undecided only close to the focal boundary, where the margin is small.
            ++undecided;
            AK_CHECK(margin < 0.1);
        }
    }
    std::printf("random patches: %d valid, %d invalid, %d undecided\n", valid, invalid, undecided);
    AK_CHECK(valid > 0 && invalid > 0);
    AK_CHECK(undecided * 4 <= valid + invalid);

    return akt::finish("offset_validity_test");
}
