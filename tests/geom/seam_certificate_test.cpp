// Certified seam bounds (geom/seam_certificate.hpp): every certified bound
// must be at least the jump measured at any sample point (soundness), close
// to it where the jump is known (tightness), and tiny where the patches
// agree to third order. The variational fill is then certified against a
// polynomial neighbour.
#include <akernel/geom/continuity.hpp>
#include <akernel/geom/nurbs_surface.hpp>
#include <akernel/geom/seam_certificate.hpp>
#include <akernel/var/fill.hpp>

#include <cmath>
#include <cstdio>
#include <vector>

#include "check.hpp"

using namespace ak::geom;

namespace {

double binom(int n, int k) {
    double c = 1.0;
    for (int i = 1; i <= k; ++i) c = c * (n - k + i) / i;
    return c;
}

/// The graph z = sum c[a][b] x^a y^b over [0,1]^2 as a Bezier patch of
/// degree (n, n), with x = u, y = v.
PolyPatchD graph(int n, const std::vector<std::vector<double>>& c) {
    PolyPatchD P(n, n);
    for (int i = 0; i <= n; ++i)
        for (int j = 0; j <= n; ++j) {
            double z = 0.0;
            for (int a = 0; a <= i && a < static_cast<int>(c.size()); ++a)
                for (int b = 0; b <= j && b < static_cast<int>(c[static_cast<std::size_t>(a)].size()); ++b)
                    z += binom(i, a) * binom(j, b) / (binom(n, a) * binom(n, b)) * c[static_cast<std::size_t>(a)][static_cast<std::size_t>(b)];
            P.at(i, j) = {static_cast<double>(i) / n, static_cast<double>(j) / n, z, 1.0};
        }
    return P;
}

struct Measured {
    double gap = 0, angle = 0, db[3] = {0, 0, 0}, dC[4] = {0, 0, 0, 0};
};

/// Jumps at 401 points of the seam, in doubles (geom/continuity.hpp).
Measured sample(const PolyPatchD& A, const PolyPatchD& B, int side) {
    const BezierSurface SA(A), SB(B);
    Measured m;
    for (int i = 0; i <= 400; ++i) {
        const double t = i / 400.0;
        const double u = side == kSeamV0 || side == kSeamV1 ? t : (side == kSeamU0 ? 0.0 : 1.0);
        const double v = side == kSeamU0 || side == kSeamU1 ? t : (side == kSeamV0 ? 0.0 : 1.0);
        const Jet3 JA = SA.jet(u, v), JB = SB.jet(u, v);
        const SeamJump j = seam_jump(JA, JB, side == kSeamV0 || side == kSeamV1 ? JA[1] : JA[2]);
        m.gap = std::max(m.gap, norm(JB[0] - JA[0]));
        m.angle = std::max(m.angle, j.normal_angle);
        const double db[3] = {j.db_tt, j.db_tn, j.db_nn}, dC[4] = {j.dC_ttt, j.dC_ttn, j.dC_tnn, j.dC_nnn};
        for (int k = 0; k < 3; ++k) m.db[k] = std::max(m.db[k], std::fabs(db[k]));
        for (int k = 0; k < 4; ++k) m.dC[k] = std::max(m.dC[k], std::fabs(dC[k]));
    }
    return m;
}

/// Every certified bound covers the sampled jump (up to the rounding of the
/// sampling itself, `slack`). Returns the largest ratio bound / sample over
/// the jumps that are clearly nonzero.
double sound(const SeamCertificate& c, const Measured& m, double slack, bool print, const char* name) {
    const double cb[9] = {c.gap, c.normal_angle, c.db_tt, c.db_tn, c.db_nn, c.dC_ttt, c.dC_ttn, c.dC_tnn, c.dC_nnn};
    const double mb[9] = {m.gap, m.angle, m.db[0], m.db[1], m.db[2], m.dC[0], m.dC[1], m.dC[2], m.dC[3]};
    double worst = 0.0;
    for (int k = 0; k < 9; ++k) {
        AK_CHECK(mb[k] <= cb[k] + slack);
        if (mb[k] > 1e3 * slack) worst = std::max(worst, cb[k] / mb[k]);
    }
    if (print)
        std::printf("  %s: certified angle %.3e (sampled %.3e), db_nn %.3e (%.3e), dC_nnn %.3e (%.3e); worst ratio %.3f\n",
                    name, c.normal_angle, m.angle, c.db_nn, m.db[2], c.dC_nnn, m.dC[3], worst);
    return worst;
}

const std::vector<std::vector<double>> kBase = {
    {0.05, -0.2, 0.2, -0.1, 0.03}, {0.1, -0.3, 0.4, 0.02}, {0.5, -0.2, 0.1}, {0.3, 0.05}, {-0.04}};

} // namespace

int main() {
    akt::Rng rng(20261115);

    // 1. Patches that agree to third order along v = 0 (rows 0..3 equal,
    //    rows 4..6 different): every jump vanishes, and the bounds are
    //    rounding-level.
    {
        const PolyPatchD A = graph(6, kBase);
        PolyPatchD B = A;
        for (int i = 0; i <= 6; ++i)
            for (int j = 4; j <= 6; ++j) B.at(i, j).z += 0.3 * std::sin(1.0 + i + 2.0 * j);
        const SeamCertificate c = certify_seam(A, B, kSeamV0, 8);
        AK_CHECK(c.ok);
        const double all[9] = {c.gap, c.normal_angle, c.db_tt, c.db_tn, c.db_nn, c.dC_ttt, c.dC_ttn, c.dC_tnn, c.dC_nnn};
        double worst = 0.0;
        for (double x : all) worst = std::max(worst, x);
        std::printf("  third-order contact: largest certified jump %.2e\n", worst);
        AK_CHECK(worst < 1e-10);
        // The same seam seen from the other side of B (v = 1 of A against B
        // continued below it) is a different seam: not checked here.
    }

    // 2. A known jump of each order: row j of B differs from A's by delta in
    //    z, so the join is C^(j-1) and jumps first in order j.
    {
        const PolyPatchD A = graph(6, kBase);
        const char* names[3] = {"G1 jump (row 1 moved)", "G2 jump (row 2 moved)", "G3 jump (row 3 moved)"};
        for (int row = 1; row <= 3; ++row) {
            PolyPatchD B = A;
            for (int i = 0; i <= 6; ++i) B.at(i, row).z += 1e-4 * (1.0 + 0.5 * std::cos(1.3 * i));
            const SeamCertificate c = certify_seam(A, B, kSeamV0, 64);
            AK_CHECK(c.ok);
            const Measured m = sample(A, B, kSeamV0);
            const double ratio = sound(c, m, 1e-12, true, names[row - 1]);
            // The jump that is first nonzero is bounded tightly.
            const double first_c = row == 1 ? c.normal_angle : row == 2 ? c.db_nn : c.dC_nnn;
            const double first_m = row == 1 ? m.angle : row == 2 ? m.db[2] : m.dC[3];
            AK_CHECK(first_c <= 1.25 * first_m);
            // Lower orders stay at rounding level.
            if (row >= 2) { AK_CHECK(c.normal_angle < 1e-11); }
            if (row >= 3) { AK_CHECK(c.db_nn < 1e-10 && c.db_tt < 1e-10 && c.db_tn < 1e-10); }
            AK_CHECK(ratio < 50.0);
        }
    }

    // 3. Soundness on random small differences in the first four rows, on
    //    all four sides.
    {
        int cases = 0;
        double worst_ratio = 0.0;
        for (int trial = 0; trial < 24; ++trial) {
            const PolyPatchD A = graph(5, kBase);
            PolyPatchD B = A;
            const double scale = std::pow(10.0, -2.0 - 4.0 * rng.unit());
            for (auto& cp : B.c) {
                cp.x += scale * rng.range(-1, 1);
                cp.y += scale * rng.range(-1, 1);
                cp.z += scale * rng.range(-1, 1);
            }
            const int side = trial % 4;
            const SeamCertificate c = certify_seam(A, B, side, 8);
            AK_CHECK(c.ok);
            if (!c.ok) continue;
            const Measured m = sample(A, B, side);
            worst_ratio = std::max(worst_ratio, sound(c, m, 1e-12, false, ""));
            ++cases;
        }
        std::printf("  random differences: %d seams, every bound covers its samples; worst bound / sample %.2f\n", cases,
                    worst_ratio);
    }

    // 4. Refusals: a rational patch, and a seam where the normals disagree in
    //    orientation (B's u reversed).
    {
        const PolyPatchD A = graph(4, kBase);
        PolyPatchD R = A;
        R.at(2, 2) = R.at(2, 2) * 1.5; // weight 1.5
        AK_CHECK(!certify_seam(A, R, kSeamV0).ok);
        PolyPatchD F = A;
        for (int i = 0; i <= 4; ++i)
            for (int j = 0; j <= 4; ++j) F.at(i, j).x = 1.0 - A.at(i, j).x; // mirrored: normal flips
        const SeamCertificate c = certify_seam(A, F, kSeamV0);
        AK_CHECK(!c.ok);
        std::printf("  refused: rational patch; mirrored neighbour (%s)\n", c.why.c_str());
    }

    // 5. The G3 fill of a hole in a degree-9 polynomial surface, certified on
    //    every boundary piece of the fill against the surface itself. The
    //    fill reproduces only cubics, so the jumps are not zero; they must
    //    fall under refinement, and the certificates must cover the samples.
    {
        const std::vector<std::vector<double>> c9 = {
            {0.05, -0.2, 0.2, -0.1, 0.0, 0.0, 0.0, 0.0, 0.0, -0.25},
            {0.1, -0.3, 0.4},
            {0.5, -0.2},
            {0.3},
            {0.0, 0.0, 0.0, 0.0, 0.0, 0.4},
            {0.0},
            {0.0},
            {0.0, 0.0, -0.2},
            {0.0},
            {0.3}};
        const PolyPatchD Q = graph(9, c9);
        const BezierSurface QS(Q);
        const ak::var::JetFn jets = [&QS](int side, int j, double s) {
            const bool u_side = side == ak::var::U0 || side == ak::var::U1;
            const double u = u_side ? (side == ak::var::U1 ? 1.0 : 0.0) : s;
            const double v = u_side ? s : (side == ak::var::V1 ? 1.0 : 0.0);
            return u_side ? QS.jet(u, v)[jet_index(j, 0)] : QS.jet(u, v)[jet_index(0, j)];
        };
        SeamCertificate prev;
        bool first = true;
        for (int n : {12, 16, 24}) {
            ak::var::FillSpec spec;
            spec.ncp = n;
            const ak::var::FillResult r = ak::var::solve_fill(spec, jets);
            AK_CHECK(r.ok);
            SeamCertificate all;
            Measured sm;
            int seams = 0;
            bool all_ok = true;
            double worst_ratio = 0.0;
            for (const auto& piece : r.surface.to_bezier()) {
                const int sides[4] = {kSeamU0, kSeamU1, kSeamV0, kSeamV1};
                const bool on[4] = {piece.u0 == 0.0, piece.u1 == 1.0, piece.v0 == 0.0, piece.v1 == 1.0};
                for (int s = 0; s < 4; ++s) {
                    if (!on[s]) continue;
                    const PolyPatchI Qr = restricted_polynomial(Q, piece.u0, piece.u1, piece.v0, piece.v1);
                    // The fill is a polynomial B-spline, so its Bezier pieces are
                    // polynomial; knot insertion in doubles leaves their weights
                    // within rounding of 1. The carrier is the polynomial piece.
                    PolyPatchD F = piece.patch;
                    for (auto& cp : F.c) cp.w = 1.0;
                    const SeamCertificate c = certify_seam(Qr, to_interval(F), sides[s], 32);
                    all_ok = all_ok && c.ok;
                    if (!c.ok) {
                        std::printf("  fill seam not certified: %s\n", c.why.c_str());
                        continue;
                    }
                    ++seams;
                    // Samples of the same pair (the restricted neighbour in doubles).
                    PolyPatchD Qd(Q.p, Q.q);
                    for (std::size_t k = 0; k < Qd.c.size(); ++k)
                        Qd.c[k] = {ak::num::mid(Qr.c[k].x), ak::num::mid(Qr.c[k].y), ak::num::mid(Qr.c[k].z), 1.0};
                    const Measured m = sample(Qd, F, sides[s]);
                    worst_ratio = std::max(worst_ratio, sound(c, m, 1e-12, false, ""));
                    all.gap = std::max(all.gap, c.gap);
                    all.normal_angle = std::max(all.normal_angle, c.normal_angle);
                    all.db_nn = std::max(all.db_nn, c.db_nn);
                    all.dC_nnn = std::max(all.dC_nnn, c.dC_nnn);
                    sm.angle = std::max(sm.angle, m.angle);
                    sm.db[2] = std::max(sm.db[2], m.db[2]);
                    sm.dC[3] = std::max(sm.dC[3], m.dC[3]);
                }
            }
            AK_CHECK(all_ok);
            std::printf("  fill %2d x %2d, %d boundary seams: certified G1 %.2e rad (sampled %.2e), G2 %.2e (%.2e), G3 %.2e "
                        "(%.2e), gap %.2e; worst bound / sample %.2f\n",
                        n, n, seams, all.normal_angle, sm.angle, all.db_nn, sm.db[2], all.dC_nnn, sm.dC[3], all.gap,
                        worst_ratio);
            if (!first) { AK_CHECK(all.normal_angle < prev.normal_angle && all.db_nn < prev.db_nn && all.dC_nnn < prev.dC_nnn); }
            prev = all;
            first = false;
        }
    }

    return akt::finish("seam_certificate_test");
}
