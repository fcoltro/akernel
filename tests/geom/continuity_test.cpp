// Differential invariants and the G2 / G3 seam reductions of blueprint §1.7.
//
// The claim under test: across a seam with unit tangent tau and conormal nu,
//   - given G1, the only component of the second fundamental form b that can
//     jump is b(nu, nu);
//   - given G2, the only component of C = nabla b that can jump is C(nu,nu,nu).
// The second statement was derived in the blueprint from Codazzi symmetry and
// flagged as needing verification. It is checked here on surfaces that are
// reparameterised nonlinearly across the seam, so nothing is aligned.
#include <akernel/geom/continuity.hpp>
#include <akernel/geom/primitives.hpp>

#include "check.hpp"

using namespace ak::geom;

namespace {

double f(double x, double y) {
    return 0.3 * x * x - 0.2 * x * y + 0.25 * y * y + 0.1 * x * x * x - 0.15 * x * x * y +
           0.05 * y * y * y + 0.2 * x + 0.1 * y;
}

// Side A: the graph of f over x in [-1, 0]; the seam is u = 1 (x = 0).
BezierSurface side_a() {
    return BezierSurface(fit_patch([](double u, double v) { return Vec3d{u - 1, v, f(u - 1, v)}; }, 3, 3));
}

// Side B: the graph of f + extra(x, y) over x >= 0, in a nonlinear
// parameterisation that keeps the seam at u = 0 with y = v there.
BezierSurface side_b(const std::function<double(double, double)>& extra) {
    return BezierSurface(fit_patch(
        [extra](double u, double v) {
            const double x = u * (1 + 0.3 * v + 0.2 * u);
            const double y = v + 0.25 * u - 0.1 * u * u;
            return Vec3d{x, y, f(x, y) + extra(x, y)};
        },
        8, 4));
}

} // namespace

int main() {
    akt::Rng rng(3);

    // Codazzi: nabla b is totally symmetric on any surface.
    {
        PolyPatchD p(3, 3);
        for (auto& c : p.c) {
            const double w = rng.range(0.7, 1.5);
            c = {w * rng.range(-1, 1), w * rng.range(-1, 1), w * rng.range(-1, 1), w};
        }
        // Make it regular: add a dominant affine part.
        for (int i = 0; i <= 3; ++i)
            for (int j = 0; j <= 3; ++j) {
                auto& c = p.at(i, j);
                c.x += c.w * 3.0 * i;
                c.y += c.w * 3.0 * j;
            }
        // The closed formula used for C is symmetric by construction, so
        // comparing its own components proves nothing. Check it against the
        // definition instead: (nabla_c b)_ab = d_c b_ab - G^d_ca b_db - G^d_cb b_ad,
        // with d_c b_ab from central differences of b.
        const BezierSurface s(p);
        const double h = 1e-5;
        double worst = 0.0;
        for (int n = 0; n < 100; ++n) {
            const double u = rng.range(0.1, 0.9), v = rng.range(0.1, 0.9);
            const Jet3 J = s.jet(u, v);
            const SurfaceFrame fr = surface_frame(J);
            const SurfaceFrame fp[2] = {surface_frame(s.jet(u + h, v)), surface_frame(s.jet(u, v + h))};
            const SurfaceFrame fm[2] = {surface_frame(s.jet(u - h, v)), surface_frame(s.jet(u, v - h))};
            const Vec3d S1[2] = {J[1], J[2]};
            double G[2][2][2];
            for (int d = 0; d < 2; ++d)
                for (int a = 0; a < 2; ++a)
                    for (int c = 0; c < 2; ++c)
                        G[d][a][c] = fr.ginv[d][0] * dot(S1[0], J[3 + a + c]) +
                                     fr.ginv[d][1] * dot(S1[1], J[3 + a + c]);
            for (int a = 0; a < 2; ++a)
                for (int b = 0; b < 2; ++b)
                    for (int c = 0; c < 2; ++c) {
                        double ref = (fp[c].b[a][b] - fm[c].b[a][b]) / (2 * h);
                        for (int d = 0; d < 2; ++d)
                            ref -= G[d][c][a] * fr.b[d][b] + G[d][c][b] * fr.b[a][d];
                        worst = std::max(worst, std::fabs(fr.C[a][b][c] - ref));
                    }
        }
        std::printf("  nabla b: closed formula vs finite-difference definition, max diff %.2e\n", worst);
        AK_CHECK(worst < 1e-6);
    }

    // Cylinder of radius r: K = 0, |H| = 1/(2r), nabla b = 0.
    {
        const double r = 1.5;
        const BezierSurface cyl(make_cylinder_quarter(r, 0.0, 2.0));
        for (int n = 0; n < 100; ++n) {
            const SurfaceFrame fr = surface_frame(cyl.jet(rng.unit(), rng.unit()));
            AK_CHECK_NEAR(gauss_curvature(fr), 0.0, 1e-12);
            AK_CHECK_NEAR(std::fabs(mean_curvature(fr)), 0.5 / r, 1e-12);
            for (int a = 0; a < 2; ++a)
                for (int b = 0; b < 2; ++b)
                    for (int c = 0; c < 2; ++c) AK_CHECK_NEAR(fr.C[a][b][c], 0.0, 1e-10);
        }
    }

    const BezierSurface A = side_a();
    const double tol = 1e-7;

    // Case 1: the same surface on both sides, reparameterised. Every jump is
    // zero although the parametric third derivatives differ.
    {
        const BezierSurface B = side_b([](double, double) { return 0.0; });
        double worst = 0.0, param_diff = 0.0;
        for (int n = 0; n <= 20; ++n) {
            const double v = n / 20.0;
            const Jet3 JA = A.jet(1.0, v), JB = B.jet(0.0, v);
            const SeamJump j = seam_jump(JA, JB, JA[2]);
            for (double x : {j.normal_angle, j.db_tt, j.db_tn, j.db_nn, j.dC_ttt, j.dC_ttn, j.dC_tnn, j.dC_nnn})
                worst = std::max(worst, std::fabs(x));
            param_diff = std::max(param_diff, norm(JA[jet_index(3, 0)] - JB[jet_index(3, 0)]));
        }
        std::printf("  same surface, reparameterised: max invariant jump %.2e, "
                    "parametric S_uuu differs by %.2f\n", worst, param_diff);
        AK_CHECK(worst < tol);
        AK_CHECK(param_diff > 0.1);
    }

    // Case 2: G1 but not G2. Extra term 0.4 x^2: only b(nu, nu) jumps, by
    // 2 * 0.4 * n_z * (nu^x)^2.
    {
        const BezierSurface B = side_b([](double x, double) { return 0.4 * x * x; });
        for (int n = 0; n <= 20; ++n) {
            const double v = n / 20.0;
            const Jet3 JA = A.jet(1.0, v), JB = B.jet(0.0, v);
            const SeamJump j = seam_jump(JA, JB, JA[2]);
            const SurfaceFrame fa = surface_frame(JA);
            const Vec3d tau = normalized(JA[2] - fa.n * dot(fa.n, JA[2]));
            const double nux = to_param(fa, cross(fa.n, tau))[0];
            AK_CHECK(j.normal_angle < tol);
            AK_CHECK_NEAR(j.db_tt, 0.0, tol);
            AK_CHECK_NEAR(j.db_tn, 0.0, tol);
            AK_CHECK_NEAR(j.db_nn, 0.8 * fa.n.z * nux * nux, tol);
            AK_CHECK(std::fabs(j.db_nn) > 0.1);
        }
    }

    // Case 3: G2 but not G3. Extra term x^3 (c0 + c1 y), varying along the
    // seam: only C(nu, nu, nu) jumps, by 6 (c0 + c1 y) n_z (nu^x)^3.
    {
        const double c0 = 0.5, c1 = -0.8;
        const BezierSurface B = side_b([=](double x, double y) { return x * x * x * (c0 + c1 * y); });
        double worst_other = 0.0, worst_pred = 0.0, largest = 0.0;
        for (int n = 0; n <= 20; ++n) {
            const double v = n / 20.0;
            const Jet3 JA = A.jet(1.0, v), JB = B.jet(0.0, v);
            const SeamJump j = seam_jump(JA, JB, JA[2]);
            const SurfaceFrame fa = surface_frame(JA);
            const Vec3d tau = normalized(JA[2] - fa.n * dot(fa.n, JA[2]));
            const double nux = to_param(fa, cross(fa.n, tau))[0];
            const double predicted = 6.0 * (c0 + c1 * v) * fa.n.z * nux * nux * nux;
            for (double x : {j.normal_angle, j.db_tt, j.db_tn, j.db_nn, j.dC_ttt, j.dC_ttn, j.dC_tnn})
                worst_other = std::max(worst_other, std::fabs(x));
            worst_pred = std::max(worst_pred, std::fabs(j.dC_nnn - predicted));
            largest = std::max(largest, std::fabs(j.dC_nnn));
        }
        std::printf("  G2-not-G3 seam: C(nu,nu,nu) jump up to %.3f (matches prediction to %.2e); "
                    "all other jumps <= %.2e\n", largest, worst_pred, worst_other);
        AK_CHECK(worst_other < tol);
        AK_CHECK(worst_pred < tol);
        AK_CHECK(largest > 0.5);
    }

    return akt::finish("continuity_test");
}
