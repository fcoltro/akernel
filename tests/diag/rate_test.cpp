// Convergence-rate estimators, and their use as an early warning of
// tangential intersection.
#include <akernel/diag/rate.hpp>
#include <akernel/geom/primitives.hpp>
#include <akernel/solve/ssi.hpp>

#include "check.hpp"

using namespace ak;

namespace {

// Residual history of minimum-norm Newton steps towards the intersection.
std::vector<double> newton_history(const solve::SsiSystem& sys, solve::Pt4 x) {
    std::vector<double> hist;
    for (int it = 0; it < 40; ++it) {
        double F[3], J[3][4];
        sys.eval(x, F, J);
        hist.push_back(std::sqrt(F[0] * F[0] + F[1] * F[1] + F[2] * F[2]));
        // One Gauss-Newton step through the public projection helper would
        // hide the history, so take the step explicitly: dx = -J^T (J J^T)^-1 F.
        double M[3][3];
        for (int r = 0; r < 3; ++r)
            for (int c = 0; c < 3; ++c) {
                M[r][c] = 0.0;
                for (int k = 0; k < 4; ++k) M[r][c] += J[r][k] * J[c][k];
            }
        const double d = solve::det3(M);
        if (std::fabs(d) < 1e-300) break;
        double y[3];
        for (int c = 0; c < 3; ++c) {
            double Mc[3][3];
            for (int r = 0; r < 3; ++r)
                for (int k = 0; k < 3; ++k) Mc[r][k] = k == c ? F[r] : M[r][k];
            y[c] = solve::det3(Mc) / d;
        }
        for (int k = 0; k < 4; ++k)
            for (int r = 0; r < 3; ++r) x[k] -= J[r][k] * y[r];
    }
    return hist;
}

} // namespace

int main() {
    // fit_rate recovers a known order.
    {
        std::vector<double> h, e;
        for (double x : {0.1, 0.05, 0.025, 0.0125}) { h.push_back(x); e.push_back(3.0 * x * x * x * x); }
        const diag::RateFit f = diag::fit_rate(h, e);
        AK_CHECK(f.ok);
        AK_CHECK_NEAR(f.order, 4.0, 1e-10);
        AK_CHECK_NEAR(f.constant, 3.0, 1e-8);
        AK_CHECK(!diag::fit_rate({0.1}, {1.0}).ok);
        AK_CHECK(!diag::fit_rate({0.1, 0.2}, {1.0, 0.0}).ok);
    }

    // iteration_order recovers quadratic and linear convergence.
    {
        AK_CHECK_NEAR(diag::iteration_order({1e-1, 1e-2, 1e-4, 1e-8}), 2.0, 1e-9);
        AK_CHECK_NEAR(diag::iteration_order({1e-1, 5e-2, 2.5e-2, 1.25e-2}), 1.0, 1e-9);
    }

    // Early warning: the same Newton iteration converges quadratically onto a
    // transversal intersection and only linearly onto a tangential one.
    {
        const geom::PolyPatchD hill = geom::make_bump(1.0);
        const solve::SsiSystem cut(hill, geom::make_plane({-0.5, -0.5, 0.5}, {2, 0, 0}, {0, 2, 0}));
        const solve::SsiSystem touch(hill, geom::make_plane({-0.5, -0.5, 1.0}, {2, 0, 0}, {0, 2, 0}));
        const double q_cut = diag::iteration_order(newton_history(cut, {0.3, 0.3, 0.4, 0.4}));
        const double q_touch = diag::iteration_order(newton_history(touch, {0.4, 0.42, 0.45, 0.46}));
        std::printf("  observed Newton order: transversal %.2f, tangential %.2f\n", q_cut, q_touch);
        AK_CHECK(q_cut > 1.7);
        AK_CHECK(q_touch < 1.3 && q_touch > 0.7);
    }

    return akt::finish("rate_test");
}
