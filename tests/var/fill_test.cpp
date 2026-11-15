// Variational fill: exact reproduction, G3 joins measured by invariants,
// and the energy-order rule of blueprint §1.7.
#include <akernel/geom/continuity.hpp>
#include <akernel/var/fill.hpp>

#include "check.hpp"

using namespace ak::geom;
using namespace ak::var;

namespace {

// A surface given as a graph z(x, y) with all partials dz(k, l).
struct Graph {
    std::function<double(int, int, double, double)> dz;

    Jet3 jet(double x, double y) const {
        Jet3 J{};
        for (int k = 0; k <= 3; ++k)
            for (int l = 0; l <= 3 - k; ++l) J[jet_index(k, l)] = {0, 0, dz(k, l, x, y)};
        J[jet_index(0, 0)].x = x;
        J[jet_index(0, 0)].y = y;
        J[jet_index(1, 0)].x = 1.0;
        J[jet_index(0, 1)].y = 1.0;
        return J;
    }
    JetFn jets() const {
        return [this](int side, int j, double s) {
            const bool u_side = side == U0 || side == U1;
            const double x = u_side ? (side == U1 ? 1.0 : 0.0) : s;
            const double y = u_side ? s : (side == V1 ? 1.0 : 0.0);
            return u_side ? jet(x, y)[jet_index(j, 0)] : jet(x, y)[jet_index(0, j)];
        };
    }
};

Graph cubic() {
    return {[](int k, int l, double x, double y) {
        // z = 0.3x^3 - 0.2x^2y + 0.4xy^2 - 0.1y^3 + 0.5x^2 - 0.3xy + 0.2y^2 + 0.1x - 0.2y + 0.05
        const double c[4][4] = {{0.05, -0.2, 0.2, -0.1}, {0.1, -0.3, 0.4, 0}, {0.5, -0.2, 0, 0}, {0.3, 0, 0, 0}};
        auto falling = [](int n, int d) { double f = 1; for (int i = 0; i < d; ++i) f *= (n - i); return f; };
        double s = 0.0;
        for (int i = k; i <= 3; ++i)
            for (int j = l; j <= 3 - i; ++j)
                s += c[i][j] * falling(i, k) * falling(j, l) * std::pow(x, i - k) * std::pow(y, j - l);
        return s;
    }};
}

Graph wave() {
    return {[](int k, int l, double x, double y) {
        const double hp = 1.57079632679489661923;
        double s = 0.2 * std::pow(2.5, k) * std::sin(2.5 * x + 0.3 + k * hp) * std::pow(2.0, l) *
                   std::cos(2.0 * y + l * hp);
        if (k == 0 && l == 0) s += 0.1 * x * y;
        if (k == 1 && l == 0) s += 0.1 * y;
        if (k == 0 && l == 1) s += 0.1 * x;
        if (k == 1 && l == 1) s += 0.1;
        return s;
    }};
}

struct Seam { double g1 = 0, g2 = 0, g3 = 0, other = 0; };

// Largest invariant jumps between the fill and the surrounding surface.
Seam seam_error(const NurbsSurface& s, const Graph& q) {
    Seam e;
    for (int side = 0; side < 4; ++side)
        for (int i = 1; i < 40; ++i) {
            const double t = i / 40.0;
            const bool u_side = side == U0 || side == U1;
            const double x = u_side ? (side == U1 ? 1.0 : 0.0) : t;
            const double y = u_side ? t : (side == V1 ? 1.0 : 0.0);
            const Jet3 JQ = q.jet(x, y), JS = s.jet(x, y);
            const SeamJump j = seam_jump(JQ, JS, u_side ? JQ[2] : JQ[1]);
            e.g1 = std::max(e.g1, j.normal_angle);
            e.g2 = std::max(e.g2, std::fabs(j.db_nn));
            e.g3 = std::max(e.g3, std::fabs(j.dC_nnn));
            for (double v : {j.db_tt, j.db_tn, j.dC_ttt, j.dC_ttn, j.dC_tnn})
                e.other = std::max(e.other, std::fabs(v));
        }
    return e;
}

} // namespace

int main() {
    // 1. A cubic surface has zero fourth-order energy, so with its own jets
    //    prescribed it is the unique minimiser: the fill must reproduce it.
    {
        const Graph q = cubic();
        FillSpec spec; // degree 7, 14 x 14, G3, energy order 4
        const FillResult r = solve_fill(spec, q.jets());
        AK_CHECK(r.ok);
        AK_CHECK(r.free_dofs == 36);
        AK_CHECK(r.corner_mismatch < 1e-11);
        double worst = 0.0;
        for (int i = 0; i <= 20; ++i)
            for (int j = 0; j <= 20; ++j) {
                const double x = i / 20.0, y = j / 20.0;
                worst = std::max(worst, std::fabs(r.surface.point(x, y).z - q.dz(0, 0, x, y)));
            }
        std::printf("  cubic reproduction: max deviation %.2e, E4 = %.2e\n", worst, r.energy[4]);
        AK_CHECK(worst < 1e-11);
        AK_CHECK(r.energy[4] < 1e-12);
    }

    // 2. G3 fill of a hole in a non-polynomial surface. The joins are judged
    //    by the parameterisation-independent jumps, which must shrink under
    //    refinement.
    {
        const Graph q = wave();
        Seam prev;
        bool first = true;
        for (int n : {12, 16, 24}) {
            FillSpec spec;
            spec.ncp = n;
            const FillResult r = solve_fill(spec, q.jets());
            AK_CHECK(r.ok);
            const Seam e = seam_error(r.surface, q);
            const double centre = std::fabs(r.surface.point(0.5, 0.5).z - q.dz(0, 0, 0.5, 0.5));
            std::printf("  wave, %2d x %2d: G1 %.2e rad, G2 %.2e, G3 %.2e, other %.2e; "
                        "centre differs from original by %.2e\n", n, n, e.g1, e.g2, e.g3, e.other, centre);
            if (!first) { AK_CHECK(e.g1 < prev.g1 && e.g2 < prev.g2 && e.g3 < prev.g3); }
            prev = e;
            first = false;
        }
        AK_CHECK(prev.g1 < 1e-8);
        AK_CHECK(prev.g2 < 1e-6);
        AK_CHECK(prev.g3 < 1e-4);
    }

    // 3. Order rule. G2 jets (k = 2) with the thin-plate energy (m = 2) is
    //    over-determined in the continuum: the third-order energy of its
    //    solutions should grow under refinement, while the matched energy
    //    (m = 3) converges.
    {
        const Graph q = wave();
        double e3_plate[3], e3_matched[3];
        const int sizes[3] = {12, 18, 27};
        for (int i = 0; i < 3; ++i) {
            FillSpec plate;
            plate.ncp = sizes[i];
            plate.continuity = 2;
            plate.energy_weights = {0, 0, 1, 0, 0};
            FillSpec matched = plate;
            matched.energy_weights = {0, 0, 0, 1, 0};
            const FillResult a = solve_fill(plate, q.jets()), b = solve_fill(matched, q.jets());
            AK_CHECK(a.ok && b.ok);
            e3_plate[i] = a.energy[3];
            e3_matched[i] = b.energy[3];
            std::printf("  G2 jets, %2d x %2d: E3 with thin-plate energy %.4f, with matched energy %.4f\n",
                        sizes[i], sizes[i], e3_plate[i], e3_matched[i]);
        }
        AK_CHECK(e3_plate[2] > 1.5 * e3_plate[0]);                        // grows
        AK_CHECK(std::fabs(e3_matched[2] - e3_matched[1]) < 0.02 * e3_matched[1]); // settles
    }

    // 4. Bad input is refused rather than half-solved.
    {
        FillSpec spec;
        spec.ncp = 6; // fewer than 2 (k + 1) rows
        AK_CHECK(!solve_fill(spec, cubic().jets()).ok);
    }

    return akt::finish("fill_test");
}
