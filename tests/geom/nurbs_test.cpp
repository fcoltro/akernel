// NURBS surfaces: evaluation, jets, and conversion to Bezier pieces.
#include <akernel/geom/nurbs_surface.hpp>

#include "check.hpp"

using namespace ak::geom;

int main() {
    akt::Rng rng(5);

    // Random rational surface, degree (3, 2), 6 x 5 control points.
    NurbsSurface s(3, 2, 6, 5, open_uniform_knots(6, 3), open_uniform_knots(5, 2));
    for (auto& c : s.cp) {
        const double w = rng.range(0.5, 2.0);
        c = {w * rng.range(-1, 1), w * rng.range(-1, 1), w * rng.range(-1, 1), w};
    }

    const auto pieces = s.to_bezier();
    AK_CHECK(pieces.size() == 9u); // (6 - 3) spans in u times (5 - 2) spans in v

    for (const auto& piece : pieces) {
        const BezierSurface bz(piece.patch);
        const double du = piece.u1 - piece.u0, dv = piece.v1 - piece.v0;
        for (int n = 0; n < 40; ++n) {
            const double a = rng.range(0.02, 0.98), b = rng.range(0.02, 0.98);
            const double u = piece.u0 + du * a, v = piece.v0 + dv * b;
            const Jet3 JN = s.jet(u, v), JB = bz.jet(a, b);
            // Chain rule: d/du = (1/du) d/da.
            for (int k = 0; k <= 3; ++k)
                for (int l = 0; l <= 3 - k; ++l) {
                    const double scale = std::pow(du, k) * std::pow(dv, l);
                    const Vec3d d = JN[jet_index(k, l)] * scale - JB[jet_index(k, l)];
                    AK_CHECK(norm(d) < 1e-9);
                }
        }
    }

    // A full cylinder as a 9-point rational quadratic circle swept along z.
    {
        const double r = 2.0, w = std::sqrt(0.5);
        const double cx[9] = {1, 1, 0, -1, -1, -1, 0, 1, 1};
        const double cy[9] = {0, 1, 1, 1, 0, -1, -1, -1, 0};
        const double cw[9] = {1, w, 1, w, 1, w, 1, w, 1};
        NurbsSurface cyl(2, 1, 9, 2, {0, 0, 0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1, 1, 1},
                         {0, 0, 1, 1});
        for (int i = 0; i < 9; ++i)
            for (int j = 0; j < 2; ++j)
                cyl.at(i, j) = {cw[i] * r * cx[i], cw[i] * r * cy[i], cw[i] * (j ? 3.0 : -1.0), cw[i]};
        for (int n = 0; n < 200; ++n) {
            const Vec3d pt = cyl.point(rng.unit(), rng.unit());
            AK_CHECK_NEAR(pt.x * pt.x + pt.y * pt.y, r * r, 1e-13);
        }
        AK_CHECK(cyl.to_bezier().size() == 4u);
        // Every Bezier piece is still on the cylinder.
        for (const auto& piece : cyl.to_bezier()) {
            const BezierSurface bz(piece.patch);
            for (int n = 0; n < 50; ++n) {
                const Vec3d pt = bz.point(rng.unit(), rng.unit());
                AK_CHECK_NEAR(pt.x * pt.x + pt.y * pt.y, r * r, 1e-13);
            }
        }
    }

    // Basis functions form a partition of unity; derivatives sum to zero.
    {
        const auto U = open_uniform_knots(9, 4);
        for (int n = 0; n < 100; ++n) {
            const double u = rng.unit();
            const auto d = ders_basis_funs(find_span(8, 4, u, U), u, 4, 3, U);
            double s0 = 0, s1 = 0, s2 = 0, s3 = 0;
            for (int j = 0; j <= 4; ++j) { s0 += d[0][j]; s1 += d[1][j]; s2 += d[2][j]; s3 += d[3][j]; }
            AK_CHECK_NEAR(s0, 1.0, 1e-13);
            AK_CHECK_NEAR(s1, 0.0, 1e-10);
            AK_CHECK_NEAR(s2, 0.0, 1e-8);
            AK_CHECK_NEAR(s3, 0.0, 1e-6);
        }
    }

    // Gauss-Legendre integrates polynomials of degree 2n-1 exactly.
    {
        std::vector<double> x, w;
        gauss_legendre(5, x, w);
        double i8 = 0, i9 = 0, i0 = 0;
        for (int i = 0; i < 5; ++i) { i0 += w[i]; i8 += w[i] * std::pow(x[i], 8); i9 += w[i] * std::pow(x[i], 9); }
        AK_CHECK_NEAR(i0, 2.0, 1e-14);
        AK_CHECK_NEAR(i8, 2.0 / 9.0, 1e-14);
        AK_CHECK_NEAR(i9, 0.0, 1e-14);
    }

    return akt::finish("nurbs_test");
}
