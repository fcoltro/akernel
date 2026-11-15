// Bezier patches: evaluation, derivatives, subdivision, rigorous enclosures.
#include <akernel/geom/bernstein.hpp>
#include <akernel/geom/primitives.hpp>

#include "check.hpp"

using namespace ak::geom;
using ak::num::Interval;

namespace {

double binom(int n, int k) {
    double c = 1.0;
    for (int i = 1; i <= k; ++i) c = c * (n - k + i) / i;
    return c;
}
double bern(int n, int i, double t) { return binom(n, i) * std::pow(t, i) * std::pow(1 - t, n - i); }

PolyPatchD random_rational(akt::Rng& rng, int p, int q) {
    PolyPatchD s(p, q);
    for (auto& c : s.c) {
        const double w = rng.range(0.5, 2.0);
        c = {w * rng.range(-1, 1), w * rng.range(-1, 1), w * rng.range(-1, 1), w};
    }
    return s;
}

double dist4(const Vec4d& a, const Vec4d& b) {
    return std::fabs(a.x - b.x) + std::fabs(a.y - b.y) + std::fabs(a.z - b.z) + std::fabs(a.w - b.w);
}

} // namespace

int main() {
    akt::Rng rng(11);
    const PolyPatchD s = random_rational(rng, 3, 2);

    // de Casteljau agrees with the explicit Bernstein sum.
    for (int n = 0; n < 200; ++n) {
        const double u = rng.unit(), v = rng.unit();
        Vec4d ref{};
        for (int i = 0; i <= s.p; ++i)
            for (int j = 0; j <= s.q; ++j) ref = ref + s.at(i, j) * (bern(s.p, i, u) * bern(s.q, j, v));
        AK_CHECK(dist4(s.eval(u, v), ref) < 1e-13);
    }

    // Subdivision reproduces the same polynomial on each half.
    {
        PolyPatchD lo, hi;
        const double t = 0.3;
        s.split_u(t, lo, hi);
        for (int n = 0; n < 100; ++n) {
            const double a = rng.unit(), v = rng.unit();
            AK_CHECK(dist4(lo.eval(a, v), s.eval(t * a, v)) < 1e-13);
            AK_CHECK(dist4(hi.eval(a, v), s.eval(t + (1 - t) * a, v)) < 1e-13);
        }
        s.split_v(t, lo, hi);
        for (int n = 0; n < 100; ++n) {
            const double u = rng.unit(), a = rng.unit();
            AK_CHECK(dist4(lo.eval(u, a), s.eval(u, t * a)) < 1e-13);
            AK_CHECK(dist4(hi.eval(u, a), s.eval(u, t + (1 - t) * a)) < 1e-13);
        }
    }

    // Restriction to a box, including boxes that leave [0,1]^2.
    for (int n = 0; n < 50; ++n) {
        const double u0 = rng.range(-0.2, 0.9), u1 = u0 + rng.range(0.01, 0.3);
        const double v0 = rng.range(-0.2, 0.9), v1 = v0 + rng.range(0.01, 0.3);
        const PolyPatchD r = s.restricted(u0, u1, v0, v1);
        for (int k = 0; k < 10; ++k) {
            const double a = rng.unit(), b = rng.unit();
            AK_CHECK(dist4(r.eval(a, b), s.eval(u0 + (u1 - u0) * a, v0 + (v1 - v0) * b)) < 1e-12);
        }
    }

    // Interval restriction: the projected bounds enclose every surface point
    // of the box, and shrink with the box.
    {
        const PolyPatchI si = to_interval(s);
        for (int n = 0; n < 50; ++n) {
            const double u0 = rng.range(0.0, 0.8), u1 = u0 + 0.1;
            const double v0 = rng.range(0.0, 0.8), v1 = v0 + 0.1;
            const auto box = projected_bounds(si.restricted(u0, u1, v0, v1));
            for (int k = 0; k < 20; ++k) {
                const double u = rng.range(u0, u1), v = rng.range(v0, v1);
                const Vec4d h = s.eval(u, v);
                AK_CHECK(contains(box.x, h.x / h.w) && contains(box.y, h.y / h.w) &&
                         contains(box.z, h.z / h.w));
            }
        }
        const auto big = projected_bounds(si.restricted(0.4, 0.6, 0.4, 0.6));
        const auto small = projected_bounds(si.restricted(0.49, 0.51, 0.49, 0.51));
        AK_CHECK(width(small.x) < 0.2 * width(big.x));

        // Interval point evaluation encloses the double evaluation.
        const auto hi = si.eval(Interval(0.37), Interval(0.61));
        const Vec4d hd = s.eval(0.37, 0.61);
        AK_CHECK(contains(hi.x, hd.x) && contains(hi.w, hd.w));
        AK_CHECK(width(hi.x) < 1e-13);
    }

    // Jets of a known polynomial surface: z = u^3 - 2u^2 v + v^3 + uv.
    {
        const BezierSurface surf(fit_patch(
            [](double u, double v) {
                return Vec3d{u, v, u * u * u - 2 * u * u * v + v * v * v + u * v};
            },
            3, 3));
        for (int n = 0; n < 50; ++n) {
            const double u = rng.unit(), v = rng.unit();
            const Jet3 J = surf.jet(u, v);
            AK_CHECK_NEAR(J[jet_index(0, 0)].z, u * u * u - 2 * u * u * v + v * v * v + u * v, 1e-12);
            AK_CHECK_NEAR(J[jet_index(1, 0)].z, 3 * u * u - 4 * u * v + v, 1e-11);
            AK_CHECK_NEAR(J[jet_index(0, 1)].z, -2 * u * u + 3 * v * v + u, 1e-11);
            AK_CHECK_NEAR(J[jet_index(2, 0)].z, 6 * u - 4 * v, 1e-10);
            AK_CHECK_NEAR(J[jet_index(1, 1)].z, -4 * u + 1, 1e-10);
            AK_CHECK_NEAR(J[jet_index(0, 2)].z, 6 * v, 1e-10);
            AK_CHECK_NEAR(J[jet_index(3, 0)].z, 6.0, 1e-9);
            AK_CHECK_NEAR(J[jet_index(2, 1)].z, -4.0, 1e-9);
            AK_CHECK_NEAR(J[jet_index(1, 2)].z, 0.0, 1e-9);
            AK_CHECK_NEAR(J[jet_index(0, 3)].z, 6.0, 1e-9);
            AK_CHECK_NEAR(J[jet_index(1, 0)].x, 1.0, 1e-12);
            AK_CHECK_NEAR(J[jet_index(2, 0)].x, 0.0, 1e-10);
        }
    }

    // Rational jets against central differences of the point evaluation.
    {
        const BezierSurface surf(s);
        const double h = 1e-4;
        for (int n = 0; n < 20; ++n) {
            const double u = rng.range(0.1, 0.9), v = rng.range(0.1, 0.9);
            const Jet3 J = surf.jet(u, v);
            const Vec3d du = (surf.point(u + h, v) - surf.point(u - h, v)) / (2 * h);
            const Vec3d dv = (surf.point(u, v + h) - surf.point(u, v - h)) / (2 * h);
            AK_CHECK(norm(J[1] - du) < 1e-6 && norm(J[2] - dv) < 1e-6);
            // Second and third order: differences of the first-order jet.
            const Jet3 Jp = surf.jet(u + h, v), Jm = surf.jet(u - h, v);
            const Jet3 Kp = surf.jet(u, v + h), Km = surf.jet(u, v - h);
            AK_CHECK(norm(J[jet_index(2, 0)] - (Jp[1] - Jm[1]) / (2 * h)) < 1e-5);
            AK_CHECK(norm(J[jet_index(1, 1)] - (Kp[1] - Km[1]) / (2 * h)) < 1e-5);
            AK_CHECK(norm(J[jet_index(0, 2)] - (Kp[2] - Km[2]) / (2 * h)) < 1e-5);
            AK_CHECK(norm(J[jet_index(3, 0)] - (Jp[3] - Jm[3]) / (2 * h)) < 1e-4);
            AK_CHECK(norm(J[jet_index(2, 1)] - (Kp[3] - Km[3]) / (2 * h)) < 1e-4);
            AK_CHECK(norm(J[jet_index(1, 2)] - (Jp[5] - Jm[5]) / (2 * h)) < 1e-4);
            AK_CHECK(norm(J[jet_index(0, 3)] - (Kp[5] - Km[5]) / (2 * h)) < 1e-4);
        }
    }

    // The cylinder patch lies on x^2 + y^2 = r^2 exactly: the interval
    // residual must contain zero and be tight.
    {
        const double r = 1.5;
        const PolyPatchI cyl = to_interval(make_cylinder_quarter(r, -1.0, 2.0));
        for (int n = 0; n < 200; ++n) {
            const double u = rng.unit(), v = rng.unit();
            const auto h = cyl.eval(Interval(u), Interval(v));
            const Interval x = h.x / h.w, y = h.y / h.w, z = h.z / h.w;
            const Interval res = sqr(x) + sqr(y) - Interval(r * r);
            AK_CHECK(contains_zero(res));
            AK_CHECK(width(res) < 1e-13);
            AK_CHECK(z.lo >= -1.0 - 1e-13 && z.hi <= 2.0 + 1e-13);
            AK_CHECK(x.hi >= 0.0 && y.hi >= 0.0);
        }
    }

    return akt::finish("bernstein_test");
}
