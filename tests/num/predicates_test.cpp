// Orientation predicates must return the exact sign, including on input
// where plain double evaluation is wrong or inconsistent.
#include <akernel/num/predicates.hpp>

#include "check.hpp"

using namespace ak::num;

namespace {

Sign negate(Sign s) {
    if (s == Sign::Positive) return Sign::Negative;
    if (s == Sign::Negative) return Sign::Positive;
    return s;
}

Sign exact2d(const double* a, const double* b, const double* c) {
    return detail::sign_from_int(
        detail::orient2d_expr<Expansion>(a[0], a[1], b[0], b[1], c[0], c[1]).sign());
}

} // namespace

int main() {
    // Basic orientation.
    {
        const double a[2] = {0, 0}, b[2] = {1, 0}, c[2] = {0, 1};
        AK_CHECK(orient2d(a, b, c) == Sign::Positive);
        AK_CHECK(orient2d(a, c, b) == Sign::Negative);
        AK_CHECK(orient2d(a, b, b) == Sign::Zero);
    }

    // The classroom example of Kettner et al.: p moves on a grid of spacing
    // one ulp near (0.5, 0.5); q and r lie on the line y = x. The exact sign
    // is the sign of (j - i).
    {
        predicate_stats() = {};
        const double ulp = 0x1.0p-53;
        const double q[2] = {12.0, 12.0}, r[2] = {24.0, 24.0};
        int naive_wrong = 0;
        for (int i = 0; i < 64; ++i) {
            for (int j = 0; j < 64; ++j) {
                const double p[2] = {0.5 + i * ulp, 0.5 + j * ulp};
                const Sign expected =
                    j > i ? Sign::Positive : (j < i ? Sign::Negative : Sign::Zero);
                AK_CHECK(orient2d(p, q, r) == expected);
                if (orient2d_naive(p, q, r) != expected) ++naive_wrong;
            }
        }
        std::printf("  grid 64x64: naive double wrong on %d of 4096 points\n", naive_wrong);
        AK_CHECK(naive_wrong > 0); // the point of having the cascade
        const PredicateStats& st = predicate_stats();
        std::printf("  grid tiers: double %llu, interval %llu, exact %llu, zeros %llu\n",
                    static_cast<unsigned long long>(st.decided_at[0]),
                    static_cast<unsigned long long>(st.decided_at[1]),
                    static_cast<unsigned long long>(st.decided_at[2]),
                    static_cast<unsigned long long>(st.zero_results));
        AK_CHECK(st.zero_results == 64);
    }

    // Consistency under permutation, on near-degenerate random input.
    akt::Rng rng(99);
    for (int i = 0; i < 20000; ++i) {
        const double t = rng.unit();
        const double a[2] = {rng.range(-1, 1), rng.range(-1, 1)};
        const double b[2] = {rng.range(-1, 1), rng.range(-1, 1)};
        // c is within a few ulps of the line through a and b.
        double c[2] = {a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1])};
        c[0] += (static_cast<int>(rng.next() % 5) - 2) * 0x1.0p-53;
        const Sign s = orient2d(a, b, c);
        AK_CHECK(s == exact2d(a, b, c));
        AK_CHECK(orient2d(b, c, a) == s);
        AK_CHECK(orient2d(c, a, b) == s);
        AK_CHECK(orient2d(a, c, b) == negate(s));
    }

    // Well-separated random input is decided by the first tier.
    {
        predicate_stats() = {};
        for (int i = 0; i < 20000; ++i) {
            const double a[2] = {rng.range(-1, 1), rng.range(-1, 1)};
            const double b[2] = {rng.range(-1, 1), rng.range(-1, 1)};
            const double c[2] = {rng.range(-1, 1), rng.range(-1, 1)};
            AK_CHECK(orient2d(a, b, c) == exact2d(a, b, c));
        }
        const PredicateStats& st = predicate_stats();
        AK_CHECK(st.decided_at[0] > 19900);
    }

    // orient3d: convention and exact coplanarity.
    {
        const double a[3] = {0, 0, 0}, b[3] = {1, 0, 0}, c[3] = {0, 1, 0};
        const double below[3] = {0, 0, -1}, above[3] = {0, 0, 1}, on[3] = {0.25, 0.25, 0};
        AK_CHECK(orient3d(a, b, c, below) == Sign::Positive);
        AK_CHECK(orient3d(a, b, c, above) == Sign::Negative);
        AK_CHECK(orient3d(a, b, c, on) == Sign::Zero);
    }
    for (int i = 0; i < 5000; ++i) {
        // Four points on the plane z = x + y with dyadic coordinates: exactly coplanar.
        double p[4][3];
        for (auto& pt : p) {
            pt[0] = std::ldexp(static_cast<double>(rng.next() % 2000001) - 1000000.0, -20);
            pt[1] = std::ldexp(static_cast<double>(rng.next() % 2000001) - 1000000.0, -20);
            pt[2] = pt[0] + pt[1];
        }
        AK_CHECK(orient3d(p[0], p[1], p[2], p[3]) == Sign::Zero);

        // Lift the last point by one ulp: the sign must flip with the lift
        // direction and with the orientation of the base triangle.
        double up[3] = {p[3][0], p[3][1], std::nextafter(p[3][2], 1e9)};
        double dn[3] = {p[3][0], p[3][1], std::nextafter(p[3][2], -1e9)};
        const Sign s_up = orient3d(p[0], p[1], p[2], up);
        const Sign s_dn = orient3d(p[0], p[1], p[2], dn);
        const double a2[2] = {p[0][0], p[0][1]}, b2[2] = {p[1][0], p[1][1]}, c2[2] = {p[2][0], p[2][1]};
        const Sign tri = orient2d(a2, b2, c2);
        if (tri == Sign::Zero) {
            AK_CHECK(s_up == Sign::Zero && s_dn == Sign::Zero);
        } else {
            AK_CHECK(s_up == negate(tri)); // above a counter-clockwise triangle is Negative
            AK_CHECK(s_dn == tri);
        }
    }

    // incircle: exact zero for cocircular points, correct sign one ulp away.
    {
        const double a[2] = {1, 0}, b[2] = {0, 1}, c[2] = {-1, 0};
        const double on[2] = {0, -1}, in[2] = {0, 0}, out[2] = {2, 0};
        AK_CHECK(incircle(a, b, c, on) == Sign::Zero);
        AK_CHECK(incircle(a, b, c, in) == Sign::Positive);
        AK_CHECK(incircle(a, b, c, out) == Sign::Negative);
        AK_CHECK(incircle(a, c, b, in) == Sign::Negative); // clockwise triangle flips the sign
        const double sq[4][2] = {{0, 0}, {1, 0}, {1, 1}, {0, 1}};
        AK_CHECK(incircle(sq[0], sq[1], sq[2], sq[3]) == Sign::Zero);

        // Pythagorean points on a circle of radius 5, shifted far from the
        // origin so that the double tier cannot decide.
        const double pts[6][2] = {{5, 0}, {3, 4}, {-4, 3}, {0, -5}, {-3, -4}, {4, -3}};
        const double off = 1048576.0 * 3.0;
        double s[6][2];
        for (int i = 0; i < 6; ++i) { s[i][0] = pts[i][0] + off; s[i][1] = pts[i][1] - off; }
        for (int i = 3; i < 6; ++i) {
            AK_CHECK(incircle(s[0], s[1], s[2], s[i]) == Sign::Zero);
            // Move the point one ulp towards the centre and one ulp away.
            const double cx = off, cy = -off;
            double inward[2] = {s[i][0], s[i][1]}, outward[2] = {s[i][0], s[i][1]};
            const int k = std::fabs(s[i][0] - cx) >= std::fabs(s[i][1] - cy) ? 0 : 1;
            const double centre = k == 0 ? cx : cy;
            inward[k] = std::nextafter(s[i][k], centre);
            outward[k] = std::nextafter(s[i][k], s[i][k] > centre ? 1e300 : -1e300);
            AK_CHECK(incircle(s[0], s[1], s[2], inward) == Sign::Positive);
            AK_CHECK(incircle(s[0], s[1], s[2], outward) == Sign::Negative);
        }
    }

    return akt::finish("predicates_test");
}
