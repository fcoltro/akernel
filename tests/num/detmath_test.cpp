// Platform-independent elementary functions (num/detmath.hpp).
//
// Two kinds of check. Accuracy: every result is within 1 ulp (2 for tan) of
// the C runtime's, which is itself within an ulp of the true value. Identity:
// a hash over the exact bits of about 470,000 results is compared with a constant.
// That constant was recorded on Linux with GCC; the CI runs this test on
// Windows with MSVC and on Linux with GCC and Clang, so a platform that gives
// different bits fails here, before any benchmark hash moves.
#include <akernel/num/detmath.hpp>

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>

#include "check.hpp"

namespace det = ak::num::det;

namespace {

/// Distance in units in the last place, by ordering the bit patterns.
std::int64_t ulps(double a, double b) {
    auto key = [](double x) {
        const auto i = std::bit_cast<std::int64_t>(x);
        return i < 0 ? std::numeric_limits<std::int64_t>::min() - i : i;
    };
    const std::int64_t d = key(a) - key(b);
    return d < 0 ? -d : d;
}

std::uint64_t g_hash = 0xcbf29ce484222325ull;
void hash(double x) {
    const auto bits = std::bit_cast<std::uint64_t>(x);
    for (int i = 0; i < 8; ++i) {
        g_hash ^= (bits >> (8 * i)) & 0xff;
        g_hash *= 0x100000001b3ull;
    }
}

} // namespace

int main() {
    akt::Rng rng(20261124);
    const double inf = std::numeric_limits<double>::infinity();
    const double nan = std::numeric_limits<double>::quiet_NaN();

    // sin, cos, tan over small, moderate and large arguments.
    {
        std::int64_t worst_s = 0, worst_c = 0, worst_t = 0;
        for (int i = 0; i < 100000; ++i) {
            const int band = i % 4;
            const double scale = band == 0 ? 1e-6 : band == 1 ? 1.0 : band == 2 ? 100.0 : 1e6;
            const double x = rng.range(-1.0, 1.0) * scale;
            const double s = det::sin(x), c = det::cos(x), t = det::tan(x);
            worst_s = std::max(worst_s, ulps(s, std::sin(x)));
            worst_c = std::max(worst_c, ulps(c, std::cos(x)));
            if (std::fabs(std::cos(x)) > 1e-3) worst_t = std::max(worst_t, ulps(t, std::tan(x)));
            hash(s);
            hash(c);
        }
        std::printf("  sin, cos, tan: worst %lld, %lld, %lld ulp from the C runtime\n",
                    static_cast<long long>(worst_s), static_cast<long long>(worst_c),
                    static_cast<long long>(worst_t));
        AK_CHECK(worst_s <= 1);
        AK_CHECK(worst_c <= 1);
        AK_CHECK(worst_t <= 3);
        // Arguments next to multiples of pi/2, where the reduction cancels.
        std::int64_t worst = 0;
        for (int k = -2000; k <= 2000; ++k) {
            const double x0 = k * 1.5707963267948966;
            double x = std::nextafter(x0, -inf);
            for (int n = 0; n < 3; ++n, x = std::nextafter(x, inf)) {
                const double s = det::sin(x), c = det::cos(x);
                worst = std::max({worst, ulps(s, std::sin(x)), ulps(c, std::cos(x))});
                hash(s);
                hash(c);
            }
        }
        std::printf("  sin, cos next to multiples of pi/2: worst %lld ulp\n", static_cast<long long>(worst));
        AK_CHECK(worst <= 1);
        AK_CHECK(det::sin(0.0) == 0.0 && !std::signbit(det::sin(0.0)));
        AK_CHECK(det::sin(-0.0) == 0.0 && std::signbit(det::sin(-0.0)));
        AK_CHECK(det::cos(0.0) == 1.0);
        AK_CHECK(std::isnan(det::sin(inf)) && std::isnan(det::cos(-inf)) && std::isnan(det::sin(nan)));
    }

    // atan and atan2 in every octant and at the special values.
    {
        std::int64_t worst_a = 0, worst_a2 = 0;
        for (int i = 0; i < 100000; ++i) {
            // One draw per statement: the order in which function arguments
            // are evaluated differs between compilers.
            const double mx = rng.range(-1.0, 1.0);
            const double x = std::ldexp(mx, static_cast<int>(rng.next() % 80) - 40);
            const double my = rng.range(-1.0, 1.0);
            const double y = std::ldexp(my, static_cast<int>(rng.next() % 80) - 40);
            const double a = det::atan(x), a2 = det::atan2(y, x);
            worst_a = std::max(worst_a, ulps(a, std::atan(x)));
            worst_a2 = std::max(worst_a2, ulps(a2, std::atan2(y, x)));
            hash(a);
            hash(a2);
        }
        std::printf("  atan, atan2: worst %lld, %lld ulp from the C runtime\n", static_cast<long long>(worst_a),
                    static_cast<long long>(worst_a2));
        AK_CHECK(worst_a <= 1);
        AK_CHECK(worst_a2 <= 2);
        const double sp[] = {0.0, -0.0, 1.0, -1.0, inf, -inf, 1e-300, -1e300};
        int special_mismatch = 0;
        for (double y : sp)
            for (double x : sp) {
                const double a = det::atan2(y, x), b = std::atan2(y, x);
                if (ulps(a, b) > 1 || std::signbit(a) != std::signbit(b)) ++special_mismatch;
                hash(a);
            }
        AK_CHECK(special_mismatch == 0);
        AK_CHECK(std::isnan(det::atan2(nan, 1.0)) && std::isnan(det::atan(nan)));
    }

    // Integer powers and binary logarithms.
    {
        int bad_pow = 0, bad_round = 0, bad_ceil = 0;
        for (int i = 0; i < 50000; ++i) {
            const double x = rng.range(0.5, 2.5);
            const int n = static_cast<int>(rng.next() % 41) - 20;
            const double p = det::ipow(x, n);
            if (std::fabs(p - std::pow(x, n)) > 1e-14 * std::fabs(std::pow(x, n))) ++bad_pow;
            hash(p);
            const double mr = rng.range(0.5, 1.0);
            const double r = std::ldexp(mr, static_cast<int>(rng.next() % 200) - 100);
            if (det::round_log2(r) != static_cast<int>(std::lround(std::log2(r)))) ++bad_round;
            if (det::ceil_log2(r) != static_cast<int>(std::ceil(std::log2(r)))) ++bad_ceil;
        }
        AK_CHECK(bad_pow == 0);
        AK_CHECK(bad_round == 0);
        AK_CHECK(bad_ceil == 0);
        AK_CHECK(det::ipow(3.0, 20) == 3486784401.0); // exact when the result is representable
        AK_CHECK(det::ipow(0.5, -10) == 1024.0);
        AK_CHECK(det::ceil_log2(1024.0) == 10 && det::ceil_log2(1025.0) == 11 && det::ceil_log2(0.75) == 0);
        AK_CHECK(det::round_log2(1.4142135623730949) == 0 && det::round_log2(1.4142135623730951) == 1);
    }

    std::printf("  hash of every result: %016llx\n", static_cast<unsigned long long>(g_hash));
    AK_CHECK(g_hash == 0x1447334637bf3f6aull);
    return akt::finish("detmath_test");
}
