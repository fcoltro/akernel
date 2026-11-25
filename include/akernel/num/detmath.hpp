// akernel — elementary functions that give the same bits on every platform.
//
// std::sin, std::cos, std::atan2, std::pow and std::log2 come from the C
// runtime, and the C runtimes of MSVC, glibc and Apple differ in the last bit
// of many results. One last bit is enough to move a random box, a Gauss node
// or a quadrature weight, and then every hash of the benchmark differs between
// Windows and Linux although the kernel itself does the same thing.
//
// The functions here use only +, -, *, / and integer operations, in a fixed
// order. IEEE-754 makes each of those correctly rounded, so with contraction
// off (/fp:strict, -ffp-contract=off) the result is the same double on every
// conforming platform. The algorithms and coefficients are those of fdlibm
// (Sun Microsystems, 1993; freely redistributable), as kept in FreeBSD's libm:
// sin and cos are within 1 ulp, atan and atan2 within 1 ulp.
//
// None of this is part of a certified computation. These functions place
// geometry (rotations, Chebyshev and Gauss nodes) and order things (angles
// around an edge); the certified layers then work on whatever doubles result.
#pragma once

#include <bit>
#include <cmath>
#include <cstdint>

namespace ak::num::det {

namespace detail {

inline std::uint32_t high_word(double x) noexcept {
    return static_cast<std::uint32_t>(std::bit_cast<std::uint64_t>(x) >> 32);
}

inline int biased_exponent(double x) noexcept {
    return static_cast<int>((std::bit_cast<std::uint64_t>(x) >> 52) & 0x7ff);
}

/// fdlibm __kernel_sin on [-pi/4, pi/4]; x + y is the argument, y its tail.
inline double kernel_sin(double x, double y) noexcept {
    constexpr double S1 = -1.66666666666666324348e-01, S2 = 8.33333333332248946124e-03,
                     S3 = -1.98412698298579493134e-04, S4 = 2.75573137070700676789e-06,
                     S5 = -2.50507602534068634195e-08, S6 = 1.58969099521155010221e-10;
    const double z = x * x, v = z * x;
    const double r = S2 + z * (S3 + z * (S4 + z * (S5 + z * S6)));
    return x - ((z * (0.5 * y - v * r) - y) - v * S1);
}

/// fdlibm __kernel_cos on [-pi/4, pi/4]; x + y is the argument, y its tail.
inline double kernel_cos(double x, double y) noexcept {
    constexpr double C1 = 4.16666666666666019037e-02, C2 = -1.38888888888741095749e-03,
                     C3 = 2.48015872894767294178e-05, C4 = -2.75573143513906633035e-07,
                     C5 = 2.08757232129817482790e-09, C6 = -1.13596475577881948265e-11;
    const double z = x * x, w = z * z;
    const double r = z * (C1 + z * (C2 + z * C3)) + w * w * (C4 + z * (C5 + z * C6));
    const double hz = 0.5 * z, v = 1.0 - hz;
    return v + (((1.0 - v) - hz) + (z * r - x * y));
}

/// Largest |x| reduced by rem_pio2: 2^20 * pi/2 (fdlibm's "medium" range).
inline constexpr double kReduceLimit = 1647099.3291652855;

/// x = n pi/2 + (y0 + y1) with |y0 + y1| <= pi/4, for |x| <= kReduceLimit.
/// Cody-Waite reduction with pi/2 split into three parts (fdlibm medium case).
inline int rem_pio2(double x, double& y0, double& y1) noexcept {
    constexpr double invpio2 = 6.36619772367581382433e-01;
    constexpr double pio2_1 = 1.57079632673412561417e+00, pio2_1t = 6.07710050650619224932e-11;
    constexpr double pio2_2 = 6.07710050630396597660e-11, pio2_2t = 2.02226624879595063154e-21;
    constexpr double pio2_3 = 2.02226624871116645580e-21, pio2_3t = 8.47842766036889956997e-32;
    const double t = x * invpio2;
    const double fn = static_cast<double>(static_cast<std::int64_t>(t + (t >= 0.0 ? 0.5 : -0.5)));
    const int n = static_cast<int>(static_cast<std::int64_t>(fn));
    double r = x - fn * pio2_1;
    double w = fn * pio2_1t;
    y0 = r - w;
    const int j = biased_exponent(x);
    if (j - biased_exponent(y0) > 16) { // cancellation: take the next 33 bits of pi/2
        double s = r;
        w = fn * pio2_2;
        r = s - w;
        w = fn * pio2_2t - ((s - r) - w);
        y0 = r - w;
        if (j - biased_exponent(y0) > 49) { // and the last 33
            s = r;
            w = fn * pio2_3;
            r = s - w;
            w = fn * pio2_3t - ((s - r) - w);
            y0 = r - w;
        }
    }
    y1 = (r - y0) - w;
    return n;
}

} // namespace detail

/// Sine; the same double on every platform for |x| <= 1.6e6.
/// Beyond that the C runtime is used (and may differ in the last bit).
inline double sin(double x) noexcept {
    if (!(std::fabs(x) <= detail::kReduceLimit)) return std::isfinite(x) ? std::sin(x) : x - x;
    if (std::fabs(x) <= 0.78539816339744828) return detail::kernel_sin(x, 0.0);
    double y0, y1;
    switch (detail::rem_pio2(x, y0, y1) & 3) {
    case 0: return detail::kernel_sin(y0, y1);
    case 1: return detail::kernel_cos(y0, y1);
    case 2: return -detail::kernel_sin(y0, y1);
    default: return -detail::kernel_cos(y0, y1);
    }
}

/// Cosine; the same double on every platform for |x| <= 1.6e6.
inline double cos(double x) noexcept {
    if (!(std::fabs(x) <= detail::kReduceLimit)) return std::isfinite(x) ? std::cos(x) : x - x;
    if (std::fabs(x) <= 0.78539816339744828) return detail::kernel_cos(x, 0.0);
    double y0, y1;
    switch (detail::rem_pio2(x, y0, y1) & 3) {
    case 0: return detail::kernel_cos(y0, y1);
    case 1: return -detail::kernel_sin(y0, y1);
    case 2: return -detail::kernel_cos(y0, y1);
    default: return detail::kernel_sin(y0, y1);
    }
}

/// Tangent as sin / cos (within 2 ulp; same bits everywhere in sin's range).
inline double tan(double x) noexcept { return det::sin(x) / det::cos(x); }

/// Arc tangent, fdlibm s_atan.c. Same double on every platform.
inline double atan(double x) noexcept {
    constexpr double atanhi[] = {4.63647609000806093515e-01, 7.85398163397448278999e-01,
                                 9.82793723247329054082e-01, 1.57079632679489655800e+00};
    constexpr double atanlo[] = {2.26987774529616870924e-17, 3.06161699786838301793e-17,
                                 1.39033110312309984516e-17, 6.12323399573676603587e-17};
    constexpr double aT[] = {3.33333333333329318027e-01,  -1.99999999998764832476e-01,
                             1.42857142725034663711e-01,  -1.11111104054623557880e-01,
                             9.09088713343650656196e-02,  -7.69187620504482999495e-02,
                             6.66107313738753120669e-02,  -5.83357013379057348645e-02,
                             4.97687799461593236017e-02,  -3.65315727442169155270e-02,
                             1.62858201153657823623e-02};
    const std::uint32_t hx = detail::high_word(x), ix = hx & 0x7fffffffu;
    const bool neg = (hx >> 31) != 0;
    if (ix >= 0x44100000u) { // |x| >= 2^66, or NaN
        if (std::isnan(x)) return x + x;
        return neg ? -(atanhi[3] + atanlo[3]) : atanhi[3] + atanlo[3];
    }
    int id;
    if (ix < 0x3fdc0000u) { // |x| < 0.4375
        if (ix < 0x3e400000u) return x; // |x| < 2^-27
        id = -1;
    } else {
        x = std::fabs(x);
        if (ix < 0x3ff30000u) {        // |x| < 1.1875
            if (ix < 0x3fe60000u) {    // 7/16 <= |x| < 11/16
                id = 0;
                x = (2.0 * x - 1.0) / (2.0 + x);
            } else {                   // 11/16 <= |x| < 19/16
                id = 1;
                x = (x - 1.0) / (x + 1.0);
            }
        } else if (ix < 0x40038000u) { // |x| < 2.4375
            id = 2;
            x = (x - 1.5) / (1.0 + 1.5 * x);
        } else {                       // 2.4375 <= |x| < 2^66
            id = 3;
            x = -1.0 / x;
        }
    }
    const double z = x * x, w = z * z;
    const double s1 = z * (aT[0] + w * (aT[2] + w * (aT[4] + w * (aT[6] + w * (aT[8] + w * aT[10])))));
    const double s2 = w * (aT[1] + w * (aT[3] + w * (aT[5] + w * (aT[7] + w * aT[9]))));
    if (id < 0) return x - x * (s1 + s2);
    const double r = atanhi[id] - ((x * (s1 + s2) - atanlo[id]) - x);
    return neg ? -r : r;
}

/// Two-argument arc tangent, fdlibm e_atan2.c. Same double on every platform.
inline double atan2(double y, double x) noexcept {
    constexpr double pi = 3.1415926535897931160e+00, pi_lo = 1.2246467991473531772e-16;
    constexpr double pi_o_2 = 1.5707963267948965580e+00, pi_o_4 = 7.8539816339744827900e-01;
    if (std::isnan(x) || std::isnan(y)) return x + y;
    if (x == 1.0) return det::atan(y);
    const std::uint32_t hx = detail::high_word(x), hy = detail::high_word(y);
    const int m = static_cast<int>(((hy >> 31) & 1u) | ((hx >> 30) & 2u)); // 2 * sign(x) + sign(y)
    if (y == 0.0) {
        switch (m) {
        case 0:
        case 1: return y; // atan(+-0, +anything) = +-0
        case 2: return pi;
        default: return -pi;
        }
    }
    if (x == 0.0) return (hy >> 31) ? -pi_o_2 : pi_o_2;
    if (std::isinf(x)) {
        if (std::isinf(y)) {
            switch (m) {
            case 0: return pi_o_4;
            case 1: return -pi_o_4;
            case 2: return 3.0 * pi_o_4;
            default: return -3.0 * pi_o_4;
            }
        }
        switch (m) {
        case 0: return 0.0;
        case 1: return -0.0;
        case 2: return pi;
        default: return -pi;
        }
    }
    if (std::isinf(y)) return (hy >> 31) ? -pi_o_2 : pi_o_2;
    const int k = (static_cast<int>((hy & 0x7fffffffu) >> 20)) - static_cast<int>((hx & 0x7fffffffu) >> 20);
    double z;
    int mm = m;
    if (k > 60) { // |y / x| > 2^60
        z = pi_o_2 + 0.5 * pi_lo;
        mm &= 1;
    } else if ((hx >> 31) && k < -60) {
        z = 0.0; // |y| / x < -2^-60
    } else {
        z = det::atan(std::fabs(y / x));
    }
    switch (mm) {
    case 0: return z;
    case 1: return -z;
    case 2: return pi - (z - pi_lo);
    default: return (z - pi_lo) - pi;
    }
}

/// x^n for an integer n by repeated squaring: only multiplications (and one
/// division for n < 0), so the same double everywhere. Not the correctly
/// rounded power; the error grows like log2(n) roundings.
inline double ipow(double x, int n) noexcept {
    if (n < 0) return 1.0 / ipow(x, -n);
    double r = 1.0;
    while (n > 0) {
        if (n & 1) r *= x;
        n >>= 1;
        if (n > 0) x *= x;
    }
    return r;
}

/// The integer nearest to log2(x) for finite x > 0, decided exactly: the
/// boundary between 2^e and 2^(e+1) is 2^e sqrt(2), and no double equals it.
inline int round_log2(double x) noexcept {
    int e;
    const double m = std::frexp(x, &e); // x = m 2^e, m in [1/2, 1)
    // log2(x) = e + log2(m), log2(m) in [-1, 0); nearest integer is e when
    // m >= 1/sqrt(2), else e - 1. 1/sqrt(2) is irrational, so m never ties;
    // the double below is the largest one under it.
    return m > 0.70710678118654746 ? e : e - 1;
}

/// The smallest integer k with 2^k >= x, for finite x > 0, decided exactly.
inline int ceil_log2(double x) noexcept {
    int e;
    const double m = std::frexp(x, &e);
    return m == 0.5 ? e - 1 : e;
}

} // namespace ak::num::det
