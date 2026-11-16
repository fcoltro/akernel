// akernel — rigorous interval arithmetic on doubles.
//
// Every operation returns an interval that is guaranteed to contain the exact
// real result for every choice of operands inside the input intervals.
//
// Rigour argument: under IEEE-754 round-to-nearest, +, -, *, / and sqrt are
// correct to half an ulp. Stepping the computed bound one representable
// double outward therefore encloses the exact value. This requires that the
// compiler neither contracts (fused multiply-add) nor reorders floating-point
// expressions: build with /fp:strict (MSVC) or -ffp-contract=off.
#pragma once

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstdint>
#include <limits>

#include <akernel/num/escalation_policy.hpp>

namespace ak::num {

class Interval {
public:
    double lo = 0.0;
    double hi = 0.0;

    constexpr Interval() noexcept = default;
    constexpr Interval(double x) noexcept : lo(x), hi(x) {} // exact point
    constexpr Interval(double l, double h) noexcept : lo(l), hi(h) {}

    static constexpr Interval entire() noexcept {
        return {-std::numeric_limits<double>::infinity(),
                std::numeric_limits<double>::infinity()};
    }
};

namespace detail {
// One representable double outward, by stepping the bit pattern: the same
// result as std::nextafter towards -inf / +inf for every input (NaN, both
// zeros, subnormals, the largest finite values and both infinities), without
// a library call. IEEE-754 doubles of one sign are ordered like their bit
// patterns read as integers.
inline double step_down(double x) noexcept {
    if (!(x == x) || x == -std::numeric_limits<double>::infinity()) return x;
    if (x == 0.0) return -std::numeric_limits<double>::denorm_min();
    const std::uint64_t b = std::bit_cast<std::uint64_t>(x);
    return std::bit_cast<double>(x > 0.0 ? b - 1 : b + 1);
}
inline double step_up(double x) noexcept {
    if (!(x == x) || x == std::numeric_limits<double>::infinity()) return x;
    if (x == 0.0) return std::numeric_limits<double>::denorm_min();
    const std::uint64_t b = std::bit_cast<std::uint64_t>(x);
    return std::bit_cast<double>(x > 0.0 ? b + 1 : b - 1);
}
/// Widen computed bounds by one double each side; NaN becomes "anything".
inline Interval outward(double l, double h) noexcept {
    if (std::isnan(l) || std::isnan(h)) return Interval::entire();
    return {step_down(l), step_up(h)};
}
} // namespace detail

// --- queries -------------------------------------------------------------

inline double lower(const Interval& a) noexcept { return a.lo; }
inline double upper(const Interval& a) noexcept { return a.hi; }
inline double mid(const Interval& a) noexcept { return 0.5 * a.lo + 0.5 * a.hi; }
/// Upper bound on the width.
inline double width(const Interval& a) noexcept { return detail::step_up(a.hi - a.lo); }
/// Largest absolute value in the interval.
inline double mag(const Interval& a) noexcept { return std::max(std::fabs(a.lo), std::fabs(a.hi)); }
inline bool contains(const Interval& a, double x) noexcept { return a.lo <= x && x <= a.hi; }
inline bool contains_zero(const Interval& a) noexcept { return a.lo <= 0.0 && 0.0 <= a.hi; }
/// a is a subset of b.
inline bool subset(const Interval& a, const Interval& b) noexcept { return b.lo <= a.lo && a.hi <= b.hi; }
/// a lies in the interior of b.
inline bool strictly_inside(const Interval& a, const Interval& b) noexcept { return b.lo < a.lo && a.hi < b.hi; }
inline bool disjoint(const Interval& a, const Interval& b) noexcept { return a.hi < b.lo || b.hi < a.lo; }
inline Interval hull(const Interval& a, const Interval& b) noexcept {
    return {std::min(a.lo, b.lo), std::max(a.hi, b.hi)};
}
inline bool is_finite(const Interval& a) noexcept { return std::isfinite(a.lo) && std::isfinite(a.hi); }

/// Proven sign of every value in the interval, or Undecided.
inline Sign sign_of(const Interval& a) noexcept {
    if (a.lo > 0.0) return Sign::Positive;
    if (a.hi < 0.0) return Sign::Negative;
    if (a.lo == 0.0 && a.hi == 0.0) return Sign::Zero;
    return Sign::Undecided;
}

// The same queries on plain doubles, so templated code can use either type.
inline double lower(double x) noexcept { return x; }
inline double upper(double x) noexcept { return x; }
inline double mid(double x) noexcept { return x; }
inline double width(double) noexcept { return 0.0; }
inline double mag(double x) noexcept { return std::fabs(x); }

// --- arithmetic ------------------------------------------------------------

inline Interval operator-(const Interval& a) noexcept { return {-a.hi, -a.lo}; }

inline Interval operator+(const Interval& a, const Interval& b) noexcept {
    return detail::outward(a.lo + b.lo, a.hi + b.hi);
}
inline Interval operator-(const Interval& a, const Interval& b) noexcept {
    return detail::outward(a.lo - b.hi, a.hi - b.lo);
}
inline Interval operator*(const Interval& a, const Interval& b) noexcept {
    const double p1 = a.lo * b.lo, p2 = a.lo * b.hi, p3 = a.hi * b.lo, p4 = a.hi * b.hi;
    if (std::isnan(p1) || std::isnan(p2) || std::isnan(p3) || std::isnan(p4))
        return Interval::entire();
    return detail::outward(std::min(std::min(p1, p2), std::min(p3, p4)),
                           std::max(std::max(p1, p2), std::max(p3, p4)));
}
inline Interval operator/(const Interval& a, const Interval& b) noexcept {
    if (contains_zero(b)) return Interval::entire();
    const double q1 = a.lo / b.lo, q2 = a.lo / b.hi, q3 = a.hi / b.lo, q4 = a.hi / b.hi;
    if (std::isnan(q1) || std::isnan(q2) || std::isnan(q3) || std::isnan(q4))
        return Interval::entire();
    return detail::outward(std::min(std::min(q1, q2), std::min(q3, q4)),
                           std::max(std::max(q1, q2), std::max(q3, q4)));
}

inline Interval& operator+=(Interval& a, const Interval& b) noexcept { return a = a + b; }
inline Interval& operator-=(Interval& a, const Interval& b) noexcept { return a = a - b; }
inline Interval& operator*=(Interval& a, const Interval& b) noexcept { return a = a * b; }
inline Interval& operator/=(Interval& a, const Interval& b) noexcept { return a = a / b; }

/// x^2, tighter than x * x because it knows both factors are the same value.
inline Interval sqr(const Interval& a) noexcept {
    const double l = a.lo * a.lo, h = a.hi * a.hi;
    if (contains_zero(a)) return {0.0, detail::step_up(std::max(l, h))};
    return detail::outward(std::min(l, h), std::max(l, h));
}
inline double sqr(double x) noexcept { return x * x; }

/// Square root of the non-negative part of the interval.
inline Interval sqrt(const Interval& a) noexcept {
    if (a.hi < 0.0) return Interval::entire(); // empty in exact terms; refuse to decide
    const double l = a.lo <= 0.0 ? 0.0 : std::max(0.0, detail::step_down(std::sqrt(a.lo)));
    return {l, detail::step_up(std::sqrt(a.hi))};
}

inline Interval abs(const Interval& a) noexcept {
    if (a.lo >= 0.0) return a;
    if (a.hi <= 0.0) return -a;
    return {0.0, std::max(-a.lo, a.hi)};
}

} // namespace ak::num
