// akernel — filtered exact orientation predicates.
//
// Each predicate runs the filter cascade of docs/BLUEPRINT.md §3.2:
//   tier Double    : plain evaluation with a static error bound
//   tier Interval  : interval evaluation
//   tier Expansion : exact evaluation (always decides for these predicates)
// Between tiers the escalation policy is consulted. The result is a proven
// sign, including an exact Zero for degenerate input.
#pragma once

#include <array>
#include <cmath>
#include <cstdint>

#include <akernel/num/escalation_policy.hpp>
#include <akernel/num/expansion.hpp>
#include <akernel/num/interval.hpp>

namespace ak::num {

/// Counters for the diagnostics loop (blueprint §2.2).
struct PredicateStats {
    std::uint64_t calls = 0;
    std::array<std::uint64_t, 3> decided_at{}; ///< [Double, Interval, Expansion]
    std::uint64_t zero_results = 0;            ///< exact degeneracies found
};
inline PredicateStats& predicate_stats() noexcept {
    thread_local PredicateStats s;
    return s;
}

namespace detail {

template <class T>
T orient2d_expr(const T& ax, const T& ay, const T& bx, const T& by, const T& cx, const T& cy) {
    return (ax - cx) * (by - cy) - (ay - cy) * (bx - cx);
}

template <class T>
T orient3d_expr(const T& ax, const T& ay, const T& az, const T& bx, const T& by, const T& bz,
                const T& cx, const T& cy, const T& cz, const T& dx, const T& dy, const T& dz) {
    const T adx = ax - dx, ady = ay - dy, adz = az - dz;
    const T bdx = bx - dx, bdy = by - dy, bdz = bz - dz;
    const T cdx = cx - dx, cdy = cy - dy, cdz = cz - dz;
    return adz * (bdx * cdy - cdx * bdy) + bdz * (cdx * ady - adx * cdy) +
           cdz * (adx * bdy - bdx * ady);
}

inline Sign sign_from_int(int s) noexcept {
    return s > 0 ? Sign::Positive : (s < 0 ? Sign::Negative : Sign::Zero);
}

/// Run the cascade. Each tier is a callable returning a Sign; the first two
/// also report the width of their undecided enclosure through `width`.
template <class DoubleTier, class IntervalTier, class ExactTier>
Sign cascade(DoubleTier&& tier_double, IntervalTier&& tier_interval, ExactTier&& tier_exact) {
    PredicateStats& st = predicate_stats();
    ++st.calls;
    double width = 0.0;

    Sign s = tier_double(width);
    if (s != Sign::Undecided) { ++st.decided_at[0]; return s; }
    if (on_undecided({Tier::Double, 0, 0, width, 0.0, false, false}) != Action::Escalate)
        return Sign::Undecided;

    s = tier_interval(width);
    if (s != Sign::Undecided) { ++st.decided_at[1]; return s; }
    if (on_undecided({Tier::Interval, 0, 0, width, 0.0, false, false}) != Action::Escalate)
        return Sign::Undecided;

    s = tier_exact();
    ++st.decided_at[2];
    if (s == Sign::Zero) ++st.zero_results;
    return s;
}

inline constexpr double kEps = 0x1.0p-53; // half an ulp of 1.0

} // namespace detail

/// Sign of the signed area of triangle (a, b, c): Positive if counter-clockwise.
inline Sign orient2d(const double* a, const double* b, const double* c) {
    using namespace detail;
    return cascade(
        [&](double& width) {
            const double l = (a[0] - c[0]) * (b[1] - c[1]);
            const double r = (a[1] - c[1]) * (b[0] - c[0]);
            const double det = l - r;
            const double bound = (3.0 + 16.0 * kEps) * kEps * (std::fabs(l) + std::fabs(r));
            width = 2.0 * bound;
            if (det > bound) return Sign::Positive;
            if (-det > bound) return Sign::Negative;
            return Sign::Undecided;
        },
        [&](double& width) {
            const Interval d = orient2d_expr<Interval>(a[0], a[1], b[0], b[1], c[0], c[1]);
            width = ak::num::width(d);
            const Sign s = sign_of(d);
            return s == Sign::Zero ? Sign::Undecided : s; // zero is only proven exactly
        },
        [&] {
            return sign_from_int(
                orient2d_expr<Expansion>(a[0], a[1], b[0], b[1], c[0], c[1]).sign());
        });
}

/// Sign of the signed volume of tetrahedron (a, b, c, d): Positive if d sees
/// a, b, c in clockwise order (Shewchuk's convention).
inline Sign orient3d(const double* a, const double* b, const double* c, const double* d) {
    using namespace detail;
    return cascade(
        [&](double& width) {
            const double adx = a[0] - d[0], ady = a[1] - d[1], adz = a[2] - d[2];
            const double bdx = b[0] - d[0], bdy = b[1] - d[1], bdz = b[2] - d[2];
            const double cdx = c[0] - d[0], cdy = c[1] - d[1], cdz = c[2] - d[2];
            const double bdxcdy = bdx * cdy, cdxbdy = cdx * bdy;
            const double cdxady = cdx * ady, adxcdy = adx * cdy;
            const double adxbdy = adx * bdy, bdxady = bdx * ady;
            const double det = adz * (bdxcdy - cdxbdy) + bdz * (cdxady - adxcdy) +
                               cdz * (adxbdy - bdxady);
            const double permanent =
                (std::fabs(bdxcdy) + std::fabs(cdxbdy)) * std::fabs(adz) +
                (std::fabs(cdxady) + std::fabs(adxcdy)) * std::fabs(bdz) +
                (std::fabs(adxbdy) + std::fabs(bdxady)) * std::fabs(cdz);
            const double bound = (7.0 + 56.0 * kEps) * kEps * permanent;
            width = 2.0 * bound;
            if (det > bound) return Sign::Positive;
            if (-det > bound) return Sign::Negative;
            return Sign::Undecided;
        },
        [&](double& width) {
            const Interval v = orient3d_expr<Interval>(a[0], a[1], a[2], b[0], b[1], b[2],
                                                      c[0], c[1], c[2], d[0], d[1], d[2]);
            width = ak::num::width(v);
            const Sign s = sign_of(v);
            return s == Sign::Zero ? Sign::Undecided : s;
        },
        [&] {
            return sign_from_int(orient3d_expr<Expansion>(a[0], a[1], a[2], b[0], b[1], b[2],
                                                          c[0], c[1], c[2], d[0], d[1], d[2])
                                     .sign());
        });
}

namespace detail {
template <class T>
T incircle_expr(const T& ax, const T& ay, const T& bx, const T& by, const T& cx, const T& cy, const T& dx,
                const T& dy) {
    const T adx = ax - dx, ady = ay - dy, bdx = bx - dx, bdy = by - dy, cdx = cx - dx, cdy = cy - dy;
    return (adx * adx + ady * ady) * (bdx * cdy - cdx * bdy) + (bdx * bdx + bdy * bdy) * (cdx * ady - adx * cdy) +
           (cdx * cdx + cdy * cdy) * (adx * bdy - bdx * ady);
}
} // namespace detail

/// Positive if d lies inside the circle through a, b, c (taken
/// counter-clockwise), Negative if outside, Zero if on it.
inline Sign incircle(const double* a, const double* b, const double* c, const double* d) {
    using namespace detail;
    return cascade(
        [&](double& width) {
            const double adx = a[0] - d[0], ady = a[1] - d[1], bdx = b[0] - d[0], bdy = b[1] - d[1];
            const double cdx = c[0] - d[0], cdy = c[1] - d[1];
            const double bdxcdy = bdx * cdy, cdxbdy = cdx * bdy, cdxady = cdx * ady, adxcdy = adx * cdy;
            const double adxbdy = adx * bdy, bdxady = bdx * ady;
            const double alift = adx * adx + ady * ady, blift = bdx * bdx + bdy * bdy, clift = cdx * cdx + cdy * cdy;
            const double det = alift * (bdxcdy - cdxbdy) + blift * (cdxady - adxcdy) + clift * (adxbdy - bdxady);
            const double permanent = (std::fabs(bdxcdy) + std::fabs(cdxbdy)) * alift +
                                     (std::fabs(cdxady) + std::fabs(adxcdy)) * blift +
                                     (std::fabs(adxbdy) + std::fabs(bdxady)) * clift;
            const double bound = (10.0 + 96.0 * kEps) * kEps * permanent;
            width = 2.0 * bound;
            if (det > bound) return Sign::Positive;
            if (-det > bound) return Sign::Negative;
            return Sign::Undecided;
        },
        [&](double& width) {
            const Interval v = incircle_expr<Interval>(a[0], a[1], b[0], b[1], c[0], c[1], d[0], d[1]);
            width = ak::num::width(v);
            const Sign s = sign_of(v);
            return s == Sign::Zero ? Sign::Undecided : s;
        },
        [&] {
            return sign_from_int(
                incircle_expr<Expansion>(a[0], a[1], b[0], b[1], c[0], c[1], d[0], d[1]).sign());
        });
}

/// Unfiltered double evaluation, for comparison in tests and benchmarks only.
inline Sign orient2d_naive(const double* a, const double* b, const double* c) noexcept {
    const double det = (a[0] - c[0]) * (b[1] - c[1]) - (a[1] - c[1]) * (b[0] - c[0]);
    return det > 0.0 ? Sign::Positive : (det < 0.0 ? Sign::Negative : Sign::Zero);
}

} // namespace ak::num
