// Interval arithmetic must enclose the exact result. The oracle is the exact
// Expansion arithmetic, which is tested independently in expansion_test.cpp.
#include <akernel/num/expansion.hpp>
#include <akernel/num/interval.hpp>

#include <bit>
#include <cstdint>
#include <limits>
#include <vector>

#include "check.hpp"

using namespace ak::num;

namespace {

// exact value of e lies in [r.lo, r.hi]
bool encloses(const Interval& r, const Expansion& e) {
    return (e - Expansion(r.lo)).sign() >= 0 && (Expansion(r.hi) - e).sign() >= 0;
}

double random_double(akt::Rng& rng) {
    // Wide dynamic range, both signs, never zero.
    const double mant = rng.range(1.0, 2.0);
    const int exp = static_cast<int>(rng.next() % 61) - 30;
    return (rng.next() & 1 ? 1.0 : -1.0) * std::ldexp(mant, exp);
}

} // namespace

int main() {
    akt::Rng rng(20261110);

    // Point operations enclose the exact result.
    for (int i = 0; i < 20000; ++i) {
        const double a = random_double(rng), b = random_double(rng);
        AK_CHECK(encloses(Interval(a) + Interval(b), Expansion(a) + Expansion(b)));
        AK_CHECK(encloses(Interval(a) - Interval(b), Expansion(a) - Expansion(b)));
        AK_CHECK(encloses(Interval(a) * Interval(b), Expansion(a) * Expansion(b)));

        // Division: q in [lo, hi]  <=>  lo*b <= a <= hi*b (b > 0), reversed for b < 0.
        const Interval q = Interval(a) / Interval(b);
        const int lo_side = (Expansion(q.lo) * Expansion(b) - Expansion(a)).sign();
        const int hi_side = (Expansion(q.hi) * Expansion(b) - Expansion(a)).sign();
        if (b > 0) { AK_CHECK(lo_side <= 0 && hi_side >= 0); }
        else       { AK_CHECK(lo_side >= 0 && hi_side <= 0); }

        // Square root: lo^2 <= x <= hi^2.
        const double x = std::fabs(a);
        const Interval s = sqrt(Interval(x));
        AK_CHECK(s.lo >= 0.0);
        AK_CHECK((Expansion(s.lo) * Expansion(s.lo) - Expansion(x)).sign() <= 0);
        AK_CHECK((Expansion(s.hi) * Expansion(s.hi) - Expansion(x)).sign() >= 0);
    }

    // Range operations enclose every pointwise result.
    for (int i = 0; i < 5000; ++i) {
        double a0 = random_double(rng), a1 = random_double(rng);
        double b0 = random_double(rng), b1 = random_double(rng);
        if (a0 > a1) std::swap(a0, a1);
        if (b0 > b1) std::swap(b0, b1);
        const Interval A(a0, a1), B(b0, b1);
        const double x = rng.range(a0, a1), y = rng.range(b0, b1);
        AK_CHECK(contains(A, x) && contains(B, y));
        AK_CHECK(encloses(A + B, Expansion(x) + Expansion(y)));
        AK_CHECK(encloses(A - B, Expansion(x) - Expansion(y)));
        AK_CHECK(encloses(A * B, Expansion(x) * Expansion(y)));
        AK_CHECK(encloses(sqr(A), Expansion(x) * Expansion(x)));
        if (!contains_zero(B)) {
            const Interval q = A / B;
            const int lo_side = (Expansion(q.lo) * Expansion(y) - Expansion(x)).sign();
            const int hi_side = (Expansion(q.hi) * Expansion(y) - Expansion(x)).sign();
            if (y > 0) { AK_CHECK(lo_side <= 0 && hi_side >= 0); }
            else       { AK_CHECK(lo_side >= 0 && hi_side <= 0); }
        }
    }

    // Division by an interval containing zero refuses to decide.
    {
        const Interval q = Interval(1.0) / Interval(-1.0, 1.0);
        AK_CHECK(std::isinf(q.lo) && std::isinf(q.hi));
    }

    // sqr is never negative, even for an interval straddling zero.
    AK_CHECK(sqr(Interval(-2.0, 3.0)).lo == 0.0);
    AK_CHECK(sqr(Interval(-2.0, 3.0)).hi >= 9.0);

    // Signs are only reported when proven.
    AK_CHECK(sign_of(Interval(1.0, 2.0)) == Sign::Positive);
    AK_CHECK(sign_of(Interval(-2.0, -1.0)) == Sign::Negative);
    AK_CHECK(sign_of(Interval(-1.0, 1.0)) == Sign::Undecided);
    AK_CHECK(sign_of(Interval(0.0, 1.0)) == Sign::Undecided);
    AK_CHECK(sign_of(Interval(0.0)) == Sign::Zero);

    // Set relations.
    AK_CHECK(strictly_inside(Interval(1.0, 2.0), Interval(0.5, 2.5)));
    AK_CHECK(!strictly_inside(Interval(1.0, 2.0), Interval(1.0, 2.5)));
    AK_CHECK(subset(Interval(1.0, 2.0), Interval(1.0, 2.5)));
    AK_CHECK(disjoint(Interval(0.0, 1.0), Interval(1.5, 2.0)));
    AK_CHECK(!disjoint(Interval(0.0, 1.0), Interval(1.0, 2.0)));

    // 0.1 + 0.2 is not 0.3 in doubles; the interval still contains 0.3's
    // double and the true sum of the two stored doubles.
    {
        const Interval s = Interval(0.1) + Interval(0.2);
        AK_CHECK(contains(s, 0.3));
        AK_CHECK(width(s) < 1e-15);
    }

    // The outward step is std::nextafter towards -inf / +inf, bit for bit:
    // special values, and random doubles over the whole exponent range.
    {
        const double inf = std::numeric_limits<double>::infinity();
        auto same = [](double a, double b) {
            return (std::isnan(a) && std::isnan(b)) || std::bit_cast<std::uint64_t>(a) == std::bit_cast<std::uint64_t>(b);
        };
        std::vector<double> xs = {0.0, -0.0, inf, -inf, std::numeric_limits<double>::quiet_NaN(),
                                  std::numeric_limits<double>::denorm_min(), -std::numeric_limits<double>::denorm_min(),
                                  std::numeric_limits<double>::min(), -std::numeric_limits<double>::min(),
                                  std::numeric_limits<double>::max(), -std::numeric_limits<double>::max(), 1.0, -1.0,
                                  0.1, -0.1};
        for (int i = 0; i < 100000; ++i) {
            std::uint64_t bits = rng.next();
            const double x = std::bit_cast<double>(bits);
            xs.push_back(x);
        }
        int mismatches = 0;
        for (double x : xs) {
            if (!same(detail::step_up(x), std::nextafter(x, inf))) ++mismatches;
            if (!same(detail::step_down(x), std::nextafter(x, -inf))) ++mismatches;
        }
        AK_CHECK(mismatches == 0);
        std::printf("  outward step against std::nextafter: %d values, %d mismatches\n", static_cast<int>(xs.size()),
                    mismatches);
    }

    return akt::finish("interval_test");
}
