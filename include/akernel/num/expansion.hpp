// akernel — exact arithmetic on doubles by floating-point expansions.
//
// An Expansion stores a real number as a sum of doubles whose bit ranges do
// not overlap, ordered by increasing magnitude. Sums, differences and
// products of expansions are computed with no rounding error at all, so the
// sign of any polynomial expression in double inputs is decided exactly.
//
// Algorithms follow Shewchuk, "Adaptive precision floating-point arithmetic
// and fast robust geometric predicates" (1997): Two-Sum, Two-Product,
// Grow-Expansion and Scale-Expansion. They are exact provided no
// intermediate overflows or underflows, and provided the compiler does not
// contract or reorder floating-point operations (see interval.hpp).
#pragma once

#include <cmath>
#include <vector>

namespace ak::num {

namespace eft {
/// a + b = x + y exactly, |y| <= ulp(x)/2.
inline void two_sum(double a, double b, double& x, double& y) noexcept {
    x = a + b;
    const double bv = x - a;
    const double av = x - bv;
    y = (a - av) + (b - bv);
}
/// As two_sum, but requires |a| >= |b|.
inline void fast_two_sum(double a, double b, double& x, double& y) noexcept {
    x = a + b;
    y = b - (x - a);
}
/// a * b = x + y exactly.
inline void two_prod(double a, double b, double& x, double& y) noexcept {
    x = a * b;
    y = std::fma(a, b, -x);
}
} // namespace eft

class Expansion {
    // Non-overlapping, increasing magnitude, no zero components.
    // The empty vector represents zero.
    std::vector<double> c_;

    static Expansion grow(const Expansion& e, double b) {
        Expansion h;
        h.c_.reserve(e.c_.size() + 1);
        double q = b;
        for (double ei : e.c_) {
            double sum, err;
            eft::two_sum(q, ei, sum, err);
            if (err != 0.0) h.c_.push_back(err);
            q = sum;
        }
        if (q != 0.0) h.c_.push_back(q);
        return h;
    }

    static Expansion scale(const Expansion& e, double b) {
        Expansion h;
        if (e.c_.empty() || b == 0.0) return h;
        h.c_.reserve(2 * e.c_.size());
        double q, err;
        eft::two_prod(e.c_[0], b, q, err);
        if (err != 0.0) h.c_.push_back(err);
        for (std::size_t i = 1; i < e.c_.size(); ++i) {
            double t_hi, t_lo, sum;
            eft::two_prod(e.c_[i], b, t_hi, t_lo);
            eft::two_sum(q, t_lo, sum, err);
            if (err != 0.0) h.c_.push_back(err);
            eft::fast_two_sum(t_hi, sum, q, err);
            if (err != 0.0) h.c_.push_back(err);
        }
        if (q != 0.0) h.c_.push_back(q);
        return h;
    }

public:
    Expansion() = default;
    Expansion(double x) { if (x != 0.0) c_.push_back(x); } // exact

    /// -1, 0 or +1. Exact: the largest component carries the sign.
    int sign() const noexcept {
        if (c_.empty()) return 0;
        return c_.back() > 0.0 ? 1 : -1;
    }
    /// Nearest-double style approximation of the value.
    double estimate() const noexcept {
        double s = 0.0;
        for (double x : c_) s += x;
        return s;
    }
    std::size_t size() const noexcept { return c_.size(); }
    /// The components, smallest first; their exact sum is the value.
    const std::vector<double>& components() const noexcept { return c_; }

    friend Expansion operator-(const Expansion& a) {
        Expansion r = a;
        for (double& x : r.c_) x = -x;
        return r;
    }
    friend Expansion operator+(const Expansion& a, const Expansion& b) {
        Expansion r = a;
        for (double x : b.c_) r = grow(r, x);
        return r;
    }
    friend Expansion operator-(const Expansion& a, const Expansion& b) {
        Expansion r = a;
        for (double x : b.c_) r = grow(r, -x);
        return r;
    }
    friend Expansion operator*(const Expansion& a, const Expansion& b) {
        Expansion r;
        for (double x : b.c_) r = r + scale(a, x);
        return r;
    }
};

} // namespace ak::num
