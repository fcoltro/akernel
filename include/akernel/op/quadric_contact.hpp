// akernel — two curved surfaces touching along a line, decided exactly.
//
// A tangency between two curved surfaces is a root of an over-determined
// system, like a tangency with a plane (op/contact.hpp): interval methods can
// prove that two surfaces miss each other, never that they touch. This file
// decides the contact of a patch X with a curved patch Y that lies on a
// quadric (a cylinder, cone or sphere sector) exactly, from the input data,
// with no tolerance. Two equal cylinders side by side, touching along a
// ruling, are the commonest case.
//
// The implicit form of Y. A quadric is x^T A x = 0 for a symmetric 4x4 A in
// homogeneous coordinates (op/surfaces.hpp). Candidate entries of A come from
// the null space of the linear conditions "the Bernstein coefficients of
// A(Y, Y) vanish", solved in doubles, scaled to integers and snapped; that
// Y lies on the candidate is then proven exactly (on_quadric). A candidate
// that is slightly off fails the proof, and nothing is claimed.
//
// The contact function. f(u, v) = A(X(u, v), X(u, v)) is a polynomial of
// bidegree (2p, 2q) with exact Bernstein coefficients (expansion arithmetic,
// scaled by one positive integer). It vanishes exactly where X meets the
// quadric. Each row of f in u (one Bernstein coefficient in v) is a
// polynomial f_l(u); if they share a factor D(u)^2, then
//     f(u, v) = D(u)^2 h(u, v),
// with D found as the greatest common divisor of the rows and their
// derivatives (pseudo-remainder sequences, exact), and h by exact
// pseudo-division (a zero remainder is checked, not assumed). If h has one
// strict sign over the parameter square (op/contact.hpp, strict_sign), then X
// meets the quadric exactly on the lines u = r, r a root of D in [0, 1], and
// lies strictly on one side of it everywhere else: X touches Y along those
// lines at most, and never crosses it. The roots of D are isolated by exact
// signs at dyadic points and interval bounds on D' (one root per interval
// where D is monotone); a dyadic root is found exactly.
//
// The same with the roles of u and v exchanged gives lines v = r.
#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include <akernel/num/expansion.hpp>
#include <akernel/num/interval.hpp>
#include <akernel/op/contact.hpp>
#include <akernel/op/surfaces.hpp>

namespace ak::op {

namespace quadric_contact_detail {

using E = num::Expansion;
using Poly = std::vector<E>; ///< power basis, index = degree; may carry zero leading terms

inline int degree(const Poly& a) {
    for (int i = static_cast<int>(a.size()) - 1; i >= 0; --i)
        if (a[static_cast<std::size_t>(i)].sign() != 0) return i;
    return -1;
}

inline Poly trimmed(Poly a) {
    a.resize(static_cast<std::size_t>(degree(a) + 1));
    return a;
}

inline Poly derivative(const Poly& a) {
    Poly d;
    for (std::size_t i = 1; i < a.size(); ++i) d.push_back(a[i] * E(static_cast<double>(i)));
    return d;
}

/// lc(g)^(n - m + 1) a = q g + r, a taken with its formal degree n =
/// a.size() - 1 (leading zeros allowed), so that every row of one size is
/// multiplied by the same constant. g must have a nonzero leading term.
inline void pseudo_divide(const Poly& a, const Poly& g, Poly& q, Poly& r) {
    const int n = static_cast<int>(a.size()) - 1, m = static_cast<int>(g.size()) - 1;
    r = a;
    q.assign(static_cast<std::size_t>(std::max(0, n - m + 1)), E());
    if (n < m) return;
    const E& lc = g.back();
    for (int i = n; i >= m; --i) {
        const E c = r[static_cast<std::size_t>(i)];
        for (E& x : r) x = lc * x;
        for (E& x : q) x = lc * x;
        q[static_cast<std::size_t>(i - m)] = q[static_cast<std::size_t>(i - m)] + c;
        if (c.sign() != 0)
            for (int k = 0; k <= m; ++k) r[static_cast<std::size_t>(i - m + k)] = r[static_cast<std::size_t>(i - m + k)] - c * g[static_cast<std::size_t>(k)];
    }
    r.resize(static_cast<std::size_t>(m));
}

/// A greatest common divisor, up to a constant: the last nonzero term of a
/// pseudo-remainder sequence. Empty if both are zero.
inline Poly gcd(Poly a, Poly b) {
    a = trimmed(std::move(a));
    b = trimmed(std::move(b));
    if (a.size() < b.size()) std::swap(a, b);
    while (!b.empty()) {
        Poly q, r;
        pseudo_divide(a, b, q, r);
        a = std::move(b);
        b = trimmed(std::move(r));
    }
    return a;
}

inline double factorial(int n) {
    double f = 1.0;
    for (int i = 2; i <= n; ++i) f *= i;
    return f;
}

inline double binom(int n, int k) { return factorial(n) / (factorial(k) * factorial(n - k)); }

/// Bernstein coefficients b_0..b_n (all scaled by one constant) to power
/// coefficients (the same constant): a_m = sum_k b_k C(n,k) C(n-k,m-k) (-1)^(m-k).
inline Poly to_power(const std::vector<E>& b) {
    const int n = static_cast<int>(b.size()) - 1;
    Poly a(b.size());
    for (int m = 0; m <= n; ++m)
        for (int k = 0; k <= m; ++k) {
            const double w = binom(n, k) * binom(n - k, m - k) * ((m - k) % 2 ? -1.0 : 1.0);
            a[static_cast<std::size_t>(m)] = a[static_cast<std::size_t>(m)] + b[static_cast<std::size_t>(k)] * E(w);
        }
    return a;
}

/// Power coefficients (degree n = a.size() - 1) to Bernstein coefficients
/// times n!: b_k n! = sum_{i <= k} a_i k! (n - i)! / (k - i)!, integers.
inline std::vector<E> to_bernstein(const Poly& a) {
    const int n = static_cast<int>(a.size()) - 1;
    std::vector<E> b(a.size());
    for (int k = 0; k <= n; ++k)
        for (int i = 0; i <= k; ++i)
            b[static_cast<std::size_t>(k)] =
                b[static_cast<std::size_t>(k)] + a[static_cast<std::size_t>(i)] * E(factorial(k) * factorial(n - i) / factorial(k - i));
    return b;
}

/// The exact value at a double x (exact for the short dyadic x used here).
inline E eval(const Poly& a, double x) {
    E s;
    for (int i = static_cast<int>(a.size()) - 1; i >= 0; --i) s = s * E(x) + a[static_cast<std::size_t>(i)];
    return s;
}

inline num::Interval enclose(const E& e) {
    num::Interval s(0.0);
    for (double x : e.components()) s = s + num::Interval(x);
    return s;
}

inline num::Interval eval(const std::vector<num::Interval>& a, const num::Interval& x) {
    num::Interval s(0.0);
    for (int i = static_cast<int>(a.size()) - 1; i >= 0; --i) s = s * x + a[static_cast<std::size_t>(i)];
    return s;
}

/// Roots of D in the open interval (0, 1), each enclosed (a point interval
/// when it is dyadic and found exactly). False if they could not be
/// isolated, or if D vanishes at 0 or 1 (contact along a side, not decided
/// here).
inline bool interior_roots(const Poly& D, std::vector<num::Interval>& roots) {
    roots.clear();
    if (eval(D, 0.0).sign() == 0 || eval(D, 1.0).sign() == 0) return false;
    const Poly dD = derivative(D);
    std::vector<num::Interval> Di, dDi;
    for (const E& c : D) Di.push_back(enclose(c));
    for (const E& c : dD) dDi.push_back(enclose(c));
    bool ok = true;
    // [a, b] with exact nonzero signs sa, sb at its ends.
    auto rec = [&](auto&& self, double a, double b, int sa, int sb, int depth) -> void {
        if (!ok) return;
        const num::Interval x(a, b);
        if (!num::contains_zero(eval(Di, x))) return;
        if (!num::contains_zero(eval(dDi, x))) {
            if (sa == sb) return; // monotone, no sign change: no root
            // One simple root; narrow it down.
            double lo = a, hi = b;
            int slo = sa;
            for (int it = 0; it < 200 && hi - lo > 0.0; ++it) {
                const double mid = 0.5 * (lo + hi);
                if (!(lo < mid && mid < hi)) break;
                const int sm = eval(D, mid).sign();
                if (sm == 0) { roots.push_back(num::Interval(mid)); return; }
                if (sm == slo) { lo = mid; } else { hi = mid; }
                if (hi - lo < 1e-15) break;
            }
            roots.push_back(num::Interval(lo, hi));
            return;
        }
        if (depth <= 0) { ok = false; return; }
        const double mid = 0.5 * (a + b);
        const int sm = eval(D, mid).sign();
        if (sm == 0) {
            // A dyadic root, exactly; D' must not vanish there for it to be
            // simple, which the caller's factorisation does not need.
            roots.push_back(num::Interval(mid));
            const double eps = std::ldexp(1.0, -40);
            const int sl = eval(D, mid - eps).sign(), sr = eval(D, mid + eps).sign();
            if (sl == 0 || sr == 0) { ok = false; return; }
            self(self, a, mid - eps, sa, sl, depth - 1);
            self(self, mid + eps, b, sr, sb, depth - 1);
            // Roots within eps of mid would be lost: none, if D' excludes 0
            // around it.
            if (num::contains_zero(eval(dDi, num::Interval(mid - eps, mid + eps)))) ok = false;
            return;
        }
        self(self, a, mid, sa, sm, depth - 1);
        self(self, mid, b, sm, sb, depth - 1);
    };
    rec(rec, 0.0, 1.0, eval(D, 0.0).sign(), eval(D, 1.0).sign(), 40);
    std::sort(roots.begin(), roots.end(), [](const num::Interval& p, const num::Interval& q) { return p.lo < q.lo; });
    return ok;
}

/// A null vector of the 10 x rows system "A(Y, Y) has zero Bernstein
/// coefficients" in the upper-triangle entries of A, in doubles, if the
/// system has rank 9 (Y lies on one quadric, not a plane or a line).
inline bool quadric_candidate(const PolyPatchD& Y, double x[10]) {
    using surfaces_detail::binom;
    const int p = Y.p, q = Y.q;
    static const int ia[10] = {0, 0, 0, 0, 1, 1, 1, 2, 2, 3}, ib[10] = {0, 1, 2, 3, 1, 2, 3, 2, 3, 3};
    std::vector<std::array<double, 10>> M;
    for (int k = 0; k <= 2 * p; ++k)
        for (int l = 0; l <= 2 * q; ++l) {
            std::array<double, 10> row{};
            for (int i = std::max(0, k - p); i <= std::min(p, k); ++i)
                for (int j = std::max(0, l - q); j <= std::min(q, l); ++j) {
                    const double w = binom(p, i) * binom(p, k - i) * binom(q, j) * binom(q, l - j);
                    const Vec4d &P = Y.at(i, j), &Q = Y.at(k - i, l - j);
                    const double pv[4] = {P.x, P.y, P.z, P.w}, qv[4] = {Q.x, Q.y, Q.z, Q.w};
                    for (int c = 0; c < 10; ++c) {
                        const int a = ia[c], b = ib[c];
                        row[static_cast<std::size_t>(c)] += w * (a == b ? pv[a] * qv[b] : pv[a] * qv[b] + pv[b] * qv[a]);
                    }
                }
            double s = 0.0;
            for (double v : row) s = std::max(s, std::fabs(v));
            if (s == 0.0) continue;
            for (double& v : row) v /= s;
            M.push_back(row);
        }
    // Gaussian elimination with full pivoting.
    int col[10];
    for (int c = 0; c < 10; ++c) col[c] = c;
    const int R = static_cast<int>(M.size());
    int rank = 0;
    double first = 0.0;
    for (; rank < 10 && rank < R; ++rank) {
        int br = -1, bc = -1;
        double best = 0.0;
        for (int r = rank; r < R; ++r)
            for (int c = rank; c < 10; ++c)
                if (std::fabs(M[static_cast<std::size_t>(r)][static_cast<std::size_t>(c)]) > best) {
                    best = std::fabs(M[static_cast<std::size_t>(r)][static_cast<std::size_t>(c)]);
                    br = r;
                    bc = c;
                }
        if (rank == 0) first = best;
        if (!(best > 1e-9 * first) || best == 0.0) break;
        std::swap(M[static_cast<std::size_t>(rank)], M[static_cast<std::size_t>(br)]);
        if (bc != rank) {
            for (auto& row : M) std::swap(row[static_cast<std::size_t>(rank)], row[static_cast<std::size_t>(bc)]);
            std::swap(col[rank], col[bc]);
        }
        for (int r = rank + 1; r < R; ++r) {
            const double f = M[static_cast<std::size_t>(r)][static_cast<std::size_t>(rank)] / M[static_cast<std::size_t>(rank)][static_cast<std::size_t>(rank)];
            if (f == 0.0) continue;
            for (int c = rank; c < 10; ++c)
                M[static_cast<std::size_t>(r)][static_cast<std::size_t>(c)] -= f * M[static_cast<std::size_t>(rank)][static_cast<std::size_t>(c)];
        }
    }
    if (rank != 9) return false;
    double y[10] = {};
    y[9] = 1.0;
    for (int r = 8; r >= 0; --r) {
        double s = 0.0;
        for (int c = r + 1; c < 10; ++c) s += M[static_cast<std::size_t>(r)][static_cast<std::size_t>(c)] * y[c];
        y[r] = -s / M[static_cast<std::size_t>(r)][static_cast<std::size_t>(r)];
    }
    for (int c = 0; c < 10; ++c) x[col[c]] = y[c];
    return true;
}

} // namespace quadric_contact_detail

/// The quadric a curved patch lies on, proven exactly: candidate integer
/// entries from the null space (doubles), checked by on_quadric. False if Y
/// lies on no single quadric (a plane lies on many) or the candidate is not
/// exact.
inline bool implicit_quadric(const PolyPatchD& Y, Quadric& out) {
    using namespace quadric_contact_detail;
    for (const Vec4d& c : Y.c)
        if (!(c.w > 0.0)) return false;
    double x[10];
    if (!quadric_candidate(Y, x)) return false;
    double big = 0.0;
    for (double v : x) big = std::max(big, std::fabs(v));
    if (!(big > 0.0) || !std::isfinite(big)) return false;
    // Rational candidates with small denominators, made integers by their
    // common denominator.
    double nu[10], den[10], l = 1.0;
    for (int c = 0; c < 10; ++c) {
        const double v = x[c] / big;
        if (std::fabs(v) < 1e-12) { nu[c] = 0.0; den[c] = 1.0; continue; }
        if (!rational_near(v, nu[c], den[c])) return false;
        // l = lcm(l, den)
        double a = l, b = den[c];
        while (b != 0.0) { const double t = std::fmod(a, b); a = b; b = t; }
        l = l / a * den[c];
        if (!(l <= 16777216.0)) return false;
    }
    double e[10];
    for (int c = 0; c < 10; ++c) e[c] = nu[c] * (l / den[c]); // integers below 2^37, exact
    out = Quadric::symmetric(e[0], e[1], e[2], e[3], e[4], e[5], e[6], e[7], e[8], e[9]);
    return on_quadric(Y, out);
}

/// f(u, v) = A(X(u, v), X(u, v)) for the quadric A: exact Bernstein
/// coefficients of bidegree (2p, 2q), all multiplied by (2p)! (2q)!.
inline ExactBernstein quadric_function(const PolyPatchD& X, const Quadric& Q) {
    using namespace surfaces_detail;
    using quadric_contact_detail::factorial;
    const int p = X.p, q = X.q;
    ExactBernstein f(2 * p, 2 * q);
    for (int k = 0; k <= 2 * p; ++k)
        for (int l = 0; l <= 2 * q; ++l) {
            Expansion c;
            for (int i = std::max(0, k - p); i <= std::min(p, k); ++i)
                for (int j = std::max(0, l - q); j <= std::min(q, l); ++j) {
                    const double w = binom(p, i) * binom(p, k - i) * binom(q, j) * binom(q, l - j);
                    c = c + Expansion(w) * form(Q, X.at(i, j), X.at(k - i, l - j));
                }
            // c = C(2p, k) C(2q, l) times the Bernstein coefficient.
            f.at(k, l) = c * Expansion(factorial(k) * factorial(2 * p - k) * factorial(l) * factorial(2 * q - l));
        }
    return f;
}

/// X touching a curved surface Y that lies on a quadric, along lines of
/// constant parameter inside X's square (see the header).
struct QuadricLineContact {
    bool found = false;
    bool fixed_u = true;             ///< the lines are u = r (else v = r)
    std::vector<num::Interval> r;    ///< their parameters, enclosed, increasing
    int sign = 0;                    ///< sign of A(X, X) off the lines (for the A found)
    Quadric Q{};                     ///< Y's quadric
};

/// Decide exactly whether X meets the quadric of Y only along lines of
/// constant parameter inside X's square, without crossing it there.
inline QuadricLineContact quadric_line_contact(const PolyPatchD& X, const PolyPatchD& Y, int depth = 12) {
    using namespace quadric_contact_detail;
    QuadricLineContact out;
    for (const Vec4d& c : X.c)
        if (!(c.w > 0.0)) return out;
    Quadric Q;
    if (!implicit_quadric(Y, Q)) return out;
    const ExactBernstein f = quadric_function(X, Q);
    bool zero = true;
    for (const E& c : f.c)
        if (c.sign() != 0) { zero = false; break; }
    if (zero) return out; // X lies on Y's quadric: coincident surfaces, not decided here
    for (int dir = 0; dir < 2; ++dir) { // 0: lines u = r (rows in u), 1: lines v = r
        const int n = dir == 0 ? f.p : f.q, rows = dir == 0 ? f.q + 1 : f.p + 1;
        if (n < 2) continue;
        std::vector<Poly> row(static_cast<std::size_t>(rows));
        for (int l = 0; l < rows; ++l) {
            std::vector<E> b(static_cast<std::size_t>(n + 1));
            for (int k = 0; k <= n; ++k) b[static_cast<std::size_t>(k)] = dir == 0 ? f.at(k, l) : f.at(l, k);
            row[static_cast<std::size_t>(l)] = to_power(b);
        }
        // D: the common factor of every row and its derivative.
        Poly D;
        for (const Poly& a : row) {
            if (degree(a) < 0) continue;
            D = D.empty() ? gcd(a, derivative(a)) : gcd(gcd(D, a), derivative(a));
            if (degree(D) < 1) break;
        }
        if (degree(D) < 1) continue;
        const Poly D2 = trimmed([&] {
            Poly s(D.size() * 2 - 1);
            for (std::size_t i = 0; i < D.size(); ++i)
                for (std::size_t j = 0; j < D.size(); ++j) s[i + j] = s[i + j] + D[i] * D[j];
            return s;
        }());
        const int m = static_cast<int>(D2.size()) - 1;
        if (m > n) continue;
        // h = f / D^2 row by row, every row times the same lc^(n - m + 1).
        ExactBernstein h(dir == 0 ? n - m : f.p, dir == 0 ? f.q : n - m);
        bool divides = true;
        for (int l = 0; l < rows && divides; ++l) {
            Poly qt, rm;
            pseudo_divide(row[static_cast<std::size_t>(l)], D2, qt, rm);
            for (const E& c : rm)
                if (c.sign() != 0) { divides = false; break; }
            if (!divides) break;
            const std::vector<E> b = to_bernstein(qt);
            for (int k = 0; k <= n - m; ++k) (dir == 0 ? h.at(k, l) : h.at(l, k)) = b[static_cast<std::size_t>(k)];
        }
        if (!divides) continue;
        const int s = contact_detail::strict_sign(h, depth);
        if (s == 0) continue;
        std::vector<num::Interval> roots;
        if (!interior_roots(D, roots) || roots.empty()) continue;
        out.found = true;
        out.fixed_u = dir == 0;
        out.r = std::move(roots);
        // The sign of f off the lines: sign of h times sign of the pseudo-
        // division constant lc(D^2)^(n - m + 1).
        const int lc = D2.back().sign();
        out.sign = ((n - m + 1) % 2 == 1 && lc < 0) ? -s : s;
        out.Q = Q;
        return out;
    }
    return out;
}

} // namespace ak::op
