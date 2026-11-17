// akernel — tangential contact with a plane, decided exactly.
//
// A tangency is a root of an over-determined system: interval methods can
// prove that a surface misses a plane, but never that it touches one, because
// any perturbation destroys the contact (docs/BLUEPRINT.md §3.3). This file
// decides the commonest such contact exactly, from the input data, with no
// tolerance: a curved patch that touches a plane along one of its sides (a
// cylinder against a plane parallel to its axis, at the cylinder's edge).
//
// The plane function. For an affine parallelogram Y with corners c00, c10, c01
// (homogeneous) and any rational patch X, the polynomial
//     g(u, v) = det(c00, c10, c01, X(u, v))
// is linear in X, so its Bernstein coefficients are det(c00, c10, c01, X_ij),
// computed exactly. It vanishes exactly where X(u, v) lies in Y's plane, and,
// the weights being positive, its sign says on which side of the plane the
// point lies.
//
// Contact along a side. If every coefficient on side s of the parameter
// square is zero, g vanishes on that side and is divisible by the side's
// linear factor (v, 1 - v, u or 1 - u); the quotient's coefficients follow
// from the degree identity, e.g. B^q_j(v) = (1 - v) (q / (q - j)) B^(q-1)_j(v),
// and are computed exactly up to one positive constant. The division is
// repeated while the side row stays zero; the number of divisions is the order
// of contact (1: the patch crosses the plane's line transversally there, 2 or
// more: it touches). If the final quotient has one strict sign over the whole
// square, then g vanishes nowhere but on side s: the patch meets the plane in
// exactly that side and lies strictly on one side of it everywhere else. The
// strict sign is proven from the Bernstein coefficients (all of one sign),
// after exact subdivision at 1/2 where needed.
//
// The same test for a curve (a patch of degree (n, 0)) proves that an edge
// meets a plane only at its ends.
#pragma once

#include <cstddef>
#include <vector>

#include <akernel/num/expansion.hpp>
#include <akernel/op/planar.hpp>

namespace ak::op {

/// An exact polynomial in Bernstein form over [0, 1]^2, of degree (p, q).
/// A curve is the case q = 0.
struct ExactBernstein {
    int p = 0, q = 0;
    std::vector<num::Expansion> c; ///< (p + 1)(q + 1) coefficients, index i (q + 1) + j

    ExactBernstein() = default;
    ExactBernstein(int p_, int q_) : p(p_), q(q_), c(static_cast<std::size_t>((p_ + 1) * (q_ + 1))) {}
    num::Expansion& at(int i, int j) { return c[static_cast<std::size_t>(i * (q + 1) + j)]; }
    const num::Expansion& at(int i, int j) const { return c[static_cast<std::size_t>(i * (q + 1) + j)]; }
};

namespace contact_detail {

/// Exact 4x4 determinant of the homogeneous coordinates (x, y, z, w) of four
/// points, expanded along the last one.
inline num::Expansion det4(const Vec4d& a, const Vec4d& b, const Vec4d& c, const Vec4d& d) {
    using E = num::Expansion;
    using planar_detail::det3;
    return E(d.w) * det3(a, b, c) - E(c.w) * det3(a, b, d) + E(b.w) * det3(a, c, d) - E(a.w) * det3(b, c, d);
}

/// Is every coefficient on side s zero? (Side 0: v = 0, 1: u = 1, 2: v = 1, 3: u = 0.)
inline bool side_zero(const ExactBernstein& f, int side) {
    const bool along_u = side == 0 || side == 2;
    const int n = along_u ? f.p : f.q;
    for (int k = 0; k <= n; ++k) {
        const num::Expansion& x = side == 0 ? f.at(k, 0) : side == 1 ? f.at(f.p, k) : side == 2 ? f.at(k, f.q) : f.at(0, k);
        if (x.sign() != 0) return false;
    }
    return true;
}

/// n! / m for 1 <= m <= n, an integer, as a double (exact for the small
/// degrees used here).
inline double fact_over(int n, int m) {
    double f = 1.0;
    for (int i = 2; i <= n; ++i) f *= i;
    return f / m;
}

/// f = (side factor) * h, for f vanishing on side s; h is returned times the
/// positive constant (n - 1)! (n the degree across the side), which keeps
/// every coefficient an exact integer multiple of one of f's.
inline ExactBernstein divide_side(const ExactBernstein& f, int side) {
    const bool across_v = side == 0 || side == 2; // the factor is in v
    ExactBernstein h(across_v ? f.p : f.p - 1, across_v ? f.q - 1 : f.q);
    for (int i = 0; i <= h.p; ++i)
        for (int j = 0; j <= h.q; ++j) {
            num::Expansion x;
            if (side == 0) x = f.at(i, j + 1) * num::Expansion(fact_over(f.q, j + 1));          // v B^(q-1)_j = (j+1)/q B^q_(j+1)
            else if (side == 2) x = f.at(i, j) * num::Expansion(fact_over(f.q, f.q - j));      // (1-v) B^(q-1)_j = (q-j)/q B^q_j
            else if (side == 3) x = f.at(i + 1, j) * num::Expansion(fact_over(f.p, i + 1));    // u
            else x = f.at(i, j) * num::Expansion(fact_over(f.p, f.p - i));                     // 1 - u
            h.at(i, j) = x;
        }
    return h;
}

/// The two halves of f at u = 1/2 (in_u) or v = 1/2, by de Casteljau in exact
/// arithmetic (halving is exact).
inline void split_half(const ExactBernstein& f, bool in_u, ExactBernstein& lo, ExactBernstein& hi) {
    lo = f;
    hi = f;
    const int n = in_u ? f.p : f.q, m = in_u ? f.q : f.p;
    const num::Expansion half(0.5);
    for (int l = 0; l <= m; ++l) {
        std::vector<num::Expansion> row(static_cast<std::size_t>(n + 1));
        for (int k = 0; k <= n; ++k) row[static_cast<std::size_t>(k)] = in_u ? f.at(k, l) : f.at(l, k);
        std::vector<num::Expansion> a(static_cast<std::size_t>(n + 1)), b(static_cast<std::size_t>(n + 1));
        for (int r = 0; r <= n; ++r) {
            a[static_cast<std::size_t>(r)] = row[0];
            b[static_cast<std::size_t>(n - r)] = row[static_cast<std::size_t>(n - r)];
            for (int k = 0; k + r < n; ++k)
                row[static_cast<std::size_t>(k)] = (row[static_cast<std::size_t>(k)] + row[static_cast<std::size_t>(k + 1)]) * half;
        }
        for (int k = 0; k <= n; ++k) {
            (in_u ? lo.at(k, l) : lo.at(l, k)) = a[static_cast<std::size_t>(k)];
            (in_u ? hi.at(k, l) : hi.at(l, k)) = b[static_cast<std::size_t>(k)];
        }
    }
}

/// f over [a, b] in u (in_u) or v, re-expressed on [0, 1]: coefficient k is
/// the blossom of f at (a, ..., a, b, ..., b) with k copies of b. Exact for
/// dyadic a and b (each de Casteljau step multiplies by 1 - x and x).
inline ExactBernstein restrict_dir(const ExactBernstein& f, bool in_u, double a, double b) {
    ExactBernstein r = f;
    const int n = in_u ? f.p : f.q, m = in_u ? f.q : f.p;
    for (int l = 0; l <= m; ++l) {
        for (int k = 0; k <= n; ++k) {
            std::vector<num::Expansion> pts(static_cast<std::size_t>(n + 1));
            for (int i = 0; i <= n; ++i) pts[static_cast<std::size_t>(i)] = in_u ? f.at(i, l) : f.at(l, i);
            for (int s = 1; s <= n; ++s) {
                const double x = s <= k ? b : a;
                const num::Expansion one_minus = num::Expansion(1.0) - num::Expansion(x), xx(x);
                for (int i = 0; i + s <= n; ++i)
                    pts[static_cast<std::size_t>(i)] = one_minus * pts[static_cast<std::size_t>(i)] + xx * pts[static_cast<std::size_t>(i + 1)];
            }
            (in_u ? r.at(k, l) : r.at(l, k)) = pts[0];
        }
    }
    return r;
}

/// +1 or -1 if f has that strict sign on all of [0, 1]^2 (proven), else 0:
/// a zero, a sign change, or not decided within `depth` subdivisions.
inline int strict_sign(const ExactBernstein& f, int depth) {
    int lo_sign = 2, hi_sign = -2;
    for (const num::Expansion& x : f.c) {
        const int s = x.sign();
        lo_sign = s < lo_sign ? s : lo_sign;
        hi_sign = s > hi_sign ? s : hi_sign;
    }
    if (lo_sign == hi_sign && lo_sign != 0) return lo_sign;
    // A corner coefficient is the value there: a zero or opposite corner
    // signs settle it.
    const int c00 = f.at(0, 0).sign(), c10 = f.at(f.p, 0).sign(), c01 = f.at(0, f.q).sign(), c11 = f.at(f.p, f.q).sign();
    if (c00 == 0 || c10 == 0 || c01 == 0 || c11 == 0) return 0;
    if (c10 != c00 || c01 != c00 || c11 != c00) return 0;
    if (depth <= 0) return 0;
    // Split along the parameter of higher degree first; a curve only in u.
    const bool in_u = f.q == 0 || (f.p > 0 && depth % 2 == 0);
    ExactBernstein a, b;
    split_half(f, in_u, a, b);
    const int sa = strict_sign(a, depth - 1);
    if (sa == 0) return 0;
    return strict_sign(b, depth - 1) == sa ? sa : 0;
}

} // namespace contact_detail

/// g(u, v) = det(c00, c10, c01, X(u, v)) for the plane of the affine
/// parallelogram Y: exact Bernstein coefficients of degree (X.p, X.q).
inline ExactBernstein plane_function(const PolyPatchD& X, const PolyPatchD& Y) {
    using planar_detail::corner;
    const Vec4d &a = corner(Y, 0, 0), &b = corner(Y, 1, 0), &c = corner(Y, 0, 1);
    ExactBernstein g(X.p, X.q);
    for (int i = 0; i <= X.p; ++i)
        for (int j = 0; j <= X.q; ++j) g.at(i, j) = contact_detail::det4(a, b, c, X.at(i, j));
    return g;
}

/// How patch X meets the plane of the affine parallelogram Y along side s.
struct SideContact {
    bool only_side = false; ///< proven: X meets the plane exactly in side s, nowhere else
    int order = 0;          ///< order of contact along s (1 transversal, 2 or more tangential)
    int sign = 0;           ///< sign of g on the rest of X: the side of the plane X lies on
};

/// Decide exactly whether X meets Y's plane only along its side s. Requires
/// every control point of the side in the plane (side_in_plane); `depth`
/// bounds the subdivisions of the sign proof. With a dyadic `margin` > 0 the
/// proof covers the square enlarged by it on every side, [-m, 1 + m]^2: then
/// a point of X's extension that close to the square and in the plane lies on
/// the side's line as well (used to make enclosures of parameters exact).
inline SideContact side_contact(const PolyPatchD& X, int side, const PolyPatchD& Y, int depth = 12, double margin = 0.0) {
    SideContact r;
    if (!planar_affine(Y)) return r;
    for (const Vec4d& c : X.c)
        if (!(c.w > 0.0)) return r;
    ExactBernstein g = plane_function(X, Y);
    const bool across_v = side == 0 || side == 2;
    while (contact_detail::side_zero(g, side) && (across_v ? g.q : g.p) > 0) {
        g = contact_detail::divide_side(g, side);
        ++r.order;
    }
    if (r.order == 0) return r;
    if (margin > 0.0) {
        g = contact_detail::restrict_dir(g, true, -margin, 1.0 + margin);
        g = contact_detail::restrict_dir(g, false, -margin, 1.0 + margin);
    }
    r.sign = contact_detail::strict_sign(g, depth);
    r.only_side = r.sign != 0;
    return r;
}

/// How patch X meets the plane of Y along several of its sides at once (a
/// hill on a flat base touches the base's plane along all four).
struct SidesContact {
    bool only_sides = false; ///< proven: X meets the plane exactly in the sides of `mask`
    int mask = 0;            ///< bit s set: side s lies in the plane
    int sign = 0;            ///< sign of g on the rest of X
};

/// Divide out every side whose coefficients all vanish, as often as they do,
/// then prove the quotient's strict sign: the plane function then vanishes
/// on those sides and nowhere else.
inline SidesContact sides_contact(const PolyPatchD& X, const PolyPatchD& Y, int depth = 12) {
    SidesContact r;
    if (!planar_affine(Y)) return r;
    for (const Vec4d& c : X.c)
        if (!(c.w > 0.0)) return r;
    ExactBernstein g = plane_function(X, Y);
    for (bool again = true; again;) {
        again = false;
        for (int s = 0; s < 4; ++s) {
            const bool across_v = s == 0 || s == 2;
            if ((across_v ? g.q : g.p) == 0 || !contact_detail::side_zero(g, s)) continue;
            g = contact_detail::divide_side(g, s);
            r.mask |= 1 << s;
            again = true;
        }
    }
    if (r.mask == 0) return r;
    r.sign = contact_detail::strict_sign(g, depth);
    r.only_sides = r.sign != 0;
    return r;
}

/// sides_contact for a flat Y that need not be an affine parallelogram (a
/// sector's cap on its rational carrier): the caller has proven every control
/// point of Y in one plane. det(c00, c10, c01, X) in homogeneous coordinates
/// is that plane's function whatever the weights; if the three corners do not
/// span the plane it vanishes identically, no strict sign is found, and
/// nothing is claimed.
inline SidesContact sides_contact_flat(const PolyPatchD& X, const PolyPatchD& Y, int depth = 12) {
    SidesContact r;
    for (const Vec4d& c : X.c)
        if (!(c.w > 0.0)) return r;
    for (const Vec4d& c : Y.c)
        if (!(c.w > 0.0)) return r;
    ExactBernstein g = plane_function(X, Y);
    for (bool again = true; again;) {
        again = false;
        for (int s = 0; s < 4; ++s) {
            const bool across_v = s == 0 || s == 2;
            if ((across_v ? g.q : g.p) == 0 || !contact_detail::side_zero(g, s)) continue;
            g = contact_detail::divide_side(g, s);
            r.mask |= 1 << s;
            again = true;
        }
    }
    if (r.mask == 0) return r;
    r.sign = contact_detail::strict_sign(g, depth);
    r.only_sides = r.sign != 0;
    return r;
}

/// Contact along a line inside the patch. For a patch X of degree 2 in one
/// parameter (a cylinder or cone sector in the angle), the plane function on
/// each control row in that parameter, c0 (1 - u)^2 + 2 c1 u (1 - u) +
/// c2 u^2, has a double root exactly when c1^2 = c0 c2; it is then
/// K (u - r)^2 with K = c0 - 2 c1 + c2 and r = (c0 - c1) / K. If every row
/// has a double root at the same r (decided exactly, by cross-multiplication),
/// then g(u, v) = (u - r)^2 h(v) with h(v) = sum_j K_j B_j(v), and if every
/// K_j has one strict sign, h does too: X meets the plane exactly on the line
/// u = r of its parameters (for every v, and for u beyond the square as well)
/// and lies strictly on one side of the plane everywhere else. A plane
/// tangent to a cylinder along one of its rulings, inside the patch.
struct PlaneLineContact {
    bool found = false;
    bool fixed_u = true;  ///< the line is u = r (else v = r)
    num::Interval r;      ///< its parameter, enclosed
    int sign = 0;         ///< sign of g off the line
};

inline PlaneLineContact line_contact(const PolyPatchD& X, const PolyPatchD& Y) {
    PlaneLineContact out;
    if (!planar_affine(Y)) return out;
    for (const Vec4d& c : X.c)
        if (!(c.w > 0.0)) return out;
    const ExactBernstein g = plane_function(X, Y);
    using E = num::Expansion;
    for (int dir = 0; dir < 2; ++dir) { // 0: along u (rows j), 1: along v (rows i)
        if ((dir == 0 ? g.p : g.q) != 2) continue;
        const int rows = dir == 0 ? g.q + 1 : g.p + 1;
        bool ok = true;
        int sign = 0, ref = -1;
        E ref_num, ref_den;
        for (int j = 0; j < rows && ok; ++j) {
            const E& c0 = dir == 0 ? g.at(0, j) : g.at(j, 0);
            const E& c1 = dir == 0 ? g.at(1, j) : g.at(j, 1);
            const E& c2 = dir == 0 ? g.at(2, j) : g.at(j, 2);
            if ((c1 * c1 - c0 * c2).sign() != 0) { ok = false; break; }
            const E K = c0 - E(2.0) * c1 + c2;
            const int s = K.sign();
            if (s == 0 || (sign != 0 && s != sign)) { ok = false; break; } // h must be strictly signed
            sign = s;
            const E num = c0 - c1;
            if (ref < 0) {
                ref = j;
                ref_num = num;
                ref_den = K;
            } else if ((num * ref_den - ref_num * K).sign() != 0) {
                ok = false; // a different double root
            }
        }
        if (!ok || ref < 0) continue;
        // r = ref_num / ref_den, enclosed: each expansion is the exact sum of
        // its components, summed here in interval arithmetic.
        auto enclose = [](const E& e) {
            num::Interval s(0.0);
            for (double x : e.components()) s = s + num::Interval(x);
            return s;
        };
        out.found = true;
        out.fixed_u = dir == 0;
        out.r = enclose(ref_num) / enclose(ref_den);
        out.sign = sign;
        return out;
    }
    return out;
}

/// Where a curve (a patch of degree (n, 0)) meets the plane of Y, decided
/// exactly when it meets it at most at its two ends.
struct EndContact {
    bool decided = false; ///< proven: the curve meets the plane at most at t = 0 and t = 1
    int order0 = 0, order1 = 0; ///< order of contact at t = 0 and t = 1 (0: not on the plane)
};

inline EndContact end_contact(const PolyPatchD& C, const PolyPatchD& Y, int depth = 16) {
    EndContact r;
    if (C.q != 0 || !planar_affine(Y)) return r;
    for (const Vec4d& c : C.c)
        if (!(c.w > 0.0)) return r;
    ExactBernstein g = plane_function(C, Y);
    while (g.p > 0 && g.at(0, 0).sign() == 0) { g = contact_detail::divide_side(g, 3); ++r.order0; }
    while (g.p > 0 && g.at(g.p, 0).sign() == 0) { g = contact_detail::divide_side(g, 1); ++r.order1; }
    r.decided = contact_detail::strict_sign(g, depth) != 0;
    return r;
}

} // namespace ak::op
