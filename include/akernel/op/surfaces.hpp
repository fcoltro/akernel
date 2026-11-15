// akernel — exact quadrics, and coincident curved surfaces decided exactly.
//
// Quadrics. A quadric is the zero set of x^T A x for a symmetric 4x4 matrix A
// in homogeneous coordinates. A rational patch with homogeneous control
// points P_ij lies on it exactly when the polynomial A(P(u,v), P(u,v)), of
// bidegree (2p, 2q), vanishes identically; its Bernstein coefficients are
//   sum over i+i'=k, j+j'=l of C(p,i)C(p,i')C(q,j)C(q,j') P_ij^T A P_i'j'
// (up to a positive factor), and they are computed here with exact expansion
// arithmetic. A patch whose weights were rounded (the usual quarter circle
// with weight sqrt(2)/2) lies on a nearby conic, not on the circle, and the
// test says so.
//
// Exact circles. The circle has the rational parameterisation
// (1 - t^2, 2t, 1 + t^2) / (1 + t^2) for all real t, whose Bezier data on any
// interval [t0, t1] with dyadic ends are exact in doubles. Every cylinder
// sector made by make_cylinder_sector is a piece of that one parameterisation
// of its cylinder, so two sectors of the same cylinder are related by an
// affine change of parameters.
//
// Coincident surfaces. Two patches P and Q are the same surface, with
// Q(u, v) = P(a u + b, c v + d), when P restricted to [b, a+b] x [d, c+d]
// equals Q up to a constant homogeneous factor. Candidate a, b, c, d come
// from locating Q's corners on P in doubles and snapping them to a dyadic
// grid; the equality is then checked exactly (blossoms in expansion
// arithmetic, cross-multiplied). A candidate that is slightly off, or a
// relation that is not dyadic, fails the check: the surfaces are then not
// treated as one, which is safe. Planes keep their own exact test
// (op/planar.hpp), which needs no common parameterisation.
#pragma once

#include <array>
#include <cmath>
#include <vector>

#include <akernel/geom/bernstein.hpp>
#include <akernel/num/expansion.hpp>
#include <akernel/op/patch_solid.hpp>
#include <akernel/op/planar.hpp>
#include <akernel/op/trim.hpp>

namespace ak::op {

using num::Expansion;

// ------------------------------------------------------------------ quadrics

/// The quadric x^T A x = 0 in homogeneous coordinates (x, y, z, w).
///
/// Every factory writes all sixteen entries from one table. An earlier form
/// that zero-initialised the matrix (a default member initialiser) and then
/// assigned the nonzero entries was miscompiled by MSVC 19.51 at /O2 with
/// /fp:strict when inlined into a large function: one entry kept a stale 1.0
/// (docs/RESULTS.md, "A compiler fault"). tests/op/quadric_test.cpp checks
/// the matrices the factories return.
struct Quadric {
    double A[4][4];
    /// The symmetric matrix with these upper-triangle entries.
    static Quadric symmetric(double a00, double a01, double a02, double a03, double a11, double a12, double a13,
                             double a22, double a23, double a33) {
        const double t[4][4] = {{a00, a01, a02, a03}, {a01, a11, a12, a13}, {a02, a12, a22, a23}, {a03, a13, a23, a33}};
        Quadric q;
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j) q.A[i][j] = t[i][j];
        return q;
    }
    /// x^2 + y^2 - r^2 w^2 = 0 shifted to the axis through (cx, cy), parallel to z.
    static Quadric cylinder_z(double cx, double cy, double r) {
        return symmetric(1.0, 0.0, 0.0, -cx, 1.0, 0.0, -cy, 0.0, 0.0, cx * cx + cy * cy - r * r);
    }
    /// |x - c|^2 = r^2.
    static Quadric sphere(double cx, double cy, double cz, double r) {
        return symmetric(1.0, 0.0, 0.0, -cx, 1.0, 0.0, -cy, 1.0, -cz, cx * cx + cy * cy + cz * cz - r * r);
    }
    /// x^2 + y^2 = (k z)^2: a cone with apex at the origin around the z axis.
    static Quadric cone_z(double k) { return symmetric(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, -k * k, 0.0, 0.0); }
};

namespace surfaces_detail {

using E4 = std::array<Expansion, 4>;

inline E4 to_e4(const Vec4d& c) { return {Expansion(c.x), Expansion(c.y), Expansion(c.z), Expansion(c.w)}; }

/// a^T A b, exactly.
inline Expansion form(const Quadric& Q, const Vec4d& a, const Vec4d& b) {
    const double av[4] = {a.x, a.y, a.z, a.w}, bv[4] = {b.x, b.y, b.z, b.w};
    Expansion s;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            if (Q.A[i][j] != 0.0) s = s + Expansion(Q.A[i][j]) * Expansion(av[i]) * Expansion(bv[j]);
    return s;
}

inline double binom(int n, int k) {
    double c = 1.0;
    for (int i = 1; i <= k; ++i) c = c * (n - k + i) / i;
    return c;
}

/// Blossom of a Bezier row (degree n) at (x1..xn): (n-k) copies of x0 and k of x1.
inline E4 blossom_row(const std::vector<E4>& row, double x0, double x1, int k) {
    std::vector<E4> pts = row;
    const int n = static_cast<int>(row.size()) - 1;
    for (int r = 1; r <= n; ++r) {
        const double x = r <= k ? x1 : x0;
        const Expansion a = Expansion(1.0) - Expansion(x), b = Expansion(x);
        for (int i = 0; i + r <= n; ++i)
            for (int c = 0; c < 4; ++c) pts[i][c] = a * pts[i][c] + b * pts[i + 1][c];
    }
    return pts[0];
}

/// P restricted to [u0, u1] x [v0, v1], exactly (all four values dyadic).
inline std::vector<E4> restrict_exact(const PolyPatchD& P, double u0, double u1, double v0, double v1) {
    const int p = P.p, q = P.q;
    std::vector<E4> mid(static_cast<std::size_t>((p + 1) * (q + 1)));
    for (int j = 0; j <= q; ++j) {
        std::vector<E4> row;
        for (int i = 0; i <= p; ++i) row.push_back(to_e4(P.at(i, j)));
        for (int k = 0; k <= p; ++k) mid[static_cast<std::size_t>(k * (q + 1) + j)] = blossom_row(row, u0, u1, k);
    }
    std::vector<E4> out(mid.size());
    for (int i = 0; i <= p; ++i) {
        std::vector<E4> col;
        for (int j = 0; j <= q; ++j) col.push_back(mid[static_cast<std::size_t>(i * (q + 1) + j)]);
        for (int l = 0; l <= q; ++l) out[static_cast<std::size_t>(i * (q + 1) + l)] = blossom_row(col, v0, v1, l);
    }
    return out;
}

/// Parameters of a point on P, by Gauss-Newton from a grid of starts, over
/// the parameter range [-2, 3]^2. A rational patch extended beyond its square
/// need not be one-to-one (an annulus reaches the same points again with a
/// negative radius), so among the starts that converge onto the point the
/// one nearest the parameter square is taken. False if none converges.
inline bool locate_on_patch(const geom::BezierSurface& S, const Vec3d& x, double scale, double& u, double& v) {
    const std::size_t iu = static_cast<std::size_t>(geom::jet_index(1, 0)), iv = static_cast<std::size_t>(geom::jet_index(0, 1));
    const double tol = 1e-9 * scale;
    double best_out = 1e300, best_res = 1e300;
    bool found = false;
    auto outside = [](double a, double b) {
        return std::max(std::max(0.0, -a), std::max(0.0, a - 1.0)) + std::max(std::max(0.0, -b), std::max(0.0, b - 1.0));
    };
    for (int a = 0; a <= 10; ++a)
        for (int b = 0; b <= 10; ++b) {
            double uu = -2.0 + 0.5 * a, vv = -2.0 + 0.5 * b;
            for (int it = 0; it < 40; ++it) {
                const geom::Jet3 J = S.jet(uu, vv);
                const Vec3d r = x - J[0];
                const double A = dot(J[iu], J[iu]), B = dot(J[iu], J[iv]), C = dot(J[iv], J[iv]);
                const double g1 = dot(J[iu], r), g2 = dot(J[iv], r), det = A * C - B * B;
                if (!(det > 0.0) || !std::isfinite(det)) break;
                const double du = (C * g1 - B * g2) / det, dv = (A * g2 - B * g1) / det;
                uu += du;
                vv += dv;
                if (std::fabs(uu) > 10.0 || std::fabs(vv) > 10.0) break;
                if (std::fabs(du) + std::fabs(dv) < 1e-15) break;
            }
            const double res = norm(S.point(uu, vv) - x);
            if (!std::isfinite(res) || res > tol) continue;
            const double out = outside(uu, vv);
            if (out < best_out - 1e-12 || (out <= best_out + 1e-12 && res < best_res)) {
                best_out = out;
                best_res = res;
                u = uu;
                v = vv;
                found = true;
            }
        }
    return found;
}

/// x snapped to the dyadic grid of step 2^-24, if it is that close to it.
inline bool snap(double x, double& s) {
    const double g = std::ldexp(1.0, 24);
    s = std::nearbyint(x * g) / g;
    return std::fabs(x - s) <= 1e-9;
}

} // namespace surfaces_detail

inline bool exact_point_at(const PolyPatchD& P, double au, double bu, double av, double bv, Vec4d& out);

/// The point of P at the parameters (u, v), exactly: homogeneous coordinates
/// computed by de Casteljau steps in expansion arithmetic. False if one of
/// them is not a double (the point is then exact but not representable here)
/// or the weight is not positive. Parameters that are doubles with few bits
/// (dyadic numbers) on control points with few bits give doubles.
inline bool exact_point_at(const PolyPatchD& P, double u, double v, Vec4d& out) {
    return exact_point_at(P, u, 1.0, v, 1.0, out);
}

/// A fraction a/b with 0 < b <= 4096 within 1e-10 of x, found by continued
/// fractions: the candidate for an exact value that is rational but not
/// dyadic (a third, say). Whether it is the exact value is for the caller
/// to prove.
inline bool rational_near(double x, double& a, double& b) {
    double h0 = 0.0, h1 = 1.0, k0 = 1.0, k1 = 0.0, y = x;
    for (int i = 0; i < 24; ++i) {
        const double q = std::floor(y);
        const double h2 = q * h1 + h0, k2 = q * k1 + k0;
        if (!(k2 <= 4096.0)) break;
        h0 = h1; h1 = h2; k0 = k1; k1 = k2;
        if (std::fabs(x - h1 / k1) <= 1e-10) { a = h1; b = k1; return true; }
        const double f = y - q;
        if (f < 1e-12) break;
        y = 1.0 / f;
    }
    return false;
}

/// The same at rational parameters (au / bu, av / bv) with bu, bv > 0: each
/// de Casteljau step is (b - a) x_i + a x_{i+1}, so the result is the point's
/// homogeneous coordinates times bu^p bv^q, which is the same point.
inline bool exact_point_at(const PolyPatchD& P, double au, double bu, double av, double bv, Vec4d& out) {
    using namespace surfaces_detail;
    if (!(bu > 0.0) || !(bv > 0.0)) return false;
    auto at = [](std::vector<E4> pts, double a, double b) {
        const int n = static_cast<int>(pts.size()) - 1;
        const Expansion ea(a), eb = Expansion(b) - Expansion(a);
        for (int r = 1; r <= n; ++r)
            for (int i = 0; i + r <= n; ++i)
                for (std::size_t c = 0; c < 4; ++c)
                    pts[static_cast<std::size_t>(i)][c] = eb * pts[static_cast<std::size_t>(i)][c] + ea * pts[static_cast<std::size_t>(i + 1)][c];
        return pts[0];
    };
    std::vector<E4> col;
    for (int j = 0; j <= P.q; ++j) {
        std::vector<E4> row;
        for (int i = 0; i <= P.p; ++i) row.push_back(to_e4(P.at(i, j)));
        col.push_back(at(std::move(row), au, bu));
    }
    const E4 h = at(std::move(col), av, bv);
    double d[4];
    for (int c = 0; c < 4; ++c) {
        d[c] = h[static_cast<std::size_t>(c)].estimate();
        if (!std::isfinite(d[c]) || (h[static_cast<std::size_t>(c)] - Expansion(d[c])).sign() != 0) return false;
    }
    out = {d[0], d[1], d[2], d[3]};
    return out.w > 0.0;
}

/// x rounded to the nearest multiple of 2^-30: the candidate for an exact
/// value where x is the midpoint of an enclosure. Enclosures of certified
/// roots are far narrower than that grid, so at most one candidate lies in
/// one. Whether the candidate is the exact value is for the caller to prove.
inline double snap_dyadic(double x) { return std::ldexp(std::nearbyint(std::ldexp(x, 30)), -30); }

/// Does the patch lie exactly on the quadric?
inline bool on_quadric(const PolyPatchD& P, const Quadric& Q) {
    using namespace surfaces_detail;
    const int p = P.p, q = P.q;
    for (int k = 0; k <= 2 * p; ++k)
        for (int l = 0; l <= 2 * q; ++l) {
            Expansion c;
            for (int i = std::max(0, k - p); i <= std::min(p, k); ++i)
                for (int j = std::max(0, l - q); j <= std::min(q, l); ++j) {
                    const double w = binom(p, i) * binom(p, k - i) * binom(q, j) * binom(q, l - j);
                    c = c + Expansion(w) * form(Q, P.at(i, j), P.at(k - i, l - j));
                }
            if (c.sign() != 0) return false;
        }
    return true;
}

/// Is Q the patch P under an affine change of parameters, Q(u,v) = P(a u + b,
/// c v + d)? Proven exactly (see the header); `q_to_p` maps Q's parameters to
/// P's, with exact (dyadic) coefficients.
inline bool affine_reparam(const PolyPatchD& P, const PolyPatchD& Q, Affine& q_to_p) {
    using namespace surfaces_detail;
    if (P.p != Q.p || P.q != Q.q || P.c.empty()) return false;
    const geom::BezierSurface S(P);
    auto corner = [&](int i, int j) {
        const Vec4d h = Q.at(i * Q.p, j * Q.q);
        return Vec3d{h.x / h.w, h.y / h.w, h.z / h.w};
    };
    double scale = 1.0;
    for (const Vec4d& c : P.c) scale = std::max(scale, std::fabs(c.x / c.w) + std::fabs(c.y / c.w) + std::fabs(c.z / c.w));
    double u00, v00, u11, v11;
    if (!locate_on_patch(S, corner(0, 0), scale, u00, v00) || !locate_on_patch(S, corner(1, 1), scale, u11, v11))
        return false;
    double b, e, d, f; // [b, e] x [d, f] in P's parameters
    if (!snap(u00, b) || !snap(u11, e) || !snap(v00, d) || !snap(v11, f)) return false;
    if (b == e || d == f) return false;
    // The scales must be powers of two, so that the map and its inverse are
    // both exact and exact parameters (vertices at corners and on sides) stay
    // exact in either direction. A relation such as 2/3 is not accepted: the
    // surfaces are then not treated as one, which is reported, not guessed.
    auto pow2 = [](double s) {
        s = std::fabs(s);
        int ex = 0;
        return s > 0.0 && std::isfinite(s) && std::frexp(s, &ex) == 0.5;
    };
    if (!pow2(e - b) || !pow2(f - d)) return false;
    const std::vector<E4> R = restrict_exact(P, b, e, d, f);
    const E4 r0 = R[0], q0 = to_e4(Q.c[0]);
    if (r0[3].sign() == 0 || q0[3].sign() == 0) return false;
    for (std::size_t n = 0; n < Q.c.size(); ++n) {
        const E4 qn = to_e4(Q.c[n]);
        for (int c = 0; c < 4; ++c)
            if ((qn[c] * r0[3] - R[n][c] * q0[3]).sign() != 0) return false;
    }
    // Coefficients as exact points where the differences are doubles (they
    // are, for dyadic ranges), so that exact parameters stay exact.
    auto diff = [](double x, double y) {
        double r = 0.0;
        return exact_value(Expansion(x) - Expansion(y), r) ? Interval(r) : Interval(x) - Interval(y);
    };
    q_to_p = Affine{};
    q_to_p.a[0][0] = Interval(b);
    q_to_p.a[0][1] = diff(e, b);
    q_to_p.a[0][2] = Interval(0.0);
    q_to_p.a[1][0] = Interval(d);
    q_to_p.a[1][1] = Interval(0.0);
    q_to_p.a[1][2] = diff(f, d);
    return true;
}

namespace surfaces_detail {
/// 1 / s, as a point when that is exact.
inline Interval reciprocal(const Interval& s) {
    if (is_point(s) && s.lo != 0.0) {
        const double r = 1.0 / s.lo;
        if ((Expansion(r) * Expansion(s.lo) - Expansion(1.0)).sign() == 0) return Interval(r);
    }
    return Interval(1.0) / s;
}
/// -c * r, as a point when that is exact.
inline Interval neg_product(const Interval& c, const Interval& r) {
    double x = 0.0;
    if (is_point(c) && is_point(r) && exact_value(-(Expansion(c.lo) * Expansion(r.lo)), x)) return Interval(x);
    return -c * r;
}
} // namespace surfaces_detail

/// The inverse of an axis-aligned affine map (u' = a + b u, v' = c + d v), or
/// of one that also exchanges the parameters (u' = a + b v, v' = c + d u).
inline Affine invert_axis_affine(const Affine& m) {
    using namespace surfaces_detail;
    Affine r;
    const bool swapped = m.a[0][1].lo == 0.0 && m.a[0][1].hi == 0.0 && m.a[1][2].lo == 0.0 && m.a[1][2].hi == 0.0;
    if (swapped) {
        // u = (v' - c) / d, v = (u' - a) / b.
        r.a[0][2] = reciprocal(m.a[1][1]);
        r.a[0][0] = neg_product(m.a[1][0], r.a[0][2]);
        r.a[0][1] = Interval(0.0);
        r.a[1][1] = reciprocal(m.a[0][2]);
        r.a[1][0] = neg_product(m.a[0][0], r.a[1][1]);
        r.a[1][2] = Interval(0.0);
        return r;
    }
    r.a[0][1] = reciprocal(m.a[0][1]);
    r.a[0][0] = neg_product(m.a[0][0], r.a[0][1]);
    r.a[0][2] = Interval(0.0);
    r.a[1][2] = reciprocal(m.a[1][2]);
    r.a[1][0] = neg_product(m.a[1][0], r.a[1][2]);
    r.a[1][1] = Interval(0.0);
    return r;
}

/// The patch with its parameters exchanged: T(u, v) = Q(v, u).
inline PolyPatchD transposed(const PolyPatchD& Q) {
    PolyPatchD T(Q.q, Q.p);
    for (int i = 0; i <= Q.p; ++i)
        for (int j = 0; j <= Q.q; ++j) T.at(j, i) = Q.at(i, j);
    return T;
}

/// affine_reparam, also with Q's parameters exchanged: Q(u, v) = P(a v + b,
/// c u + d). Two faces of one cylinder made as the inner and the outer face
/// of different solids are parameterised that way round.
inline bool affine_reparam_any(const PolyPatchD& P, const PolyPatchD& Q, Affine& q_to_p) {
    if (affine_reparam(P, Q, q_to_p)) return true;
    Affine t;
    if (!affine_reparam(P, transposed(Q), t)) return false;
    // t maps T's (u, v) = Q's (v, u): Q's (u, v) goes to t applied to (v, u).
    q_to_p.a[0][0] = t.a[0][0];
    q_to_p.a[0][1] = t.a[0][2];
    q_to_p.a[0][2] = t.a[0][1];
    q_to_p.a[1][0] = t.a[1][0];
    q_to_p.a[1][1] = t.a[1][2];
    q_to_p.a[1][2] = t.a[1][1];
    return true;
}

/// Are `from` and `to` the same surface (exactly), and if so the map from
/// `from`'s parameters to `to`'s: flat parallelograms in one plane, or one
/// patch an affine reparameterisation of the other (possibly with the
/// parameters exchanged).
inline bool same_surface_map(const PolyPatchD& from, const PolyPatchD& to, Affine& m) {
    if (planar_affine(from) && planar_affine(to)) return coplanar(from, to) && chart_map(from, to, m);
    Affine q_to_p;
    if (affine_reparam(to, from, q_to_p)) { m = q_to_p; return true; }
    if (affine_reparam(from, to, q_to_p)) { m = invert_axis_affine(q_to_p); return true; }
    if (affine_reparam_any(to, from, q_to_p)) { m = q_to_p; return true; }
    if (affine_reparam_any(from, to, q_to_p)) { m = invert_axis_affine(q_to_p); return true; }
    return false;
}

// ------------------------------------------------------- exact cylinder sectors

/// The sector r0 <= rho <= r1, z0 <= z <= z1 of the cylinders around the z
/// axis, over the circle parameter range [t0, t1] of (1 - t^2, 2t, 1 + t^2)
/// (t = 0 at angle 0, t = 1 at 90 degrees, t = -1 at -90 degrees). With dyadic
/// t0, t1, r0, r1, z0, z1 of moderate size every control point is exact, so
/// the faces lie exactly on their cylinders and planes, and sectors of one
/// cylinder are affine reparameterisations of each other.
inline PatchSolid make_cylinder_sector(double r0, double r1, double t0, double t1, double z0, double z1) {
    Lattice L(1, 2, 1);
    const double rad[2] = {r0, r1}, z[2] = {z0, z1};
    // Blossom of the circle's homogeneous form at (x, y): (1 - xy, x + y, 1 + xy).
    auto circ = [](double x, double y, double rho, double h) {
        const double w = 1.0 + x * y;
        return Vec4d{rho * (1.0 - x * y), rho * (x + y), h * w, w};
    };
    const double pairs[3][2] = {{t0, t0}, {t0, t1}, {t1, t1}};
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 2; ++k) L.at(i, j, k) = circ(pairs[j][0], pairs[j][1], rad[i], z[k]);
    return make_hex(L);
}

} // namespace ak::op
