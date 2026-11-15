// akernel — exact facts about planar faces.
//
// Two surfaces can coincide without being the same by construction: two boxes
// placed face to face both have a face in the plane they touch along, each
// with its own parameterisation and extent. Coincidence of that kind cannot
// be decided by intervals (a coincidence is destroyed by any perturbation);
// it needs exact arithmetic on the input data. This file decides it exactly
// for the simplest and most common case, parallelogram patches:
//   - a patch is an affine parallelogram if its degree is (1,1), its weights
//     are 1, and c00 + c11 = c10 + c01 exactly;
//   - two such patches lie in the same plane if the corners of one are
//     exactly coplanar with three corners of the other (exact orient3d);
//   - the map between their parameter planes is affine; it is computed with
//     interval coefficients that enclose the exact map.
// A patch whose corners were rounded off a plane (a rotated box built in
// doubles) is, correctly, not planar, and is not merged with anything.
#pragma once

#include <array>
#include <cmath>
#include <vector>

#include <akernel/geom/bernstein.hpp>
#include <akernel/num/expansion.hpp>
#include <akernel/num/predicates.hpp>
#include <akernel/op/patch_solid.hpp>
#include <akernel/op/trim.hpp>

namespace ak::op {

namespace planar_detail {

inline std::array<double, 3> xyz(const Vec4d& c) { return {c.x, c.y, c.z}; }

/// Exact: a - b - c + d == 0 for one coordinate.
inline bool exact_zero4(double a, double b, double c, double d) {
    return (num::Expansion(a) - num::Expansion(b) - num::Expansion(c) + num::Expansion(d)).sign() == 0;
}

/// Exact 3x3 determinant of the xyz parts of three homogeneous points.
inline num::Expansion det3(const Vec4d& a, const Vec4d& b, const Vec4d& c) {
    using E = num::Expansion;
    return E(a.x) * (E(b.y) * E(c.z) - E(b.z) * E(c.y)) - E(a.y) * (E(b.x) * E(c.z) - E(b.z) * E(c.x)) +
           E(a.z) * (E(b.x) * E(c.y) - E(b.y) * E(c.x));
}

/// Are four homogeneous points in one plane? The 4x4 determinant of their
/// coordinates (x, y, z, w) vanishes, exactly; no division by the weights.
inline bool coplanar4(const Vec4d& a, const Vec4d& b, const Vec4d& c, const Vec4d& d) {
    using E = num::Expansion;
    const E det = E(d.w) * det3(a, b, c) - E(c.w) * det3(a, b, d) + E(b.w) * det3(a, c, d) - E(a.w) * det3(b, c, d);
    return det.sign() == 0;
}

/// The point of a homogeneous control point, in doubles (approximate unless
/// the weight is a power of two). For estimates only.
inline Vec4d affine(const Vec4d& c) { return {c.x / c.w, c.y / c.w, c.z / c.w, 1.0}; }

/// Corner control point (i, j in {0, 1}) of a patch of any degree.
inline const Vec4d& corner(const PolyPatchD& P, int i, int j) { return P.at(i ? P.p : 0, j ? P.q : 0); }

inline bool on_plane(const PolyPatchD& P, const Vec4d& q) {
    return coplanar4(corner(P, 0, 0), corner(P, 1, 0), corner(P, 0, 1), q);
}

/// Do two homogeneous points denote the same point? Exact cross-multiplication.
inline bool same_point_h(const Vec4d& a, const Vec4d& b) {
    using E = num::Expansion;
    if (!(a.w > 0.0) || !(b.w > 0.0)) return false;
    return (E(a.x) * E(b.w) - E(b.x) * E(a.w)).sign() == 0 && (E(a.y) * E(b.w) - E(b.y) * E(a.w)).sign() == 0 &&
           (E(a.z) * E(b.w) - E(b.z) * E(a.w)).sign() == 0;
}

} // namespace planar_detail

/// Corner (i, j) of a parallelogram patch as an affine point, in doubles
/// (exact when the weight is a power of two). For estimates only.
inline Vec4d corner_point(const PolyPatchD& P, int i, int j) { return planar_detail::affine(planar_detail::corner(P, i, j)); }

/// An affine map of the square onto a parallelogram: one positive weight on
/// every control point, the corners an exact parallelogram (c00 + c11 = c10 +
/// c01 in homogeneous coordinates, which with a common weight is the same as
/// for the points), and, for degree above (1,1), every control point exactly
/// where degree elevation of that parallelogram puts it:
///   p q c_ij = (p-i)(q-j) c00 + i(q-j) c_p0 + (p-i) j c_0q + i j c_pq,
/// checked in exact arithmetic (the factors are small integers). Such a patch
/// is its corners' parallelogram; everything planar uses only the corners.
inline bool planar_affine(const PolyPatchD& P) {
    if (P.p < 1 || P.q < 1 || P.c.empty()) return false;
    const double w = P.c[0].w;
    if (!(w > 0.0) || !std::isfinite(w)) return false;
    for (const Vec4d& c : P.c)
        if (c.w != w) return false;
    using planar_detail::corner;
    const Vec4d &c00 = corner(P, 0, 0), &c10 = corner(P, 1, 0), &c01 = corner(P, 0, 1), &c11 = corner(P, 1, 1);
    if (!(planar_detail::exact_zero4(c11.x, c10.x, c01.x, c00.x) && planar_detail::exact_zero4(c11.y, c10.y, c01.y, c00.y) &&
          planar_detail::exact_zero4(c11.z, c10.z, c01.z, c00.z)))
        return false;
    if (P.p == 1 && P.q == 1) return true;
    using E = num::Expansion;
    const int p = P.p, q = P.q;
    for (int i = 0; i <= p; ++i)
        for (int j = 0; j <= q; ++j) {
            const Vec4d& c = P.at(i, j);
            const double a00 = (p - i) * (q - j), a10 = i * (q - j), a01 = (p - i) * j, a11 = i * j, pq = p * q;
            const double cv[3] = {c.x, c.y, c.z}, v00[3] = {c00.x, c00.y, c00.z}, v10[3] = {c10.x, c10.y, c10.z},
                         v01[3] = {c01.x, c01.y, c01.z}, v11[3] = {c11.x, c11.y, c11.z};
            for (int k = 0; k < 3; ++k) {
                const E lhs = E(pq) * E(cv[k]);
                const E rhs = E(a00) * E(v00[k]) + E(a10) * E(v10[k]) + E(a01) * E(v01[k]) + E(a11) * E(v11[k]);
                if ((lhs - rhs).sign() != 0) return false;
            }
        }
    return true;
}

/// Two affine parallelograms in the same plane (exact).
inline bool coplanar(const PolyPatchD& A, const PolyPatchD& B) {
    if (!planar_affine(A) || !planar_affine(B)) return false;
    for (const Vec4d& q : B.c)
        if (!planar_detail::on_plane(A, q)) return false;
    return true;
}

/// Does side `side` of patch X lie exactly in the plane of the affine
/// parallelogram Y? (Every control point of the side, a positive weight, in
/// the plane: then the whole rational side curve is.)
inline bool side_in_plane(const PolyPatchD& X, int side, const PolyPatchD& Y) {
    if (!planar_affine(Y)) return false;
    for (const Vec4d& c : side_control(X, side))
        if (!(c.w > 0.0) || !planar_detail::on_plane(Y, c)) return false;
    return true;
}

/// Does all of patch X lie exactly in a plane that contains Y, where Y is any
/// patch whose control points are coplanar? (Every control point of X, with
/// a positive weight, in the plane of Y's first three non-collinear control
/// points, and all of Y's in it as well.) Used to name the case of flat faces
/// in one plane whose parameterisations are not affinely related.
inline bool same_plane_patches(const PolyPatchD& X, const PolyPatchD& Y) {
    const std::vector<Vec4d>& c = Y.c;
    // Three homogeneous points span a plane when they are linearly
    // independent: some coordinate unit vector is not coplanar with them.
    const Vec4d units[4] = {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};
    for (std::size_t i = 0; i < c.size(); ++i)
        for (std::size_t j = i + 1; j < c.size(); ++j)
            for (std::size_t k = j + 1; k < c.size(); ++k) {
                bool spans = false;
                for (const Vec4d& e : units) spans = spans || !planar_detail::coplanar4(c[i], c[j], c[k], e);
                if (!spans) continue;
                for (const Vec4d& p : Y.c)
                    if (!planar_detail::coplanar4(c[i], c[j], c[k], p)) return false;
                for (const Vec4d& p : X.c)
                    if (!(p.w > 0.0) || !planar_detail::coplanar4(c[i], c[j], c[k], p)) return false;
                return true;
            }
    return false;
}

/// Three control points of a flat patch Y that span its plane, with every
/// control point of Y in that plane (exact). False if Y is not flat.
inline bool flat_triple(const PolyPatchD& Y, Vec4d t[3]) {
    const std::vector<Vec4d>& c = Y.c;
    const Vec4d units[4] = {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};
    for (std::size_t i = 0; i < c.size(); ++i)
        for (std::size_t j = i + 1; j < c.size(); ++j)
            for (std::size_t k = j + 1; k < c.size(); ++k) {
                bool spans = false;
                for (const Vec4d& e : units) spans = spans || !planar_detail::coplanar4(c[i], c[j], c[k], e);
                if (!spans) continue;
                for (const Vec4d& p : c)
                    if (!planar_detail::coplanar4(c[i], c[j], c[k], p)) return false;
                t[0] = c[i];
                t[1] = c[j];
                t[2] = c[k];
                return true;
            }
    return false;
}

/// Are the points exactly on one line? (3D collinearity by three exact
/// projections.)
inline bool collinear(const Vec4d& a, const Vec4d& b, const Vec4d& c) {
    const double ax[2] = {a.x, a.y}, bx[2] = {b.x, b.y}, cx[2] = {c.x, c.y};
    const double ay[2] = {a.y, a.z}, by[2] = {b.y, b.z}, cy[2] = {c.y, c.z};
    const double az[2] = {a.x, a.z}, bz[2] = {b.x, b.z}, cz[2] = {c.x, c.z};
    return num::orient2d(ax, bx, cx) == num::Sign::Zero && num::orient2d(ay, by, cy) == num::Sign::Zero &&
           num::orient2d(az, bz, cz) == num::Sign::Zero;
}

/// The map from the parameters of P to those of Q, for affine parallelograms
/// in the same plane: Q^-1 o P, with interval coefficients. Returns false if
/// Q is degenerate.
inline bool chart_map(const PolyPatchD& P, const PolyPatchD& Q, Affine& m) {
    using I = Interval;
    // Corners as points: exact for unit weights, an enclosure otherwise.
    auto pt = [](const Vec4d& c) {
        if (c.w == 1.0) return std::array<I, 3>{I(c.x), I(c.y), I(c.z)};
        return std::array<I, 3>{I(c.x) / I(c.w), I(c.y) / I(c.w), I(c.z) / I(c.w)};
    };
    auto vec = [&](const Vec4d& a, const Vec4d& b) {
        const auto pa = pt(a), pb = pt(b);
        return std::array<I, 3>{pa[0] - pb[0], pa[1] - pb[1], pa[2] - pb[2]};
    };
    using planar_detail::corner;
    const std::array<I, 3> q1 = vec(corner(Q, 1, 0), corner(Q, 0, 0)), q2 = vec(corner(Q, 0, 1), corner(Q, 0, 0));
    const std::array<I, 3> p1 = vec(corner(P, 1, 0), corner(P, 0, 0)), p2 = vec(corner(P, 0, 1), corner(P, 0, 0));
    const std::array<I, 3> d = vec(corner(P, 0, 0), corner(Q, 0, 0));
    // Solve q1 s + q2 t = r in the two coordinates where [q1 q2] is best conditioned.
    int bi = 0, bj = 1;
    double best = -1.0;
    for (int i = 0; i < 3; ++i)
        for (int j = i + 1; j < 3; ++j) {
            const double det = std::fabs(num::mid(q1[i]) * num::mid(q2[j]) - num::mid(q1[j]) * num::mid(q2[i]));
            if (det > best) { best = det; bi = i; bj = j; }
        }
    auto cross = [&](const std::array<I, 3>& a, const std::array<I, 3>& b) { return a[bi] * b[bj] - a[bj] * b[bi]; };
    const I D = cross(q1, q2);
    if (num::contains_zero(D)) return false;
    m.a[0][0] = cross(d, q2) / D;
    m.a[0][1] = cross(p1, q2) / D;
    m.a[0][2] = cross(p2, q2) / D;
    m.a[1][0] = cross(q1, d) / D;
    m.a[1][1] = cross(q1, p1) / D;
    m.a[1][2] = cross(q1, p2) / D;

    // The same in exact arithmetic where every corner is an exact point (its
    // weight a power of two): each coefficient whose exact value is a double
    // becomes that point, so that exact parameters stay exact through the map
    // (a corner of one face on a side of the other, say).
    auto pow2 = [](double w) {
        int e = 0;
        return w > 0.0 && std::isfinite(w) && std::frexp(w, &e) == 0.5;
    };
    const Vec4d cs[6] = {corner(Q, 0, 0), corner(Q, 1, 0), corner(Q, 0, 1), corner(P, 0, 0), corner(P, 1, 0), corner(P, 0, 1)};
    bool exact = true;
    for (const Vec4d& c : cs) exact = exact && pow2(c.w);
    if (exact) {
        using E = num::Expansion;
        auto ept = [](const Vec4d& c) { return std::array<double, 3>{c.x / c.w, c.y / c.w, c.z / c.w}; };
        auto evec = [&](const Vec4d& a, const Vec4d& b) {
            const auto pa = ept(a), pb = ept(b);
            return std::array<E, 3>{E(pa[0]) - E(pb[0]), E(pa[1]) - E(pb[1]), E(pa[2]) - E(pb[2])};
        };
        const auto eq1 = evec(cs[1], cs[0]), eq2 = evec(cs[2], cs[0]);
        const auto ep1 = evec(cs[4], cs[3]), ep2 = evec(cs[5], cs[3]), ed = evec(cs[3], cs[0]);
        auto ecross = [&](const std::array<E, 3>& a, const std::array<E, 3>& b) { return a[bi] * b[bj] - a[bj] * b[bi]; };
        const E eD = ecross(eq1, eq2);
        // n / eD as an exact double, if it is one.
        auto quotient = [&](const E& n, I& out) {
            const double q = n.estimate() / eD.estimate();
            if (std::isfinite(q) && (E(q) * eD - n).sign() == 0) out = I(q);
        };
        quotient(ecross(ed, eq2), m.a[0][0]);
        quotient(ecross(ep1, eq2), m.a[0][1]);
        quotient(ecross(ep2, eq2), m.a[0][2]);
        quotient(ecross(eq1, ed), m.a[1][0]);
        quotient(ecross(eq1, ep1), m.a[1][1]);
        quotient(ecross(eq1, ep2), m.a[1][2]);
    }
    return true;
}

/// Parameters (s, t) of a point P that lies in the plane of the affine
/// parallelogram Q, as enclosures.
inline bool plane_params(const PolyPatchD& Q, const geom::Vec3<Interval>& P, Interval& s, Interval& t) {
    using I = Interval;
    // Corners as points: exact for unit weights, an enclosure otherwise.
    auto pt = [](const Vec4d& c) {
        if (c.w == 1.0) return std::array<I, 3>{I(c.x), I(c.y), I(c.z)};
        return std::array<I, 3>{I(c.x) / I(c.w), I(c.y) / I(c.w), I(c.z) / I(c.w)};
    };
    auto vec = [&](const Vec4d& a, const Vec4d& b) {
        const auto pa = pt(a), pb = pt(b);
        return std::array<I, 3>{pa[0] - pb[0], pa[1] - pb[1], pa[2] - pb[2]};
    };
    using planar_detail::corner;
    const auto q0 = pt(corner(Q, 0, 0));
    const std::array<I, 3> q1 = vec(corner(Q, 1, 0), corner(Q, 0, 0)), q2 = vec(corner(Q, 0, 1), corner(Q, 0, 0));
    const std::array<I, 3> r = {P.x - q0[0], P.y - q0[1], P.z - q0[2]};
    int bi = 0, bj = 1;
    double best = -1.0;
    for (int i = 0; i < 3; ++i)
        for (int j = i + 1; j < 3; ++j) {
            const double det = std::fabs(num::mid(q1[i]) * num::mid(q2[j]) - num::mid(q1[j]) * num::mid(q2[i]));
            if (det > best) { best = det; bi = i; bj = j; }
        }
    auto cross = [&](const std::array<I, 3>& a, const std::array<I, 3>& b) { return a[bi] * b[bj] - a[bj] * b[bi]; };
    const I D = cross(q1, q2);
    if (num::contains_zero(D)) return false;
    s = cross(r, q2) / D;
    t = cross(q1, r) / D;
    return true;
}

/// If the chart point is exactly a corner of its affine parallelogram, that
/// corner's homogeneous coordinates: an exact point, comparable across
/// surfaces with planar_detail::same_point_h.
inline bool exact_corner(const PolyPatchD& P, const Interval& u, const Interval& v, Vec4d& p) {
    if (u.lo != u.hi || v.lo != v.hi) return false;
    if ((u.lo != 0.0 && u.lo != 1.0) || (v.lo != 0.0 && v.lo != 1.0)) return false;
    p = planar_detail::corner(P, u.lo == 1.0 ? 1 : 0, v.lo == 1.0 ? 1 : 0);
    return p.w > 0.0;
}

} // namespace ak::op
