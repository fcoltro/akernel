// akernel — certified classification of an isolated tangent contact between
// two surface patches (docs/SURVEY.md, roadmap phase 2; RESULTS.md §3,
// "Tangent points").
//
// Where two surfaces touch at one point, their intersection is a singular
// root of A(u,v) - B(s,t) = 0: interval methods can prove that no root
// exists, but never that one does (an arbitrarily small perturbation
// removes it). Deflation turns the question into a regular one. The square
// system in x = (u, v, s, t, mu)
//
//     A(u,v) - B(s,t) + mu N_A(u,v) = 0      (3 equations, N_A = A_u x A_v)
//     N_A . B_s = 0,  N_A . B_t = 0           (parallel tangent planes)
//
// has, at a point where the surfaces have parallel tangent planes and a
// definite relative curvature, a regular root: an interval Krawczyk test
// proves it unique in a box. mu is the signed distance from A's point to
// B's along A's normal (in units of |N_A|): positive when B lies above. Then:
//  - mu proven nonzero: the surfaces have a pair of facing points apart.
//    With the relative curvature D = II_B - II_A definite, the height of B
//    over A near there is mu |N_A| + x^T D x / 2 + ...: if mu and D have the
//    same sign the surfaces are separated there (`Separated`), otherwise
//    they cross in a small closed loop around it (`Loop`).
//  - mu's enclosure contains 0: an exact candidate (dyadic parameters near
//    the root) is checked in expansion arithmetic, in homogeneous form so
//    rational patches are covered: the points are proportional and B_s, B_t
//    lie in A's projective tangent plane span(P, P_u, P_v), i.e. two 4x4
//    determinants vanish. The candidate is then a root with mu = 0, the
//    only one in the box: the surfaces touch exactly there. With D definite
//    the contact is isolated (`Touch`); with D indefinite two intersection
//    curves cross at the point (`SingularCrossing`).
//  - Otherwise `Undecided`, with the reason.
// "Separated" and "Loop" are local statements about the neighbourhood of
// the certified pair; deciding that no other intersection exists in a whole
// region is the job of the caller (the surface/surface search).
#pragma once

#include <array>
#include <cmath>
#include <string>
#include <vector>

#include <akernel/geom/bernstein.hpp>
#include <akernel/geom/seam_certificate.hpp>
#include <akernel/num/expansion.hpp>
#include <akernel/num/interval.hpp>

namespace ak::solve {

using geom::PolyPatchD;
using geom::PolyPatchI;
using geom::Vec3;
using geom::Vec4;
using num::Interval;

enum class ContactKind { NotFound, Separated, Loop, Touch, SingularCrossing, Undecided };

inline const char* to_string(ContactKind k) {
    switch (k) {
    case ContactKind::NotFound: return "not found";
    case ContactKind::Separated: return "separated";
    case ContactKind::Loop: return "small loop";
    case ContactKind::Touch: return "isolated touch";
    case ContactKind::SingularCrossing: return "singular crossing";
    case ContactKind::Undecided: return "undecided";
    }
    return "?";
}

struct ContactResult {
    ContactKind kind = ContactKind::NotFound;
    std::string why;                ///< for NotFound / Undecided
    std::array<Interval, 5> enc{};  ///< (u, v, s, t, mu): the only root in its box
    int d_sign = 0;                 ///< sign of the relative curvature if definite (+1, -1), else 0
    bool exact = false;             ///< `point` is the contact, exactly (Touch, SingularCrossing)
    std::array<double, 4> point{};  ///< exact (u, v, s, t) of the contact
};

namespace contact_detail {

using Box5 = std::array<Interval, 5>;
using Mat5 = std::array<std::array<double, 5>, 5>;

/// Homogeneous derivative patches up to order 3, in intervals and doubles.
struct Derivs {
    std::array<PolyPatchI, geom::kJetSize> hi;
    geom::BezierSurface d; // double evaluation
    explicit Derivs(const PolyPatchD& P) : d(P) {
        hi[0] = geom::to_interval(P);
        for (int n = 1; n <= 3; ++n)
            for (int l = 0; l <= n; ++l) {
                const int k = n - l;
                hi[static_cast<std::size_t>(geom::jet_index(k, l))] =
                    k > 0 ? hi[static_cast<std::size_t>(geom::jet_index(k - 1, l))].deriv_u()
                          : hi[static_cast<std::size_t>(geom::jet_index(k, l - 1))].deriv_v();
            }
    }
    /// Projected jet up to order 3 over the parameter box.
    std::array<Vec3<Interval>, geom::kJetSize> jet(const Interval& u, const Interval& v) const {
        std::array<Vec4<Interval>, geom::kJetSize> H{};
        for (int e = 0; e < geom::kJetSize; ++e) H[static_cast<std::size_t>(e)] = hi[static_cast<std::size_t>(e)].eval(u, v);
        return geom::seam_detail::project_jet(H);
    }
};

template <class T> Vec3<T> J1(const std::array<Vec3<T>, geom::kJetSize>& J, int k, int l) {
    return J[static_cast<std::size_t>(geom::jet_index(k, l))];
}

/// F and its Jacobian, generic over the scalar (double or Interval), from
/// the projected jets of A at (u, v) and B at (s, t).
template <class T>
void system(const std::array<Vec3<T>, geom::kJetSize>& a, const std::array<Vec3<T>, geom::kJetSize>& b, const T& mu,
            std::array<T, 5>& F, std::array<std::array<T, 5>, 5>& J) {
    const Vec3<T> Au = J1(a, 1, 0), Av = J1(a, 0, 1), Auu = J1(a, 2, 0), Auv = J1(a, 1, 1), Avv = J1(a, 0, 2);
    const Vec3<T> Bs = J1(b, 1, 0), Bt = J1(b, 0, 1), Bss = J1(b, 2, 0), Bst = J1(b, 1, 1), Btt = J1(b, 0, 2);
    const Vec3<T> N = geom::cross(Au, Av);
    const Vec3<T> Nu = geom::cross(Auu, Av) + geom::cross(Au, Auv);
    const Vec3<T> Nv = geom::cross(Auv, Av) + geom::cross(Au, Avv);
    const Vec3<T> R = a[0] - b[0] + N * mu;
    F = {R.x, R.y, R.z, geom::dot(N, Bs), geom::dot(N, Bt)};
    const Vec3<T> cu = Au + Nu * mu, cv = Av + Nv * mu;
    const T cols[3][5] = {{cu.x, cv.x, -Bs.x, -Bt.x, N.x}, {cu.y, cv.y, -Bs.y, -Bt.y, N.y}, {cu.z, cv.z, -Bs.z, -Bt.z, N.z}};
    for (int r = 0; r < 3; ++r)
        for (int c = 0; c < 5; ++c) J[static_cast<std::size_t>(r)][static_cast<std::size_t>(c)] = cols[r][c];
    J[3] = {geom::dot(Nu, Bs), geom::dot(Nv, Bs), geom::dot(N, Bss), geom::dot(N, Bst), T(0.0)};
    J[4] = {geom::dot(Nu, Bt), geom::dot(Nv, Bt), geom::dot(N, Bst), geom::dot(N, Btt), T(0.0)};
}

inline std::array<geom::Vec3d, geom::kJetSize> jet_d(const geom::BezierSurface& S, double u, double v) {
    const geom::Jet3 j = S.jet(u, v);
    std::array<geom::Vec3d, geom::kJetSize> r;
    for (int e = 0; e < geom::kJetSize; ++e) r[static_cast<std::size_t>(e)] = j[static_cast<std::size_t>(e)];
    return r;
}

/// Inverse of a 5x5 matrix by Gauss-Jordan with partial pivoting (doubles; an
/// approximate inverse is all the Krawczyk test needs).
inline bool invert(Mat5 M, Mat5& inv) {
    for (int i = 0; i < 5; ++i)
        for (int j = 0; j < 5; ++j) inv[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] = i == j ? 1.0 : 0.0;
    for (int c = 0; c < 5; ++c) {
        int p = c;
        for (int r = c + 1; r < 5; ++r)
            if (std::fabs(M[static_cast<std::size_t>(r)][static_cast<std::size_t>(c)]) > std::fabs(M[static_cast<std::size_t>(p)][static_cast<std::size_t>(c)])) p = r;
        if (!(std::fabs(M[static_cast<std::size_t>(p)][static_cast<std::size_t>(c)]) > 0.0)) return false;
        std::swap(M[static_cast<std::size_t>(p)], M[static_cast<std::size_t>(c)]);
        std::swap(inv[static_cast<std::size_t>(p)], inv[static_cast<std::size_t>(c)]);
        const double d = M[static_cast<std::size_t>(c)][static_cast<std::size_t>(c)];
        for (int j = 0; j < 5; ++j) {
            M[static_cast<std::size_t>(c)][static_cast<std::size_t>(j)] /= d;
            inv[static_cast<std::size_t>(c)][static_cast<std::size_t>(j)] /= d;
        }
        for (int r = 0; r < 5; ++r) {
            if (r == c) continue;
            const double f = M[static_cast<std::size_t>(r)][static_cast<std::size_t>(c)];
            if (f == 0.0) continue;
            for (int j = 0; j < 5; ++j) {
                M[static_cast<std::size_t>(r)][static_cast<std::size_t>(j)] -= f * M[static_cast<std::size_t>(c)][static_cast<std::size_t>(j)];
                inv[static_cast<std::size_t>(r)][static_cast<std::size_t>(j)] -= f * inv[static_cast<std::size_t>(c)][static_cast<std::size_t>(j)];
            }
        }
    }
    return true;
}

// --- exact verification ----------------------------------------------------

using E = num::Expansion;
using E4 = std::array<E, 4>;

/// Homogeneous point and first derivatives of P at dyadic (u, v), exactly.
inline void exact_first(const PolyPatchD& P, double u, double v, E4& pt, E4& du, E4& dv) {
    const int p = P.p, q = P.q;
    // (1 - t) as an expansion: exact, whatever t is.
    auto lerpx = [](const E& a, const E& b, const E& one_minus_t, double t) { return a * one_minus_t + b * E(t); };
    const E su = E(1.0) - E(u), sv = E(1.0) - E(v);
    // Rows in v first: for each i, the curve in v at v, and its v-derivative.
    std::vector<E4> col(static_cast<std::size_t>(p + 1)), dcol(static_cast<std::size_t>(p + 1));
    for (int i = 0; i <= p; ++i) {
        std::vector<E4> b(static_cast<std::size_t>(q + 1));
        for (int j = 0; j <= q; ++j) {
            const auto& c = P.at(i, j);
            b[static_cast<std::size_t>(j)] = {E(c.x), E(c.y), E(c.z), E(c.w)};
        }
        // de Casteljau to level q - 1, keeping the last two points for the derivative.
        for (int r = 1; r < q; ++r)
            for (int j = 0; j + r <= q; ++j)
                for (int k = 0; k < 4; ++k)
                    b[static_cast<std::size_t>(j)][static_cast<std::size_t>(k)] =
                        lerpx(b[static_cast<std::size_t>(j)][static_cast<std::size_t>(k)], b[static_cast<std::size_t>(j + 1)][static_cast<std::size_t>(k)], sv, v);
        E4 val, der;
        if (q == 0) {
            val = b[0];
            der = {E(0.0), E(0.0), E(0.0), E(0.0)};
        } else {
            for (int k = 0; k < 4; ++k) {
                val[static_cast<std::size_t>(k)] = lerpx(b[0][static_cast<std::size_t>(k)], b[1][static_cast<std::size_t>(k)], sv, v);
                der[static_cast<std::size_t>(k)] = (b[1][static_cast<std::size_t>(k)] - b[0][static_cast<std::size_t>(k)]) * E(static_cast<double>(q));
            }
        }
        col[static_cast<std::size_t>(i)] = val;
        dcol[static_cast<std::size_t>(i)] = der;
    }
    // In u: the value and u-derivative of `col`, the value of `dcol`.
    auto casteljau = [&](std::vector<E4> b, bool want_der, E4& value, E4& deriv) {
        const int n = static_cast<int>(b.size()) - 1;
        for (int r = 1; r < n; ++r)
            for (int j = 0; j + r <= n; ++j)
                for (int k = 0; k < 4; ++k)
                    b[static_cast<std::size_t>(j)][static_cast<std::size_t>(k)] =
                        lerpx(b[static_cast<std::size_t>(j)][static_cast<std::size_t>(k)], b[static_cast<std::size_t>(j + 1)][static_cast<std::size_t>(k)], su, u);
        if (n == 0) {
            value = b[0];
            deriv = {E(0.0), E(0.0), E(0.0), E(0.0)};
            return;
        }
        for (int k = 0; k < 4; ++k) {
            value[static_cast<std::size_t>(k)] = lerpx(b[0][static_cast<std::size_t>(k)], b[1][static_cast<std::size_t>(k)], su, u);
            if (want_der) deriv[static_cast<std::size_t>(k)] = (b[1][static_cast<std::size_t>(k)] - b[0][static_cast<std::size_t>(k)]) * E(static_cast<double>(n));
        }
    };
    E4 unused;
    casteljau(col, true, pt, du);
    casteljau(dcol, false, dv, unused);
}

inline E det3(const E& a, const E& b, const E& c, const E& d, const E& e, const E& f, const E& g, const E& h, const E& i) {
    return a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g);
}

/// det of the 4x4 matrix with rows r0..r3.
inline E det4(const E4& r0, const E4& r1, const E4& r2, const E4& r3) {
    const E m0 = det3(r1[1], r1[2], r1[3], r2[1], r2[2], r2[3], r3[1], r3[2], r3[3]);
    const E m1 = det3(r1[0], r1[2], r1[3], r2[0], r2[2], r2[3], r3[0], r3[2], r3[3]);
    const E m2 = det3(r1[0], r1[1], r1[3], r2[0], r2[1], r2[3], r3[0], r3[1], r3[3]);
    const E m3 = det3(r1[0], r1[1], r1[2], r2[0], r2[1], r2[2], r3[0], r3[1], r3[2]);
    return r0[0] * m0 - r0[1] * m1 + r0[2] * m2 - r0[3] * m3;
}

/// Is (u, v, s, t) exactly a point where A and B touch (same point, B's
/// tangents in A's tangent plane)?
inline bool exact_touch(const PolyPatchD& A, const PolyPatchD& B, double u, double v, double s, double t) {
    E4 a, au, av, b, bs, bt;
    exact_first(A, u, v, a, au, av);
    exact_first(B, s, t, b, bs, bt);
    for (int k = 0; k < 3; ++k)
        if ((a[static_cast<std::size_t>(k)] * b[3] - b[static_cast<std::size_t>(k)] * a[3]).sign() != 0) return false;
    return det4(a, au, av, bs).sign() == 0 && det4(a, au, av, bt).sign() == 0;
}

} // namespace contact_detail

/// Look for a pair of points with parallel tangent planes near `seed`
/// (u, v, s, t, mu), prove it the only one in a box, and classify it.
inline ContactResult classify_contact(const PolyPatchD& A, const PolyPatchD& B, std::array<double, 5> x) {
    using namespace contact_detail;
    ContactResult out;
    const Derivs DA(A), DB(B);
    // Newton in doubles.
    for (int it = 0; it < 60; ++it) {
        const auto a = jet_d(DA.d, x[0], x[1]), b = jet_d(DB.d, x[2], x[3]);
        std::array<double, 5> F;
        Mat5 J;
        system(a, b, x[4], F, J);
        Mat5 Y;
        if (!invert(J, Y)) {
            out.why = "the contact system is singular at the seed";
            return out;
        }
        double step = 0.0;
        std::array<double, 5> dx{};
        for (int i = 0; i < 5; ++i) {
            for (int j = 0; j < 5; ++j) dx[static_cast<std::size_t>(i)] += Y[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] * F[static_cast<std::size_t>(j)];
            x[static_cast<std::size_t>(i)] -= dx[static_cast<std::size_t>(i)];
            step = std::max(step, std::fabs(dx[static_cast<std::size_t>(i)]));
        }
        if (step < 1e-15) break;
    }
    for (int i = 0; i < 4; ++i)
        if (!(x[static_cast<std::size_t>(i)] > -0.5 && x[static_cast<std::size_t>(i)] < 1.5)) {
            out.why = "Newton left the patches";
            return out;
        }
    // Krawczyk on boxes of growing radius around x.
    const auto evalJ = [&](const Box5& X, std::array<Interval, 5>& F, std::array<std::array<Interval, 5>, 5>& J) {
        const auto a = DA.jet(X[0], X[1]), b = DB.jet(X[2], X[3]);
        system(a, b, X[4], F, J);
    };
    Mat5 Y;
    {
        const auto a = jet_d(DA.d, x[0], x[1]), b = jet_d(DB.d, x[2], x[3]);
        std::array<double, 5> F;
        Mat5 J;
        system(a, b, x[4], F, J);
        if (!invert(J, Y)) {
            out.why = "the contact system is singular at the root (not a regular contact)";
            return out;
        }
    }
    bool proven = false;
    Box5 K{};
    for (int k = 0; k < 24 && !proven; ++k) {
        const double r = 1e-13 * std::pow(4.0, k);
        Box5 X, P;
        for (int i = 0; i < 5; ++i) {
            const double w = r * (1.0 + std::fabs(x[static_cast<std::size_t>(i)]));
            X[static_cast<std::size_t>(i)] = Interval(x[static_cast<std::size_t>(i)] - w, x[static_cast<std::size_t>(i)] + w);
            P[static_cast<std::size_t>(i)] = Interval(x[static_cast<std::size_t>(i)]);
        }
        std::array<Interval, 5> Fp, FX;
        std::array<std::array<Interval, 5>, 5> Jp, JX;
        evalJ(P, Fp, Jp);
        evalJ(X, FX, JX);
        bool inside = true;
        for (int i = 0; i < 5; ++i) {
            Interval s = Interval(x[static_cast<std::size_t>(i)]);
            for (int j = 0; j < 5; ++j) s = s - Interval(Y[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)]) * Fp[static_cast<std::size_t>(j)];
            for (int j = 0; j < 5; ++j) {
                Interval m = Interval(i == j ? 1.0 : 0.0);
                for (int l = 0; l < 5; ++l) m = m - Interval(Y[static_cast<std::size_t>(i)][static_cast<std::size_t>(l)]) * JX[static_cast<std::size_t>(l)][static_cast<std::size_t>(j)];
                s = s + m * (X[static_cast<std::size_t>(j)] - Interval(x[static_cast<std::size_t>(j)]));
            }
            K[static_cast<std::size_t>(i)] = s;
            inside = inside && num::strictly_inside(s, X[static_cast<std::size_t>(i)]);
        }
        proven = inside;
    }
    if (!proven) {
        out.why = "no box around the contact passed the Krawczyk test";
        return out;
    }
    out.enc = K;
    // Relative curvature over K: second forms of A and of B in A's tangent
    // coordinates (u, v), both with A's normal.
    {
        const auto a = DA.jet(K[0], K[1]), b = DB.jet(K[2], K[3]);
        const Vec3<Interval> Au = J1(a, 1, 0), Av = J1(a, 0, 1), Bs = J1(b, 1, 0), Bt = J1(b, 0, 1);
        const Vec3<Interval> N = geom::cross(Au, Av);
        const Interval nN = num::sqrt(geom::dot(N, N));
        const Vec3<Interval> n = N / nN;
        // A_u, A_v in B's tangent basis: solve [Bs Bt] M = [Au Av] by B's metric.
        const Interval g11 = geom::dot(Bs, Bs), g12 = geom::dot(Bs, Bt), g22 = geom::dot(Bt, Bt);
        const Interval det = g11 * g22 - g12 * g12;
        auto coords = [&](const Vec3<Interval>& w, Interval& c1, Interval& c2) {
            const Interval r1 = geom::dot(Bs, w), r2 = geom::dot(Bt, w);
            c1 = (g22 * r1 - g12 * r2) / det;
            c2 = (g11 * r2 - g12 * r1) / det;
        };
        Interval m11, m21, m12, m22;
        coords(Au, m11, m21);
        coords(Av, m12, m22);
        const Interval bB11 = geom::dot(n, J1(b, 2, 0)), bB12 = geom::dot(n, J1(b, 1, 1)), bB22 = geom::dot(n, J1(b, 0, 2));
        // H_B = M^T b_B M
        const Interval h11 = m11 * m11 * bB11 + Interval(2.0) * m11 * m21 * bB12 + m21 * m21 * bB22;
        const Interval h12 = m11 * m12 * bB11 + (m11 * m22 + m21 * m12) * bB12 + m21 * m22 * bB22;
        const Interval h22 = m12 * m12 * bB11 + Interval(2.0) * m12 * m22 * bB12 + m22 * m22 * bB22;
        const Interval d11 = h11 - geom::dot(n, J1(a, 2, 0)), d12 = h12 - geom::dot(n, J1(a, 1, 1)),
                       d22 = h22 - geom::dot(n, J1(a, 0, 2));
        const Interval dd = d11 * d22 - d12 * d12;
        if (dd.lo > 0.0 && d11.lo > 0.0) out.d_sign = 1;
        else if (dd.lo > 0.0 && d11.hi < 0.0) out.d_sign = -1;
        else if (dd.hi < 0.0) out.d_sign = 2; // indefinite, proven
        else out.d_sign = 0;
    }
    const Interval mu = K[4];
    if (mu.lo > 0.0 || mu.hi < 0.0) {
        const int ms = mu.lo > 0.0 ? 1 : -1;
        if (out.d_sign == 1 || out.d_sign == -1) {
            out.kind = ms == out.d_sign ? ContactKind::Separated : ContactKind::Loop;
            return out;
        }
        out.kind = ContactKind::Undecided;
        out.why = "the facing points are apart, but the relative curvature is not proven definite";
        return out;
    }
    // mu may be zero: try an exact contact at dyadic parameters near x.
    for (int bits : {30, 40, 50}) {
        std::array<double, 4> c;
        bool in = true;
        for (int i = 0; i < 4; ++i) {
            c[static_cast<std::size_t>(i)] = std::ldexp(std::nearbyint(std::ldexp(x[static_cast<std::size_t>(i)], bits)), -bits);
            in = in && num::contains(K[static_cast<std::size_t>(i)], c[static_cast<std::size_t>(i)]);
        }
        if (!in) continue;
        if (exact_touch(A, B, c[0], c[1], c[2], c[3])) {
            out.exact = true;
            out.point = c;
            if (out.d_sign == 1 || out.d_sign == -1) out.kind = ContactKind::Touch;
            else if (out.d_sign == 2) out.kind = ContactKind::SingularCrossing;
            else {
                out.kind = ContactKind::Undecided;
                out.why = "the surfaces touch exactly here, but the relative curvature is not proven definite or indefinite";
            }
            return out;
        }
    }
    out.kind = ContactKind::Undecided;
    out.why = "the gap between the facing points is not proven nonzero, and no exact contact was found";
    return out;
}

namespace contact_detail {

/// The normal offset from A to B over a box of A's parameters: for every
/// (u, v) in U x V, the normal line A(u,v) + lambda N_A(u,v) meets B exactly
/// once in the box Y (s, t, lambda) (a Krawczyk test with (u, v) as interval
/// parameters, centred at yc), and H encloses the Hessian of lambda(u, v)
/// over U x V (implicit differentiation, A's jet to order 3 and B's to
/// order 2 enclosed over the boxes).
struct Offset {
    bool ok = false;
    std::string why;
    std::array<Interval, 3> Y{};
    Interval H[2][2];
};

inline bool inverse3(const Vec3<Interval>& c0, const Vec3<Interval>& c1, const Vec3<Interval>& c2, std::array<std::array<Interval, 3>, 3>& inv) {
    const Vec3<Interval> r0 = geom::cross(c1, c2), r1 = geom::cross(c2, c0), r2 = geom::cross(c0, c1);
    const Interval det = geom::dot(c0, r0);
    if (num::contains_zero(det)) return false;
    const Vec3<Interval> rows[3] = {r0, r1, r2};
    for (int i = 0; i < 3; ++i) inv[static_cast<std::size_t>(i)] = {rows[i].x / det, rows[i].y / det, rows[i].z / det};
    return true;
}

inline std::array<Interval, 3> apply3(const std::array<std::array<Interval, 3>, 3>& M, const Vec3<Interval>& x) {
    return {M[0][0] * x.x + M[0][1] * x.y + M[0][2] * x.z, M[1][0] * x.x + M[1][1] * x.y + M[1][2] * x.z,
            M[2][0] * x.x + M[2][1] * x.y + M[2][2] * x.z};
}

/// The solution box Y over U x V (see Offset), centred at yc with at least
/// the radii `rad`. False if no inflation passed.
inline bool offset_solution(const Derivs& DA, const Derivs& DB, const Interval& U, const Interval& V, const double yc[3], double rad[3],
                            std::array<Interval, 3>& Y, std::string& why) {
    const auto a = DA.jet(U, V);
    const Vec3<Interval> N = geom::cross(J1(a, 1, 0), J1(a, 0, 1));
    // An approximate inverse: the Jacobian at the centre, in doubles.
    std::array<std::array<double, 3>, 3> Yd{};
    {
        const auto ad = jet_d(DA.d, num::mid(U), num::mid(V)), bd = jet_d(DB.d, yc[0], yc[1]);
        const geom::Vec3d n = geom::cross(ad[1], ad[2]);
        const geom::Vec3d c0 = -bd[1], c1 = -bd[2], c2 = n;
        const geom::Vec3d r0 = geom::cross(c1, c2), r1 = geom::cross(c2, c0), r2 = geom::cross(c0, c1);
        const double det = geom::dot(c0, r0);
        if (!(std::fabs(det) > 0.0)) {
            why = "the normal lines of A are tangent to B";
            return false;
        }
        const geom::Vec3d rows[3] = {r0, r1, r2};
        for (int i = 0; i < 3; ++i) Yd[static_cast<std::size_t>(i)] = {rows[i].x / det, rows[i].y / det, rows[i].z / det};
    }
    const auto bc = DB.jet(Interval(yc[0]), Interval(yc[1]));
    const Vec3<Interval> G0 = a[0] + N * Interval(yc[2]) - bc[0]; // G at the centre, for all (u, v)
    for (int it = 0; it < 12; ++it) {
        for (int i = 0; i < 3; ++i) Y[static_cast<std::size_t>(i)] = Interval(yc[i] - rad[i], yc[i] + rad[i]);
        const auto b = DB.jet(Y[0], Y[1]);
        const Vec3<Interval> cols[3] = {-J1(b, 1, 0), -J1(b, 0, 1), N};
        bool inside = true;
        double grow[3] = {0, 0, 0};
        for (int i = 0; i < 3; ++i) {
            const auto& yr = Yd[static_cast<std::size_t>(i)];
            Interval k = Interval(yc[i]) - (Interval(yr[0]) * G0.x + Interval(yr[1]) * G0.y + Interval(yr[2]) * G0.z);
            for (int j = 0; j < 3; ++j) {
                const Interval m = Interval(i == j ? 1.0 : 0.0) - (Interval(yr[0]) * cols[j].x + Interval(yr[1]) * cols[j].y + Interval(yr[2]) * cols[j].z);
                k = k + m * (Y[static_cast<std::size_t>(j)] - Interval(yc[j]));
            }
            inside = inside && num::strictly_inside(k, Y[static_cast<std::size_t>(i)]);
            grow[i] = std::max(yc[i] - k.lo, k.hi - yc[i]);
        }
        if (inside) return true;
        for (int i = 0; i < 3; ++i) rad[i] = std::max(rad[i], 2.0 * grow[i]);
    }
    why = "the normal lines of A could not be proven to meet B once each in a box";
    return false;
}

/// First derivatives of y = (s, t, lambda) in (u, v), and the Hessian of
/// lambda, from A's jet over (U, V) and B's over Y.
inline bool offset_derivatives(const Derivs& DA, const Derivs& DB, const Interval& U, const Interval& V, const std::array<Interval, 3>& Y,
                               std::array<Interval, 3> yx[2], Interval H[2][2], std::string& why) {
    const auto a = DA.jet(U, V);
    const auto b = DB.jet(Y[0], Y[1]);
    const Interval lam = Y[2];
    const Vec3<Interval> Au = J1(a, 1, 0), Av = J1(a, 0, 1), Auu = J1(a, 2, 0), Auv = J1(a, 1, 1), Avv = J1(a, 0, 2);
    const Vec3<Interval> N = geom::cross(Au, Av);
    const Vec3<Interval> Nu = geom::cross(Auu, Av) + geom::cross(Au, Auv);
    const Vec3<Interval> Nv = geom::cross(Auv, Av) + geom::cross(Au, Avv);
    const Vec3<Interval> Bs = J1(b, 1, 0), Bt = J1(b, 0, 1), Bss = J1(b, 2, 0), Bst = J1(b, 1, 1), Btt = J1(b, 0, 2);
    std::array<std::array<Interval, 3>, 3> Gi;
    if (!inverse3(-Bs, -Bt, N, Gi)) {
        why = "the normal lines of A are not proven transversal to B over the box";
        return false;
    }
    const Vec3<Interval> Gx[2] = {Au + Nu * lam, Av + Nv * lam};
    for (int i = 0; i < 2; ++i) {
        const auto r = apply3(Gi, Gx[i]);
        yx[i] = {-r[0], -r[1], -r[2]};
    }
    const Vec3<Interval> Auuu = J1(a, 3, 0), Auuv = J1(a, 2, 1), Auvv = J1(a, 1, 2), Avvv = J1(a, 0, 3);
    const Vec3<Interval> Nuu = geom::cross(Auuu, Av) + geom::cross(Auu, Auv) * Interval(2.0) + geom::cross(Au, Auuv);
    const Vec3<Interval> Nuv = geom::cross(Auuv, Av) + geom::cross(Auu, Avv) + geom::cross(Au, Auvv);
    const Vec3<Interval> Nvv = geom::cross(Auvv, Av) + geom::cross(Auv, Avv) * Interval(2.0) + geom::cross(Au, Avvv);
    const Vec3<Interval> A2[2][2] = {{Auu, Auv}, {Auv, Avv}}, N1[2] = {Nu, Nv}, N2[2][2] = {{Nuu, Nuv}, {Nuv, Nvv}};
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j) {
            const Interval si = yx[i][0], ti = yx[i][1], li = yx[i][2], sj = yx[j][0], tj = yx[j][1], lj = yx[j][2];
            const Vec3<Interval> r = A2[i][j] + N2[i][j] * lam + N1[i] * lj + N1[j] * li -
                                     (Bss * (si * sj) + Bst * (si * tj + sj * ti) + Btt * (ti * tj));
            H[i][j] = -apply3(Gi, r)[2];
        }
    return true;
}

/// The offset over U x V around the centre (u, v) -> (s, t, lambda) = yc,
/// with Y covering S x T.
inline Offset offset_over(const Derivs& DA, const Derivs& DB, const Interval& U, const Interval& V, const Interval& S, const Interval& T,
                          const double yc[3]) {
    Offset o;
    double rad[3] = {std::max(yc[0] - S.lo, S.hi - yc[0]), std::max(yc[1] - T.lo, T.hi - yc[1]), std::max(U.hi - U.lo, V.hi - V.lo)};
    if (!offset_solution(DA, DB, U, V, yc, rad, o.Y, o.why)) return o;
    std::array<Interval, 3> yx[2];
    if (!offset_derivatives(DA, DB, U, V, o.Y, yx, o.H, o.why)) return o;
    o.ok = true;
    return o;
}

/// +1 or -1 if the symmetric interval matrix H is proven positive or negative
/// definite, with a lower bound m > 0 on the magnitude of its eigenvalues
/// (all of the same sign); else 0.
inline int definite(const Interval H[2][2], double& m) {
    const Interval det = H[0][0] * H[1][1] - H[0][1] * H[1][0];
    if (!(det.lo > 0.0)) return 0;
    const int sg = H[0][0].lo > 0.0 ? 1 : H[0][0].hi < 0.0 ? -1 : 0;
    if (sg == 0) return 0;
    // Smallest |eigenvalue| = (|a + c| - sqrt((a - c)^2 + 4 b^2)) / 2.
    const Interval tr = sg > 0 ? H[0][0] + H[1][1] : -(H[0][0] + H[1][1]);
    const Interval dis = num::sqr(H[0][0] - H[1][1]) + Interval(4.0) * num::sqr(H[0][1]);
    const Interval e = (tr - num::sqrt(dis)) * Interval(0.5);
    if (!(e.lo > 0.0)) return 0;
    m = e.lo;
    return sg;
}

} // namespace contact_detail

/// The pair touches only at `p` (exact parameters (u, v, s, t) of a contact,
/// e.g. from classify_contact) within the parameter boxes U x V of A and
/// S x T of B? Proven by the normal offset function. For (u, v) in U x V
/// the equation A(u,v) + lambda N_A(u,v) = B(s,t) has exactly one solution
/// (s, t, lambda) in a box Y containing S x T (a Krawczyk test with (u, v)
/// as interval parameters), so every intersection point in the boxes has
/// lambda(u, v) = 0. lambda vanishes at p with zero gradient (the tangent
/// planes are parallel there), and its Hessian, by implicit differentiation
/// enclosed over U x V x Y, is definite: lambda is strictly convex or concave
/// on the convex box U x V, and p is its only zero. False, with the reason,
/// if a step is not proven.
inline bool unique_touch_in(const PolyPatchD& A, const PolyPatchD& B, const std::array<double, 4>& p, const Interval& U,
                            const Interval& V, const Interval& S, const Interval& T, std::string& why) {
    using namespace contact_detail;
    if (!(num::contains(U, p[0]) && num::contains(V, p[1]) && num::contains(S, p[2]) && num::contains(T, p[3]))) {
        why = "the contact is not inside the boxes";
        return false;
    }
    const Derivs DA(A), DB(B);
    const double yc[3] = {p[2], p[3], 0.0};
    const Offset o = offset_over(DA, DB, U, V, S, T, yc);
    if (!o.ok) {
        why = o.why;
        return false;
    }
    double m = 0.0;
    if (definite(o.H, m) == 0) {
        why = "the offset between the surfaces is not proven strictly convex over the box";
        return false;
    }
    return true;
}

/// A near miss: do A and B not meet at all within the parameter boxes U x V
/// of A and S x T of B? `c` is the facing pair (u, v, s, t, mu) found by
/// classify_contact (mu proven nonzero there). The normal offset lambda is
/// proven defined over U x V as in unique_touch_in, with a Hessian of one
/// sign whose eigenvalues are at least m in magnitude; lambda and its
/// gradient g are enclosed at the single point (u, v) of c. By Taylor's
/// theorem lambda(c + d) >= lambda(c) + g.d + m |d|^2 / 2 >= lambda(c) -
/// |g|^2 / (2 m) for a convex offset (mirrored for a concave one): if that
/// stays on the side of lambda(c), lambda has no zero in U x V, so the
/// surfaces do not meet there. Returns +1 (B on the side of A's normal) or
/// -1, or 0 with the reason.
inline int separated_in(const PolyPatchD& A, const PolyPatchD& B, const std::array<double, 5>& c, const Interval& U, const Interval& V,
                        const Interval& S, const Interval& T, std::string& why) {
    using namespace contact_detail;
    if (!(num::contains(U, c[0]) && num::contains(V, c[1]))) {
        why = "the facing pair is not inside the boxes";
        return 0;
    }
    const Derivs DA(A), DB(B);
    const double yc[3] = {c[2], c[3], c[4]};
    const Offset o = offset_over(DA, DB, U, V, S, T, yc);
    if (!o.ok) {
        why = o.why;
        return 0;
    }
    double m = 0.0;
    const int sg = definite(o.H, m);
    if (sg == 0) {
        why = "the offset between the surfaces is not proven strictly convex over the box";
        return 0;
    }
    // lambda and its gradient at the point (u, v) of c.
    const Interval Uc(c[0]), Vc(c[1]);
    double rad[3] = {1e-12 * (1.0 + std::fabs(c[2])), 1e-12 * (1.0 + std::fabs(c[3])), 1e-12 * (1.0 + std::fabs(c[4]))};
    std::array<Interval, 3> Yc;
    if (!offset_solution(DA, DB, Uc, Vc, yc, rad, Yc, why)) return 0;
    std::array<Interval, 3> yx[2];
    Interval Hc[2][2];
    if (!offset_derivatives(DA, DB, Uc, Vc, Yc, yx, Hc, why)) return 0;
    const Interval g2 = num::sqr(yx[0][2]) + num::sqr(yx[1][2]);
    const Interval drop = g2 / (Interval(2.0) * Interval(m));
    if (sg > 0 && (Yc[2] - drop).lo > 0.0) return 1;
    if (sg < 0 && (Yc[2] + drop).hi < 0.0) return -1;
    why = "the gap between the surfaces is not proven to keep its sign over the box";
    return 0;
}

} // namespace ak::solve
