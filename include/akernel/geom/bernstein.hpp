// akernel — tensor-product Bernstein (Bezier) patches in homogeneous form.
//
// PolyPatch<T> is a polynomial map [0,1]^2 -> T^4, (u,v) -> (w*x, w*y, w*z, w).
// A rational Bezier surface is its projection (x, y, z). All subdivision and
// evaluation is by de Casteljau, which is exact arithmetic on the control
// points; with T = Interval every result is a rigorous enclosure.
//
// The polynomial is defined for all real (u,v), so restriction to a box that
// extends outside [0,1]^2 is valid, and the convex-hull property holds on
// whatever box the patch was restricted to.
#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <cassert>
#include <vector>

#include <akernel/geom/vec.hpp>
#include <akernel/num/interval.hpp>

namespace ak::geom {

inline constexpr int kMaxDegree = 15;

/// Index of the mixed partial d^(k+l)/du^k dv^l in a jet of total order <= 3:
/// (0,0) (1,0) (0,1) (2,0) (1,1) (0,2) (3,0) (2,1) (1,2) (0,3).
constexpr int jet_index(int k, int l) noexcept { return (k + l) * (k + l + 1) / 2 + l; }
inline constexpr int kJetSize = 10;
using Jet3 = std::array<Vec3d, kJetSize>;

namespace detail {

/// In-place de Casteljau on b[0..n]; returns the value at t.
template <class P, class T>
P decasteljau(P* b, int n, const T& t) {
    const T s = T(1.0) - t;
    for (int r = 1; r <= n; ++r)
        for (int i = 0; i <= n - r; ++i) b[i] = b[i] * s + b[i + 1] * t;
    return b[0];
}

/// Split b[0..n] at t into the control points of [0,t] (lo) and [t,1] (hi).
template <class P, class T>
void split1d(const P* b, int n, const T& t, P* lo, P* hi) {
    std::array<P, kMaxDegree + 1> w;
    for (int i = 0; i <= n; ++i) w[i] = b[i];
    const T s = T(1.0) - t;
    lo[0] = w[0];
    hi[n] = w[n];
    for (int r = 1; r <= n; ++r) {
        for (int i = 0; i <= n - r; ++i) w[i] = w[i] * s + w[i + 1] * t;
        lo[r] = w[0];
        hi[n - r] = w[n - r];
    }
}

/// Control points of b restricted to [a, c].
template <class P, class T>
void restrict1d(P* b, int n, double a, double c) {
    std::array<P, kMaxDegree + 1> lo, hi;
    if (std::fabs(1.0 - a) >= std::fabs(c)) {
        // keep [a,1], then cut at the image of c
        split1d(b, n, T(a), lo.data(), hi.data());
        const T t = (T(c) - T(a)) / (T(1.0) - T(a));
        split1d(hi.data(), n, t, lo.data(), b);
        for (int i = 0; i <= n; ++i) b[i] = lo[i];
    } else {
        // keep [0,c], then cut at the image of a
        split1d(b, n, T(c), lo.data(), hi.data());
        const T t = T(a) / T(c);
        split1d(lo.data(), n, t, hi.data(), b);
    }
}

} // namespace detail

template <class T>
class PolyPatch {
public:
    int p = 0; ///< degree in u
    int q = 0; ///< degree in v
    std::vector<Vec4<T>> c; ///< (p+1)*(q+1) control points, row-major in u

    PolyPatch() = default;
    PolyPatch(int p_, int q_) : p(p_), q(q_), c(static_cast<std::size_t>((p_ + 1) * (q_ + 1))) {
        assert(p_ >= 0 && q_ >= 0 && p_ <= kMaxDegree && q_ <= kMaxDegree);
    }

    Vec4<T>& at(int i, int j) { return c[static_cast<std::size_t>(i * (q + 1) + j)]; }
    const Vec4<T>& at(int i, int j) const { return c[static_cast<std::size_t>(i * (q + 1) + j)]; }

    /// Homogeneous value at (u, v).
    Vec4<T> eval(const T& u, const T& v) const {
        std::array<Vec4<T>, kMaxDegree + 1> col, row;
        for (int i = 0; i <= p; ++i) {
            for (int j = 0; j <= q; ++j) row[j] = at(i, j);
            col[i] = detail::decasteljau(row.data(), q, v);
        }
        return detail::decasteljau(col.data(), p, u);
    }

    /// Partial derivative in u as a patch of degree (p-1, q).
    PolyPatch deriv_u() const {
        if (p == 0) return PolyPatch(0, q); // identically zero
        PolyPatch d(p - 1, q);
        const T scale(static_cast<double>(p));
        for (int i = 0; i < p; ++i)
            for (int j = 0; j <= q; ++j) d.at(i, j) = (at(i + 1, j) - at(i, j)) * scale;
        return d;
    }
    PolyPatch deriv_v() const {
        if (q == 0) return PolyPatch(p, 0);
        PolyPatch d(p, q - 1);
        const T scale(static_cast<double>(q));
        for (int i = 0; i <= p; ++i)
            for (int j = 0; j < q; ++j) d.at(i, j) = (at(i, j + 1) - at(i, j)) * scale;
        return d;
    }

    void split_u(const T& t, PolyPatch& lo, PolyPatch& hi) const {
        lo = PolyPatch(p, q);
        hi = PolyPatch(p, q);
        std::array<Vec4<T>, kMaxDegree + 1> b, l, h;
        for (int j = 0; j <= q; ++j) {
            for (int i = 0; i <= p; ++i) b[i] = at(i, j);
            detail::split1d(b.data(), p, t, l.data(), h.data());
            for (int i = 0; i <= p; ++i) { lo.at(i, j) = l[i]; hi.at(i, j) = h[i]; }
        }
    }
    void split_v(const T& t, PolyPatch& lo, PolyPatch& hi) const {
        lo = PolyPatch(p, q);
        hi = PolyPatch(p, q);
        std::array<Vec4<T>, kMaxDegree + 1> b, l, h;
        for (int i = 0; i <= p; ++i) {
            for (int j = 0; j <= q; ++j) b[j] = at(i, j);
            detail::split1d(b.data(), q, t, l.data(), h.data());
            for (int j = 0; j <= q; ++j) { lo.at(i, j) = l[j]; hi.at(i, j) = h[j]; }
        }
    }

    /// The same polynomial reparameterised so that [u0,u1] x [v0,v1] maps to [0,1]^2.
    PolyPatch restricted(double u0, double u1, double v0, double v1) const {
        PolyPatch r = *this;
        std::array<Vec4<T>, kMaxDegree + 1> b;
        for (int j = 0; j <= q; ++j) {
            for (int i = 0; i <= p; ++i) b[i] = r.at(i, j);
            detail::restrict1d<Vec4<T>, T>(b.data(), p, u0, u1);
            for (int i = 0; i <= p; ++i) r.at(i, j) = b[i];
        }
        for (int i = 0; i <= p; ++i) {
            for (int j = 0; j <= q; ++j) b[j] = r.at(i, j);
            detail::restrict1d<Vec4<T>, T>(b.data(), q, v0, v1);
            for (int j = 0; j <= q; ++j) r.at(i, j) = b[j];
        }
        return r;
    }
};

using PolyPatchD = PolyPatch<double>;
using PolyPatchI = PolyPatch<num::Interval>;

/// Exact conversion of double control points to point intervals.
inline PolyPatchI to_interval(const PolyPatchD& a) {
    PolyPatchI r(a.p, a.q);
    for (std::size_t i = 0; i < a.c.size(); ++i)
        r.c[i] = {num::Interval(a.c[i].x), num::Interval(a.c[i].y), num::Interval(a.c[i].z),
                  num::Interval(a.c[i].w)};
    return r;
}

/// Componentwise hull of the homogeneous control points. By the convex-hull
/// property this encloses the homogeneous patch over [0,1]^2.
inline Vec4<num::Interval> hull(const PolyPatchI& a) {
    Vec4<num::Interval> h = a.c[0];
    for (const auto& cp : a.c) {
        h.x = num::hull(h.x, cp.x);
        h.y = num::hull(h.y, cp.y);
        h.z = num::hull(h.z, cp.z);
        h.w = num::hull(h.w, cp.w);
    }
    return h;
}

namespace detail {

/// Grid used by hull_restricted: on multiples of 2^-53 in [0,1], 1 - t is exact.
inline double grid_down(double t) { return std::max(0.0, std::floor(t * 0x1p53) * 0x1p-53); }
inline double grid_up(double t) { return std::min(1.0, std::ceil(t * 0x1p53) * 0x1p-53); }

/// One de Casteljau split in doubles at t (on the grid), as split1d.
inline void split1d_grid(const Vec4d* b, int n, double t, Vec4d* lo, Vec4d* hi) {
    std::array<Vec4d, kMaxDegree + 1> w;
    for (int i = 0; i <= n; ++i) w[i] = b[i];
    const double s = 1.0 - t;
    lo[0] = w[0];
    hi[n] = w[n];
    for (int r = 1; r <= n; ++r) {
        for (int i = 0; i <= n - r; ++i)
            w[i] = {w[i].x * s + w[i + 1].x * t, w[i].y * s + w[i + 1].y * t, w[i].z * s + w[i + 1].z * t,
                    w[i].w * s + w[i + 1].w * t};
        lo[r] = w[0];
        hi[n - r] = w[n - r];
    }
}

/// b restricted, in doubles, to a parameter range that contains [a, c]
/// (0 <= a <= c <= 1); both cuts are at grid points, so every de Casteljau
/// step is a convex combination with exact weights.
inline void restrict1d_grid(Vec4d* b, int n, double a, double c) {
    std::array<Vec4d, kMaxDegree + 1> lo, hi;
    const double a0 = grid_down(a), c0 = grid_up(c);
    if (1.0 - a0 >= c0) {
        split1d_grid(b, n, a0, lo.data(), hi.data());
        // [a0, 1] -> [0, 1]; c0 maps to (c0 - a0) / (1 - a0), rounded up past it.
        const double t = a0 < 1.0 ? std::min(1.0, std::ceil(((c0 - a0) / (1.0 - a0)) * 0x1p53 + 1.0) * 0x1p-53) : 1.0;
        split1d_grid(hi.data(), n, t, b, hi.data());
    } else {
        split1d_grid(b, n, c0, lo.data(), hi.data());
        const double t = std::max(0.0, std::floor((a0 / c0) * 0x1p53 - 1.0) * 0x1p-53);
        split1d_grid(lo.data(), n, t, lo.data(), b);
    }
}

} // namespace detail

/// Componentwise hull of the control points of a restricted to
/// [u0,u1] x [v0,v1], the same enclosure as hull(a.restricted(...)) but
/// computed in doubles, several times faster.
///
/// Rigour: the patch is split at the midpoints m of its control points, with
/// radii r. Both cuts in each direction are moved outward to multiples of
/// 2^-53, so the restriction is to a box containing the asked one (still a
/// valid enclosure) and every step is a convex combination with exact
/// weights s + t = 1. A convex combination moves no point further than r
/// from its exact image, and each rounded step a s + b t adds at most
/// 3u max(|a|, |b|) (u = 2^-53), with at most 2(p + q) steps on any path. The
/// hull is widened by R + 4u (2p + 2q) M (M the largest |mid|), plus an
/// absolute term for underflow. Boxes outside [0,1]^2 fall back to intervals.
inline Vec4<num::Interval> hull_restricted(const PolyPatchI& a, double u0, double u1, double v0, double v1) {
    using num::Interval;
    if (!(0.0 <= u0 && u0 <= u1 && u1 <= 1.0 && 0.0 <= v0 && v0 <= v1 && v1 <= 1.0))
        return hull(a.restricted(u0, u1, v0, v1));
    const int p = a.p, q = a.q;
    std::vector<Vec4d> m(a.c.size());
    double R[4] = {0, 0, 0, 0}, M[4] = {0, 0, 0, 0};
    for (std::size_t k = 0; k < a.c.size(); ++k) {
        const Interval* iv[4] = {&a.c[k].x, &a.c[k].y, &a.c[k].z, &a.c[k].w};
        double mv[4];
        for (int j = 0; j < 4; ++j) {
            const Interval& x = *iv[j];
            if (!num::is_finite(x)) return hull(a.restricted(u0, u1, v0, v1));
            mv[j] = num::mid(x);
            R[j] = std::max(R[j], num::detail::step_up(std::max(x.hi - mv[j], mv[j] - x.lo)));
            M[j] = std::max(M[j], std::fabs(mv[j]));
        }
        m[k] = {mv[0], mv[1], mv[2], mv[3]};
    }
    std::array<Vec4d, kMaxDegree + 1> b;
    for (int j = 0; j <= q; ++j) {
        for (int i = 0; i <= p; ++i) b[i] = m[static_cast<std::size_t>(i * (q + 1) + j)];
        detail::restrict1d_grid(b.data(), p, u0, u1);
        for (int i = 0; i <= p; ++i) m[static_cast<std::size_t>(i * (q + 1) + j)] = b[i];
    }
    double lo[4], hi[4];
    for (int j = 0; j < 4; ++j) { lo[j] = std::numeric_limits<double>::infinity(); hi[j] = -lo[j]; }
    for (int i = 0; i <= p; ++i) {
        for (int j = 0; j <= q; ++j) b[j] = m[static_cast<std::size_t>(i * (q + 1) + j)];
        detail::restrict1d_grid(b.data(), q, v0, v1);
        for (int j = 0; j <= q; ++j) {
            const double v[4] = {b[j].x, b[j].y, b[j].z, b[j].w};
            for (int k = 0; k < 4; ++k) { lo[k] = std::min(lo[k], v[k]); hi[k] = std::max(hi[k], v[k]); }
        }
    }
    const double steps = 2.0 * (p + q) + 1.0;
    Interval out[4];
    for (int k = 0; k < 4; ++k) {
        if (!std::isfinite(lo[k]) || !std::isfinite(hi[k])) return hull(a.restricted(u0, u1, v0, v1));
        // Upper bound of R + 4u * steps * M + tiny, rounded up.
        const Interval pad = Interval(R[k]) + Interval(4.0 * 0x1p-53 * steps) * Interval(M[k]) + Interval(0x1p-1000);
        out[k] = Interval(lo[k], hi[k]) + Interval(-pad.hi, pad.hi);
    }
    return {out[0], out[1], out[2], out[3]};
}

/// Axis-aligned box enclosing the projected (x, y, z) surface over [0,1]^2.
/// Valid when every weight is positive; otherwise returns the whole space.
inline Vec3<num::Interval> projected_bounds(const PolyPatchI& a) {
    using num::Interval;
    Vec3<Interval> box;
    bool first = true;
    for (const auto& cp : a.c) {
        if (!(cp.w.lo > 0.0)) return {Interval::entire(), Interval::entire(), Interval::entire()};
        const Vec3<Interval> pt{cp.x / cp.w, cp.y / cp.w, cp.z / cp.w};
        if (first) { box = pt; first = false; }
        else {
            box.x = num::hull(box.x, pt.x);
            box.y = num::hull(box.y, pt.y);
            box.z = num::hull(box.z, pt.z);
        }
    }
    return box;
}

/// Non-rigorous double version, used only for seeding.
inline void projected_bounds(const PolyPatchD& a, Vec3d& lo, Vec3d& hi) {
    bool first = true;
    for (const auto& cp : a.c) {
        const Vec3d pt{cp.x / cp.w, cp.y / cp.w, cp.z / cp.w};
        if (first) { lo = hi = pt; first = false; }
        else {
            lo = {std::min(lo.x, pt.x), std::min(lo.y, pt.y), std::min(lo.z, pt.z)};
            hi = {std::max(hi.x, pt.x), std::max(hi.y, pt.y), std::max(hi.z, pt.z)};
        }
    }
}

/// Derivatives of a rational surface from derivatives of its homogeneous form
/// (Piegl & Tiller, eq. 4.20). A[jet_index(k,l)] is the (k,l) partial of the
/// homogeneous map; the result is the (k,l) partial of the projected surface.
inline Jet3 rational_jet(const std::array<Vec4d, kJetSize>& A) {
    static constexpr double binom[4][4] = {{1, 0, 0, 0}, {1, 1, 0, 0}, {1, 2, 1, 0}, {1, 3, 3, 1}};
    Jet3 S{};
    const double w = A[0].w;
    for (int k = 0; k <= 3; ++k) {
        for (int l = 0; l <= 3 - k; ++l) {
            Vec3d v = xyz(A[jet_index(k, l)]);
            for (int j = 1; j <= l; ++j)
                v -= S[jet_index(k, l - j)] * (binom[l][j] * A[jet_index(0, j)].w);
            for (int i = 1; i <= k; ++i) {
                v -= S[jet_index(k - i, l)] * (binom[k][i] * A[jet_index(i, 0)].w);
                Vec3d v2{};
                for (int j = 1; j <= l; ++j)
                    v2 += S[jet_index(k - i, l - j)] * (binom[l][j] * A[jet_index(i, j)].w);
                v -= v2 * binom[k][i];
            }
            S[jet_index(k, l)] = v / w;
        }
    }
    return S;
}

/// A rational Bezier patch with cached derivative patches, for point and
/// derivative evaluation in doubles up to total order 3.
class BezierSurface {
    std::array<PolyPatchD, kJetSize> d_; // d_[jet_index(k,l)] = (k,l) homogeneous partial

public:
    BezierSurface() = default;
    explicit BezierSurface(const PolyPatchD& patch) {
        d_[0] = patch;
        for (int n = 1; n <= 3; ++n)
            for (int l = 0; l <= n; ++l) {
                const int k = n - l;
                d_[jet_index(k, l)] = k > 0 ? d_[jet_index(k - 1, l)].deriv_u()
                                            : d_[jet_index(k, l - 1)].deriv_v();
            }
    }
    const PolyPatchD& patch() const { return d_[0]; }

    Vec3d point(double u, double v) const {
        const Vec4d h = d_[0].eval(u, v);
        return {h.x / h.w, h.y / h.w, h.z / h.w};
    }
    /// All partials of total order <= 3, indexed by jet_index.
    Jet3 jet(double u, double v) const {
        std::array<Vec4d, kJetSize> A;
        for (int i = 0; i < kJetSize; ++i) A[i] = d_[i].eval(u, v);
        return rational_jet(A);
    }
};

} // namespace ak::geom
