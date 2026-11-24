// akernel — generative surfaces: extrusion, revolution, translational sweep,
// ruled loft and skinned loft (docs/FEATURES.md, step F1).
//
// Every constructor returns a NurbsSurface; NurbsSurface::to_bezier() turns it
// into the rational Bezier patches the rest of the kernel works on.
//
// Exactness. Extrusion, revolution by quarter turns and translational sweep
// are closed-form tensor products of the input control data: the result lies
// exactly on the swept point set when the arithmetic on the control data is
// exact (it is for revolution, which only multiplies by 0, +-1, +-2, 2^k, and
// for inputs on round numbers in general). Otherwise the rounded control data
// become the definition, as everywhere in the kernel. Lofts interpolate their
// sections by solving a linear system in doubles; the sections are met up to
// that solve's rounding, which the tests measure.
//
// Not here: sweeps of a profile along a space curve with a moving frame.
// The rotation-minimising frame of a general curve is not rational, so such
// a sweep is either restricted to curves with rational RMF (Farouki's RRMF
// curves) or approximated with a certified deviation bound (FEATURES.md, F1b).
#pragma once

#include <algorithm>
#include <cmath>
#include <optional>
#include <vector>

#include <akernel/geom/bspline.hpp>
#include <akernel/geom/nurbs_curve.hpp>
#include <akernel/geom/nurbs_surface.hpp>
#include <akernel/num/dense.hpp>

namespace ak::geom {

/// Extrusion of `profile` along the vector d: S(u, v) = C(u) + v d, v in [0, 1].
/// Degree (p, 1). A line of constant u is a straight ruling.
inline NurbsSurface extrude(const NurbsCurve& profile, const Vec3d& d) {
    const int n = static_cast<int>(profile.cp.size());
    NurbsSurface s(profile.p, 1, n, 2, profile.U, {0, 0, 1, 1});
    for (int i = 0; i < n; ++i) {
        const Vec4d c = profile.cp[static_cast<std::size_t>(i)];
        s.at(i, 0) = c;
        s.at(i, 1) = {c.x + c.w * d.x, c.y + c.w * d.y, c.z + c.w * d.z, c.w};
    }
    return s;
}

/// Revolution of `profile` about the z axis by `quarters` quarter turns
/// (1..4, counter-clockwise seen from +z). u runs around the axis on the
/// exact rational circle (exact_circle), v along the profile.
///
/// With the circle's homogeneous controls (c, s, h) (cos = c/h, sin = s/h)
/// and a profile control (X, Y, Z, W), the rotated homogeneous point is
/// (X c - Y s, X s + Y c, Z h, W h): a tensor product of degree (2, p). The
/// multiplications are by 0, +-1, +-2 and powers of two, so the surface
/// control data are exact whenever the profile's are not near overflow:
/// every point of the result lies on the surface of revolution exactly.
/// Pieces of a full turn share their boundary rows bitwise.
///
/// A profile control on the axis gives a degenerate (pole) row; the patch is
/// then singular there, which callers must handle (a sphere's poles).
inline NurbsSurface revolve(const NurbsCurve& profile, int quarters) {
    const NurbsCurve circ = exact_circle(1.0, 1.0, quarters); // (c, s, h, h)
    const int nu = static_cast<int>(circ.cp.size());
    const int nv = static_cast<int>(profile.cp.size());
    NurbsSurface s(2, profile.p, nu, nv, circ.U, profile.U);
    for (int i = 0; i < nu; ++i) {
        const Vec4d a = circ.cp[static_cast<std::size_t>(i)];
        const double c = a.x, sn = a.y, h = a.w;
        for (int j = 0; j < nv; ++j) {
            const Vec4d P = profile.cp[static_cast<std::size_t>(j)];
            s.at(i, j) = {P.x * c - P.y * sn, P.x * sn + P.y * c, P.z * h, P.w * h};
        }
    }
    return s;
}

/// Translational sweep: S(u, v) = P(u) + C(v) - C(t0), the profile P moved
/// along the path C without rotation. Degree (p_P, p_C). For rational inputs
/// the homogeneous form is the product
///   (P_i w_j + C_j w_i - C0 w_i w_j, w_i w_j),
/// which divides out to P + C - C0 exactly in real arithmetic.
inline NurbsSurface translational_sweep(const NurbsCurve& profile, const NurbsCurve& path) {
    const Vec3d C0 = path.point(path.t0());
    const int nu = static_cast<int>(profile.cp.size()), nv = static_cast<int>(path.cp.size());
    NurbsSurface s(profile.p, path.p, nu, nv, profile.U, path.U);
    for (int i = 0; i < nu; ++i)
        for (int j = 0; j < nv; ++j) {
            const Vec4d P = profile.cp[static_cast<std::size_t>(i)];
            const Vec4d C = path.cp[static_cast<std::size_t>(j)];
            const double w = P.w * C.w;
            s.at(i, j) = {P.x * C.w + C.x * P.w - C0.x * w, P.y * C.w + C.y * P.w - C0.y * w,
                          P.z * C.w + C.z * P.w - C0.z * w, w};
        }
    return s;
}

/// Ruled surface between two curves: S(u, v) = (1 - v) A(u) + v B(u) in
/// homogeneous form, after making the curves compatible. Degree (p, 1).
/// Also the carrier of a distance-distance chamfer between two rails.
inline NurbsSurface ruled(const NurbsCurve& a, const NurbsCurve& b) {
    const auto cs = make_compatible({a, b});
    const int n = static_cast<int>(cs[0].cp.size());
    NurbsSurface s(cs[0].p, 1, n, 2, cs[0].U, {0, 0, 1, 1});
    for (int i = 0; i < n; ++i) {
        s.at(i, 0) = cs[0].cp[static_cast<std::size_t>(i)];
        s.at(i, 1) = cs[1].cp[static_cast<std::size_t>(i)];
    }
    return s;
}

struct SkinResult {
    NurbsSurface surface;
    std::vector<double> section_params; ///< v at which each section is interpolated
};

/// Skinned loft through `sections` (Piegl & Tiller §10.3): the sections are
/// made compatible, section parameters v_k come from averaged chord lengths
/// between corresponding control points, the v knot vector from averaging
/// them, and each column of homogeneous control points is interpolated by a
/// B-spline of degree q (lowered to sections - 1 if needed).
///
/// Interpolating homogeneous points means the surface passes through each
/// section as a rational curve, up to the solve's rounding. Returns nullopt
/// when the system is singular or an interpolated weight is not positive
/// (a surface with a weight crossing zero has a pole; it is refused, not
/// returned).
inline std::optional<SkinResult> skin(const std::vector<NurbsCurve>& sections, int q) {
    const int K = static_cast<int>(sections.size());
    if (K < 2) return std::nullopt;
    q = std::clamp(q, 1, K - 1);
    const auto cs = make_compatible(sections);
    const int n = static_cast<int>(cs[0].cp.size());

    // Averaged chord-length section parameters.
    std::vector<double> v(static_cast<std::size_t>(K), 0.0);
    bool apart = false;
    for (int i = 0; i < n; ++i) {
        std::vector<double> d(static_cast<std::size_t>(K), 0.0);
        double total = 0.0;
        for (int k = 1; k < K; ++k) {
            const Vec4d a = cs[static_cast<std::size_t>(k - 1)].cp[static_cast<std::size_t>(i)];
            const Vec4d b = cs[static_cast<std::size_t>(k)].cp[static_cast<std::size_t>(i)];
            const Vec3d pa{a.x / a.w, a.y / a.w, a.z / a.w}, pb{b.x / b.w, b.y / b.w, b.z / b.w};
            d[static_cast<std::size_t>(k)] = norm(pb - pa);
            total += d[static_cast<std::size_t>(k)];
        }
        apart = apart || total > 0.0;
        double acc = 0.0;
        for (int k = 1; k < K; ++k) {
            acc += total > 0.0 ? d[static_cast<std::size_t>(k)] / total : 1.0 / (K - 1);
            v[static_cast<std::size_t>(k)] += acc / n;
        }
    }
    if (!apart) return std::nullopt; // every section the same curve
    v.front() = 0.0;
    v.back() = 1.0;
    for (int k = 1; k < K; ++k)
        if (!(v[static_cast<std::size_t>(k)] > v[static_cast<std::size_t>(k - 1)]))
            return std::nullopt; // coincident sections

    // Knots by averaging (P&T eq. 9.8).
    std::vector<double> V(static_cast<std::size_t>(K + q + 1), 0.0);
    for (int j = 1; j <= K - 1 - q; ++j) {
        double s = 0.0;
        for (int i = j; i < j + q; ++i) s += v[static_cast<std::size_t>(i)];
        V[static_cast<std::size_t>(j + q)] = s / q;
    }
    for (int i = K; i <= K + q; ++i) V[static_cast<std::size_t>(i)] = 1.0;

    // Collocation matrix A[k][j] = N_j(v_k), shared by every column.
    std::vector<double> A(static_cast<std::size_t>(K * K), 0.0);
    for (int k = 0; k < K; ++k) {
        const double t = v[static_cast<std::size_t>(k)];
        const int span = find_span(K - 1, q, t, V);
        const auto N = ders_basis_funs(span, t, q, 0, V);
        for (int j = 0; j <= q; ++j) A[static_cast<std::size_t>(k * K + span - q + j)] = N[0][j];
    }
    std::vector<double> B(static_cast<std::size_t>(K * n * 4));
    for (int k = 0; k < K; ++k)
        for (int i = 0; i < n; ++i) {
            const Vec4d c = cs[static_cast<std::size_t>(k)].cp[static_cast<std::size_t>(i)];
            const std::size_t o = static_cast<std::size_t>(k * n * 4 + i * 4);
            B[o] = c.x; B[o + 1] = c.y; B[o + 2] = c.z; B[o + 3] = c.w;
        }
    if (!num::solve_linear(A, B, K, n * 4)) return std::nullopt;

    SkinResult r{NurbsSurface(cs[0].p, q, n, K, cs[0].U, V), v};
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < K; ++j) {
            const std::size_t o = static_cast<std::size_t>(j * n * 4 + i * 4);
            if (!(B[o + 3] > 0.0)) return std::nullopt;
            r.surface.at(i, j) = {B[o], B[o + 1], B[o + 2], B[o + 3]};
        }
    return r;
}

} // namespace ak::geom
