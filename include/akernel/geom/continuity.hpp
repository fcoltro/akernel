// akernel — differential invariants of a surface and geometric-continuity
// jumps across a seam (docs/BLUEPRINT.md §1.2 and §1.7).
//
// Everything here is computed from a parametric jet but is independent of the
// parameterisation: the second fundamental form b and its covariant
// derivative C = nabla b are tensors on the surface.
#pragma once

#include <array>
#include <cmath>

#include <akernel/num/detmath.hpp>
#include <akernel/geom/bernstein.hpp>

namespace ak::geom {

/// First, second and third-order invariants at one surface point.
struct SurfaceFrame {
    Vec3d Su, Sv, n;
    double g[2][2];    ///< first fundamental form
    double ginv[2][2];
    double b[2][2];    ///< second fundamental form, b_ab = n . S_ab
    double C[2][2][2]; ///< C[a][b][c] = (nabla_c b)_ab, the Codazzi tensor
};

inline SurfaceFrame surface_frame(const Jet3& J) {
    SurfaceFrame f{};
    f.Su = J[1];
    f.Sv = J[2];
    f.n = normalized(cross(f.Su, f.Sv));
    const Vec3d S1[2] = {f.Su, f.Sv};
    auto S2 = [&](int a, int b2) -> const Vec3d& { return J[3 + a + b2]; };
    auto S3 = [&](int a, int b2, int c) -> const Vec3d& { return J[6 + a + b2 + c]; };

    for (int a = 0; a < 2; ++a)
        for (int c = 0; c < 2; ++c) {
            f.g[a][c] = dot(S1[a], S1[c]);
            f.b[a][c] = dot(f.n, S2(a, c));
        }
    const double det = f.g[0][0] * f.g[1][1] - f.g[0][1] * f.g[1][0];
    f.ginv[0][0] = f.g[1][1] / det;
    f.ginv[1][1] = f.g[0][0] / det;
    f.ginv[0][1] = f.ginv[1][0] = -f.g[0][1] / det;

    // Christoffel symbols Gamma^d_ab = g^{de} S_e . S_ab
    double G[2][2][2];
    for (int d = 0; d < 2; ++d)
        for (int a = 0; a < 2; ++a)
            for (int c = 0; c < 2; ++c) {
                double s = 0.0;
                for (int e = 0; e < 2; ++e) s += f.ginv[d][e] * dot(S1[e], S2(a, c));
                G[d][a][c] = s;
            }

    // (nabla_c b)_ab = n.S_abc - b_cd G^d_ab - b_db G^d_ca - b_ad G^d_cb
    for (int a = 0; a < 2; ++a)
        for (int b2 = 0; b2 < 2; ++b2)
            for (int c = 0; c < 2; ++c) {
                double s = dot(f.n, S3(a, b2, c));
                for (int d = 0; d < 2; ++d)
                    s -= f.b[c][d] * G[d][a][b2] + f.b[d][b2] * G[d][c][a] + f.b[a][d] * G[d][c][b2];
                f.C[a][b2][c] = s;
            }
    return f;
}

/// Parameter-space components of the tangential part of a 3D vector.
inline std::array<double, 2> to_param(const SurfaceFrame& f, const Vec3d& d) {
    const double r0 = dot(f.Su, d), r1 = dot(f.Sv, d);
    return {f.ginv[0][0] * r0 + f.ginv[0][1] * r1, f.ginv[1][0] * r0 + f.ginv[1][1] * r1};
}

inline double second_form(const SurfaceFrame& f, const std::array<double, 2>& x,
                          const std::array<double, 2>& y) {
    double s = 0.0;
    for (int a = 0; a < 2; ++a)
        for (int c = 0; c < 2; ++c) s += f.b[a][c] * x[a] * y[c];
    return s;
}

inline double cubic_form(const SurfaceFrame& f, const std::array<double, 2>& x,
                         const std::array<double, 2>& y, const std::array<double, 2>& z) {
    double s = 0.0;
    for (int a = 0; a < 2; ++a)
        for (int c = 0; c < 2; ++c)
            for (int e = 0; e < 2; ++e) s += f.C[a][c][e] * x[a] * y[c] * z[e];
    return s;
}

inline double mean_curvature(const SurfaceFrame& f) {
    double s = 0.0;
    for (int a = 0; a < 2; ++a)
        for (int c = 0; c < 2; ++c) s += f.ginv[a][c] * f.b[a][c];
    return 0.5 * s;
}

inline double gauss_curvature(const SurfaceFrame& f) {
    return (f.b[0][0] * f.b[1][1] - f.b[0][1] * f.b[1][0]) /
           (f.g[0][0] * f.g[1][1] - f.g[0][1] * f.g[1][0]);
}

/// Jumps of the continuity invariants across a seam at one point.
/// tau is the unit seam tangent, nu the unit conormal (tangent to the
/// surface, normal to the seam). G1 <=> normal_angle = 0. Given G1,
/// G2 <=> db_nn = 0. Given G2, G3 <=> dC_nnn = 0 (blueprint §1.7).
struct SeamJump {
    double normal_angle;             ///< angle between the two normals, radians
    double db_tt, db_tn, db_nn;      ///< jump of b in the (tau, nu) frame
    double dC_ttt, dC_ttn, dC_tnn, dC_nnn; ///< jump of nabla b in the same frame
};

inline SeamJump seam_jump(const Jet3& JA, const Jet3& JB, const Vec3d& seam_tangent) {
    const SurfaceFrame A = surface_frame(JA), B = surface_frame(JB);
    const double c = dot(A.n, B.n);
    const double sgn = c < 0.0 ? -1.0 : 1.0; // opposite parameter orientation flips b and C

    // Frame from side A: tau along the seam, nu = n x tau.
    Vec3d tau = seam_tangent - A.n * dot(A.n, seam_tangent);
    tau = normalized(tau);
    const Vec3d nu = cross(A.n, tau);

    const auto tA = to_param(A, tau), nA = to_param(A, nu);
    const auto tB = to_param(B, tau), nB = to_param(B, nu);

    SeamJump j{};
    j.normal_angle = num::det::atan2(norm(cross(A.n, B.n)), std::fabs(c));
    j.db_tt = sgn * second_form(B, tB, tB) - second_form(A, tA, tA);
    j.db_tn = sgn * second_form(B, tB, nB) - second_form(A, tA, nA);
    j.db_nn = sgn * second_form(B, nB, nB) - second_form(A, nA, nA);
    j.dC_ttt = sgn * cubic_form(B, tB, tB, tB) - cubic_form(A, tA, tA, tA);
    j.dC_ttn = sgn * cubic_form(B, tB, tB, nB) - cubic_form(A, tA, tA, nA);
    j.dC_tnn = sgn * cubic_form(B, tB, nB, nB) - cubic_form(A, tA, nA, nA);
    j.dC_nnn = sgn * cubic_form(B, nB, nB, nB) - cubic_form(A, nA, nA, nA);
    return j;
}

} // namespace ak::geom
