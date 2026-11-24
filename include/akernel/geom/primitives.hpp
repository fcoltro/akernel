// akernel — constructors for common patches, used by tests and benchmarks.
#pragma once

#include <cmath>
#include <functional>
#include <vector>

#include <akernel/num/detmath.hpp>
#include <akernel/geom/bernstein.hpp>
#include <akernel/geom/nurbs_surface.hpp>
#include <akernel/num/dense.hpp>

namespace ak::geom {

/// Planar parallelogram patch: origin + u*du + v*dv.
inline PolyPatchD make_plane(const Vec3d& o, const Vec3d& du, const Vec3d& dv) {
    PolyPatchD s(1, 1);
    auto h = [](const Vec3d& x) { return Vec4d{x.x, x.y, x.z, 1.0}; };
    s.at(0, 0) = h(o);
    s.at(1, 0) = h(o + du);
    s.at(0, 1) = h(o + dv);
    s.at(1, 1) = h(o + du + dv);
    return s;
}

/// Quarter of the cylinder x^2 + y^2 = r^2 for z in [z0, z1], first quadrant.
/// u runs around the arc, v along the axis. The arc uses the homogeneous
/// form (1 - t^2, 2t, 1 + t^2), whose control data are exact in doubles when
/// r is, so the patch lies on the cylinder exactly.
inline PolyPatchD make_cylinder_quarter(double r, double z0, double z1) {
    PolyPatchD s(2, 1);
    const double z[2] = {z0, z1};
    for (int j = 0; j < 2; ++j) {
        s.at(0, j) = {r, 0.0, z[j], 1.0};
        s.at(1, j) = {r, r, z[j], 1.0};
        s.at(2, j) = {0.0, 2.0 * r, 2.0 * z[j], 2.0};
    }
    return s;
}

/// Biquadratic hill over the unit square: z = 16 h u(1-u) v(1-v), peak h.
inline PolyPatchD make_bump(double h) {
    PolyPatchD s(2, 2);
    for (int i = 0; i <= 2; ++i)
        for (int j = 0; j <= 2; ++j)
            s.at(i, j) = {0.5 * i, 0.5 * j, (i == 1 && j == 1) ? 4.0 * h : 0.0, 1.0};
    return s;
}

/// Apply x -> M x + t to a patch (M row-major 3x3). Rounding in the transform
/// becomes part of the new patch's definition.
inline PolyPatchD transformed(const PolyPatchD& a, const double M[9], const Vec3d& t) {
    PolyPatchD r = a;
    for (auto& c : r.c) {
        const Vec4d h = c;
        c.x = M[0] * h.x + M[1] * h.y + M[2] * h.z + t.x * h.w;
        c.y = M[3] * h.x + M[4] * h.y + M[5] * h.z + t.y * h.w;
        c.z = M[6] * h.x + M[7] * h.y + M[8] * h.z + t.z * h.w;
    }
    return r;
}

/// Polynomial patch of bidegree (p, q) interpolating f on a tensor grid.
/// Exact up to rounding when f is a polynomial of bidegree <= (p, q).
inline PolyPatchD fit_patch(const std::function<Vec3d(double, double)>& f, int p, int q) {
    auto nodes = [](int n) {
        std::vector<double> x(n + 1);
        const double pi = 3.14159265358979323846;
        for (int i = 0; i <= n; ++i) x[i] = n == 0 ? 0.5 : 0.5 - 0.5 * num::det::cos(pi * i / n);
        return x;
    };
    auto bern = [](int n, int i, double t) {
        double c = 1.0;
        for (int k = 1; k <= i; ++k) c = c * (n - i + k) / k;
        return c * num::det::ipow(t, i) * num::det::ipow(1.0 - t, n - i);
    };
    const auto un = nodes(p), vn = nodes(q);
    const int nu = p + 1, nv = q + 1;

    // Solve in u for each v-node, then in v for each u-index.
    std::vector<double> Au(static_cast<std::size_t>(nu * nu));
    for (int r = 0; r < nu; ++r)
        for (int i = 0; i < nu; ++i) Au[static_cast<std::size_t>(r * nu + i)] = bern(p, i, un[r]);
    std::vector<double> Bu(static_cast<std::size_t>(nu * nv * 3));
    for (int r = 0; r < nu; ++r)
        for (int s = 0; s < nv; ++s) {
            const Vec3d val = f(un[r], vn[s]);
            Bu[static_cast<std::size_t>(r * nv * 3 + s * 3 + 0)] = val.x;
            Bu[static_cast<std::size_t>(r * nv * 3 + s * 3 + 1)] = val.y;
            Bu[static_cast<std::size_t>(r * nv * 3 + s * 3 + 2)] = val.z;
        }
    num::solve_linear(Au, Bu, nu, nv * 3); // Bu[i][s] = coefficient i in u at v-node s

    std::vector<double> Av(static_cast<std::size_t>(nv * nv));
    for (int s = 0; s < nv; ++s)
        for (int j = 0; j < nv; ++j) Av[static_cast<std::size_t>(s * nv + j)] = bern(q, j, vn[s]);
    std::vector<double> Bv(static_cast<std::size_t>(nv * nu * 3));
    for (int s = 0; s < nv; ++s)
        for (int i = 0; i < nu; ++i)
            for (int k = 0; k < 3; ++k)
                Bv[static_cast<std::size_t>(s * nu * 3 + i * 3 + k)] =
                    Bu[static_cast<std::size_t>(i * nv * 3 + s * 3 + k)];
    num::solve_linear(Av, Bv, nv, nu * 3); // Bv[j][i] = control point (i, j)

    PolyPatchD patch(p, q);
    for (int i = 0; i < nu; ++i)
        for (int j = 0; j < nv; ++j) {
            const std::size_t o = static_cast<std::size_t>(j * nu * 3 + i * 3);
            patch.at(i, j) = {Bv[o], Bv[o + 1], Bv[o + 2], 1.0};
        }
    return patch;
}

} // namespace ak::geom
