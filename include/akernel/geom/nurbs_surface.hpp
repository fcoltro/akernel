// akernel — NURBS surfaces: evaluation with derivatives and exact-structure
// conversion to rational Bezier patches (the kernel's carriers).
#pragma once

#include <array>
#include <vector>

#include <akernel/geom/bernstein.hpp>
#include <akernel/geom/bspline.hpp>

namespace ak::geom {

class NurbsSurface {
public:
    int p = 0, q = 0;        ///< degrees in u and v
    int nu = 0, nv = 0;      ///< control-point counts
    std::vector<double> U, V;
    std::vector<Vec4d> cp;   ///< homogeneous (w*x, w*y, w*z, w), row-major in u

    NurbsSurface() = default;
    NurbsSurface(int p_, int q_, int nu_, int nv_, std::vector<double> U_, std::vector<double> V_)
        : p(p_), q(q_), nu(nu_), nv(nv_), U(std::move(U_)), V(std::move(V_)),
          cp(static_cast<std::size_t>(nu_ * nv_)) {}

    Vec4d& at(int i, int j) { return cp[static_cast<std::size_t>(i * nv + j)]; }
    const Vec4d& at(int i, int j) const { return cp[static_cast<std::size_t>(i * nv + j)]; }

    /// All partials of total order <= 3 at (u, v), indexed by jet_index.
    Jet3 jet(double u, double v) const {
        const int su = find_span(nu - 1, p, u, U), sv = find_span(nv - 1, q, v, V);
        const auto Nu = ders_basis_funs(su, u, p, 3, U);
        const auto Nv = ders_basis_funs(sv, v, q, 3, V);
        std::array<Vec4d, kJetSize> A{};
        for (int k = 0; k <= 3; ++k)
            for (int l = 0; l <= 3 - k; ++l) {
                Vec4d s{};
                for (int i = 0; i <= p; ++i)
                    for (int j = 0; j <= q; ++j)
                        s = s + at(su - p + i, sv - q + j) * (Nu[k][i] * Nv[l][j]);
                A[jet_index(k, l)] = s;
            }
        return rational_jet(A);
    }
    Vec3d point(double u, double v) const { return jet(u, v)[0]; }

    struct BezierPiece {
        PolyPatchD patch;
        double u0, u1, v0, v1; ///< parameter rectangle of this piece in the NURBS domain
    };

    /// Split into rational Bezier patches by inserting every interior knot up
    /// to full multiplicity. Done in doubles: the pieces differ from the
    /// input by rounding only, and from here on the pieces are the definition.
    std::vector<BezierPiece> to_bezier() const;
};

namespace detail {

/// Boehm knot insertion of `u` into one row of homogeneous control points.
inline void insert_knot(std::vector<double>& U, std::vector<Vec4d>& P, int p, double u) {
    const int n = static_cast<int>(P.size()) - 1;
    const int k = find_span(n, p, u, U);
    std::vector<Vec4d> Q(P.size() + 1);
    for (int i = 0; i <= n + 1; ++i) {
        if (i <= k - p) Q[i] = P[i];
        else if (i >= k + 1) Q[i] = P[i - 1];
        else {
            const double alpha = (u - U[i]) / (U[i + p] - U[i]);
            Q[i] = P[i] * alpha + P[i - 1] * (1.0 - alpha);
        }
    }
    P = std::move(Q);
    U.insert(U.begin() + k + 1, u);
}

/// Knots to insert so that every interior knot reaches multiplicity p.
inline std::vector<double> bezier_refinement(const std::vector<double>& U, int p) {
    std::vector<double> ins;
    const std::size_t m = U.size();
    std::size_t i = static_cast<std::size_t>(p) + 1;
    while (i < m - p - 1) {
        std::size_t j = i;
        while (j < m - p - 1 && U[j] == U[i]) ++j;
        const int mult = static_cast<int>(j - i);
        for (int r = mult; r < p; ++r) ins.push_back(U[i]);
        i = j;
    }
    return ins;
}

inline std::vector<double> distinct_knots(const std::vector<double>& U) {
    std::vector<double> d;
    for (double u : U)
        if (d.empty() || u != d.back()) d.push_back(u);
    return d;
}

} // namespace detail

inline std::vector<NurbsSurface::BezierPiece> NurbsSurface::to_bezier() const {
    // Refine in u: each column of control points is a curve in u.
    std::vector<std::vector<Vec4d>> cols(nv, std::vector<Vec4d>(nu));
    for (int j = 0; j < nv; ++j)
        for (int i = 0; i < nu; ++i) cols[j][i] = at(i, j);
    const auto ins_u = detail::bezier_refinement(U, p);
    for (int j = 0; j < nv; ++j) {
        std::vector<double> Uc = U;
        for (double u : ins_u) detail::insert_knot(Uc, cols[j], p, u);
    }
    const int nu2 = static_cast<int>(cols[0].size());

    // Refine in v: each row of the refined net is a curve in v.
    std::vector<std::vector<Vec4d>> rows(nu2, std::vector<Vec4d>(nv));
    for (int i = 0; i < nu2; ++i)
        for (int j = 0; j < nv; ++j) rows[i][j] = cols[j][i];
    const auto ins_v = detail::bezier_refinement(V, q);
    for (int i = 0; i < nu2; ++i) {
        std::vector<double> Vc = V;
        for (double v : ins_v) detail::insert_knot(Vc, rows[i], q, v);
    }

    const auto du = detail::distinct_knots(U), dv = detail::distinct_knots(V);
    std::vector<BezierPiece> pieces;
    for (std::size_t a = 0; a + 1 < du.size(); ++a)
        for (std::size_t b = 0; b + 1 < dv.size(); ++b) {
            BezierPiece piece{PolyPatchD(p, q), du[a], du[a + 1], dv[b], dv[b + 1]};
            for (int i = 0; i <= p; ++i)
                for (int j = 0; j <= q; ++j)
                    piece.patch.at(i, j) = rows[a * p + i][b * q + j];
            pieces.push_back(std::move(piece));
        }
    return pieces;
}

} // namespace ak::geom
