// akernel — small dense linear algebra in doubles (row-major storage).
// Used for collocation, constraint elimination and the variational solver.
#pragma once

#include <cmath>
#include <cstddef>
#include <vector>

namespace ak::num {

/// Solve A X = B in place by Gaussian elimination with partial pivoting.
/// A is n x n, B is n x m (both row-major). Returns false if A is singular.
inline bool solve_linear(std::vector<double>& A, std::vector<double>& B, int n, int m) {
    auto a = [&](int i, int j) -> double& { return A[static_cast<std::size_t>(i * n + j)]; };
    auto b = [&](int i, int j) -> double& { return B[static_cast<std::size_t>(i * m + j)]; };
    for (int k = 0; k < n; ++k) {
        int piv = k;
        for (int i = k + 1; i < n; ++i)
            if (std::fabs(a(i, k)) > std::fabs(a(piv, k))) piv = i;
        if (a(piv, k) == 0.0) return false;
        if (piv != k) {
            for (int j = 0; j < n; ++j) std::swap(a(k, j), a(piv, j));
            for (int j = 0; j < m; ++j) std::swap(b(k, j), b(piv, j));
        }
        for (int i = k + 1; i < n; ++i) {
            const double f = a(i, k) / a(k, k);
            if (f == 0.0) continue;
            for (int j = k; j < n; ++j) a(i, j) -= f * a(k, j);
            for (int j = 0; j < m; ++j) b(i, j) -= f * b(k, j);
        }
    }
    for (int k = n - 1; k >= 0; --k) {
        for (int j = 0; j < m; ++j) {
            double s = b(k, j);
            for (int i = k + 1; i < n; ++i) s -= a(k, i) * b(i, j);
            b(k, j) = s / a(k, k);
        }
    }
    return true;
}

/// Cholesky factorisation A = L L^T in place (lower triangle).
/// Returns false if A is not positive definite.
inline bool cholesky(std::vector<double>& A, int n) {
    auto a = [&](int i, int j) -> double& { return A[static_cast<std::size_t>(i * n + j)]; };
    for (int j = 0; j < n; ++j) {
        double d = a(j, j);
        for (int k = 0; k < j; ++k) d -= a(j, k) * a(j, k);
        if (!(d > 0.0)) return false;
        const double ljj = std::sqrt(d);
        a(j, j) = ljj;
        for (int i = j + 1; i < n; ++i) {
            double s = a(i, j);
            for (int k = 0; k < j; ++k) s -= a(i, k) * a(j, k);
            a(i, j) = s / ljj;
        }
    }
    return true;
}

/// Solve L L^T X = B in place, given the factor from cholesky(). B is n x m.
inline void cholesky_solve(const std::vector<double>& L, std::vector<double>& B, int n, int m) {
    auto l = [&](int i, int j) { return L[static_cast<std::size_t>(i * n + j)]; };
    auto b = [&](int i, int j) -> double& { return B[static_cast<std::size_t>(i * m + j)]; };
    for (int c = 0; c < m; ++c) {
        for (int i = 0; i < n; ++i) {
            double s = b(i, c);
            for (int k = 0; k < i; ++k) s -= l(i, k) * b(k, c);
            b(i, c) = s / l(i, i);
        }
        for (int i = n - 1; i >= 0; --i) {
            double s = b(i, c);
            for (int k = i + 1; k < n; ++k) s -= l(k, i) * b(k, c);
            b(i, c) = s / l(i, i);
        }
    }
}

} // namespace ak::num
