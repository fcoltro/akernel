// akernel — where two rolling-ball blends on a shared face meet: the
// touching point and the local branches of their intersection, proven.
//
// Two blends tangent to one face G along contact curves g1, g2 that cross
// at A touch each other at A (both tangent planes are G's). Their
// intersection curve is singular there: two branches cross, and only one
// lies inside the first blend's arc. intersect_carriers needs transversal
// crossings, so this file proves the local picture instead
// (docs/FEATURES.md §5.1, "Fillet corners"):
//
//  1. A: the contact curves cross in G's parameters, w1(sigma1) = w2(sigma2),
//     a transversal 2-D problem (Krawczyk).
//  2. The branches: on the second blend, h = (distance to the first spine)
//     - r1 and its gradient vanish exactly at A (the blends share their
//     normal there), so along a ray A + rho d in the second blend's scaled
//     parameters, h / rho^2 lies in d^T H d / 2 with H the Hessian of h over
//     the ray's box (Taylor's integral form). Directions d cover the closed
//     half-plane in three charts; where the enclosure excludes zero there is
//     no intersection for any rho in (0, rho_max]; near each root a 1-D
//     Krawczyk proves one branch, a graph over rho.
//  3. Which branch lies on the first blend's arc: k = (P - p1) . (q1 - p1)
//     (p1 the first blend's contact on G, q1 its other contact) vanishes at
//     A, and k / rho lies in grad k . d; inside the arc k >= 0.
#pragma once

#include <algorithm>
#include <cmath>
#include <memory>
#include <numbers>
#include <optional>
#include <vector>

#include <akernel/solve/carrier_ssi.hpp>
#include <akernel/num/detmath.hpp>

namespace ak::solve {

/// A carrier over a rectangle of another's parameters, [u0, u1] x [v0, v1]
/// taken to [0, 1]^2 (u1 < u0 reverses). Enclosures go through the interval
/// image of each box, so they hold for the exact affine map.
class RectCarrier final : public geom::Carrier {
    std::shared_ptr<const geom::Carrier> base_;
    double u0_, u1_, v0_, v1_;
    Interval map(double a, double b, double lo, double hi) const {
        return Interval(lo) + Interval(a, b) * (Interval(hi) - Interval(lo));
    }

public:
    RectCarrier(std::shared_ptr<const geom::Carrier> base, double u0, double u1, double v0, double v1)
        : base_(std::move(base)), u0_(u0), u1_(u1), v0_(v0), v1_(v1) {}
    const geom::Carrier& base() const { return *base_; }
    void to_base(double u, double v, double& bu, double& bv) const {
        bu = u0_ + u * (u1_ - u0_);
        bv = v0_ + v * (v1_ - v0_);
    }
    void eval(double u, double v, Vec3d& P, Vec3d& Pu, Vec3d& Pv) const override {
        double bu, bv;
        to_base(u, v, bu, bv);
        base_->eval(bu, bv, P, Pu, Pv);
        Pu = Pu * (u1_ - u0_);
        Pv = Pv * (v1_ - v0_);
    }
    geom::IVec3 point_at(double u, double v) const override {
        geom::IVec3 P, Pu, Pv;
        if (!jet1_over(u, u, v, v, P, Pu, Pv)) return {Interval::entire(), Interval::entire(), Interval::entire()};
        return P;
    }
    bool jet1_over(double a0, double a1, double b0, double b1, geom::IVec3& P, geom::IVec3& Pu,
                   geom::IVec3& Pv) const override {
        Interval U = map(a0, a1, u0_, u1_), V = map(b0, b1, v0_, v1_);
        if (!base_->jet1_over(U.lo, U.hi, V.lo, V.hi, P, Pu, Pv)) return false;
        const Interval du = Interval(u1_) - Interval(u0_), dv = Interval(v1_) - Interval(v0_);
        Pu = {Pu.x * du, Pu.y * du, Pu.z * du};
        Pv = {Pv.x * dv, Pv.y * dv, Pv.z * dv};
        return true;
    }
    std::vector<double> breaks_u() const override {
        std::vector<double> out{0.0};
        const double lo = std::min(u0_, u1_), hi = std::max(u0_, u1_);
        for (double b : base_->breaks_u())
            if (b > lo && b < hi) out.push_back((b - u0_) / (u1_ - u0_));
        out.push_back(1.0);
        std::sort(out.begin(), out.end());
        return out;
    }
};

/// A blend's contact on the shared face: face coordinate k = alpha[k] x[coord
/// + k] + beta[k] of the spine coordinates x, at theta = theta_c.
struct ContactOnFace {
    const BlendCarrier* blend = nullptr;
    int coord = 0;
    double alpha[2] = {1.0, 1.0}, beta[2] = {0.0, 0.0};
    double theta_c = 0.0;
};

struct CornerPoint {
    Interval s1, s2; ///< sigma of A on each blend
};

namespace corner_detail {

using geom::IVec3;
using num::Interval;

inline IVec3 ivec(const Vec3d& v) { return {Interval(v.x), Interval(v.y), Interval(v.z)}; }
inline Interval idot(const IVec3& a, const IVec3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline IVec3 iscale(const IVec3& a, const Interval& s) { return {a.x * s, a.y * s, a.z * s}; }
inline IVec3 iadd(const IVec3& a, const IVec3& b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }
inline IVec3 isub(const IVec3& a, const IVec3& b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
inline bool inside(const Interval& a, const Interval& b) { return b.lo < a.lo && a.hi < b.hi; }
inline bool finite(const IVec3& v) { return num::is_finite(v.x) && num::is_finite(v.y) && num::is_finite(v.z); }

/// Face coordinates of the contact and their sigma derivatives over a range.
inline bool face_w(const ContactOnFace& f, double s0, double s1, Interval w[2], Interval wp[2]) {
    BlendCarrier::SpineBox b;
    if (!f.blend->spine_box(s0, s1, b)) return false;
    for (int k = 0; k < 2; ++k) {
        w[k] = Interval(f.alpha[k]) * b.x[static_cast<std::size_t>(f.coord + k)] + Interval(f.beta[k]);
        wp[k] = Interval(f.alpha[k]) * b.xp[static_cast<std::size_t>(f.coord + k)];
    }
    return true;
}
inline bool face_w(const ContactOnFace& f, double s, double w[2]) {
    Pt4 x;
    if (!f.blend->spine_at(s, x)) return false;
    for (int k = 0; k < 2; ++k) w[k] = f.alpha[k] * x[static_cast<std::size_t>(f.coord + k)] + f.beta[k];
    return true;
}

} // namespace corner_detail

/// The crossing of the two contact curves on the shared face, from guesses
/// g1, g2 of sigma on each blend: Newton in doubles, then Krawczyk.
inline std::optional<CornerPoint> corner_point(const ContactOnFace& f1, const ContactOnFace& f2, double g1, double g2) {
    using namespace corner_detail;
    double s[2] = {g1, g2};
    for (int it = 0; it < 60; ++it) {
        double w1[2], w2[2], a1[2], a2[2];
        const double h = 1e-7;
        if (!face_w(f1, s[0], w1) || !face_w(f2, s[1], w2)) return std::nullopt;
        double p1[2], p2[2];
        if (!face_w(f1, std::min(1.0, s[0] + h), p1) || !face_w(f2, std::min(1.0, s[1] + h), p2)) return std::nullopt;
        for (int k = 0; k < 2; ++k) {
            a1[k] = (p1[k] - w1[k]) / (std::min(1.0, s[0] + h) - s[0]);
            a2[k] = (p2[k] - w2[k]) / (std::min(1.0, s[1] + h) - s[1]);
        }
        // [a1, -a2] ds = -(w1 - w2)
        const double r0 = w2[0] - w1[0], r1 = w2[1] - w1[1];
        const double det = -a1[0] * a2[1] + a2[0] * a1[1];
        if (!(std::fabs(det) > 0.0)) return std::nullopt;
        const double d0 = (r0 * -a2[1] + a2[0] * r1) / det, d1 = (a1[0] * r1 - a1[1] * r0) / det;
        s[0] = std::clamp(s[0] + d0, 0.0, 1.0);
        s[1] = std::clamp(s[1] + d1, 0.0, 1.0);
        if (std::fabs(d0) + std::fabs(d1) < 1e-16) break;
    }
    // Krawczyk on Z = s +- delta: K = s - Y F(s) + (I - Y J(Z)) (Z - s).
    for (double delta : {1e-12, 1e-10, 1e-8, 1e-6}) {
        const double z0[2] = {s[0] - delta, s[0] + delta}, z1[2] = {s[1] - delta, s[1] + delta};
        Interval W1[2], W1p[2], W2[2], W2p[2], F1[2], F1p[2], F2[2], F2p[2];
        if (!face_w(f1, s[0], s[0], F1, F1p) || !face_w(f2, s[1], s[1], F2, F2p)) return std::nullopt;
        if (!face_w(f1, std::max(0.0, z0[0]), std::min(1.0, z0[1]), W1, W1p) ||
            !face_w(f2, std::max(0.0, z1[0]), std::min(1.0, z1[1]), W2, W2p))
            continue;
        const Interval F[2] = {F1[0] - F2[0], F1[1] - F2[1]};
        const Interval J[2][2] = {{W1p[0], -W2p[0]}, {W1p[1], -W2p[1]}};
        const double m[2][2] = {{num::mid(J[0][0]), num::mid(J[0][1])}, {num::mid(J[1][0]), num::mid(J[1][1])}};
        const double det = m[0][0] * m[1][1] - m[0][1] * m[1][0];
        if (!(std::fabs(det) > 0.0)) return std::nullopt;
        const double Y[2][2] = {{m[1][1] / det, -m[0][1] / det}, {-m[1][0] / det, m[0][0] / det}};
        const Interval Z[2] = {Interval(z0[0], z0[1]), Interval(z1[0], z1[1])};
        const Interval D[2] = {Z[0] - Interval(s[0]), Z[1] - Interval(s[1])};
        Interval K[2];
        for (int i = 0; i < 2; ++i) {
            K[i] = Interval(s[i]) - (Interval(Y[i][0]) * F[0] + Interval(Y[i][1]) * F[1]);
            for (int j = 0; j < 2; ++j) {
                Interval c = Interval(i == j ? 1.0 : 0.0) - (Interval(Y[i][0]) * J[0][j] + Interval(Y[i][1]) * J[1][j]);
                K[i] = K[i] + c * D[j];
            }
        }
        if (inside(K[0], Z[0]) && inside(K[1], Z[1])) return CornerPoint{K[0], K[1]};
    }
    return std::nullopt;
}

/// One branch of the blends' intersection at A, in the second blend's
/// parameters: sigma2 = s2 + k_s rho d_s, theta2 = theta_c + k_t rho d_t with
/// d of chart_dir (chart 0 or 1) at tau = m(rho), rho in (0, rho_max] and
/// m(rho) the unique root in `m`.
struct CornerBranch {
    int chart = 0;
    Interval m;
    bool on_first = false; ///< inside the first blend's arc
};

struct CornerLocal {
    CornerPoint A;
    double k_s = 1.0, k_t = 1.0, rho_max = 0.0;
    std::vector<CornerBranch> branches;
    /// Point of the second blend's parameters on a branch at rho (doubles,
    /// by Newton on h in m from the branch's middle).
    bool point(const ContactOnFace& f1, const ContactOnFace& f2, const CornerBranch& br, double rho, double& s,
               double& t) const;
};

namespace corner_detail {

/// Enclosures of h's gradient sign data over the second blend's box
/// [s0, s1] x [t0, t1] (one spine segment): the Hessian of h and the
/// gradient of k, both in the scaled coordinates (a, b) with sigma =
/// s + k_s a, theta = theta_c + k_t b.
struct BoxData {
    Interval H[2][2];
    Interval gk[2];
};

inline bool box_data(const ContactOnFace& f1, const ContactOnFace& f2, const CornerPoint& A, double s0, double s1,
                     double t0, double t1, double ks, double kt, BoxData& out) {
    const BlendCarrier& b2 = *f2.blend;
    IVec3 P, Ps, Pt, sec[3];
    if (!b2.jet2_over(s0, s1, std::min(t0, t1), std::max(t0, t1), P, Ps, Pt, sec)) return false;
    // The footpoint on the first spine of every point of P: interval Newton
    // on phi(s) = (P - c(s)) . c'(s) from A's sigma, widened until it holds.
    const BlendCarrier& b1 = *f1.blend;
    const double sm = num::mid(A.s1);
    BlendCarrier::SpineBox cm;
    if (!b1.spine_box(sm, sm, cm)) return false;
    Interval S1;
    bool found = false;
    for (double w : {1e-6, 1e-5, 1e-4, 1e-3, 1e-2}) {
        const Interval S(std::max(0.0, A.s1.lo - w), std::min(1.0, A.s1.hi + w));
        BlendCarrier::SpineBox sb;
        if (!b1.spine_box(S.lo, S.hi, sb)) continue;
        const Interval phi_m = idot(isub(P, cm.c), cm.c1);
        const Interval dphi = -idot(sb.c1, sb.c1) + idot(isub(P, sb.c), sb.c2);
        if (!(dphi.hi < 0.0)) continue;
        const Interval N = Interval(sm) - phi_m / dphi;
        if (inside(N, S)) { S1 = N; found = true; break; }
    }
    if (!found) return false;
    BlendCarrier::SpineBox s1b;
    if (!b1.spine_box(S1.lo, S1.hi, s1b)) return false;
    const IVec3 d = isub(P, s1b.c);
    const Interval rho = num::sqrt(idot(d, d));
    if (!(rho.lo > 0.0)) return false;
    const IVec3 u = iscale(d, Interval(1.0) / rho);
    const Interval den = idot(s1b.c1, s1b.c1) - idot(d, s1b.c2);
    if (!(den.lo > 0.0)) return false;
    // Hessian of g = |P - c(s(P))| - r1: (I - u u^T - c' c'^T / den) / rho.
    auto Hg = [&](const IVec3& x, const IVec3& y) {
        return (idot(x, y) - idot(u, x) * idot(u, y) - idot(s1b.c1, x) * idot(s1b.c1, y) / den) / rho;
    };
    const IVec3 J[2] = {iscale(Ps, Interval(ks)), iscale(Pt, Interval(kt))};
    const IVec3 S2[3] = {iscale(sec[0], Interval(ks * ks)), iscale(sec[1], Interval(ks * kt)), iscale(sec[2], Interval(kt * kt))};
    out.H[0][0] = Hg(J[0], J[0]) + idot(u, S2[0]);
    out.H[0][1] = out.H[1][0] = Hg(J[0], J[1]) + idot(u, S2[1]);
    out.H[1][1] = Hg(J[1], J[1]) + idot(u, S2[2]);
    // k = (P - p)(q - p), p the first blend's contact on the face, q the other.
    const bool at_a = f1.theta_c == 0.0;
    const IVec3 p = at_a ? s1b.a : s1b.b, p1 = at_a ? s1b.a1 : s1b.b1;
    const IVec3 q = at_a ? s1b.b : s1b.a, q1 = at_a ? s1b.b1 : s1b.a1;
    const IVec3 gs = iscale(s1b.c1, Interval(1.0) / den); // grad of the footpoint sigma
    const IVec3 qp = isub(q, p), qp1 = isub(q1, p1), Pp = isub(P, p);
    for (int i = 0; i < 2; ++i) {
        const Interval ds = idot(gs, J[i]);
        out.gk[i] = idot(isub(J[i], iscale(p1, ds)), qp) + idot(Pp, qp1) * ds;
    }
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j)
            if (!num::is_finite(out.H[i][j])) return false;
    return num::is_finite(out.gk[0]) && num::is_finite(out.gk[1]);
}

/// Directions of the closed upper half-plane in two charts: chart 0,
/// d = (1 - tau^2, 2 tau) (angles 0 to about 100 degrees for tau in [0,
/// 1.25]); chart 1, d = (tau^2 - 1, 2 tau) (180 down to about 80 degrees).
inline void chart_dir(int chart, const Interval& tau, Interval& da, Interval& db) {
    const Interval t2 = num::sqr(tau);
    da = chart == 0 ? Interval(1.0) - t2 : t2 - Interval(1.0);
    db = Interval(2.0) * tau;
}
inline void chart_ddir(int chart, const Interval& tau, Interval& ea, Interval& eb) {
    ea = chart == 0 ? Interval(-2.0) * tau : Interval(2.0) * tau;
    eb = Interval(2.0);
}
/// The angle of chart direction tau (doubles, monotone in tau).
inline double chart_angle(int chart, double tau) {
    const double a = 2.0 * num::det::atan(tau);
    return chart == 0 ? a : std::numbers::pi - a;
}

} // namespace corner_detail

/// The local picture at A (see the file comment); `rho_max` in the scaled
/// coordinates, which are near arclength on the second blend. Empty
/// optional if anything fails to be proven.
inline std::optional<CornerLocal> corner_local(const ContactOnFace& f1, const ContactOnFace& f2, const CornerPoint& A,
                                               double rho_max) {
    using namespace corner_detail;
    const BlendCarrier& b2 = *f2.blend;
    CornerLocal L;
    L.A = A;
    L.rho_max = rho_max;
    const double s2 = num::mid(A.s2);
    {
        Vec3d P, Ps, Pt;
        b2.eval(s2, f2.theta_c, P, Ps, Pt);
        L.k_s = 1.0 / norm(Ps);
        L.k_t = (f2.theta_c == 0.0 ? 1.0 : -1.0) / norm(Pt);
    }
    double seg_lo, seg_hi;
    b2.segment_range(s2, seg_lo, seg_hi);
    // The box of the sector {A + rho d : rho in [0, rho_max], d in chart(M)}.
    auto data = [&](int chart, const Interval& M, BoxData& out) {
        Interval da, db;
        chart_dir(chart, M, da, db);
        const Interval ra = num::hull(Interval(0.0), Interval(rho_max) * da);
        const Interval rb = num::hull(Interval(0.0), Interval(rho_max) * db);
        const Interval S = A.s2 + Interval(L.k_s) * ra, T = Interval(f2.theta_c) + Interval(L.k_t) * rb;
        if (S.lo < seg_lo || S.hi > seg_hi) return false;
        return box_data(f1, f2, A, S.lo, S.hi, T.lo, T.hi, L.k_s, L.k_t, out);
    };
    auto quad = [](const BoxData& D, const Interval& da, const Interval& db) {
        return (D.H[0][0] * num::sqr(da) + Interval(2.0) * D.H[0][1] * (da * db) + D.H[1][1] * num::sqr(db)) / Interval(2.0);
    };
    // Exclusion over each chart's tau in pieces; the rest are candidate
    // runs, each proven to hold one root by Krawczyk on tau.
    const int n = 64;
    const double tmax = 1.25;
    struct Found {
        CornerBranch br;
        double ang_lo, ang_hi; ///< the run's angles, for the overlap of the charts
    };
    std::vector<Found> found;
    for (int chart = 0; chart < 2; ++chart) {
        std::vector<std::pair<double, double>> runs;
        for (int i = 0; i < n; ++i) {
            const double a = tmax * i / n, b = tmax * (i + 1) / n;
            BoxData D;
            if (!data(chart, Interval(a, b), D)) {
                return std::nullopt;
            }
            Interval da, db;
            chart_dir(chart, Interval(a, b), da, db);
            const Interval q = quad(D, da, db);
            if (q.lo > 0.0 || q.hi < 0.0) continue;
            if (!runs.empty() && runs.back().second == a) runs.back().second = b;
            else runs.push_back({a, b});
        }
        for (const auto& [a, b] : runs) {
            // Krawczyk on tau: every root in M lies in K(M), so M <- K(M) & M
            // keeps all of the run's roots; K(M) inside M's interior proves
            // exactly one.
            Interval M(a, b), K, da, db;
            BoxData D;
            bool proven = false, empty = false, inflated = false;
            for (int it = 0; it < 30 && !proven && !empty; ++it) {
                if (!data(chart, M, D)) return std::nullopt;
                const double mh = num::mid(M);
                Interval ea, eb;
                chart_dir(chart, Interval(mh), da, db);
                const Interval fm = quad(D, da, db);
                chart_dir(chart, M, da, db);
                chart_ddir(chart, M, ea, eb);
                const Interval fp = ea * (D.H[0][0] * da + D.H[0][1] * db) + eb * (D.H[1][0] * da + D.H[1][1] * db);
                if (fp.lo <= 0.0 && fp.hi >= 0.0) break;
                const double y = 1.0 / num::mid(fp);
                K = Interval(mh) - Interval(y) * fm + (Interval(1.0) - Interval(y) * fp) * (M - Interval(mh));
                if (K.lo > M.lo && K.hi < M.hi) { proven = true; break; }
                if (K.lo > M.hi || K.hi < M.lo) { empty = true; break; }
                const Interval Mn(std::max(M.lo, K.lo), std::min(M.hi, K.hi));
                if (Mn.hi - Mn.lo > 0.9 * (M.hi - M.lo)) {
                    // Stalled: inflate around the middle of what is left.
                    const double c = num::mid(Mn), w = 0.25 * (Mn.hi - Mn.lo);
                    M = Interval(c - w, c + w);
                    inflated = true;
                    // Inflating can lose roots: the rest of the run is then
                    // excluded below.
                } else {
                    M = Mn;
                }
            }
            if (empty) continue;
            if (!proven) {
                if (b == tmax) continue; // at the far end of a chart: the other chart's
                return std::nullopt;
            }
            // Whatever of the run lies outside M must hold no root (already
            // proven unless M was inflated).
            if (inflated) for (const auto& piece : {std::make_pair(a, M.lo), std::make_pair(M.hi, b)}) {
                if (!(piece.second > piece.first)) continue;
                const int m = 16;
                for (int i = 0; i < m; ++i) {
                    const double p0 = piece.first + (piece.second - piece.first) * i / m;
                    const double p1 = piece.first + (piece.second - piece.first) * (i + 1) / m;
                    BoxData E;
                    Interval ga, gb;
                    if (!data(chart, Interval(p0, p1), E)) return std::nullopt;
                    chart_dir(chart, Interval(p0, p1), ga, gb);
                    const Interval q = quad(E, ga, gb);
                    if (q.lo <= 0.0 && q.hi >= 0.0) return std::nullopt;
                }
            }
            CornerBranch br;
            br.chart = chart;
            br.m = K;
            const Interval kr = D.gk[0] * da + D.gk[1] * db;
            if (kr.lo > 0.0) br.on_first = true;
            else if (kr.hi < 0.0) br.on_first = false;
            else return std::nullopt;
            const double x0 = chart_angle(chart, a), x1 = chart_angle(chart, b);
            found.push_back({br, std::min(x0, x1), std::max(x0, x1)});
        }
    }
    // The same root proven in both charts: a root of one chart whose angles
    // lie inside the other's proven run (unique there).
    for (std::size_t i = 0; i < found.size(); ++i) {
        bool dup = false;
        const double r0 = chart_angle(found[i].br.chart, found[i].br.m.lo), r1 = chart_angle(found[i].br.chart, found[i].br.m.hi);
        for (std::size_t j = 0; j < i && !dup; ++j)
            dup = found[j].br.chart != found[i].br.chart && std::min(r0, r1) > found[j].ang_lo + 1e-12 &&
                  std::max(r0, r1) < found[j].ang_hi - 1e-12;
        if (!dup) L.branches.push_back(found[i].br);
    }
    return L;
}

inline bool CornerLocal::point(const ContactOnFace& f1, const ContactOnFace& f2, const CornerBranch& br, double rho,
                               double& s, double& t) const {
    const BlendCarrier& b1 = *f1.blend;
    const BlendCarrier& b2 = *f2.blend;
    auto at = [&](double m, double& ss, double& tt) {
        const double da = br.chart == 0 ? 1.0 - m * m : m * m - 1.0, db = 2.0 * m;
        ss = num::mid(A.s2) + k_s * rho * da;
        tt = f2.theta_c + k_t * rho * db;
    };
    // h by footpoint Newton on the first spine.
    auto h = [&](double m) {
        double ss, tt;
        at(m, ss, tt);
        Vec3d P, Ps, Pt;
        b2.eval(ss, tt, P, Ps, Pt);
        double sg = num::mid(A.s1);
        Vec3d c{}, a, b;
        for (int it = 0; it < 40; ++it) {
            BlendCarrier::Section S0 = b1.section(sg);
            const double dh = 1e-7;
            const double sg2 = std::min(1.0, sg + dh);
            BlendCarrier::Section S1 = b1.section(sg2);
            c = S0.centre;
            const Vec3d cp = (S1.centre - S0.centre) / (sg2 - sg);
            const double step = dot(P - c, cp) / dot(cp, cp);
            sg = std::clamp(sg + step, 0.0, 1.0);
            if (std::fabs(step) < 1e-16) break;
        }
        (void)a;
        (void)b;
        return norm(P - b1.section(sg).centre) - b1.radius();
    };
    double m = num::mid(br.m);
    for (int it = 0; it < 60; ++it) {
        const double e = 1e-7, f = h(m), g = (h(m + e) - h(m - e)) / (2 * e);
        if (!(std::fabs(g) > 0.0)) break;
        const double step = f / g;
        m -= step;
        if (std::fabs(step) < 1e-15) break;
    }
    at(m, s, t);
    return std::isfinite(s) && std::isfinite(t);
}

} // namespace ak::solve
