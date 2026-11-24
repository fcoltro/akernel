// akernel — procedural edges that are traced intersection curves of two
// carriers (docs/FEATURES.md §5.1). A blend cut by a plane, or any face
// meeting a procedural face along a curve that is not a side of its square,
// has such an edge: the curve is the certified intersection of the two
// surfaces (solve/carrier_ssi.hpp, intersect_carriers), kept as its tube
// segments, and each face reads it in its own parameters through the
// curve's coordinates. One definition, no fitted trimming curve.
#pragma once

#include <algorithm>
#include <cmath>
#include <memory>
#include <vector>

#include <akernel/op/brep.hpp>
#include <akernel/solve/carrier_ssi.hpp>

namespace ak::op {

/// A traced curve in the parameters (u, v, s, t) of two carriers,
/// parametrised by sigma in [0, 1]: each tube segment takes a share of
/// sigma proportional to its length in space, and inside a segment sigma is
/// linear in the segment's slice value.
struct CurvePath {
    std::vector<solve::TubeSegment> segs; ///< in curve order
    std::vector<double> brk;              ///< sigma at the segments' ends

    int segment_of(double sigma) const {
        const auto it = std::upper_bound(brk.begin() + 1, brk.end() - 1, sigma);
        return static_cast<int>(it - brk.begin()) - 1;
    }
    double slice(int k, double sigma) const {
        const solve::TubeSegment& s = segs[static_cast<std::size_t>(k)];
        const double a = brk[static_cast<std::size_t>(k)], b = brk[static_cast<std::size_t>(k) + 1];
        const double f = std::clamp((sigma - a) / (b - a), 0.0, 1.0);
        return (s.reversed ? 1.0 - f : f) * s.len;
    }
    template <class Sys>
    bool point(const Sys& sys, double sigma, solve::Pt4& x) const {
        const int k = segment_of(sigma);
        const solve::TubeSegment& seg = segs[static_cast<std::size_t>(k)];
        const double c = slice(k, sigma);
        for (int i = 0; i < 4; ++i) x[i] = seg.a[i] + c * seg.m[i];
        return solve::newton_slice(sys, seg.m, seg.a, c, x);
    }
    /// Shares of sigma by length in space, measured on carrier `a`.
    template <class Sys>
    void make_breaks(const Sys& sys, const geom::Carrier& a) {
        const std::size_t n = segs.size();
        brk.assign(n + 1, 0.0);
        for (std::size_t k = 0; k < n; ++k) brk[k + 1] = static_cast<double>(k + 1) / static_cast<double>(n);
        std::vector<double> w(n, 0.0);
        double total = 0.0;
        for (std::size_t k = 0; k < n; ++k) {
            const solve::TubeSegment& seg = segs[k];
            Vec3d prev{};
            for (int i = 0; i <= 4; ++i) {
                const double c = (seg.reversed ? 1.0 - i / 4.0 : i / 4.0) * seg.len;
                solve::Pt4 x;
                for (int d = 0; d < 4; ++d) x[d] = seg.a[d] + c * seg.m[d];
                solve::newton_slice(sys, seg.m, seg.a, c, x);
                Vec3d P, Pu, Pv;
                a.eval(x[0], x[1], P, Pu, Pv);
                if (i > 0) w[k] += norm(P - prev);
                prev = P;
            }
            if (!std::isfinite(w[k])) return; // keep the even shares
            total += w[k];
        }
        if (!(total > 0.0)) return;
        double acc = 0.0, sum = 0.0;
        for (double x : w) sum += std::max(x, 1e-6 * total);
        for (std::size_t k = 0; k < n; ++k) {
            acc += std::max(w[k], 1e-6 * total);
            brk[k + 1] = acc / sum;
        }
        brk.back() = 1.0;
    }
};

/// An edge on the intersection of two carriers: the part sigma in [from, to]
/// of a traced curve (to < from runs it backwards), tau = 0 at `from`.
/// Each side names the Brep surface it is: a polynomial carrier (index into
/// Brep::carriers) or a procedural one (index into Brep::procs).
class TracedEdge final : public ProcEdge {
public:
    std::shared_ptr<const geom::Carrier> a, b;
    int carrier[2] = {-1, -1}; ///< Brep::carriers index of side a / b, or -1
    int proc[2] = {-1, -1};    ///< Brep::procs index of side a / b, or -1
    CurvePath path;
    double from = 0.0, to = 1.0;

    solve::CarrierPairSystem sys() const { return solve::CarrierPairSystem(*a, *b); }
    double sigma(double tau) const { return from + tau * (to - from); }
    bool at(double tau, solve::Pt4& x) const { return path.point(sys(), sigma(tau), x); }

    Vec3d point(double tau) const override {
        solve::Pt4 x;
        at(tau, x);
        Vec3d P, Pu, Pv;
        a->eval(x[0], x[1], P, Pu, Pv);
        return P;
    }
    std::vector<double> breaks() const override {
        const double lo = std::min(from, to), hi = std::max(from, to);
        std::vector<double> t{0.0};
        for (double s : path.brk)
            if (s > lo && s < hi) t.push_back((s - from) / (to - from));
        t.push_back(1.0);
        std::sort(t.begin(), t.end());
        return t;
    }
    bool coords(int side, double tau, double& u, double& v) const {
        solve::Pt4 x;
        if (!at(tau, x)) return false;
        u = x[static_cast<std::size_t>(2 * side)];
        v = x[static_cast<std::size_t>(2 * side + 1)];
        return true;
    }
    bool on_carrier(int c, double tau, double& u, double& v, double& du, double& dv) const override {
        const int side = carrier[0] == c ? 0 : carrier[1] == c ? 1 : -1;
        if (side < 0 || !coords(side, tau, u, v)) return false;
        return piece_derivative(breaks(), tau, [&](double t, double& p, double& q) { return coords(side, t, p, q); }, du, dv);
    }
    bool on_proc(int p, double tau, double& s, double& t) const override {
        const int side = proc[0] == p ? 0 : proc[1] == p ? 1 : -1;
        return side >= 0 && coords(side, tau, s, t);
    }
};

} // namespace ak::op
