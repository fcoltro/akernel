// akernel — watertight tessellation of a Brep with trimmed faces.
//
// Crack-free by construction (docs/BLUEPRINT.md §5.1, ak::viz): every edge is
// sampled exactly once, and the faces on its two sides both use those same
// sample points as the boundary of their triangulations. No point is ever
// welded by distance, so the mesh is closed whenever the Brep is closed,
// whatever the accuracy asked for.
//
//   edges    sampled adaptively until the chord deviation is below
//            `chord`, and no longer than the interior spacing of the faces
//            that use them. Side edges are sampled on their carrier,
//            intersection edges on their certified tube (a Newton point per
//            slice), which also gives their parameters on both carriers.
//   faces    triangulated in the carrier's parameter plane: the loops become
//            the constraint polygon, interior points come from a grid whose
//            spacing follows the carrier's second derivatives, and a
//            constrained Delaunay triangulation joins them.
//
// The triangulation decides with the exact predicates orient2d and incircle,
// so its combinatorics are valid for any input that is valid. Parameters are
// scaled by a power of two (exact) to make triangles closer to isotropic in
// space. What is numerical is the mesh's geometry: its deviation from the
// surfaces is measured (mesh_deviation), not bounded.
#pragma once

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <deque>
#include <limits>
#include <map>
#include <ostream>
#include <set>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <akernel/num/predicates.hpp>
#include <akernel/op/brep.hpp>

namespace ak::op {

struct TessOptions {
    double chord = 1e-3;    ///< wanted deviation of the mesh from the surfaces (model units)
    double max_edge = 0.0;  ///< longest triangle side wanted (0: no limit)
    int max_depth = 16;     ///< subdivision depth limit when sampling an edge
    int max_grid = 200;     ///< limit on interior grid lines per direction and face
};

struct Mesh {
    std::vector<Vec3d> points;
    std::vector<std::array<int, 3>> tris;        ///< counter-clockwise seen from outside
    std::vector<int> tri_face;                   ///< Brep face of every triangle
    std::vector<std::array<double, 6>> tri_uv;   ///< carrier parameters (u0 v0 u1 v1 u2 v2)
    std::vector<std::vector<int>> edges;         ///< per Brep edge: its points from v0 to v1
};

struct TessStats {
    int edge_points = 0;      ///< points inside edges (not counting Brep vertices)
    int interior_points = 0;  ///< points inside faces
    int triangles = 0;
    long long flips = 0;
    long long constraint_flips = 0;
    int lens_splits = 0;      ///< chords two faces had in common, split (see tessellate)
    int axis_points = 0;      ///< edge samples taken on a tube segment's axis where Newton failed (see tessellate)
    double ms = 0.0;
};

struct TessResult {
    bool ok = false;
    std::string error;
    Mesh mesh;
    TessStats stats;
};

namespace tess_detail {

using P2 = std::array<double, 2>;
using num::Sign;

inline double seg_dist(const Vec3d& p, const Vec3d& a, const Vec3d& b) {
    const Vec3d ab = b - a, ap = p - a;
    const double l2 = dot(ab, ab);
    double t = l2 > 0.0 ? dot(ap, ab) / l2 : 0.0;
    t = std::min(1.0, std::max(0.0, t));
    return norm(ap - ab * t);
}

inline double seg_dist2(const P2& p, const P2& a, const P2& b) {
    const double abx = b[0] - a[0], aby = b[1] - a[1], apx = p[0] - a[0], apy = p[1] - a[1];
    const double l2 = abx * abx + aby * aby;
    double t = l2 > 0.0 ? (apx * abx + apy * aby) / l2 : 0.0;
    t = std::min(1.0, std::max(0.0, t));
    const double dx = apx - t * abx, dy = apy - t * aby;
    return std::sqrt(dx * dx + dy * dy);
}

/// Constrained Delaunay triangulation by incremental insertion (Lawson
/// flips) and constraint recovery by flips (Sloan). Topological decisions
/// use exact orient2d; Delaunay flips use exact incircle.
class Cdt {
public:
    struct Tri {
        std::array<int, 3> v{};      ///< counter-clockwise
        std::array<int, 3> n{};      ///< n[i]: neighbour across the edge opposite v[i], or -1
        std::array<char, 3> c{};     ///< c[i]: that edge is a constraint
    };
    std::vector<P2> P;
    std::vector<Tri> T;
    std::vector<int> vt;  ///< a triangle incident to each vertex
    long long flips = 0, constraint_flips = 0;
    std::string error;

    void init(const P2& lo, const P2& hi) {
        const double cx = 0.5 * (lo[0] + hi[0]), cy = 0.5 * (lo[1] + hi[1]);
        const double r = std::max(hi[0] - lo[0], hi[1] - lo[1]) + 1.0;
        P = {{cx - 30.0 * r, cy - 10.0 * r}, {cx + 30.0 * r, cy - 10.0 * r}, {cx, cy + 30.0 * r}};
        T.assign(1, Tri{{0, 1, 2}, {-1, -1, -1}, {0, 0, 0}});
        vt = {0, 0, 0};
        last_ = 0;
    }

    Sign orient(int a, int b, int c) const { return num::orient2d(P[a].data(), P[b].data(), P[c].data()); }

    /// Insert a point. Returns its vertex index, or the index of an existing
    /// vertex at the same place, or -1 on failure.
    int insert(const P2& q) {
        int t = -1, edge = -1, vertex = -1;
        if (!locate(q, t, edge, vertex)) return -1;
        if (vertex >= 0) return vertex;
        const int p = static_cast<int>(P.size());
        P.push_back(q);
        vt.push_back(t);
        if (edge >= 0) split_edge(t, edge, p);
        else split_tri(t, p);
        return p;
    }

    /// Make segment a-b an edge of the triangulation and mark it.
    bool constrain(int a, int b) {
        if (a == b) return true;
        int t, i;
        if (!find_edge(a, b, t, i)) {
            std::deque<std::pair<int, int>> crossed;
            if (!crossing_edges(a, b, crossed)) return false;
            long long guard = 0;
            const long long limit = 1000 + 50LL * static_cast<long long>(crossed.size()) *
                                               static_cast<long long>(crossed.size() + 4);
            while (!crossed.empty()) {
                if (++guard > limit) { error = "constraint recovery did not converge"; return false; }
                const auto [q, r] = crossed.front();
                crossed.pop_front();
                int te, ie;
                if (!find_edge(q, r, te, ie)) { error = "lost an edge during constraint recovery"; return false; }
                if (T[te].c[ie]) { error = "two boundary segments of a face cross in its parameters"; return false; }
                const int u = T[te].n[ie];
                const int p = T[te].v[ie];
                const int d = third(u, te);
                if (strictly_opposite(orient(p, d, q), orient(p, d, r))) {
                    flip(te, ie);
                    ++constraint_flips;
                    if (strictly_opposite(orient(a, b, p), orient(a, b, d)) &&
                        strictly_opposite(orient(p, d, a), orient(p, d, b)))
                        crossed.emplace_back(p, d);
                } else {
                    crossed.emplace_back(q, r);
                }
            }
            if (!find_edge(a, b, t, i)) { error = "constraint recovery failed"; return false; }
        }
        T[t].c[i] = 1;
        const int u = T[t].n[i];
        if (u >= 0)
            for (int j = 0; j < 3; ++j)
                if (T[u].n[j] == t) T[u].c[j] = 1;
        return true;
    }

    /// Restore the Delaunay property on edges that are not constraints.
    void delaunay_pass() {
        std::deque<std::pair<int, int>> q;
        for (int t = 0; t < static_cast<int>(T.size()); ++t)
            for (int i = 0; i < 3; ++i) q.emplace_back(t, i);
        long long guard = 0;
        const long long limit = 100LL * static_cast<long long>(T.size()) + 1000;
        while (!q.empty() && ++guard < limit) {
            const auto [t, i] = q.front();
            q.pop_front();
            if (!should_flip(t, i)) continue;
            const int u = T[t].n[i];
            const int p = T[t].v[i], d = third(u, t);
            const int a = T[t].v[(i + 1) % 3], b = T[t].v[(i + 2) % 3];
            if (!strictly_opposite(orient(p, d, a), orient(p, d, b))) continue;
            flip(t, i);
            ++flips;
            q.emplace_back(t, 0);
            q.emplace_back(t, 2);
            q.emplace_back(u, 0);
            q.emplace_back(u, 1);
        }
    }

    /// Triangles inside the constraint polygon (odd number of constraints
    /// crossed from the outside). False if the parity is inconsistent.
    bool inside(std::vector<char>& in) {
        in.assign(T.size(), 0);
        // Start from the triangle at a corner of the super triangle: outside.
        const int start = vt[0];
        std::vector<signed char> state(T.size(), -1);
        state[static_cast<std::size_t>(start)] = 0;
        std::vector<int> stack{start};
        while (!stack.empty()) {
            const int t = stack.back();
            stack.pop_back();
            for (int i = 0; i < 3; ++i) {
                const int u = T[t].n[i];
                if (u < 0) continue;
                const signed char s = static_cast<signed char>(state[static_cast<std::size_t>(t)] ^ (T[t].c[i] ? 1 : 0));
                if (state[static_cast<std::size_t>(u)] < 0) {
                    state[static_cast<std::size_t>(u)] = s;
                    stack.push_back(u);
                } else if (state[static_cast<std::size_t>(u)] != s) {
                    error = "the loops of a face do not bound a region consistently";
                    return false;
                }
            }
        }
        for (std::size_t t = 0; t < T.size(); ++t) {
            in[t] = state[t] == 1;
            if (in[t] && (T[t].v[0] < 3 || T[t].v[1] < 3 || T[t].v[2] < 3)) {
                error = "a face's loops are open";
                return false;
            }
        }
        return true;
    }

private:
    int last_ = 0;

    static bool strictly_opposite(Sign s, Sign t) {
        return (s == Sign::Positive && t == Sign::Negative) || (s == Sign::Negative && t == Sign::Positive);
    }

    int third(int u, int t) const {
        for (int j = 0; j < 3; ++j)
            if (T[u].n[j] == t) return T[u].v[j];
        return -1;
    }

    void repoint(int nb, int from, int to) {
        if (nb < 0) return;
        for (int j = 0; j < 3; ++j)
            if (T[nb].n[j] == from) T[nb].n[j] = to;
    }

    void set(int t, int a, int b, int c, int na, int nb, int nc, char ca, char cb, char cc) {
        T[t].v = {a, b, c};
        T[t].n = {na, nb, nc};
        T[t].c = {ca, cb, cc};
        vt[a] = vt[b] = vt[c] = t;
    }

    bool locate(const P2& q, int& t, int& edge, int& vertex) {
        edge = vertex = -1;
        t = last_;
        const long long guard_limit = 4LL * static_cast<long long>(T.size()) + 100;
        for (long long guard = 0; guard < guard_limit; ++guard) {
            const Tri& tr = T[t];
            int next = -2, zeros = 0, z[2] = {-1, -1};
            for (int k0 = 0; k0 < 3; ++k0) {
                const int k = static_cast<int>((k0 + guard) % 3);
                const Sign s = num::orient2d(P[tr.v[(k + 1) % 3]].data(), P[tr.v[(k + 2) % 3]].data(), q.data());
                if (s == Sign::Negative) { next = tr.n[k]; break; }
                if (s == Sign::Zero) { if (zeros < 2) z[zeros] = k; ++zeros; }
            }
            if (next == -1) { error = "a point lies outside the triangulation"; return false; }
            if (next >= 0) { t = next; continue; }
            last_ = t;
            if (zeros == 1) edge = z[0];
            else if (zeros >= 2) vertex = tr.v[3 - z[0] - z[1]];
            return true;
        }
        error = "point location did not terminate";
        return false;
    }

    void split_tri(int t, int p) {
        const Tri old = T[t];
        const int a = old.v[0], b = old.v[1], c = old.v[2];
        const int t1 = static_cast<int>(T.size()), t2 = t1 + 1;
        T.resize(T.size() + 2);
        set(t, p, b, c, old.n[0], t1, t2, old.c[0], 0, 0);
        set(t1, p, c, a, old.n[1], t2, t, old.c[1], 0, 0);
        set(t2, p, a, b, old.n[2], t, t1, old.c[2], 0, 0);
        repoint(old.n[1], t, t1);
        repoint(old.n[2], t, t2);
        legalize(t);
        legalize(t1);
        legalize(t2);
    }

    void split_edge(int t, int k, int p) {
        const Tri A0 = T[t];
        const int c = A0.v[k], a = A0.v[(k + 1) % 3], b = A0.v[(k + 2) % 3];
        const int u = A0.n[k];
        const int N_bc = A0.n[(k + 1) % 3], N_ca = A0.n[(k + 2) % 3];
        const char C_bc = A0.c[(k + 1) % 3], C_ca = A0.c[(k + 2) % 3];
        const int B = static_cast<int>(T.size());
        if (u < 0) {
            // On the outer boundary: cannot happen inside the super triangle.
            T.resize(T.size() + 1);
            set(t, p, c, a, N_ca, -1, B, C_ca, 0, 0);
            set(B, p, b, c, N_bc, t, -1, C_bc, 0, 0);
            repoint(N_bc, t, B);
            legalize(t);
            legalize(B);
            return;
        }
        const Tri U0 = T[u];
        int j = 0;
        while (U0.n[j] != t) ++j;
        const int d = U0.v[j];
        const int N_ad = U0.n[(j + 1) % 3], N_db = U0.n[(j + 2) % 3];
        const char C_ad = U0.c[(j + 1) % 3], C_db = U0.c[(j + 2) % 3];
        const char C_ab = A0.c[k];
        const int D = B + 1;
        T.resize(T.size() + 2);
        set(t, p, c, a, N_ca, D, B, C_ca, C_ab, 0);
        set(B, p, b, c, N_bc, t, u, C_bc, 0, C_ab);
        set(u, p, d, b, N_db, B, D, C_db, C_ab, 0);
        set(D, p, a, d, N_ad, u, t, C_ad, 0, C_ab);
        repoint(N_bc, t, B);
        repoint(N_ad, u, D);
        legalize(t);
        legalize(B);
        legalize(u);
        legalize(D);
    }

    bool should_flip(int t, int i) const {
        if (T[t].c[i]) return false;
        const int u = T[t].n[i];
        if (u < 0) return false;
        const int d = third(u, t);
        return num::incircle(P[T[t].v[0]].data(), P[T[t].v[1]].data(), P[T[t].v[2]].data(), P[d].data()) ==
               Sign::Positive;
    }

    /// Legalize the edge opposite vertex 0 (the new point) of t, recursively.
    void legalize(int t0) {
        std::vector<int> stack{t0};
        while (!stack.empty()) {
            const int t = stack.back();
            stack.pop_back();
            if (!should_flip(t, 0)) continue;
            const int u = T[t].n[0];
            flip(t, 0);
            ++flips;
            stack.push_back(t);
            stack.push_back(u);
        }
    }

    /// Flip the edge opposite v[i] of t. Afterwards t = (p, a, d) and
    /// u = (p, d, b), both with p = old v[i] at index 0.
    void flip(int t, int i) {
        const Tri A = T[t];
        const int p = A.v[i], a = A.v[(i + 1) % 3], b = A.v[(i + 2) % 3];
        const int u = A.n[i];
        const Tri U = T[u];
        int j = 0;
        while (U.n[j] != t) ++j;
        const int d = U.v[j];
        const int N_pa = A.n[(i + 2) % 3], N_bp = A.n[(i + 1) % 3];
        const char C_pa = A.c[(i + 2) % 3], C_bp = A.c[(i + 1) % 3];
        const int N_ad = U.n[(j + 1) % 3], N_db = U.n[(j + 2) % 3];
        const char C_ad = U.c[(j + 1) % 3], C_db = U.c[(j + 2) % 3];
        set(t, p, a, d, N_ad, u, N_pa, C_ad, 0, C_pa);
        set(u, p, d, b, N_db, N_bp, t, C_db, C_bp, 0);
        repoint(N_ad, u, t);
        repoint(N_bp, t, u);
    }

    /// The triangle and index of edge q-r (index of the vertex opposite it).
    bool find_edge(int q, int r, int& t, int& i) const {
        auto check = [&](int tt) {
            for (int k = 0; k < 3; ++k)
                if (T[tt].v[k] == q) {
                    if (T[tt].v[(k + 1) % 3] == r) { t = tt; i = (k + 2) % 3; return true; }
                    if (T[tt].v[(k + 2) % 3] == r) { t = tt; i = (k + 1) % 3; return true; }
                }
            return false;
        };
        const int start = vt[q];
        for (int dir = 0; dir < 2; ++dir) {
            int tt = start;
            for (int guard = 0; guard < 100000; ++guard) {
                if (check(tt)) return true;
                int k = 0;
                while (T[tt].v[k] != q) ++k;
                const int nt = T[tt].n[dir == 0 ? (k + 2) % 3 : (k + 1) % 3];
                if (nt < 0 || nt == start) break;
                tt = nt;
            }
        }
        return false;
    }

    /// Edges crossed by segment a-b, walking from a.
    bool crossing_edges(int a, int b, std::deque<std::pair<int, int>>& out) {
        // Find the triangle around a whose wedge contains the direction to b.
        int q = -1, r = -1;
        const int start = vt[a];
        int tt = start;
        for (int guard = 0; guard < 100000; ++guard) {
            int k = 0;
            while (T[tt].v[k] != a) ++k;
            const int qq = T[tt].v[(k + 1) % 3], rr = T[tt].v[(k + 2) % 3];
            const Sign s1 = orient(a, b, qq), s2 = orient(a, b, rr);
            if ((s1 == Sign::Zero && ahead(a, b, qq)) || (s2 == Sign::Zero && ahead(a, b, rr))) {
                error = "a point lies on a boundary segment of a face";
                return false;
            }
            if (s1 == Sign::Negative && s2 == Sign::Positive) { q = qq; r = rr; break; }
            tt = T[tt].n[(k + 2) % 3];
            if (tt < 0 || tt == start) break;
        }
        if (q < 0) { error = "could not start a boundary segment"; return false; }
        int prev = a;
        out.emplace_back(q, r);
        for (int guard = 0; guard < 1000000; ++guard) {
            int te, ie;
            if (!find_edge(q, r, te, ie)) { error = "lost an edge while walking a segment"; return false; }
            const int u = T[te].n[ie];
            const int s = T[te].v[ie] == prev ? third(u, te) : T[te].v[ie];
            if (s < 0) { error = "a boundary segment leaves the triangulation"; return false; }
            if (s == b) return true;
            const Sign o = orient(a, b, s);
            if (o == Sign::Zero) { error = "a point lies on a boundary segment of a face"; return false; }
            if (o == Sign::Negative) { prev = q; q = s; }
            else { prev = r; r = s; }
            out.emplace_back(q, r);
        }
        error = "segment walk did not terminate";
        return false;
    }

    bool ahead(int a, int b, int q) const {
        return (P[q][0] - P[a][0]) * (P[b][0] - P[a][0]) + (P[q][1] - P[a][1]) * (P[b][1] - P[a][1]) > 0.0;
    }
};

/// Samples of one edge, from v0 to v1. For a side edge, `t` is the edge
/// parameter; for an arc edge, `x` is the curve point (both carriers'
/// parameters). End samples are the Brep vertices.
struct EdgeSamples {
    std::vector<int> gid;
    std::vector<double> t;
    std::vector<Pt4> x;
};

/// Second-derivative and first-derivative sizes of a carrier over its square.
struct CarrierShape {
    double Lu = 0.0, Lv = 0.0;        ///< largest |Su|, |Sv|
    double Lu_mean = 0.0, Lv_mean = 0.0;
    double Kuu = 0.0, Kuv = 0.0, Kvv = 0.0;
};

inline CarrierShape carrier_shape(const geom::BezierSurface& S) {
    CarrierShape c;
    const int n = 9;
    int cnt = 0;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) {
            const geom::Jet3 J = S.jet(static_cast<double>(i) / (n - 1), static_cast<double>(j) / (n - 1));
            const double su = norm(J[static_cast<std::size_t>(geom::jet_index(1, 0))]);
            const double sv = norm(J[static_cast<std::size_t>(geom::jet_index(0, 1))]);
            c.Lu = std::max(c.Lu, su);
            c.Lv = std::max(c.Lv, sv);
            c.Lu_mean += su;
            c.Lv_mean += sv;
            ++cnt;
            c.Kuu = std::max(c.Kuu, norm(J[static_cast<std::size_t>(geom::jet_index(2, 0))]));
            c.Kuv = std::max(c.Kuv, norm(J[static_cast<std::size_t>(geom::jet_index(1, 1))]));
            c.Kvv = std::max(c.Kvv, norm(J[static_cast<std::size_t>(geom::jet_index(0, 2))]));
        }
    c.Lu_mean /= cnt;
    c.Lv_mean /= cnt;
    return c;
}

/// The same for a procedural carrier, from its first partials and their
/// central differences, sampled inside each of its smooth pieces in u.
inline CarrierShape carrier_shape(const geom::Carrier& C) {
    CarrierShape c;
    const std::vector<double> bs = C.breaks_u();
    const int n = 9;
    int cnt = 0;
    for (std::size_t k = 0; k + 1 < bs.size(); ++k)
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j) {
                const double h = 1e-4 * (bs[k + 1] - bs[k]);
                const double u = bs[k] + (bs[k + 1] - bs[k]) * (i + 0.5) / n;
                const double v = std::clamp(static_cast<double>(j) / (n - 1), 1e-4, 1 - 1e-4);
                Vec3d P, Pu, Pv, A, Au, Av, B, Bu, Bv, D, Du, Dv, E, Eu, Ev;
                C.eval(u, v, P, Pu, Pv);
                C.eval(u + h, v, A, Au, Av);
                C.eval(u - h, v, B, Bu, Bv);
                C.eval(u, v + 1e-4, D, Du, Dv);
                C.eval(u, v - 1e-4, E, Eu, Ev);
                const double su = norm(Pu), sv = norm(Pv);
                c.Lu = std::max(c.Lu, su);
                c.Lv = std::max(c.Lv, sv);
                c.Lu_mean += su;
                c.Lv_mean += sv;
                ++cnt;
                c.Kuu = std::max(c.Kuu, norm(Au - Bu) / (2 * h));
                c.Kuv = std::max(c.Kuv, norm(Av - Bv) / (2 * h));
                c.Kvv = std::max(c.Kvv, norm(Dv - Ev) / 2e-4);
            }
    c.Lu_mean /= cnt;
    c.Lv_mean /= cnt;
    return c;
}

/// Parameter steps for a face's interior grid: each term of the second-order
/// deviation of a triangle is kept below a share of `chord`. Infinity for a
/// direction that needs no interior lines.
inline void grid_steps(const CarrierShape& c, const TessOptions& opt, double& du, double& dv) {
    const double inf = std::numeric_limits<double>::infinity();
    // Curvature measured relative to the patch size: below this, flat.
    const double flat = 1e-12 * std::max(1.0, std::max(c.Lu, c.Lv));
    du = c.Kuu > flat ? std::sqrt(2.0 * opt.chord / c.Kuu) : inf;
    dv = c.Kvv > flat ? std::sqrt(2.0 * opt.chord / c.Kvv) : inf;
    if (c.Kuv > flat) {
        const double cap = 2.0 * opt.chord / c.Kuv; // du * dv <= cap
        if (!std::isfinite(du) && !std::isfinite(dv)) du = dv = std::sqrt(cap);
        else if (!std::isfinite(du)) du = cap / dv;
        else if (!std::isfinite(dv)) dv = cap / du;
        else if (du * dv > cap) {
            const double f = std::sqrt(cap / (du * dv));
            du *= f;
            dv *= f;
        }
    }
    if (opt.max_edge > 0.0) {
        if (c.Lu > 0.0) du = std::min(du, opt.max_edge / c.Lu);
        if (c.Lv > 0.0) dv = std::min(dv, opt.max_edge / c.Lv);
    }
}

/// Adaptive sampling of a curve f on [s0, s1] (end points p0, p1 known).
/// Interior samples are appended in order.
template <class F>
bool subdivide(const F& f, double s0, double s1, const Vec3d& p0, const Vec3d& p1, int depth, int min_depth,
               double max_len, const TessOptions& opt, std::vector<double>& out_s, std::vector<Vec3d>& out_p) {
    const double sm = 0.5 * (s0 + s1);
    Vec3d pm;
    if (!f(sm, pm)) return false;
    bool split = depth < min_depth;
    if (!split && depth < opt.max_depth) {
        Vec3d q1, q3;
        if (!f(0.75 * s0 + 0.25 * s1, q1) || !f(0.25 * s0 + 0.75 * s1, q3)) return false;
        const double dev = std::max(seg_dist(pm, p0, p1), std::max(seg_dist(q1, p0, p1), seg_dist(q3, p0, p1)));
        split = dev > opt.chord || norm(p1 - p0) > max_len;
    }
    if (!split) return true;
    if (!subdivide(f, s0, sm, p0, pm, depth + 1, min_depth, max_len, opt, out_s, out_p)) return false;
    out_s.push_back(sm);
    out_p.push_back(pm);
    return subdivide(f, sm, s1, pm, p1, depth + 1, min_depth, max_len, opt, out_s, out_p);
}

} // namespace tess_detail

inline TessResult tessellate(const Brep& b, const TessOptions& opt = {}) {
    using namespace tess_detail;
    const auto t_start = std::chrono::steady_clock::now();
    TessResult res;
    Mesh& M = res.mesh;
    auto fail = [&](std::string msg) {
        res.ok = false;
        res.error = std::move(msg);
        return res;
    };

    std::vector<geom::BezierSurface> surf;
    std::vector<CarrierShape> shape;
    surf.reserve(b.carriers.size());
    for (const PolyPatchD& c : b.carriers) {
        surf.emplace_back(c);
        shape.push_back(carrier_shape(surf.back()));
    }
    // Procedural carriers follow the polynomial ones in `shape`.
    for (const auto& pc : b.procs) shape.push_back(carrier_shape(*pc));
    auto shape_of = [&](const BFace& f) -> const CarrierShape& {
        return shape[static_cast<std::size_t>(f.proc >= 0 ? static_cast<int>(b.carriers.size()) + f.proc : f.carrier)];
    };
    // Point of a face's surface at parameters (u, v).
    auto face_point = [&](const BFace& f, double u, double v) {
        if (f.proc < 0) return surf[static_cast<std::size_t>(f.carrier)].point(u, v);
        Vec3d P, Pu, Pv;
        b.procs[static_cast<std::size_t>(f.proc)]->eval(u, v, P, Pu, Pv);
        return P;
    };

    // Spacing wanted along each face, in space: edges are sampled no coarser.
    const double inf = std::numeric_limits<double>::infinity();
    std::vector<double> face_len(b.faces.size(), inf);
    for (std::size_t f = 0; f < b.faces.size(); ++f) {
        const CarrierShape& c = shape_of(b.faces[f]);
        double du, dv;
        grid_steps(c, opt, du, dv);
        face_len[f] = std::min(du * c.Lu, dv * c.Lv);
    }
    std::vector<double> edge_len(b.edges.size(), opt.max_edge > 0.0 ? opt.max_edge : inf);
    for (std::size_t f = 0; f < b.faces.size(); ++f)
        for (const auto& loop : b.faces[f].loops)
            for (const BCoedge& ce : loop)
                edge_len[static_cast<std::size_t>(ce.edge)] =
                    std::min(edge_len[static_cast<std::size_t>(ce.edge)], face_len[f]);

    // Points: Brep vertices first.
    for (const BVertex& v : b.vertices) M.points.push_back(v.point);

    // Edges, each sampled once.
    std::vector<EdgeSamples> samples(b.edges.size());
    M.edges.resize(b.edges.size());
    for (std::size_t ei = 0; ei < b.edges.size(); ++ei) {
        const BEdge& e = b.edges[ei];
        EdgeSamples& S = samples[ei];
        const Vec3d p0 = b.vertices[static_cast<std::size_t>(e.v0)].point;
        const Vec3d p1 = b.vertices[static_cast<std::size_t>(e.v1)].point;
        std::vector<double> s;
        std::vector<Vec3d> pts;
        S.gid.push_back(e.v0);
        if (e.proc) {
            // Procedural edge: sampled in tau, piece by piece.
            auto f = [&](double tau, Vec3d& p) {
                p = e.proc->point(tau);
                return std::isfinite(p.x) && std::isfinite(p.y) && std::isfinite(p.z);
            };
            const std::vector<double> bt = e.proc->breaks();
            const int n = static_cast<int>(bt.size()) - 1;
            S.t.push_back(0.0);
            Vec3d prev = p0;
            for (int k = 0; k < n; ++k) {
                const double ta = bt[static_cast<std::size_t>(k)], tb = bt[static_cast<std::size_t>(k) + 1];
                Vec3d end = p1;
                if (k + 1 < n && !f(tb, end)) return fail("could not evaluate a procedural edge");
                s.clear();
                pts.clear();
                if (!subdivide(f, ta, tb, prev, end, 0, n == 1 ? 1 : 0, edge_len[ei], opt, s, pts))
                    return fail("could not evaluate a procedural edge");
                if (k + 1 < n) {
                    s.push_back(tb);
                    pts.push_back(end);
                }
                for (std::size_t q = 0; q < s.size(); ++q) {
                    S.t.push_back(s[q]);
                    S.gid.push_back(static_cast<int>(M.points.size()));
                    M.points.push_back(pts[q]);
                }
                prev = end;
            }
            S.t.push_back(1.0);
        } else if (!e.arc) {
            const BEdge::SideRef& ref = e.on[0];
            const geom::BezierSurface& C = surf[static_cast<std::size_t>(ref.carrier)];
            auto f = [&](double t, Vec3d& p) {
                double u, v;
                side_uv(ref.side, ref.reversed ? 1.0 - t : t, u, v);
                p = C.point(u, v);
                return std::isfinite(p.x) && std::isfinite(p.y) && std::isfinite(p.z);
            };
            const double t0 = num::mid(e.t0), t1 = num::mid(e.t1);
            if (!subdivide(f, t0, t1, p0, p1, 0, e.v0 == e.v1 ? 2 : 0, edge_len[ei], opt, s, pts))
                return fail("could not evaluate a side edge");
            S.t.push_back(t0);
            for (std::size_t k = 0; k < s.size(); ++k) {
                S.t.push_back(s[k]);
                S.gid.push_back(static_cast<int>(M.points.size()));
                M.points.push_back(pts[k]);
            }
            S.t.push_back(t1);
        } else {
            const ArcGeom& arc = *e.geom;
            const SsiSystem& sys = *arc.base->sys;
            const auto& segs = arc.base->curve->segs;
            const geom::BezierSurface& C0 = surf[static_cast<std::size_t>(e.carrier[0])];
            const geom::BezierSurface& C1 = surf[static_cast<std::size_t>(e.carrier[1])];
            S.x.push_back(Pt4{});
            const int ns = static_cast<int>(arc.spans.size());
            const bool closed = e.v0 == e.v1;
            Vec3d prev = p0;
            // Is a curve point in the given slice range of a segment (and in
            // its certified box)?
            auto in_range = [](const Pt4& x, const TubeSegment& sg, double lo, double hi) {
                for (int i = 0; i < 4; ++i)
                    if (!num::contains(sg.X[i], x[i])) return false;
                double c = 0.0;
                for (int i = 0; i < 4; ++i) c += sg.m[i] * (x[i] - sg.a[i]);
                return c >= lo && c <= hi;
            };
            // A span travelled against its segment's direction is the closing
            // span of a closed curve whose last segment overlaps the first:
            // the curve goes back over the overlap and forward again. Such a
            // piece adds nothing to a winding number, but a polyline must not
            // double back, so the samples that overshoot are dropped instead.
            int cut_seg = -1;
            double cut_lo = 0.0, cut_hi = 0.0;
            for (int j = 0; j < ns; ++j) {
                const Span& sp = arc.spans[static_cast<std::size_t>(j)];
                const TubeSegment& seg = segs[static_cast<std::size_t>(sp.seg)];
                const bool last = j + 1 == ns;
                const bool backward =
                    ns > 1 && (sp.cb - sp.ca) * (seg_end_slice(seg) - seg_start_slice(seg)) < 0.0;
                if (backward) {
                    const double lo = std::min(sp.ca, sp.cb) - kSlicePad, hi = std::max(sp.ca, sp.cb) + kSlicePad;
                    if (j == 0) {
                        // The arc starts inside the overlap: drop what the next
                        // span samples before the start.
                        cut_seg = sp.seg;
                        cut_lo = lo;
                        cut_hi = hi;
                        continue;
                    }
                    while (S.x.size() > 1 && in_range(S.x.back(), seg, lo, hi)) {
                        S.x.pop_back();
                        S.gid.pop_back();
                        M.points.pop_back();
                    }
                    if (last) break;
                    Pt4 xj;
                    for (int i = 0; i < 4; ++i) xj[i] = seg.a[i] + (sp.cb / seg.len) * (seg.b[i] - seg.a[i]);
                    const Pt4 axis = xj;
                    if (!solve::newton_slice(sys, seg.m, seg.a, sp.cb, xj)) {
                        bool in = true; // as in point_at below: the axis point, in the certified box
                        for (int i = 0; i < 4; ++i) in = in && num::contains(seg.X[static_cast<std::size_t>(i)], axis[static_cast<std::size_t>(i)]);
                        if (!in) return fail("could not evaluate an intersection edge");
                        xj = axis;
                        ++res.stats.axis_points;
                    }
                    prev = (C0.point(xj[0], xj[1]) + C1.point(xj[2], xj[3])) * 0.5;
                    S.x.push_back(xj);
                    S.gid.push_back(static_cast<int>(M.points.size()));
                    M.points.push_back(prev);
                    continue;
                }
                auto point_at = [&](double c, Pt4& x) {
                    for (int i = 0; i < 4; ++i) x[i] = seg.a[i] + (c / seg.len) * (seg.b[i] - seg.a[i]);
                    const Pt4 axis = x;
                    bool ok = solve::newton_slice(sys, seg.m, seg.a, c, x);
                    for (int i = 0; i < 4 && ok; ++i) ok = std::isfinite(x[i]);
                    if (ok) return true;
                    // Newton can fail where the surfaces meet at a very small angle
                    // (a tiny loop near a tangency): the sliced system is then too
                    // badly conditioned for its tolerance. The point on the
                    // segment's axis lies in its certified box, so it is within the
                    // box's size of the curve; the mesh is measured, not certified.
                    for (int i = 0; i < 4; ++i)
                        if (!num::contains(seg.X[static_cast<std::size_t>(i)], axis[static_cast<std::size_t>(i)])) return false;
                    x = axis;
                    ++res.stats.axis_points;
                    return true;
                };
                auto f = [&](double c, Vec3d& p) {
                    Pt4 x;
                    if (!point_at(c, x)) return false;
                    p = (C0.point(x[0], x[1]) + C1.point(x[2], x[3])) * 0.5;
                    return true;
                };
                Vec3d end = p1;
                Pt4 xend{};
                if (!last) {
                    if (!point_at(sp.cb, xend)) return fail("could not evaluate an intersection edge");
                    end = (C0.point(xend[0], xend[1]) + C1.point(xend[2], xend[3])) * 0.5;
                }
                // At least one interior point on an arc (two pieces at least),
                // and four pieces on a closed one.
                const int min_depth = (closed && ns < 4) ? 2 : (ns == 1 ? 1 : 0);
                s.clear();
                pts.clear();
                if (!subdivide(f, sp.ca, sp.cb, prev, end, 0, min_depth, edge_len[ei], opt, s, pts))
                    return fail("could not evaluate an intersection edge");
                for (std::size_t k = 0; k < s.size(); ++k) {
                    Pt4 x;
                    point_at(s[k], x);
                    if (cut_seg >= 0 && in_range(x, segs[static_cast<std::size_t>(cut_seg)], cut_lo, cut_hi)) continue;
                    S.x.push_back(x);
                    S.gid.push_back(static_cast<int>(M.points.size()));
                    M.points.push_back(pts[k]);
                }
                if (!last &&
                    !(cut_seg >= 0 && in_range(xend, segs[static_cast<std::size_t>(cut_seg)], cut_lo, cut_hi))) {
                    S.x.push_back(xend);
                    S.gid.push_back(static_cast<int>(M.points.size()));
                    M.points.push_back(end);
                }
                cut_seg = -1;
                prev = end;
            }
            S.x.push_back(Pt4{});
        }
        S.gid.push_back(e.v1);
        res.stats.edge_points += static_cast<int>(S.gid.size()) - 2;
        M.edges[ei] = S.gid;
    }

    // Faces.
    for (std::size_t fi = 0; fi < b.faces.size(); ++fi) {
        const BFace& face = b.faces[fi];
        const CarrierShape& cs = shape_of(face);
        struct Node {
            double u, v;
            int gid;
        };
        std::vector<std::vector<Node>> loops;
        for (const auto& loop : face.loops) {
            std::vector<Node> L;
            for (const BCoedge& ce : loop) {
                const BEdge& e = b.edges[static_cast<std::size_t>(ce.edge)];
                const EdgeSamples& S = samples[static_cast<std::size_t>(ce.edge)];
                const std::size_t n = S.gid.size();
                std::vector<Node> pts(n);
                if (face.proc >= 0) {
                    // A procedural face: its edges are sides of its square.
                    if (!e.proc) return fail("a procedural face is bounded by an edge that is not procedural");
                    for (std::size_t k = 0; k < n; ++k) {
                        double u, v;
                        if (!e.proc->on_proc(face.proc, S.t[k], u, v))
                            return fail("a procedural edge does not lie on the procedural face it bounds");
                        pts[k] = {u, v, S.gid[k]};
                    }
                    if (ce.reversed) std::reverse(pts.begin(), pts.end());
                    if (!L.empty() && L.back().gid != pts.front().gid) return fail("a face loop is not connected");
                    for (std::size_t k = L.empty() ? 0 : 1; k < pts.size(); ++k) L.push_back(pts[k]);
                    continue;
                }
                const Chart* c0 = b.vertices[static_cast<std::size_t>(e.v0)].chart(face.carrier);
                const Chart* c1 = b.vertices[static_cast<std::size_t>(e.v1)].chart(face.carrier);
                if (!c0 || !c1) return fail("a vertex has no parameters on a face it bounds");
                const double u0 = num::mid(c0->u), v0 = num::mid(c0->v), u1 = num::mid(c1->u), v1 = num::mid(c1->v);
                pts.front() = {u0, v0, S.gid.front()};
                pts.back() = {u1, v1, S.gid.back()};
                if (e.proc) {
                    for (std::size_t k = 1; k + 1 < n; ++k) {
                        double u, v, du, dv;
                        if (!e.proc->on_carrier(face.carrier, S.t[k], u, v, du, dv))
                            return fail("a procedural edge does not lie on the carrier of a face it bounds");
                        pts[k] = {u, v, S.gid[k]};
                    }
                } else if (!e.arc) {
                    const double t0 = S.t.front(), t1 = S.t.back();
                    for (std::size_t k = 1; k + 1 < n; ++k) {
                        const double w = (S.t[k] - t0) / (t1 - t0);
                        pts[k] = {u0 + w * (u1 - u0), v0 + w * (v1 - v0), S.gid[k]};
                    }
                } else {
                    int co = -1;
                    bool mapped = false;
                    Affine map;
                    if (!arc_on(b, e, face.carrier, co, mapped, map))
                        return fail("an arc edge does not lie on the carrier of a face it bounds");
                    for (std::size_t k = 1; k + 1 < n; ++k) {
                        const Pt4& x = S.x[k];
                        double u = x[static_cast<std::size_t>(co)], v = x[static_cast<std::size_t>(co) + 1];
                        if (mapped) {
                            Interval mu, mv;
                            map.apply(Interval(u), Interval(v), mu, mv);
                            u = num::mid(mu);
                            v = num::mid(mv);
                        }
                        pts[k] = {u, v, S.gid[k]};
                    }
                }
                if (ce.reversed) std::reverse(pts.begin(), pts.end());
                if (!L.empty() && L.back().gid != pts.front().gid)
                    return fail("a face loop is not connected");
                for (std::size_t k = L.empty() ? 0 : 1; k < pts.size(); ++k) L.push_back(pts[k]);
            }
            if (L.size() < 2 || L.front().gid != L.back().gid) return fail("a face loop is not closed");
            L.pop_back();
            loops.push_back(std::move(L));
        }

        // Exact scaling of v by a power of two, to make the parameter plane
        // closer to isotropic in space.
        double sv = 1.0;
        if (cs.Lu_mean > 0.0 && cs.Lv_mean > 0.0)
            sv = std::ldexp(1.0, static_cast<int>(std::lround(std::log2(cs.Lv_mean / cs.Lu_mean))));
        auto Q = [&](double u, double v) { return P2{u, v * sv}; };

        P2 lo{inf, inf}, hi{-inf, -inf};
        for (const auto& L : loops)
            for (const Node& nd : L) {
                const P2 q = Q(nd.u, nd.v);
                lo = {std::min(lo[0], q[0]), std::min(lo[1], q[1])};
                hi = {std::max(hi[0], q[0]), std::max(hi[1], q[1])};
            }

        Cdt cdt;
        cdt.init(lo, hi);
        std::vector<std::vector<int>> loop_v;  // CDT vertex of every loop node
        std::vector<int> gid_of(3, -1);         // per CDT vertex
        std::vector<std::array<double, 2>> uv_of(3, {0.0, 0.0});
        for (const auto& L : loops) {
            std::vector<int> lv;
            for (const Node& nd : L) {
                const int v = cdt.insert(Q(nd.u, nd.v));
                if (v < 0) return fail("face " + std::to_string(fi) + ": " + cdt.error);
                if (v >= static_cast<int>(gid_of.size())) {
                    gid_of.push_back(nd.gid);
                    uv_of.push_back({nd.u, nd.v});
                } else if (gid_of[static_cast<std::size_t>(v)] != nd.gid) {
                    return fail("face " + std::to_string(fi) + ": two boundary points coincide in its parameters");
                }
                lv.push_back(v);
            }
            loop_v.push_back(std::move(lv));
        }

        // Interior grid.
        double du, dv;
        grid_steps(cs, opt, du, dv);
        const bool curved = std::isfinite(du) || std::isfinite(dv);
        const double wu = hi[0] - lo[0], wv = hi[1] - lo[1];
        if (curved && wu > 0.0 && wv > 0.0) {
            // One step in both directions of the scaled plane: the Delaunay
            // triangles are then close to equilateral there, so no triangle
            // side runs far along the more curved direction.
            const double h = std::min(du, dv * sv);
            const int nu = std::min(std::max(static_cast<int>(std::ceil(wu / h)), 2), opt.max_grid);
            const int nv = std::min(std::max(static_cast<int>(std::ceil(wv / h)), 2), opt.max_grid);
            const double hu = wu / nu, hv = wv / nv; // steps in the scaled plane
            const double margin = 0.4 * std::min(hu, hv);
            std::vector<char> ok(static_cast<std::size_t>((nu - 1) * (nv - 1)), 1);
            auto at = [&](int i, int j) { return P2{lo[0] + i * hu, lo[1] + j * hv}; };
            // Reject grid points near the boundary.
            for (const auto& L : loops)
                for (std::size_t k = 0; k < L.size(); ++k) {
                    const P2 a = Q(L[k].u, L[k].v), c = Q(L[(k + 1) % L.size()].u, L[(k + 1) % L.size()].v);
                    const int i0 = std::max(1, static_cast<int>(std::floor((std::min(a[0], c[0]) - margin - lo[0]) / hu)));
                    const int i1 = std::min(nu - 1, static_cast<int>(std::ceil((std::max(a[0], c[0]) + margin - lo[0]) / hu)));
                    const int j0 = std::max(1, static_cast<int>(std::floor((std::min(a[1], c[1]) - margin - lo[1]) / hv)));
                    const int j1 = std::min(nv - 1, static_cast<int>(std::ceil((std::max(a[1], c[1]) + margin - lo[1]) / hv)));
                    for (int i = i0; i <= i1; ++i)
                        for (int j = j0; j <= j1; ++j)
                            if (seg_dist2(at(i, j), a, c) < margin)
                                ok[static_cast<std::size_t>((i - 1) * (nv - 1) + (j - 1))] = 0;
                }
            // Keep points inside the loops (scanline parity per grid row).
            std::vector<double> xs;
            for (int j = 1; j < nv; ++j) {
                const double y = lo[1] + j * hv;
                xs.clear();
                for (const auto& L : loops)
                    for (std::size_t k = 0; k < L.size(); ++k) {
                        const P2 a = Q(L[k].u, L[k].v), c = Q(L[(k + 1) % L.size()].u, L[(k + 1) % L.size()].v);
                        if ((a[1] <= y) != (c[1] <= y)) xs.push_back(a[0] + (y - a[1]) * (c[0] - a[0]) / (c[1] - a[1]));
                    }
                std::sort(xs.begin(), xs.end());
                for (int i = 1; i < nu; ++i) {
                    char& keep = ok[static_cast<std::size_t>((i - 1) * (nv - 1) + (j - 1))];
                    if (!keep) continue;
                    const double x = lo[0] + i * hu;
                    const auto crossings = std::lower_bound(xs.begin(), xs.end(), x) - xs.begin();
                    if (crossings % 2 == 0) { keep = 0; continue; }
                    const P2 q = at(i, j);
                    const int v = cdt.insert(q);
                    if (v < 0) return fail("face " + std::to_string(fi) + ": " + cdt.error);
                    if (v < static_cast<int>(gid_of.size())) continue; // already there
                    const double u = q[0], vv = q[1] / sv;
                    gid_of.push_back(static_cast<int>(M.points.size()));
                    uv_of.push_back({u, vv});
                    M.points.push_back(face_point(face, u, vv));
                    ++res.stats.interior_points;
                }
            }
        }

        for (const auto& lv : loop_v)
            for (std::size_t k = 0; k < lv.size(); ++k)
                if (!cdt.constrain(lv[k], lv[(k + 1) % lv.size()]))
                    return fail("face " + std::to_string(fi) + ": " + cdt.error);
        cdt.delaunay_pass();
        std::vector<char> in;
        if (!cdt.inside(in)) return fail("face " + std::to_string(fi) + ": " + cdt.error);
        res.stats.flips += cdt.flips;
        res.stats.constraint_flips += cdt.constraint_flips;
        for (std::size_t t = 0; t < cdt.T.size(); ++t) {
            if (!in[t]) continue;
            std::array<int, 3> v = cdt.T[t].v;
            if (face.flipped) std::swap(v[1], v[2]);
            std::array<int, 3> g;
            std::array<double, 6> uv;
            for (int k = 0; k < 3; ++k) {
                g[static_cast<std::size_t>(k)] = gid_of[static_cast<std::size_t>(v[static_cast<std::size_t>(k)])];
                uv[static_cast<std::size_t>(2 * k)] = uv_of[static_cast<std::size_t>(v[static_cast<std::size_t>(k)])][0];
                uv[static_cast<std::size_t>(2 * k + 1)] = uv_of[static_cast<std::size_t>(v[static_cast<std::size_t>(k)])][1];
            }
            M.tris.push_back(g);
            M.tri_face.push_back(static_cast<int>(fi));
            M.tri_uv.push_back(uv);
        }
    }
    // Two faces that bound a thin lens (two surfaces meeting along all of
    // its boundary, such as two tops that cross along a curve and share the
    // lines it ends on) have the same outline, and each may choose the same
    // chord between two boundary points. The mesh would then touch itself
    // along that segment although the surfaces are apart there. Such a chord
    // (used by more than one face, and not a piece of a Brep edge) is split
    // in each face at that face's own parameter midpoint.
    {
        std::set<std::pair<int, int>> on_edges;
        for (const auto& pts : M.edges)
            for (std::size_t k = 0; k + 1 < pts.size(); ++k)
                on_edges.insert({std::min(pts[k], pts[k + 1]), std::max(pts[k], pts[k + 1])});
        for (int round = 0; round < 4; ++round) {
            std::map<std::pair<int, int>, std::vector<int>> use;
            for (int t = 0; t < static_cast<int>(M.tris.size()); ++t)
                for (int k = 0; k < 3; ++k) {
                    const int a = M.tris[static_cast<std::size_t>(t)][static_cast<std::size_t>(k)];
                    const int c = M.tris[static_cast<std::size_t>(t)][static_cast<std::size_t>((k + 1) % 3)];
                    use[{std::min(a, c), std::max(a, c)}].push_back(t);
                }
            int splits = 0;
            std::vector<char> touched(M.tris.size(), 0); // a triangle is split at most once a round
            for (const auto& [key, ts] : use) {
                if (ts.size() <= 2 || on_edges.count(key)) continue;
                std::map<int, std::vector<int>> by_face;
                bool fresh = true;
                for (int t : ts) {
                    by_face[M.tri_face[static_cast<std::size_t>(t)]].push_back(t);
                    fresh = fresh && !touched[static_cast<std::size_t>(t)];
                }
                bool pairs = fresh && by_face.size() > 1;
                for (const auto& [f, list] : by_face) pairs = pairs && list.size() == 2;
                if (!pairs) continue;
                for (int t : ts) touched[static_cast<std::size_t>(t)] = 1;
                for (const auto& [f, list] : by_face) {
                    const BFace& lf = b.faces[static_cast<std::size_t>(f)];
                    int mid = -1;
                    double mu = 0.0, mv = 0.0;
                    for (int t : list) {
                        const std::array<int, 3> g = M.tris[static_cast<std::size_t>(t)];
                        const std::array<double, 6> uv = M.tri_uv[static_cast<std::size_t>(t)];
                        int i = 0;
                        for (; i < 3; ++i) {
                            const int a = g[static_cast<std::size_t>(i)], c = g[static_cast<std::size_t>((i + 1) % 3)];
                            if (std::min(a, c) == key.first && std::max(a, c) == key.second) break;
                        }
                        const int j = (i + 1) % 3;
                        if (mid < 0) {
                            mu = 0.5 * (uv[static_cast<std::size_t>(2 * i)] + uv[static_cast<std::size_t>(2 * j)]);
                            mv = 0.5 * (uv[static_cast<std::size_t>(2 * i + 1)] + uv[static_cast<std::size_t>(2 * j + 1)]);
                            mid = static_cast<int>(M.points.size());
                            M.points.push_back(face_point(lf, mu, mv));
                        }
                        // (p_i, p_j, p_l) becomes (p_i, m, p_l) and (m, p_j, p_l).
                        std::array<int, 3> g1 = g, g2 = g;
                        std::array<double, 6> uv1 = uv, uv2 = uv;
                        g1[static_cast<std::size_t>(j)] = mid;
                        uv1[static_cast<std::size_t>(2 * j)] = mu;
                        uv1[static_cast<std::size_t>(2 * j + 1)] = mv;
                        g2[static_cast<std::size_t>(i)] = mid;
                        uv2[static_cast<std::size_t>(2 * i)] = mu;
                        uv2[static_cast<std::size_t>(2 * i + 1)] = mv;
                        M.tris[static_cast<std::size_t>(t)] = g1;
                        M.tri_uv[static_cast<std::size_t>(t)] = uv1;
                        M.tris.push_back(g2);
                        M.tri_uv.push_back(uv2);
                        M.tri_face.push_back(f);
                    }
                }
                ++splits;
                ++res.stats.lens_splits;
            }
            if (splits == 0) break;
        }
    }
    res.stats.triangles = static_cast<int>(M.tris.size());
    res.stats.ms =
        std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t_start).count();
    res.ok = true;
    return res;
}

/// Integer and measured checks of a mesh.
struct MeshAudit {
    int points = 0, triangles = 0, edges = 0;
    int open_edges = 0;        ///< used by one triangle
    int nonmanifold_edges = 0; ///< used by four or more triangles (parts touching along an edge)
    int misoriented_edges = 0; ///< used more often in one direction than in the other
    int unused_points = 0;
    int degenerate = 0;        ///< triangles with a repeated point or zero area
    int components = 0;        ///< connected by shared edges
    int euler = 0;             ///< V - E + F over used points
    double volume = 0.0;       ///< divergence theorem
    double area = 0.0;
    /// Closed: every edge used as often in each direction (once each, or
    /// more where parts of the solid touch along an edge).
    bool watertight() const { return open_edges == 0 && misoriented_edges == 0; }
    bool manifold() const { return nonmanifold_edges == 0; }
};

inline MeshAudit audit(const Mesh& m) {
    MeshAudit a;
    a.points = static_cast<int>(m.points.size());
    a.triangles = static_cast<int>(m.tris.size());
    struct Use {
        int fwd = 0, bwd = 0;
    };
    std::unordered_map<std::uint64_t, Use> uses;
    uses.reserve(m.tris.size() * 2);
    std::vector<char> used(m.points.size(), 0);
    std::vector<int> parent(m.tris.size());
    for (std::size_t i = 0; i < parent.size(); ++i) parent[i] = static_cast<int>(i);
    auto find = [&](int x) {
        while (parent[static_cast<std::size_t>(x)] != x) {
            parent[static_cast<std::size_t>(x)] = parent[static_cast<std::size_t>(parent[static_cast<std::size_t>(x)])];
            x = parent[static_cast<std::size_t>(x)];
        }
        return x;
    };
    std::unordered_map<std::uint64_t, int> first_tri;
    for (std::size_t t = 0; t < m.tris.size(); ++t) {
        const auto& T = m.tris[t];
        const Vec3d& p0 = m.points[static_cast<std::size_t>(T[0])];
        const Vec3d& p1 = m.points[static_cast<std::size_t>(T[1])];
        const Vec3d& p2 = m.points[static_cast<std::size_t>(T[2])];
        const Vec3d n = cross(p1 - p0, p2 - p0);
        a.area += 0.5 * norm(n);
        a.volume += dot(p0, cross(p1, p2)) / 6.0;
        if (T[0] == T[1] || T[1] == T[2] || T[0] == T[2] || norm(n) == 0.0) ++a.degenerate;
        for (int k = 0; k < 3; ++k) {
            used[static_cast<std::size_t>(T[static_cast<std::size_t>(k)])] = 1;
            const int x = T[static_cast<std::size_t>(k)], y = T[static_cast<std::size_t>((k + 1) % 3)];
            const std::uint64_t key = (static_cast<std::uint64_t>(std::min(x, y)) << 32) |
                                      static_cast<std::uint32_t>(std::max(x, y));
            Use& u = uses[key];
            (x < y ? u.fwd : u.bwd) += 1;
            auto it = first_tri.find(key);
            if (it == first_tri.end()) first_tri.emplace(key, static_cast<int>(t));
            else parent[static_cast<std::size_t>(find(static_cast<int>(t)))] = find(it->second);
        }
    }
    a.edges = static_cast<int>(uses.size());
    for (const auto& [key, u] : uses) {
        const int total = u.fwd + u.bwd;
        if (total == 1) ++a.open_edges;
        else {
            if (total > 2) ++a.nonmanifold_edges;
            if (u.fwd != u.bwd) ++a.misoriented_edges;
        }
    }
    int used_points = 0;
    for (char c : used) used_points += c ? 1 : 0;
    a.unused_points = a.points - used_points;
    for (std::size_t t = 0; t < m.tris.size(); ++t)
        if (find(static_cast<int>(t)) == static_cast<int>(t)) ++a.components;
    a.euler = used_points - a.edges + a.triangles;
    return a;
}

/// Distance from q to a carrier: closest point by Gauss-Newton, started at
/// the parameters (u, v) and kept in the parameter square.
inline double surface_distance(const geom::BezierSurface& S, const Vec3d& q, double u, double v) {
    const std::size_t iu = static_cast<std::size_t>(geom::jet_index(1, 0));
    const std::size_t iv = static_cast<std::size_t>(geom::jet_index(0, 1));
    for (int it = 0; it < 12; ++it) {
        const geom::Jet3 J = S.jet(u, v);
        const Vec3d r = q - J[0];
        const double a = dot(J[iu], J[iu]), c = dot(J[iu], J[iv]), d = dot(J[iv], J[iv]);
        const double g1 = dot(J[iu], r), g2 = dot(J[iv], r);
        const double det = a * d - c * c;
        if (!(det > 0.0)) break;
        const double du = (d * g1 - c * g2) / det, dv = (a * g2 - c * g1) / det;
        const double nu = std::min(1.0, std::max(0.0, u + du)), nv = std::min(1.0, std::max(0.0, v + dv));
        const bool done = std::fabs(nu - u) + std::fabs(nv - v) < 1e-15;
        u = nu;
        v = nv;
        if (done) break;
    }
    return norm(q - S.point(u, v));
}

/// The same for a procedural carrier (its partials are differences, so the
/// steps stop at 1e-12 rather than at rounding).
inline double surface_distance(const geom::Carrier& S, const Vec3d& q, double u, double v) {
    Vec3d P, Pu, Pv;
    for (int it = 0; it < 12; ++it) {
        S.eval(u, v, P, Pu, Pv);
        const Vec3d r = q - P;
        const double a = dot(Pu, Pu), c = dot(Pu, Pv), d = dot(Pv, Pv);
        const double g1 = dot(Pu, r), g2 = dot(Pv, r);
        const double det = a * d - c * c;
        if (!(det > 0.0)) break;
        const double du = (d * g1 - c * g2) / det, dv = (a * g2 - c * g1) / det;
        const double nu = std::min(1.0, std::max(0.0, u + du)), nv = std::min(1.0, std::max(0.0, v + dv));
        const bool done = std::fabs(nu - u) + std::fabs(nv - v) < 1e-12;
        u = nu;
        v = nv;
        if (done) break;
    }
    S.eval(u, v, P, Pu, Pv);
    return norm(q - P);
}

/// Largest distance from the mesh to its surfaces, over the centroid and the
/// side midpoints of every triangle (closest points on the triangle's
/// carrier). A measurement of how well the mesh follows the surfaces, not a
/// bound.
inline double mesh_deviation(const Mesh& m, const Brep& b) {
    std::vector<geom::BezierSurface> surf;
    surf.reserve(b.carriers.size());
    for (const PolyPatchD& c : b.carriers) surf.emplace_back(c);
    double worst = 0.0;
    for (std::size_t t = 0; t < m.tris.size(); ++t) {
        const BFace& face = b.faces[static_cast<std::size_t>(m.tri_face[t])];
        const auto& T = m.tris[t];
        const auto& uv = m.tri_uv[t];
        const Vec3d p[3] = {m.points[static_cast<std::size_t>(T[0])], m.points[static_cast<std::size_t>(T[1])],
                            m.points[static_cast<std::size_t>(T[2])]};
        static constexpr double W[4][3] = {
            {1.0 / 3, 1.0 / 3, 1.0 / 3}, {0.5, 0.5, 0.0}, {0.0, 0.5, 0.5}, {0.5, 0.0, 0.5}};
        for (const auto& w : W) {
            const double u = w[0] * uv[0] + w[1] * uv[2] + w[2] * uv[4];
            const double v = w[0] * uv[1] + w[1] * uv[3] + w[2] * uv[5];
            const Vec3d q = p[0] * w[0] + p[1] * w[1] + p[2] * w[2];
            worst = std::max(worst, face.proc >= 0
                                        ? surface_distance(*b.procs[static_cast<std::size_t>(face.proc)], q, u, v)
                                        : surface_distance(surf[static_cast<std::size_t>(face.carrier)], q, u, v));
        }
    }
    return worst;
}

/// Wavefront OBJ: one group per Brep face, and the Brep edges as polylines.
inline void write_obj(std::ostream& out, const Mesh& m, bool with_edges = true) {
    out.precision(12);
    for (const Vec3d& p : m.points) out << "v " << p.x << " " << p.y << " " << p.z << "\n";
    int face = -1;
    for (std::size_t t = 0; t < m.tris.size(); ++t) {
        if (m.tri_face[t] != face) {
            face = m.tri_face[t];
            out << "g face_" << face << "\n";
        }
        out << "f " << m.tris[t][0] + 1 << " " << m.tris[t][1] + 1 << " " << m.tris[t][2] + 1 << "\n";
    }
    if (!with_edges) return;
    for (std::size_t e = 0; e < m.edges.size(); ++e) {
        out << "g edge_" << e << "\nl";
        for (int g : m.edges[e]) out << " " << g + 1;
        out << "\n";
    }
}

} // namespace ak::op
