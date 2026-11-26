// akernel — solids bounded by whole rational Bezier patches.
//
// This is the input model of the first Boolean operation: a closed shell in
// which every face is one untrimmed patch and the outward normal of every
// face is Su x Sv.
//
// Identity by construction (docs/BLUEPRINT.md, thesis T1): two face sides
// are the same edge only if their boundary control points are identical, in
// the same or the opposite order. Then the two boundary curves are the same
// polynomial curve exactly, and the edge has one definition. No distance is
// measured and no tolerance is involved.
#pragma once

#include <array>
#include <numeric>
#include <string>
#include <vector>

#include <akernel/geom/bernstein.hpp>
#include <akernel/geom/bspline.hpp>
#include <akernel/geom/nurbs_surface.hpp>
#include <akernel/solve/csi.hpp>

namespace ak::op {

using geom::PolyPatchD;
using geom::Vec3d;
using geom::Vec4d;

// Sides of the parameter square, numbered counter-clockwise:
//   0: v = 0 (parameter u)   1: u = 1 (parameter v)
//   2: v = 1 (parameter u)   3: u = 0 (parameter v)
// Corners: 0 = (0,0), 1 = (1,0), 2 = (1,1), 3 = (0,1).

/// Control points of one side as a curve in the side's own parameter.
inline std::vector<Vec4d> side_control(const PolyPatchD& P, int side) {
    std::vector<Vec4d> c;
    if (side == 0) for (int i = 0; i <= P.p; ++i) c.push_back(P.at(i, 0));
    else if (side == 1) for (int j = 0; j <= P.q; ++j) c.push_back(P.at(P.p, j));
    else if (side == 2) for (int i = 0; i <= P.p; ++i) c.push_back(P.at(i, P.q));
    else for (int j = 0; j <= P.q; ++j) c.push_back(P.at(0, j));
    return c;
}

/// Point of the parameter square on a side, for scalar type T (double or Interval).
template <class T>
void side_uv(int side, const T& tau, T& u, T& v) {
    if (side == 0) { u = tau; v = T(0.0); }
    else if (side == 1) { u = T(1.0); v = tau; }
    else if (side == 2) { u = tau; v = T(1.0); }
    else { u = T(0.0); v = tau; }
}

/// Does the counter-clockwise direction of the boundary run with the side parameter?
inline bool side_ccw_forward(int side) { return side < 2; }
/// Corner indices at side parameter 0 and 1.
inline int side_corner(int side, int end) {
    static constexpr int table[4][2] = {{0, 1}, {1, 2}, {3, 2}, {0, 3}};
    return table[side][end];
}

struct PatchSolid {
    struct SideRef {
        int edge = -1;
        bool reversed = false; ///< side parameter = 1 - edge parameter
    };
    struct ModelEdge {
        int face[2] = {-1, -1};
        int side[2] = {-1, -1};
        bool rev[2] = {false, false};
        PolyPatchD curve;      ///< the one definition of the edge
        int v0 = -1, v1 = -1;  ///< vertices at edge parameter 0 and 1
    };

    std::vector<PolyPatchD> faces;
    std::vector<std::array<SideRef, 4>> sides;
    std::vector<ModelEdge> edges;
    std::vector<std::array<int, 4>> corners; ///< vertex index of each face corner
    std::vector<Vec3d> vertex_points;

    /// Derive edges and vertices from the faces. Fails unless every side is
    /// shared by exactly two faces; with `open`, a side of one face only is a
    /// boundary edge (face[1] = -1): the faces of a sheet.
    bool build(std::string& error, bool open = false) {
        const int nF = static_cast<int>(faces.size());
        sides.assign(static_cast<std::size_t>(nF), {});
        edges.clear();
        auto same = [](const Vec4d& a, const Vec4d& b) {
            return a.x == b.x && a.y == b.y && a.z == b.z && a.w == b.w;
        };
        std::vector<std::vector<Vec4d>> edge_cp;
        for (int f = 0; f < nF; ++f)
            for (int s = 0; s < 4; ++s) {
                const auto cp = side_control(faces[static_cast<std::size_t>(f)], s);
                int match = -1;
                bool rev = false;
                for (std::size_t e = 0; e < edges.size() && match < 0; ++e) {
                    const auto& ref = edge_cp[e];
                    if (ref.size() != cp.size()) continue;
                    bool fwd = true, bwd = true;
                    for (std::size_t i = 0; i < cp.size(); ++i) {
                        fwd = fwd && same(cp[i], ref[i]);
                        bwd = bwd && same(cp[i], ref[cp.size() - 1 - i]);
                    }
                    if (fwd || bwd) { match = static_cast<int>(e); rev = !fwd; }
                }
                if (match < 0) {
                    ModelEdge me;
                    me.face[0] = f;
                    me.side[0] = s;
                    me.curve = solve::make_curve(cp);
                    edges.push_back(me);
                    edge_cp.push_back(cp);
                    sides[static_cast<std::size_t>(f)][static_cast<std::size_t>(s)] =
                        {static_cast<int>(edges.size()) - 1, false};
                } else {
                    ModelEdge& me = edges[static_cast<std::size_t>(match)];
                    if (me.face[1] >= 0) { error = "an edge is shared by more than two face sides"; return false; }
                    me.face[1] = f;
                    me.side[1] = s;
                    me.rev[1] = rev;
                    sides[static_cast<std::size_t>(f)][static_cast<std::size_t>(s)] = {match, rev};
                }
            }
        if (!open)
            for (const ModelEdge& me : edges)
                if (me.face[1] < 0) { error = "a face side has no matching side: the shell is open"; return false; }

        // Vertices: corner slots joined through the shared edges.
        std::vector<int> parent(static_cast<std::size_t>(4 * nF));
        std::iota(parent.begin(), parent.end(), 0);
        auto find = [&](int x) {
            while (parent[static_cast<std::size_t>(x)] != x) x = parent[static_cast<std::size_t>(x)];
            return x;
        };
        auto slot = [](int f, int corner) { return 4 * f + corner; };
        for (const ModelEdge& me : edges)
            for (int end = 0; end < 2 && me.face[1] >= 0; ++end) {
                const int a = slot(me.face[0], side_corner(me.side[0], end));
                const int b = slot(me.face[1], side_corner(me.side[1], me.rev[1] ? 1 - end : end));
                parent[static_cast<std::size_t>(find(a))] = find(b);
            }
        std::vector<int> id(static_cast<std::size_t>(4 * nF), -1);
        corners.assign(static_cast<std::size_t>(nF), {});
        vertex_points.clear();
        for (int f = 0; f < nF; ++f)
            for (int c = 0; c < 4; ++c) {
                const int root = find(slot(f, c));
                if (id[static_cast<std::size_t>(root)] < 0) {
                    id[static_cast<std::size_t>(root)] = static_cast<int>(vertex_points.size());
                    const PolyPatchD& P = faces[static_cast<std::size_t>(f)];
                    const Vec4d h = P.at(c == 1 || c == 2 ? P.p : 0, c >= 2 ? P.q : 0);
                    vertex_points.push_back({h.x / h.w, h.y / h.w, h.z / h.w});
                }
                corners[static_cast<std::size_t>(f)][static_cast<std::size_t>(c)] = id[static_cast<std::size_t>(root)];
            }
        for (ModelEdge& me : edges) {
            me.v0 = corners[static_cast<std::size_t>(me.face[0])][static_cast<std::size_t>(side_corner(me.side[0], 0))];
            me.v1 = corners[static_cast<std::size_t>(me.face[0])][static_cast<std::size_t>(side_corner(me.side[0], 1))];
        }
        return true;
    }
};

/// Integrand of the divergence-theorem volume: S . (Su x Sv) at (u, v).
inline double volume_density(const geom::BezierSurface& S, double u, double v) {
    const geom::Jet3 j = S.jet(u, v);
    return dot(j[0], cross(j[1], j[2]));
}

/// Gauss-Legendre rule mapped to [0, 1].
struct UnitRule {
    std::vector<double> x, w;
    explicit UnitRule(int n) {
        geom::gauss_legendre(n, x, w);
        for (std::size_t i = 0; i < x.size(); ++i) {
            x[i] = 0.5 * (x[i] + 1.0);
            w[i] *= 0.5;
        }
    }
};

/// Volume enclosed by the shell: (1/3) * sum over faces of the integral of
/// S . (Su x Sv). Numerical quadrature in doubles, not a certified bound.
inline double volume(const PatchSolid& s) {
    const UnitRule rule(16);
    double total = 0.0;
    for (const PolyPatchD& P : s.faces) {
        const geom::BezierSurface S(P);
        for (std::size_t i = 0; i < rule.x.size(); ++i)
            for (std::size_t j = 0; j < rule.x.size(); ++j)
                total += rule.w[i] * rule.w[j] * volume_density(S, rule.x[i], rule.x[j]);
    }
    return total / 3.0;
}

/// A trivariate Bezier lattice of degree (p, q, r) in homogeneous form.
struct Lattice {
    int p = 1, q = 1, r = 1;
    std::vector<Vec4d> c;
    Lattice(int p_, int q_, int r_)
        : p(p_), q(q_), r(r_), c(static_cast<std::size_t>((p_ + 1) * (q_ + 1) * (r_ + 1))) {}
    Vec4d& at(int i, int j, int k) { return c[static_cast<std::size_t>((i * (q + 1) + j) * (r + 1) + k)]; }
    const Vec4d& at(int i, int j, int k) const {
        return c[static_cast<std::size_t>((i * (q + 1) + j) * (r + 1) + k)];
    }
};

/// The six boundary faces of a lattice, each oriented so that Su x Sv points
/// out when the lattice map preserves orientation. Neighbouring faces take
/// their shared control points from the same lattice nodes, so edges match
/// exactly.
inline PatchSolid make_hex(const Lattice& L, std::string* error = nullptr) {
    PatchSolid s;
    const int p = L.p, q = L.q, r = L.r;
    PolyPatchD x0(r, q), x1(q, r), y0(p, r), y1(r, p), z0(q, p), z1(p, q);
    for (int a = 0; a <= r; ++a) for (int b = 0; b <= q; ++b) x0.at(a, b) = L.at(0, b, a);
    for (int a = 0; a <= q; ++a) for (int b = 0; b <= r; ++b) x1.at(a, b) = L.at(p, a, b);
    for (int a = 0; a <= p; ++a) for (int b = 0; b <= r; ++b) y0.at(a, b) = L.at(a, 0, b);
    for (int a = 0; a <= r; ++a) for (int b = 0; b <= p; ++b) y1.at(a, b) = L.at(b, q, a);
    for (int a = 0; a <= q; ++a) for (int b = 0; b <= p; ++b) z0.at(a, b) = L.at(b, a, 0);
    for (int a = 0; a <= p; ++a) for (int b = 0; b <= q; ++b) z1.at(a, b) = L.at(a, b, r);
    s.faces = {x0, x1, y0, y1, z0, z1};
    std::string err;
    if (!s.build(err)) s.faces.clear();
    if (error) *error = err;
    return s;
}

/// Parallelepiped o + x ex + y ey + z ez; (ex, ey, ez) must be right-handed.
inline PatchSolid make_box(const Vec3d& o, const Vec3d& ex, const Vec3d& ey, const Vec3d& ez) {
    Lattice L(1, 1, 1);
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j)
            for (int k = 0; k < 2; ++k) {
                const Vec3d pt = o + ex * static_cast<double>(i) + ey * static_cast<double>(j) +
                                 ez * static_cast<double>(k);
                L.at(i, j, k) = {pt.x, pt.y, pt.z, 1.0};
            }
    return make_hex(L);
}
inline PatchSolid make_box(const Vec3d& lo, const Vec3d& hi) {
    return make_box(lo, {hi.x - lo.x, 0.0, 0.0}, {0.0, hi.y - lo.y, 0.0}, {0.0, 0.0, hi.z - lo.z});
}

/// Box with footprint [o.x, o.x+sx] x [o.y, o.y+sy] whose top is the graph
/// z = o.z + sz + 16 bump x'(1-x') y'(1-y') in footprint coordinates.
inline PatchSolid make_bump_box(const Vec3d& o, double sx, double sy, double sz, double bump) {
    Lattice L(2, 2, 1);
    for (int i = 0; i <= 2; ++i)
        for (int j = 0; j <= 2; ++j) {
            const double x = o.x + 0.5 * sx * i, y = o.y + 0.5 * sy * j;
            L.at(i, j, 0) = {x, y, o.z, 1.0};
            L.at(i, j, 1) = {x, y, o.z + sz + ((i == 1 && j == 1) ? 4.0 * bump : 0.0), 1.0};
        }
    return make_hex(L);
}

/// The outer faces of a grid of lattices that fill a solid column by column:
/// cell (a, b) of an m x n grid, row-major in a. Neighbouring cells must take
/// their shared face from the same nodes; those inner faces are left out,
/// and the rest are matched into one shell by PatchSolid::build (identical
/// side control points). Faces of one carrier surface split into several
/// patches (a B-spline's Bezier pieces) are separate faces sharing their
/// sides: multi-patch faces.
inline PatchSolid make_hex_grid(const std::vector<Lattice>& cells, int m, int n, std::string* error = nullptr) {
    PatchSolid s;
    for (int a = 0; a < m; ++a)
        for (int b = 0; b < n; ++b) {
            const Lattice& L = cells[static_cast<std::size_t>(a * n + b)];
            const int p = L.p, q = L.q, r = L.r;
            if (a == 0) {
                PolyPatchD x0(r, q);
                for (int i = 0; i <= r; ++i) for (int j = 0; j <= q; ++j) x0.at(i, j) = L.at(0, j, i);
                s.faces.push_back(x0);
            }
            if (a == m - 1) {
                PolyPatchD x1(q, r);
                for (int i = 0; i <= q; ++i) for (int j = 0; j <= r; ++j) x1.at(i, j) = L.at(p, i, j);
                s.faces.push_back(x1);
            }
            if (b == 0) {
                PolyPatchD y0(p, r);
                for (int i = 0; i <= p; ++i) for (int j = 0; j <= r; ++j) y0.at(i, j) = L.at(i, 0, j);
                s.faces.push_back(y0);
            }
            if (b == n - 1) {
                PolyPatchD y1(r, p);
                for (int i = 0; i <= r; ++i) for (int j = 0; j <= p; ++j) y1.at(i, j) = L.at(j, q, i);
                s.faces.push_back(y1);
            }
            PolyPatchD z0(q, p), z1(p, q);
            for (int i = 0; i <= q; ++i) for (int j = 0; j <= p; ++j) z0.at(i, j) = L.at(j, i, 0);
            for (int i = 0; i <= p; ++i) for (int j = 0; j <= q; ++j) z1.at(i, j) = L.at(i, j, r);
            s.faces.push_back(z0);
            s.faces.push_back(z1);
        }
    std::string err;
    if (!s.build(err)) s.faces.clear();
    if (error) *error = err;
    return s;
}

/// Box with footprint [o.x, o.x+sx] x [o.y, o.y+sy] whose top is the graph
/// z = o.z + sz + f(x', y') of a quadratic B-spline f in footprint
/// coordinates, with clamped uniform knots of m and n spans and control
/// heights h (row-major in x', (m + 2) x (n + 2) values). The x and y of the
/// control points are the Greville abscissae, so x and y are affine in each
/// piece's parameters: every top piece is a graph. The top is split into its
/// m x n Bezier pieces, one face each (C1 across the seams), and so are the
/// bottom and the sides below them. With dyadic o, sx, sy, sz, h and spans a
/// power of two, every control point is exact.
inline PatchSolid make_spline_box(const Vec3d& o, double sx, double sy, double sz, int m, int n, const std::vector<double>& h) {
    const int p = 2, nu = m + p, nv = n + p;
    if (static_cast<int>(h.size()) != nu * nv) return {};
    geom::NurbsSurface top(p, p, nu, nv, geom::open_uniform_knots(nu, p), geom::open_uniform_knots(nv, p));
    const auto gu = geom::greville(nu, p, top.U), gv = geom::greville(nv, p, top.V);
    for (int i = 0; i < nu; ++i)
        for (int j = 0; j < nv; ++j)
            top.at(i, j) = {o.x + sx * gu[static_cast<std::size_t>(i)], o.y + sy * gv[static_cast<std::size_t>(j)],
                            o.z + sz + h[static_cast<std::size_t>(i * nv + j)], 1.0};
    const auto pieces = top.to_bezier();
    std::vector<Lattice> cells;
    for (const auto& piece : pieces) {
        Lattice L(p, p, 1);
        for (int i = 0; i <= p; ++i)
            for (int j = 0; j <= p; ++j) {
                const Vec4d& c = piece.patch.at(i, j);
                L.at(i, j, 1) = c;
                L.at(i, j, 0) = {c.x, c.y, o.z, 1.0};
            }
        cells.push_back(L);
    }
    // to_bezier lists pieces row-major in u: (a, b) at a * n + b.
    return make_hex_grid(cells, m, n);
}

/// Quarter of a cylindrical shell: r0 <= sqrt(x^2 + y^2) <= r1, x, y >= 0,
/// z0 <= z <= z1. The circular arcs are exact (rational quadratics).
inline PatchSolid make_annular_sector(double r0, double r1, double z0, double z1) {
    Lattice L(1, 2, 1);
    const double rad[2] = {r0, r1}, z[2] = {z0, z1};
    for (int i = 0; i < 2; ++i)
        for (int k = 0; k < 2; ++k) {
            L.at(i, 0, k) = {rad[i], 0.0, z[k], 1.0};
            L.at(i, 1, k) = {rad[i], rad[i], z[k], 1.0};
            L.at(i, 2, k) = {0.0, 2.0 * rad[i], 2.0 * z[k], 2.0};
        }
    return make_hex(L);
}

} // namespace ak::op
