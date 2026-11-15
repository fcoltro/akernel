// akernel — are two solids the same solid?
//
// Used for the algebraic identities of docs/BLUEPRINT.md §5.4 (for example
// A - (A - B) = A and B, commutativity, associativity), where two results
// computed along different routes must be the same solid with the same
// topology, not just the same volume.
//
// Two correct results can differ in ways that do not change the solid:
//   - a seam: an edge whose two faces lie on the same surface with the same
//     outward side (two pieces of one face, glued back together);
//   - a vertex of degree two where the edge simply continues between the same
//     two faces.
// Both are removed first: faces joined by seams form one face group, and
// edges joined at such vertices form one super-edge. What remains is
// compared:
//   - essential vertices: matched by their certified 3D enclosures, which
//     come from the vertices' certified parameters (interval evaluation).
//     Each vertex must meet exactly one vertex of the other solid, and the
//     matching must be one to one. Disjoint enclosures prove two vertices
//     different; meeting ones are taken as the same point;
//   - surfaces: the same surface means identical control points, or two flat
//     parallelograms in the same plane (exact test, op/planar.hpp), with the
//     outward side compared through the plane's normal;
//   - super-edges: by end vertices and the two face groups' surfaces and
//     sides;
//   - face groups: by surface, side and the super-edges that bound them.
// The solids are the same when these multisets are equal. If two super-edges
// have equal keys (two arcs of one curve between the same two vertices) they
// are interchangeable in the comparison; this is reported.
#pragma once

#include <algorithm>
#include <map>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

#include <akernel/op/brep.hpp>
#include <akernel/op/planar.hpp>
#include <akernel/op/surfaces.hpp>

namespace ak::op {

struct SolidComparison {
    bool same = false;
    std::string reason;            ///< first difference found, if any
    int vertices[2] = {0, 0};      ///< essential vertices
    int super_edges[2] = {0, 0};
    int face_groups[2] = {0, 0};
    int removed_vertices[2] = {0, 0}; ///< vertices of degree two removed
    int seams[2] = {0, 0};            ///< seam edges removed
    int ambiguous = 0;                ///< super-edges sharing a key with another
};

namespace compare_detail {

inline bool identical(const PolyPatchD& a, const PolyPatchD& b) {
    if (a.p != b.p || a.q != b.q || a.c.size() != b.c.size()) return false;
    for (std::size_t i = 0; i < a.c.size(); ++i)
        if (a.c[i].x != b.c[i].x || a.c[i].y != b.c[i].y || a.c[i].z != b.c[i].z || a.c[i].w != b.c[i].w)
            return false;
    return true;
}

inline Vec3d plane_normal(const PolyPatchD& P) {
    const Vec4d &c00 = planar_detail::corner(P, 0, 0), &c10 = planar_detail::corner(P, 1, 0), &c01 = planar_detail::corner(P, 0, 1);
    return cross(Vec3d{c10.x - c00.x, c10.y - c00.y, c10.z - c00.z}, Vec3d{c01.x - c00.x, c01.y - c00.y, c01.z - c00.z});
}

/// One solid with seams and inessential vertices removed.
struct Reduced {
    std::vector<int> group_of_face;
    std::vector<int> group_class, group_side; ///< surface class and outward side (+1 / -1) per group
    std::vector<char> seam;                   ///< per edge
    std::vector<char> essential;              ///< per vertex
    struct Super {
        int v0 = -1, v1 = -1; ///< essential end vertices, -1 for a closed chain
        std::vector<int> g;   ///< face groups along it (two, or more where parts touch)
    };
    std::vector<Super> supers;
    int groups = 0, removed = 0, seams = 0;
};

inline Reduced reduce(const Brep& b, const std::vector<int>& cls, const std::vector<int>& side_of_carrier) {
    Reduced R;
    const std::size_t nf = b.faces.size(), ne = b.edges.size(), nv = b.vertices.size();
    // Faces of every edge.
    std::vector<std::vector<int>> faces_of(ne);
    for (std::size_t f = 0; f < nf; ++f)
        for (const auto& loop : b.faces[f].loops)
            for (const BCoedge& ce : loop) faces_of[static_cast<std::size_t>(ce.edge)].push_back(static_cast<int>(f));
    auto face_class = [&](int f) { return cls[static_cast<std::size_t>(b.faces[static_cast<std::size_t>(f)].carrier)]; };
    auto face_side = [&](int f) {
        const BFace& face = b.faces[static_cast<std::size_t>(f)];
        return (face.flipped ? -1 : 1) * side_of_carrier[static_cast<std::size_t>(face.carrier)];
    };
    // Seams and face groups.
    std::vector<int> parent(nf);
    std::iota(parent.begin(), parent.end(), 0);
    auto find = [&](int x) {
        while (parent[static_cast<std::size_t>(x)] != x)
            x = parent[static_cast<std::size_t>(x)] = parent[static_cast<std::size_t>(parent[static_cast<std::size_t>(x)])];
        return x;
    };
    R.seam.assign(ne, 0);
    for (std::size_t e = 0; e < ne; ++e) {
        const auto& fs = faces_of[e];
        if (fs.size() != 2) continue;
        if (face_class(fs[0]) == face_class(fs[1]) && face_side(fs[0]) == face_side(fs[1])) {
            R.seam[e] = 1;
            ++R.seams;
            parent[static_cast<std::size_t>(find(fs[0]))] = find(fs[1]);
        }
    }
    std::map<int, int> group_id;
    R.group_of_face.resize(nf);
    for (std::size_t f = 0; f < nf; ++f) {
        const int root = find(static_cast<int>(f));
        auto it = group_id.find(root);
        if (it == group_id.end()) {
            it = group_id.emplace(root, R.groups++).first;
            R.group_class.push_back(face_class(static_cast<int>(f)));
            R.group_side.push_back(face_side(static_cast<int>(f)));
        }
        R.group_of_face[f] = it->second;
    }
    auto groups_of = [&](std::size_t e) {
        std::vector<int> g;
        for (int f : faces_of[e]) g.push_back(R.group_of_face[static_cast<std::size_t>(f)]);
        std::sort(g.begin(), g.end());
        return g;
    };
    // Vertices: an edge-end count over the edges that are not seams.
    std::vector<std::vector<int>> inc(nv);
    for (std::size_t e = 0; e < ne; ++e) {
        if (R.seam[e] || faces_of[e].size() < 2) continue;
        inc[static_cast<std::size_t>(b.edges[e].v0)].push_back(static_cast<int>(e));
        inc[static_cast<std::size_t>(b.edges[e].v1)].push_back(static_cast<int>(e));
    }
    R.essential.assign(nv, 0);
    for (std::size_t v = 0; v < nv; ++v) {
        if (inc[v].empty()) continue;
        const bool through = inc[v].size() == 2 &&
                             groups_of(static_cast<std::size_t>(inc[v][0])) == groups_of(static_cast<std::size_t>(inc[v][1]));
        R.essential[v] = through ? 0 : 1;
        if (through) ++R.removed;
    }
    // Super-edges: chains of edges through inessential vertices.
    std::vector<char> visited(ne, 0);
    auto other_end = [&](int e, int v) { return b.edges[static_cast<std::size_t>(e)].v0 == v ? b.edges[static_cast<std::size_t>(e)].v1 : b.edges[static_cast<std::size_t>(e)].v0; };
    auto other_edge = [&](int v, int e) {
        const auto& in = inc[static_cast<std::size_t>(v)];
        return in[0] == e ? in[1] : in[0];
    };
    for (std::size_t e0 = 0; e0 < ne; ++e0) {
        if (R.seam[e0] || faces_of[e0].size() < 2 || visited[e0]) continue;
        visited[e0] = 1;
        Reduced::Super s;
        s.g = groups_of(e0);
        bool cycle = false;
        // Walk forward from v1, then backward from v0.
        int ends[2] = {b.edges[e0].v1, b.edges[e0].v0};
        for (int dir = 0; dir < 2 && !cycle; ++dir) {
            int v = ends[dir], cur = static_cast<int>(e0);
            for (std::size_t guard = 0; guard <= ne; ++guard) {
                if (R.essential[static_cast<std::size_t>(v)]) break;
                const int next = other_edge(v, cur);
                if (next == static_cast<int>(e0) || visited[static_cast<std::size_t>(next)]) { cycle = true; break; }
                visited[static_cast<std::size_t>(next)] = 1;
                v = other_end(next, v);
                cur = next;
            }
            ends[dir] = v;
        }
        if (!cycle) {
            s.v0 = ends[1];
            s.v1 = ends[0];
        }
        R.supers.push_back(s);
    }
    return R;
}

inline bool vertex_box(const Brep& b, const std::vector<geom::PolyPatchI>& ip, int v, Vec3<Interval>& box) {
    bool have = false;
    for (const Chart& c : b.vertices[static_cast<std::size_t>(v)].charts) {
        const geom::Vec4<Interval> h = ip[static_cast<std::size_t>(c.carrier)].eval(c.u, c.v);
        if (!(h.w.lo > 0.0) && !(h.w.hi < 0.0)) continue;
        const Vec3<Interval> p{h.x / h.w, h.y / h.w, h.z / h.w};
        if (!have) { box = p; have = true; continue; }
        // Every chart is certified, so their enclosures all hold the point.
        box = {Interval(std::max(box.x.lo, p.x.lo), std::min(box.x.hi, p.x.hi)),
               Interval(std::max(box.y.lo, p.y.lo), std::min(box.y.hi, p.y.hi)),
               Interval(std::max(box.z.lo, p.z.lo), std::min(box.z.hi, p.z.hi))};
        if (box.x.lo > box.x.hi || box.y.lo > box.y.hi || box.z.lo > box.z.hi) return false;
    }
    return have;
}

inline bool meet(const Vec3<Interval>& a, const Vec3<Interval>& b) {
    return !(a.x.hi < b.x.lo || b.x.hi < a.x.lo || a.y.hi < b.y.lo || b.y.hi < a.y.lo || a.z.hi < b.z.lo ||
             b.z.hi < a.z.lo);
}

} // namespace compare_detail

inline SolidComparison compare_solids(const Brep& A, const Brep& B) {
    using namespace compare_detail;
    SolidComparison r;
    const Brep* S[2] = {&A, &B};
    const bool empty[2] = {A.faces.empty(), B.faces.empty()};
    if (empty[0] || empty[1]) {
        r.same = empty[0] && empty[1];
        if (!r.same) r.reason = "one solid is empty and the other is not";
        return r;
    }
    for (int k = 0; k < 2; ++k)
        if (!audit(*S[k]).closed()) {
            r.reason = k == 0 ? "the first solid does not pass the audit" : "the second solid does not pass the audit";
            return r;
        }

    // Surface classes over the carriers of both solids.
    std::vector<const PolyPatchD*> all;
    for (int k = 0; k < 2; ++k)
        for (const PolyPatchD& c : S[k]->carriers) all.push_back(&c);
    const std::size_t n = all.size(), na = A.carriers.size();
    std::vector<int> parent(n);
    std::iota(parent.begin(), parent.end(), 0);
    auto find = [&](int x) {
        while (parent[static_cast<std::size_t>(x)] != x)
            x = parent[static_cast<std::size_t>(x)] = parent[static_cast<std::size_t>(parent[static_cast<std::size_t>(x)])];
        return x;
    };
    std::vector<char> flat(n);
    for (std::size_t i = 0; i < n; ++i) flat[i] = planar_affine(*all[i]) ? 1 : 0;
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = i + 1; j < n; ++j) {
            if (find(static_cast<int>(i)) == find(static_cast<int>(j))) continue;
            Affine m;
            const bool curved_same = !(flat[i] && flat[j]) && all[i]->p == all[j]->p && all[i]->q == all[j]->q &&
                                     !(all[i]->p == 1 && all[i]->q == 1) && affine_reparam(*all[i], *all[j], m);
            if (identical(*all[i], *all[j]) || (flat[i] && flat[j] && coplanar(*all[i], *all[j])) || curved_same)
                parent[static_cast<std::size_t>(find(static_cast<int>(i)))] = find(static_cast<int>(j));
        }
    // Outward side relative to the class: flat carriers compare their normal
    // with the class representative's (parallel, so the sign is clear).
    std::vector<int> cls(n), side(n, 1);
    for (std::size_t i = 0; i < n; ++i) {
        cls[i] = find(static_cast<int>(i));
        const std::size_t c = static_cast<std::size_t>(cls[i]);
        if (flat[i] && flat[c]) {
            side[i] = dot(plane_normal(*all[i]), plane_normal(*all[c])) > 0.0 ? 1 : -1;
        } else if (c != i && !identical(*all[i], *all[c])) {
            // One curved surface under two parameterisations: the normals agree
            // where the affine map between them keeps orientation.
            Affine m;
            if (same_surface_map(*all[i], *all[c], m)) side[i] = m.det().lo > 0.0 ? 1 : -1;
        }
    }
    std::vector<int> cls_k[2], side_k[2];
    for (int k = 0; k < 2; ++k) {
        const std::size_t off = k == 0 ? 0 : na, cnt = S[k]->carriers.size();
        cls_k[k].assign(cls.begin() + static_cast<std::ptrdiff_t>(off), cls.begin() + static_cast<std::ptrdiff_t>(off + cnt));
        side_k[k].assign(side.begin() + static_cast<std::ptrdiff_t>(off), side.begin() + static_cast<std::ptrdiff_t>(off + cnt));
    }
    const Reduced R[2] = {reduce(A, cls_k[0], side_k[0]), reduce(B, cls_k[1], side_k[1])};

    // Essential vertices, matched by certified enclosures.
    std::vector<Vec3<Interval>> box[2];
    std::vector<int> ess[2];
    for (int k = 0; k < 2; ++k) {
        std::vector<geom::PolyPatchI> ip;
        for (const PolyPatchD& c : S[k]->carriers) ip.push_back(geom::to_interval(c));
        box[k].resize(S[k]->vertices.size());
        for (std::size_t v = 0; v < S[k]->vertices.size(); ++v) {
            if (!R[k].essential[v]) continue;
            if (!vertex_box(*S[k], ip, static_cast<int>(v), box[k][v])) {
                r.reason = "a vertex has no usable certified enclosure";
                return r;
            }
            ess[k].push_back(static_cast<int>(v));
        }
        r.vertices[k] = static_cast<int>(ess[k].size());
        r.super_edges[k] = static_cast<int>(R[k].supers.size());
        r.face_groups[k] = R[k].groups;
        r.removed_vertices[k] = R[k].removed;
        r.seams[k] = R[k].seams;
    }
    if (r.vertices[0] != r.vertices[1]) {
        r.reason = "different numbers of vertices (" + std::to_string(r.vertices[0]) + " and " +
                   std::to_string(r.vertices[1]) + ")";
        return r;
    }
    // id[k][v]: common id of an essential vertex (its index among A's).
    std::vector<int> id[2];
    id[0].assign(A.vertices.size(), -1);
    id[1].assign(B.vertices.size(), -1);
    for (int v : ess[0]) id[0][static_cast<std::size_t>(v)] = v;
    std::vector<char> taken(A.vertices.size(), 0);
    for (int w : ess[1]) {
        int match = -1, count = 0;
        for (int v : ess[0])
            if (meet(box[0][static_cast<std::size_t>(v)], box[1][static_cast<std::size_t>(w)])) { match = v; ++count; }
        if (count != 1) {
            r.reason = count == 0 ? "a vertex of the second solid has no counterpart in the first"
                                  : "a vertex of the second solid meets several vertices of the first";
            return r;
        }
        if (taken[static_cast<std::size_t>(match)]) {
            r.reason = "two vertices of the second solid match the same vertex of the first";
            return r;
        }
        taken[static_cast<std::size_t>(match)] = 1;
        id[1][static_cast<std::size_t>(w)] = match;
    }

    // Keys.
    auto group_key = [&](int k, int g) {
        return std::to_string(R[k].group_class[static_cast<std::size_t>(g)]) +
               (R[k].group_side[static_cast<std::size_t>(g)] > 0 ? "+" : "-");
    };
    std::vector<std::string> super_keys[2], group_keys[2];
    for (int k = 0; k < 2; ++k) {
        std::vector<std::vector<std::string>> bounds(static_cast<std::size_t>(R[k].groups));
        for (const Reduced::Super& s : R[k].supers) {
            std::string ends = "cycle";
            if (s.v0 >= 0) {
                const int a = id[k][static_cast<std::size_t>(s.v0)], c = id[k][static_cast<std::size_t>(s.v1)];
                ends = std::to_string(std::min(a, c)) + "-" + std::to_string(std::max(a, c));
            }
            std::vector<std::string> gk;
            for (int g : s.g) gk.push_back(group_key(k, g));
            std::sort(gk.begin(), gk.end());
            std::string key = ends + " between";
            for (const std::string& x : gk) key += " " + x;
            super_keys[k].push_back(key);
            for (int g : s.g) bounds[static_cast<std::size_t>(g)].push_back(key);
        }
        for (int g = 0; g < R[k].groups; ++g) {
            auto& bk = bounds[static_cast<std::size_t>(g)];
            std::sort(bk.begin(), bk.end());
            std::string key = "surface " + group_key(k, g) + ":";
            for (const std::string& s : bk) key += " [" + s + "]";
            group_keys[k].push_back(key);
        }
        std::sort(super_keys[k].begin(), super_keys[k].end());
        std::sort(group_keys[k].begin(), group_keys[k].end());
    }
    for (std::size_t i = 1; i < super_keys[0].size(); ++i)
        if (super_keys[0][i] == super_keys[0][i - 1]) ++r.ambiguous;
    if (super_keys[0] != super_keys[1]) {
        r.reason = "the edges differ (" + std::to_string(super_keys[0].size()) + " and " +
                   std::to_string(super_keys[1].size()) + " after merging)";
        for (std::size_t i = 0; i < std::min(super_keys[0].size(), super_keys[1].size()); ++i)
            if (super_keys[0][i] != super_keys[1][i]) {
                r.reason += "; first difference: " + super_keys[0][i] + " / " + super_keys[1][i];
                break;
            }
        return r;
    }
    if (group_keys[0] != group_keys[1]) {
        r.reason = "the faces differ (" + std::to_string(group_keys[0].size()) + " and " +
                   std::to_string(group_keys[1].size()) + " after merging)";
        return r;
    }
    r.same = true;
    return r;
}

} // namespace ak::op
