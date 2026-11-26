// akernel — boundary-representation cell complex and its integer audit
// (docs/BLUEPRINT.md §1.5).
//
// The audit is pure integer arithmetic on the incidence structure:
//   - every loop closes;
//   - the chain condition  d1 d2 = 0;
//   - shell closure        d2 . 1 = 0  (each edge used twice, opposite ways);
//   - Euler-Poincare per shell, giving an integer genus;
//   - vertex manifoldness: the face corners at a vertex, joined through the
//     edges at that vertex, form a single fan (no pinch point, and no two
//     shells sharing a vertex).
// Geometry is referenced by index only. Under the kernel's edge model (one
// definition shared by both faces) combinatorial closure implies geometric
// closure, so this audit is the watertightness test.
//
// The same audit classifies a body by its topology (Parasolid's body types):
//   solid   closed manifold shells (watertight());
//   sheet   open manifold oriented shells: every edge used once (boundary)
//           or twice in opposite directions, each shell an orientable
//           surface of integer genus g with b boundary cycles,
//           V - E + F - rings = 2 - 2g - b;
//   wire    edges used by no face, each vertex at the end of at most two
//           edge ends: every component is a simple path or a simple cycle;
//   acorn   vertices only.
// Anything else that is a valid chain (a closed non-manifold boundary,
// faces and wires together) is a general body.
//
// Not checked here: self-intersection (op/validity.hpp checks it for
// Breps of whole-carrier faces).
#pragma once

#include <algorithm>
#include <map>
#include <numeric>
#include <utility>
#include <vector>

namespace ak::topo {

struct Vertex {
    int point = -1;     ///< index of the defining point variety
    double tol = 0.0;   ///< tolerance-ledger entry; 0 unless a merge was recorded
};
struct Edge {
    int v0 = -1, v1 = -1;
    int curve = -1;     ///< index of the defining edge variety
    double tol = 0.0;
};
struct Coedge {
    int edge = -1;
    bool reversed = false;
};
using Loop = std::vector<Coedge>;
struct Face {
    int surface = -1;
    std::vector<Loop> loops; ///< first is the outer loop, the rest are rings
};

struct ShellReport {
    int vertices = 0, edges = 0, faces = 0, rings = 0;
    int euler = 0;        ///< V - E + F - rings
    bool closed = false;  ///< no boundary edges
    int boundary_edges = 0; ///< edges of the shell used once
    int boundary_loops = 0; ///< cycles the boundary edges form (0 for a closed shell)
    /// (2 - euler - boundary_loops) / 2 for a manifold shell whose boundary
    /// is a disjoint union of cycles (closed or open shell); -1 if that is
    /// not an integer >= 0 or the shell is not manifold.
    int genus = -1;
};

/// One connected component of the edges used by no face.
struct WireReport {
    int vertices = 0, edges = 0;
    int max_degree = 0;   ///< most edge ends at one vertex (a ring edge counts twice)
    bool closed = false;  ///< every vertex has degree 2: a simple cycle
    /// A simple path or a simple cycle: max_degree <= 2 and V - E is 1 or 0.
    bool manifold() const { return max_degree <= 2 && (vertices - edges == 1 || (closed && vertices == edges)); }
};

/// Parasolid's body types, decided on the incidence structure alone.
enum class BodyKind { Empty, Acorn, Wire, Sheet, Solid, General, Invalid };

inline const char* to_string(BodyKind k) {
    switch (k) {
    case BodyKind::Empty: return "empty";
    case BodyKind::Acorn: return "acorn";
    case BodyKind::Wire: return "wire";
    case BodyKind::Sheet: return "sheet";
    case BodyKind::Solid: return "solid";
    case BodyKind::General: return "general";
    case BodyKind::Invalid: return "invalid";
    }
    return "?";
}

struct AuditReport {
    bool loops_closed = true;
    bool chain_condition = true; ///< d1 d2 = 0
    int boundary_edges = 0;      ///< used once
    int nonmanifold_edges = 0;   ///< used more than twice
    int misoriented_edges = 0;   ///< used as often in one direction as in the other: no
    int same_way_edges = 0;      ///< used by exactly two faces, both the same way: orientations disagree
    int unused_edges = 0;
    int nonmanifold_vertices = 0; ///< corners at the vertex form more than one fan
    int open_fans = 0;            ///< fans at a vertex that do not close around it
    int euler = 0;                ///< V - E + F - rings over the whole complex (used entities)
    std::vector<ShellReport> shells;
    std::vector<WireReport> wires; ///< components of the edges used by no face
    int isolated_vertices = 0;     ///< vertices on no edge (acorns)
    int ledger_entries = 0;      ///< entities with a recorded tolerance
    double max_tolerance = 0.0;

    /// A closed, oriented boundary: every loop closes, every edge is used an
    /// even number of times, as often in each direction (the 2-chain is a
    /// cycle), and the faces at every vertex form closed fans. Edges used by
    /// four or more faces and vertices with several fans are allowed: the
    /// boundary of a solid whose parts touch along an edge or at a point.
    /// That the faces around such an edge alternate between material and
    /// empty space is geometry, checked by op::check_edge_fans.
    bool closed() const {
        if (!loops_closed || !chain_condition) return false;
        if (boundary_edges || misoriented_edges || unused_edges || open_fans) return false;
        return !shells.empty();
    }
    /// No edge with more than two faces and no vertex with more than one fan.
    bool manifold() const { return nonmanifold_edges == 0 && nonmanifold_vertices == 0; }
    /// A valid closed manifold solid boundary: watertight, oriented, integer genus.
    bool watertight() const {
        if (!closed() || !manifold()) return false;
        for (const auto& s : shells)
            if (!s.closed || s.genus < 0) return false;
        return true;
    }
    /// A valid manifold sheet: oriented, every shell open with boundary
    /// cycles and an integer genus, no wire edges and no acorns.
    bool sheet() const {
        if (!loops_closed || !chain_condition || misoriented_edges || unused_edges || isolated_vertices) return false;
        if (!manifold() || shells.empty()) return false;
        for (const auto& s : shells)
            if (s.closed || s.boundary_loops < 1 || s.genus < 0) return false;
        return true;
    }
    /// A valid manifold wire: no faces, at least one edge, every component a
    /// simple path or cycle, no acorns.
    bool wire() const {
        if (!shells.empty() || wires.empty() || isolated_vertices) return false;
        for (const auto& w : wires)
            if (!w.manifold()) return false;
        return true;
    }
    /// The body type. A body that is not a valid chain (an open loop, a
    /// broken chain condition) or whose faces disagree on orientation across
    /// an edge of two faces is Invalid; a valid one that is not of a special
    /// type is General (an edge of three faces, say, cannot balance its uses).
    BodyKind kind() const {
        if (!loops_closed || !chain_condition || same_way_edges) return BodyKind::Invalid;
        if (shells.empty()) {
            if (wires.empty()) return isolated_vertices ? BodyKind::Acorn : BodyKind::Empty;
            return wire() ? BodyKind::Wire : BodyKind::General;
        }
        if (!unused_edges && !isolated_vertices && watertight()) return BodyKind::Solid;
        if (sheet()) return BodyKind::Sheet;
        return BodyKind::General;
    }
};

class Complex {
public:
    std::vector<Vertex> vertices;
    std::vector<Edge> edges;
    std::vector<Face> faces;

    int add_vertex(int point = -1) {
        vertices.push_back({point, 0.0});
        return static_cast<int>(vertices.size()) - 1;
    }
    int add_edge(int v0, int v1, int curve = -1) {
        edges.push_back({v0, v1, curve, 0.0});
        return static_cast<int>(edges.size()) - 1;
    }
    int add_face(int surface, std::vector<Loop> loops) {
        faces.push_back({surface, std::move(loops)});
        return static_cast<int>(faces.size()) - 1;
    }

    /// Coedge running from vertex a to vertex b, creating the edge if needed.
    Coedge coedge_between(int a, int b) {
        for (std::size_t e = 0; e < edges.size(); ++e) {
            if (edges[e].v0 == a && edges[e].v1 == b) return {static_cast<int>(e), false};
            if (edges[e].v0 == b && edges[e].v1 == a) return {static_cast<int>(e), true};
        }
        return {add_edge(a, b), false};
    }
    /// Face whose loops are given as vertex cycles (outer loop first).
    int add_polygon_face(int surface, const std::vector<std::vector<int>>& cycles) {
        std::vector<Loop> loops;
        for (const auto& cyc : cycles) {
            Loop loop;
            for (std::size_t i = 0; i < cyc.size(); ++i)
                loop.push_back(coedge_between(cyc[i], cyc[(i + 1) % cyc.size()]));
            loops.push_back(std::move(loop));
        }
        return add_face(surface, std::move(loops));
    }

    int start_vertex(const Coedge& c) const { return c.reversed ? edges[c.edge].v1 : edges[c.edge].v0; }
    int end_vertex(const Coedge& c) const { return c.reversed ? edges[c.edge].v0 : edges[c.edge].v1; }

    AuditReport audit() const {
        AuditReport r;
        const int nE = static_cast<int>(edges.size()), nF = static_cast<int>(faces.size());

        std::vector<int> uses(nE, 0), orient_sum(nE, 0);
        // Union-find over faces to find shells.
        std::vector<int> parent(nF);
        std::iota(parent.begin(), parent.end(), 0);
        auto find = [&](int x) {
            while (parent[x] != x) x = parent[x] = parent[parent[x]];
            return x;
        };
        std::vector<int> first_face_of_edge(nE, -1);

        for (int f = 0; f < nF; ++f) {
            std::map<int, int> boundary_of_boundary; // vertex -> coefficient of d1 d2 (face f)
            for (const Loop& loop : faces[f].loops) {
                for (std::size_t i = 0; i < loop.size(); ++i) {
                    const Coedge& c = loop[i];
                    const Coedge& next = loop[(i + 1) % loop.size()];
                    if (end_vertex(c) != start_vertex(next)) r.loops_closed = false;
                    ++uses[c.edge];
                    const int sign = c.reversed ? -1 : 1;
                    orient_sum[c.edge] += sign;
                    // d1(edge) = v1 - v0, taken with the coedge sign.
                    boundary_of_boundary[edges[c.edge].v1] += sign;
                    boundary_of_boundary[edges[c.edge].v0] -= sign;
                    if (first_face_of_edge[c.edge] < 0) first_face_of_edge[c.edge] = f;
                    else parent[find(f)] = find(first_face_of_edge[c.edge]);
                }
            }
            for (const auto& kv : boundary_of_boundary)
                if (kv.second != 0) r.chain_condition = false;
        }

        for (int e = 0; e < nE; ++e) {
            if (uses[e] == 0) ++r.unused_edges;
            else if (uses[e] == 1) ++r.boundary_edges;
            else {
                if (uses[e] > 2) ++r.nonmanifold_edges;
                if (orient_sum[e] != 0) ++r.misoriented_edges;
                if (uses[e] == 2 && orient_sum[e] != 0) ++r.same_way_edges;
            }
            if (edges[e].tol > 0.0) { ++r.ledger_entries; r.max_tolerance = std::max(r.max_tolerance, edges[e].tol); }
        }
        for (const Vertex& v : vertices)
            if (v.tol > 0.0) { ++r.ledger_entries; r.max_tolerance = std::max(r.max_tolerance, v.tol); }

        // Vertex manifoldness. A corner is a pair of consecutive coedges of a
        // loop meeting at a vertex; two corners at the same vertex are in the
        // same fan when they use the same edge there. One fan per vertex.
        std::vector<int> pinched;
        {
            std::vector<int> cparent;
            std::vector<int> corner_vertex;
            std::map<std::pair<int, int>, int> first_corner; // (vertex, edge) -> corner
            auto cfind = [&](int x) {
                while (cparent[static_cast<std::size_t>(x)] != x)
                    x = cparent[static_cast<std::size_t>(x)] = cparent[static_cast<std::size_t>(cparent[static_cast<std::size_t>(x)])];
                return x;
            };
            auto join = [&](int corner, int vertex, int edge) {
                const auto key = std::make_pair(vertex, edge);
                const auto it = first_corner.find(key);
                if (it == first_corner.end()) first_corner.emplace(key, corner);
                else cparent[static_cast<std::size_t>(cfind(corner))] = cfind(it->second);
            };
            for (const Face& face : faces)
                for (const Loop& loop : face.loops)
                    for (std::size_t i = 0; i < loop.size(); ++i) {
                        const Coedge& in = loop[i];
                        const Coedge& out = loop[(i + 1) % loop.size()];
                        const int v = end_vertex(in);
                        const int corner = static_cast<int>(cparent.size());
                        cparent.push_back(corner);
                        corner_vertex.push_back(v);
                        join(corner, v, in.edge);
                        join(corner, v, out.edge);
                    }
            std::map<int, std::map<int, bool>> fans; // vertex -> fan roots
            for (std::size_t c = 0; c < cparent.size(); ++c)
                fans[corner_vertex[c]][cfind(static_cast<int>(c))] = true;
            for (const auto& kv : fans)
                if (kv.second.size() > 1) {
                    ++r.nonmanifold_vertices;
                    pinched.push_back(kv.first);
                }
            // A fan closes around its vertex when each edge end in it is used an
            // even number of times.
            std::map<std::pair<int, int>, int> ends;
            for (const Face& face : faces)
                for (const Loop& loop : face.loops)
                    for (std::size_t i = 0; i < loop.size(); ++i) {
                        const int v = end_vertex(loop[i]);
                        ++ends[{v, loop[i].edge}];
                        ++ends[{v, loop[(i + 1) % loop.size()].edge}];
                    }
            std::map<int, bool> open_at;
            for (const auto& kv : ends)
                if (kv.second % 2 != 0) open_at[kv.first.first] = true;
            r.open_fans = static_cast<int>(open_at.size());
        }

        // Per-shell Euler-Poincare.
        std::map<int, int> shell_of_root;
        for (int f = 0; f < nF; ++f) {
            const int root = find(f);
            if (!shell_of_root.count(root)) {
                shell_of_root[root] = static_cast<int>(r.shells.size());
                r.shells.push_back({});
                r.shells.back().closed = true;
            }
        }
        std::vector<std::map<int, bool>> shell_vertices(r.shells.size()), shell_edges(r.shells.size());
        std::vector<char> shell_nonmanifold(r.shells.size(), 0);
        std::map<int, bool> all_vertices, all_edges;
        int all_rings = 0;
        for (int f = 0; f < nF; ++f) {
            const std::size_t s = static_cast<std::size_t>(shell_of_root[find(f)]);
            ++r.shells[s].faces;
            r.shells[s].rings += static_cast<int>(faces[f].loops.size()) - 1;
            all_rings += static_cast<int>(faces[f].loops.size()) - 1;
            for (const Loop& loop : faces[f].loops)
                for (const Coedge& c : loop) {
                    shell_edges[s][c.edge] = true;
                    shell_vertices[s][edges[c.edge].v0] = true;
                    shell_vertices[s][edges[c.edge].v1] = true;
                    all_edges[c.edge] = true;
                    all_vertices[edges[c.edge].v0] = true;
                    all_vertices[edges[c.edge].v1] = true;
                    if (uses[c.edge] < 2 || uses[c.edge] % 2 != 0 || orient_sum[c.edge] != 0) r.shells[s].closed = false;
                    if (uses[c.edge] > 2) shell_nonmanifold[s] = 1;
                }
        }
        r.euler = static_cast<int>(all_vertices.size()) - static_cast<int>(all_edges.size()) + nF - all_rings;
        // A vertex pinched inside one shell makes that shell's genus undefined
        // too (two shells touching at a point are each fine).
        for (int v : pinched) {
            std::map<int, int> shells_at;
            for (int f = 0; f < nF; ++f)
                for (const Loop& loop : faces[f].loops)
                    for (const Coedge& c : loop)
                        if (end_vertex(c) == v) ++shells_at[shell_of_root[find(f)]];
            // All fans in one shell: that shell is pinched at v.
            if (shells_at.size() == 1) shell_nonmanifold[static_cast<std::size_t>(shells_at.begin()->first)] = 1;
        }
        // Boundary cycles per shell: components of the boundary edges. At a
        // vertex with one fan the boundary edges number 0 or 2 (an open fan
        // has two free ends), so on a manifold shell they form cycles.
        std::vector<int> vparent(vertices.size());
        std::iota(vparent.begin(), vparent.end(), 0);
        auto vfind = [&](int x) {
            while (vparent[static_cast<std::size_t>(x)] != x)
                x = vparent[static_cast<std::size_t>(x)] = vparent[static_cast<std::size_t>(vparent[static_cast<std::size_t>(x)])];
            return x;
        };
        auto valid_vertex = [&](int v) { return v >= 0 && v < static_cast<int>(vertices.size()); };
        for (std::size_t s = 0; s < r.shells.size(); ++s)
            for (const auto& kv : shell_edges[s])
                if (uses[kv.first] == 1) {
                    ++r.shells[s].boundary_edges;
                    const Edge& e = edges[static_cast<std::size_t>(kv.first)];
                    if (valid_vertex(e.v0) && valid_vertex(e.v1)) vparent[static_cast<std::size_t>(vfind(e.v0))] = vfind(e.v1);
                }
        for (std::size_t s = 0; s < r.shells.size(); ++s) {
            std::map<int, bool> roots;
            for (const auto& kv : shell_edges[s])
                if (uses[kv.first] == 1 && valid_vertex(edges[static_cast<std::size_t>(kv.first)].v0))
                    roots[vfind(edges[static_cast<std::size_t>(kv.first)].v0)] = true;
            r.shells[s].boundary_loops = static_cast<int>(roots.size());
        }
        for (std::size_t s = 0; s < r.shells.size(); ++s) {
            ShellReport& sh = r.shells[s];
            sh.vertices = static_cast<int>(shell_vertices[s].size());
            sh.edges = static_cast<int>(shell_edges[s].size());
            sh.euler = sh.vertices - sh.edges + sh.faces - sh.rings;
            const int chi = sh.euler + sh.boundary_loops; // 2 - 2g
            // An open shell needs manifold edges used once or twice in
            // opposite directions; a closed shell was checked above.
            bool ok = !shell_nonmanifold[s];
            if (!sh.closed)
                for (const auto& kv : shell_edges[s])
                    if (uses[kv.first] > 2 || (uses[kv.first] == 2 && orient_sum[kv.first] != 0)) ok = false;
            if (ok && chi <= 2 && (2 - chi) % 2 == 0) sh.genus = (2 - chi) / 2;
        }

        // Wires: edges used by no face, as components over their vertices.
        {
            std::vector<int> wparent(vertices.size());
            std::iota(wparent.begin(), wparent.end(), 0);
            auto wfind = [&](int x) {
                while (wparent[static_cast<std::size_t>(x)] != x)
                    x = wparent[static_cast<std::size_t>(x)] = wparent[static_cast<std::size_t>(wparent[static_cast<std::size_t>(x)])];
                return x;
            };
            std::vector<int> degree(vertices.size(), 0);
            std::vector<char> on_edge(vertices.size(), 0);
            for (int e = 0; e < nE; ++e) {
                const Edge& E = edges[static_cast<std::size_t>(e)];
                if (!valid_vertex(E.v0) || !valid_vertex(E.v1)) continue;
                on_edge[static_cast<std::size_t>(E.v0)] = on_edge[static_cast<std::size_t>(E.v1)] = 1;
                if (uses[e] != 0) continue;
                ++degree[static_cast<std::size_t>(E.v0)];
                ++degree[static_cast<std::size_t>(E.v1)];
                wparent[static_cast<std::size_t>(wfind(E.v0))] = wfind(E.v1);
            }
            std::map<int, int> wire_of_root;
            for (int e = 0; e < nE; ++e) {
                const Edge& E = edges[static_cast<std::size_t>(e)];
                if (uses[e] != 0 || !valid_vertex(E.v0) || !valid_vertex(E.v1)) continue;
                const int root = wfind(E.v0);
                if (!wire_of_root.count(root)) {
                    wire_of_root[root] = static_cast<int>(r.wires.size());
                    r.wires.push_back({});
                }
                ++r.wires[static_cast<std::size_t>(wire_of_root[root])].edges;
            }
            std::vector<char> all_two(r.wires.size(), 1);
            for (std::size_t v = 0; v < vertices.size(); ++v) {
                if (!on_edge[v]) { ++r.isolated_vertices; continue; }
                if (degree[v] == 0) continue;
                const auto it = wire_of_root.find(wfind(static_cast<int>(v)));
                if (it == wire_of_root.end()) continue;
                WireReport& w = r.wires[static_cast<std::size_t>(it->second)];
                ++w.vertices;
                w.max_degree = std::max(w.max_degree, degree[v]);
                if (degree[v] != 2) all_two[static_cast<std::size_t>(it->second)] = 0;
            }
            for (std::size_t w = 0; w < r.wires.size(); ++w) r.wires[w].closed = all_two[w] != 0;
        }
        return r;
    }
};

} // namespace ak::topo
