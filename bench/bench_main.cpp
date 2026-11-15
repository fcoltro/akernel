// akernel benchmark: prints the measurements recorded in docs/RESULTS.md.
// Deterministic: fixed seeds, no threads. The final line is a hash of every
// numeric result, used to compare runs, build types and machines.
#include <akernel/diag/rate.hpp>
#include <akernel/geom/continuity.hpp>
#include <akernel/geom/primitives.hpp>
#include <akernel/num/predicates.hpp>
#include <akernel/op/boolean.hpp>
#include <akernel/op/compare.hpp>
#include <akernel/op/tessellate.hpp>
#include <akernel/solve/ssi.hpp>
#include <akernel/var/fill.hpp>

#include <chrono>
#include <cstdio>
#include <cstring>
#include <functional>
#include <map>
#include <string>

#include "check.hpp"

using namespace ak;
using geom::PolyPatchD;
using geom::Vec3d;

namespace {

std::uint64_t g_hash = 0xcbf29ce484222325ull;
bool g_profile = false; ///< print the surface-pair profile in the chain sections
bool g_plane_lines = true; ///< BooleanOptions::plane_lines in the chain and scaling sections
bool g_exact_points = true; ///< BooleanOptions::exact_points, in every section
/// Default options, with exact points as the command line asks.
op::BooleanOptions with_exact() {
    op::BooleanOptions o;
    o.exact_points = g_exact_points;
    return o;
}
bool g_dump_faces = false; ///< print the faces of a result that could not be tessellated (curved chains with trace)
void hash_into(std::uint64_t& h, double x) {
    std::uint64_t bits;
    std::memcpy(&bits, &x, sizeof bits);
    for (int i = 0; i < 8; ++i) {
        h ^= (bits >> (8 * i)) & 0xff;
        h *= 0x100000001b3ull;
    }
}
void hash_double(double x) { hash_into(g_hash, x); }

/// Tessellation of every committed chain result (section 9), per section.
/// Kept apart from the chain hashes so those stay comparable with earlier runs.
constexpr double kTessChord = 1e-3;
struct TessTally {
    std::string name;
    int meshes = 0, failed = 0, watertight = 0, euler_ok = 0, components_ok = 0, nonmanifold = 0;
    long long triangles = 0;
    int max_triangles = 0;
    double worst_dev = 0.0, worst_vol = 0.0, total_ms = 0.0, worst_ms = 0.0;
    std::vector<std::string> failures;
};
std::vector<TessTally> g_tess;
std::uint64_t g_tess_hash = 0xcbf29ce484222325ull;

TessTally& tally_named(const std::string& name) {
    for (TessTally& t : g_tess)
        if (t.name == name) return t;
    g_tess.push_back({});
    g_tess.back().name = name;
    return g_tess.back();
}

void tally_tessellation(TessTally& t, const op::Brep& b, const std::string& where) {
    op::TessOptions opt;
    opt.chord = kTessChord;
    const op::TessResult r = op::tessellate(b, opt);
    ++t.meshes;
    if (!r.ok) {
        ++t.failed;
        t.failures.push_back(where + ": " + r.error);
        hash_into(g_tess_hash, -1.0);
        if (g_dump_faces) {
            // The loops of every face: edges, their kind and carriers, and the
            // parameters of their end vertices on the face's carrier.
            for (std::size_t f = 0; f < b.faces.size(); ++f) {
                const op::BFace& face = b.faces[f];
                std::printf("  face %d on carrier %d (%dx%d)%s\n", static_cast<int>(f), face.carrier,
                            b.carriers[static_cast<std::size_t>(face.carrier)].p, b.carriers[static_cast<std::size_t>(face.carrier)].q,
                            face.flipped ? " flipped" : "");
                for (const auto& loop : face.loops) {
                    std::printf("    loop:\n");
                    for (const op::BCoedge& ce : loop) {
                        const op::BEdge& e = b.edges[static_cast<std::size_t>(ce.edge)];
                        auto at = [&](int v) {
                            const op::Chart* c = b.vertices[static_cast<std::size_t>(v)].chart(face.carrier);
                            char buf[96];
                            if (c) std::snprintf(buf, sizeof buf, "v%d (%.6f, %.6f)", v, num::mid(c->u), num::mid(c->v));
                            else std::snprintf(buf, sizeof buf, "v%d (no chart)", v);
                            return std::string(buf);
                        };
                        if (e.arc)
                            std::printf("      edge %d arc on %d,%d%s: %s -> %s\n", ce.edge, e.carrier[0], e.carrier[1],
                                        ce.reversed ? " reversed" : "", at(e.v0).c_str(), at(e.v1).c_str());
                        else
                            std::printf("      edge %d side %d of %d / side %d of %d [%g, %g]%s: %s -> %s\n", ce.edge, e.on[0].side,
                                        e.on[0].carrier, e.on[1].side, e.on[1].carrier, num::mid(e.t0), num::mid(e.t1),
                                        ce.reversed ? " reversed" : "", at(e.v0).c_str(), at(e.v1).c_str());
                    }
                }
            }
        }
        return;
    }
    const op::MeshAudit a = op::audit(r.mesh);
    const topo::AuditReport br = op::audit(b);
    const int euler = br.euler; // of the whole complex, which the mesh must reproduce
    if (!a.manifold()) ++t.nonmanifold;
    if (a.watertight() && a.degenerate == 0) ++t.watertight;
    else t.failures.push_back(where + ": mesh not closed (open " + std::to_string(a.open_edges) + ", misoriented " +
                              std::to_string(a.misoriented_edges) + ", degenerate " + std::to_string(a.degenerate) + ")");
    if (a.euler == euler) ++t.euler_ok;
    else t.failures.push_back(where + ": mesh Euler characteristic " + std::to_string(a.euler) + ", Brep " +
                              std::to_string(euler) + " (unused points " + std::to_string(a.unused_points) + ")");
    if (a.components == static_cast<int>(br.shells.size())) ++t.components_ok;
    t.triangles += a.triangles;
    t.max_triangles = std::max(t.max_triangles, a.triangles);
    t.worst_dev = std::max(t.worst_dev, op::mesh_deviation(r.mesh, b));
    t.worst_vol = std::max(t.worst_vol, std::fabs(a.volume - op::volume(b)) / std::max(a.area, 1e-300));
    t.total_ms += r.stats.ms;
    t.worst_ms = std::max(t.worst_ms, r.stats.ms);
    hash_into(g_tess_hash, static_cast<double>(a.triangles));
    hash_into(g_tess_hash, a.volume);
}

void print_tessellation() {
    std::printf("## 9. Tessellation of every committed result\n\n");
    std::printf("Every result of section 5 (random boxes against the curved solids) and of sections 6 to 8 that "
                "passed the audit and is not empty, tessellated with a chord "
                "tolerance of %.0e. A mesh is counted as closed when every mesh edge is used as often in each "
                "direction (once each, or more where parts of a solid touch along an edge) and no triangle is "
                "degenerate; its Euler characteristic and number of components are "
                "compared with the Brep's integer audit. Deviation is the largest distance from triangle centroids "
                "and side midpoints to the surface; volume error is |mesh volume - Brep volume| / mesh area.\n\n",
                kTessChord);
    std::printf("| Section | Meshes | Not tessellated | Closed | Of those, not manifold | Euler characteristic matches | "
                "Components match | Mean triangles | Most triangles | Worst deviation / chord | Worst volume error / chord | "
                "Mean time, ms | Worst time, ms |\n|---|---|---|---|---|---|---|---|---|---|---|---|---|\n");
    for (const TessTally& t : g_tess) {
        const int ok = std::max(1, t.meshes - t.failed);
        std::printf("| %s | %d | %d | %d | %d | %d | %d | %.0f | %d | %.2f | %.2f | %.1f | %.1f |\n", t.name.c_str(),
                    t.meshes, t.failed, t.watertight, t.nonmanifold, t.euler_ok, t.components_ok,
                    static_cast<double>(t.triangles) / ok, t.max_triangles, t.worst_dev / kTessChord,
                    t.worst_vol / kTessChord, t.total_ms / ok, t.worst_ms);
    }
    std::printf("\n");
    for (const TessTally& t : g_tess)
        for (const std::string& f : t.failures) std::printf("- %s, %s\n", t.name.c_str(), f.c_str());
    std::printf("Tessellation hash: %016llx\n", static_cast<unsigned long long>(g_tess_hash));
}

double now_ms() {
    return std::chrono::duration<double, std::milli>(
               std::chrono::steady_clock::now().time_since_epoch()).count();
}

// ---------------------------------------------------------------- predicates

void bench_predicates() {
    std::printf("## 1. Orientation predicates\n\n");
    std::printf("| Input (1,000,000 calls) | decided: double | interval | exact | exact zeros | naive double wrong | naive inconsistent | ns/call (cascade) | ns/call (naive) |\n");
    std::printf("|---|---|---|---|---|---|---|---|---|\n");
    const int N = 1000000;
    for (int mode = 0; mode < 2; ++mode) {
        akt::Rng rng(1234 + mode);
        std::vector<double> pts(static_cast<std::size_t>(N) * 6);
        for (int i = 0; i < N; ++i) {
            double* p = &pts[static_cast<std::size_t>(i) * 6];
            p[0] = rng.range(-1, 1); p[1] = rng.range(-1, 1);
            p[2] = rng.range(-1, 1); p[3] = rng.range(-1, 1);
            if (mode == 0) { p[4] = rng.range(-1, 1); p[5] = rng.range(-1, 1); }
            else {
                const double t = rng.unit();
                p[4] = p[0] + t * (p[2] - p[0]) + (static_cast<int>(rng.next() % 5) - 2) * 0x1.0p-53;
                p[5] = p[1] + t * (p[3] - p[1]);
            }
        }
        num::predicate_stats() = {};
        long long sum = 0;
        double t0 = now_ms();
        for (int i = 0; i < N; ++i) {
            const double* p = &pts[static_cast<std::size_t>(i) * 6];
            sum += static_cast<int>(num::orient2d(p, p + 2, p + 4));
        }
        const double t_cascade = now_ms() - t0;
        const num::PredicateStats st = num::predicate_stats();

        long long naive_sum = 0;
        t0 = now_ms();
        for (int i = 0; i < N; ++i) {
            const double* p = &pts[static_cast<std::size_t>(i) * 6];
            naive_sum += static_cast<int>(num::orient2d_naive(p, p + 2, p + 4));
        }
        const double t_naive = now_ms() - t0;

        int wrong = 0, inconsistent = 0;
        for (int i = 0; i < N; ++i) {
            const double* p = &pts[static_cast<std::size_t>(i) * 6];
            const num::Sign truth = num::orient2d(p, p + 2, p + 4);
            const num::Sign n1 = num::orient2d_naive(p, p + 2, p + 4);
            const num::Sign n2 = num::orient2d_naive(p + 2, p + 4, p);
            const num::Sign n3 = num::orient2d_naive(p + 4, p, p + 2);
            if (n1 != truth) ++wrong;
            if (n1 != n2 || n2 != n3) ++inconsistent;
        }
        std::printf("| %s | %llu | %llu | %llu | %llu | %d | %d | %.0f | %.0f |\n",
                    mode == 0 ? "uniform random" : "within 2 ulp of collinear",
                    static_cast<unsigned long long>(st.decided_at[0]),
                    static_cast<unsigned long long>(st.decided_at[1]),
                    static_cast<unsigned long long>(st.decided_at[2]),
                    static_cast<unsigned long long>(st.zero_results), wrong, inconsistent,
                    1e6 * t_cascade / N, 1e6 * t_naive / N);
        hash_double(static_cast<double>(sum));
        hash_double(static_cast<double>(naive_sum));
        hash_double(static_cast<double>(wrong));
    }
    std::printf("\n");
}

// ------------------------------------------------------------------------ SSI

struct SsiRow {
    std::size_t curves = 0, segments = 0;
    bool complete = false;
    std::size_t unverified = 0;
    double gap = 0.0, fitted_cubic = 0.0, fitted_linear = 0.0, min_sin = 1.0, ms = 0.0;
};

double surf_gap(const solve::SsiSystem& sys, const solve::Pt4& x) {
    const auto a = sys.patch_a().eval(x[0], x[1]);
    const auto b = sys.patch_b().eval(x[2], x[3]);
    const Vec3d d{a.x / a.w - b.x / b.w, a.y / a.w - b.y / b.w, a.z / a.w - b.z / b.w};
    return geom::norm(d);
}

SsiRow run_ssi(const PolyPatchD& A, const PolyPatchD& B, solve::SsiOptions opt = {}) {
    const solve::SsiSystem sys(A, B);
    SsiRow row;
    const double t0 = now_ms();
    const solve::SsiResult res = solve::intersect(sys, opt);
    row.ms = now_ms() - t0;
    row.curves = res.curves.size();
    row.complete = res.complete;
    row.unverified = res.unverified_boxes;
    for (const auto& c : res.curves) {
        row.segments += c.segs.size();
        row.min_sin = std::min(row.min_sin, c.min_sin_angle);
        // Certified points at the segment junctions, in curve order.
        std::vector<solve::Pt4> pts;
        for (std::size_t i = 0; i < c.segs.size(); ++i) {
            for (double f : {0.0, 1.0}) {
                if (f == 0.0 && i > 0) continue;
                solve::Box4 enc;
                if (!solve::curve_point(sys, c.segs[i], f, enc)) continue;
                const auto pa = sys.point_a(enc), pb = sys.point_b(enc);
                const num::Interval d2 = sqr(pa.x - pb.x) + sqr(pa.y - pb.y) + sqr(pa.z - pb.z);
                row.gap = std::max(row.gap, std::sqrt(d2.hi));
                pts.push_back({mid(enc[0]), mid(enc[1]), mid(enc[2]), mid(enc[3])});
                for (const auto& e : enc) hash_double(mid(e));
            }
        }
        // What interpolating parameter-space curves through the same points give.
        const int n = static_cast<int>(pts.size());
        for (int i = 0; i + 1 < n; ++i) {
            const solve::Pt4& p0 = pts[std::max(i - 1, 0)];
            const solve::Pt4& p1 = pts[i];
            const solve::Pt4& p2 = pts[i + 1];
            const solve::Pt4& p3 = pts[std::min(i + 2, n - 1)];
            solve::Pt4 lin, cub;
            for (int k = 0; k < 4; ++k) {
                lin[k] = 0.5 * (p1[k] + p2[k]);
                cub[k] = (-p0[k] + 9.0 * p1[k] + 9.0 * p2[k] - p3[k]) / 16.0; // Catmull-Rom at t = 1/2
            }
            row.fitted_linear = std::max(row.fitted_linear, surf_gap(sys, lin));
            row.fitted_cubic = std::max(row.fitted_cubic, surf_gap(sys, cub));
        }
    }
    return row;
}

void bench_ssi() {
    std::printf("## 2. Certified surface/surface intersection\n\n");
    std::printf("| Case | curves | segments | proven complete | unresolved boxes | min sin(angle) | certified gap bound | gap of cubic interpolant | gap of linear interpolant | ms |\n");
    std::printf("|---|---|---|---|---|---|---|---|---|---|\n");
    const double swap_xz[9] = {0, 0, 1, 0, 1, 0, 1, 0, 0};
    struct Case { const char* name; PolyPatchD a, b; solve::SsiOptions opt; };
    solve::SsiOptions tangent_opt;
    tangent_opt.verify_depth = 24;
    tangent_opt.verify_node_limit = 200000;
    const std::vector<Case> cases = {
        {"plane x cylinder", geom::make_cylinder_quarter(1.0, 0.0, 2.0),
         geom::make_plane({-0.5, -0.5, 0.35}, {2, 0, 0.6}, {0, 2, 0}), {}},
        {"plane x hill (closed loop)", geom::make_bump(1.0),
         geom::make_plane({-0.5, -0.5, 0.5}, {2, 0, 0}, {0, 2, 0}), {}},
        {"plane above hill (empty)", geom::make_bump(1.0),
         geom::make_plane({-0.5, -0.5, 2.0}, {2, 0, 0}, {0, 2, 0}), {}},
        {"cylinder x cylinder", geom::make_cylinder_quarter(1.0, 0.0, 2.0),
         geom::transformed(geom::make_cylinder_quarter(0.8, -1.0, 2.0), swap_xz, {0, 0, 0}), {}},
        {"bicubic x bicubic",
         geom::fit_patch([](double u, double v) { return Vec3d{u, v, 0.3 * u * u * u - 0.4 * u * v + 0.2 * v * v + 0.1}; }, 3, 3),
         geom::fit_patch([](double s, double t) { return Vec3d{0.1 + 0.8 * s + 0.1 * t * t, 0.9 * t + 0.05, 0.6 * s - 0.2 * t * s * s + 0.05 * t}; }, 3, 3), {}},
        {"plane tangent to hill", geom::make_bump(1.0),
         geom::make_plane({-0.5, -0.5, 1.0}, {2, 0, 0}, {0, 2, 0}), tangent_opt},
    };
    for (const Case& c : cases) {
        const SsiRow r = run_ssi(c.a, c.b, c.opt);
        std::printf("| %s | %zu | %zu | %s | %zu | %.3f | %.1e | %.1e | %.1e | %.1f |\n", c.name, r.curves,
                    r.segments, r.complete ? "yes" : "no", r.unverified, r.min_sin, r.gap,
                    r.fitted_cubic, r.fitted_linear, r.ms);
    }

    // Random freeform pairs: a robustness sample, not a curated success list.
    const int trials = 200;
    akt::Rng rng(777);
    int complete = 0, with_curves = 0, incomplete = 0;
    std::size_t segs = 0;
    double worst_gap = 0.0, total_ms = 0.0, worst_ms = 0.0, lowest_sin = 1.0, worst_fit = 0.0;
    for (int t = 0; t < trials; ++t) {
        PolyPatchD a(3, 3), b(3, 3);
        for (int i = 0; i <= 3; ++i)
            for (int j = 0; j <= 3; ++j) {
                a.at(i, j) = {i / 3.0, j / 3.0, rng.range(-0.15, 0.15), 1.0};
                b.at(i, j) = {i / 3.0, 0.5 + rng.range(-0.15, 0.15), -0.5 + j / 3.0, 1.0};
            }
        const SsiRow r = run_ssi(a, b);
        if (r.complete) ++complete; else ++incomplete;
        if (r.curves > 0) ++with_curves;
        segs += r.segments;
        worst_gap = std::max(worst_gap, r.gap);
        worst_fit = std::max(worst_fit, r.fitted_cubic);
        total_ms += r.ms;
        worst_ms = std::max(worst_ms, r.ms);
        if (r.curves > 0) lowest_sin = std::min(lowest_sin, r.min_sin);
    }
    std::printf("\nRandom bicubic pairs (%d trials, fixed seed): proven complete %d, not proven %d, "
                "with at least one curve %d, segments %zu, worst certified gap bound %.1e, "
                "worst cubic-interpolant gap %.1e, lowest sin(angle) %.3f, mean %.1f ms, worst %.1f ms.\n\n",
                trials, complete, incomplete, with_curves, segs, worst_gap, worst_fit, lowest_sin,
                total_ms / trials, worst_ms);
    hash_double(static_cast<double>(complete));
    hash_double(static_cast<double>(segs));

    // Hard family: two random graphs z = f(x,y), z = g(x,y) over the same
    // square. The intersection is the zero set of f - g, with loops, several
    // components and near-tangencies (saddles of f - g close to zero).
    //
    // Independent cross-check of "proven complete": along fixed lines
    // u = const, count the roots of f - g by dense sampling, and compare with
    // the number of times the returned curves cross that line.
    const double lines[6] = {0.137, 0.291, 0.463, 0.618, 0.794, 0.912};
    int h_complete = 0, h_incomplete = 0, h_curves = 0, h_closed = 0, h_singular_ends = 0;
    int lines_checked = 0, lines_agree = 0, lines_disagree_complete = 0;
    double h_total_ms = 0.0, h_worst_ms = 0.0, h_lowest_sin = 1.0, h_worst_gap = 0.0;
    std::size_t h_segs = 0;
    akt::Rng rng2(4242);
    for (int t = 0; t < trials; ++t) {
        PolyPatchD a(3, 3), b(3, 3);
        for (int i = 0; i <= 3; ++i)
            for (int j = 0; j <= 3; ++j) {
                a.at(i, j) = {i / 3.0, j / 3.0, rng2.range(-0.4, 0.4), 1.0};
                b.at(i, j) = {i / 3.0, j / 3.0, rng2.range(-0.4, 0.4), 1.0};
            }
        const solve::SsiSystem sys(a, b);
        const double t0 = now_ms();
        const solve::SsiResult res = solve::intersect(sys);
        const double ms = now_ms() - t0;
        h_total_ms += ms;
        h_worst_ms = std::max(h_worst_ms, ms);
        if (res.complete) ++h_complete; else ++h_incomplete;
        h_curves += static_cast<int>(res.curves.size());

        // Certified sample points of every curve, in order.
        std::vector<std::vector<double>> curve_u;
        for (const auto& c : res.curves) {
            if (c.closed) ++h_closed;
            if (c.end_front == solve::EndKind::Singular || c.end_back == solve::EndKind::Singular) ++h_singular_ends;
            h_segs += c.segs.size();
            h_lowest_sin = std::min(h_lowest_sin, c.min_sin_angle);
            std::vector<double> us;
            for (const auto& seg : c.segs)
                for (double f : {0.0, 0.25, 0.5, 0.75, 1.0}) {
                    solve::Box4 enc;
                    if (!solve::curve_point(sys, seg, f, enc)) continue;
                    us.push_back(mid(enc[0]));
                    const auto pa = sys.point_a(enc), pb = sys.point_b(enc);
                    h_worst_gap = std::max(h_worst_gap, std::sqrt((sqr(pa.x - pb.x) + sqr(pa.y - pb.y) + sqr(pa.z - pb.z)).hi));
                }
            if (c.closed && !us.empty()) us.push_back(us.front());
            curve_u.push_back(std::move(us));
        }
        for (double u0 : lines) {
            int truth = 0;
            double prev = 0.0;
            for (int k = 0; k <= 20000; ++k) {
                const double v = k / 20000.0;
                const double d = a.eval(u0, v).z - b.eval(u0, v).z;
                if (k > 0 && ((prev < 0.0) != (d < 0.0))) ++truth;
                prev = d;
            }
            int found = 0;
            for (const auto& us : curve_u)
                for (std::size_t k = 1; k < us.size(); ++k)
                    if ((us[k - 1] < u0) != (us[k] < u0)) ++found;
            ++lines_checked;
            if (truth == found) ++lines_agree;
            else if (res.complete) ++lines_disagree_complete;
        }
        hash_double(static_cast<double>(res.curves.size()));
    }
    std::printf("Random graph pairs with loops and near-tangencies (%d trials, fixed seed): proven complete %d, "
                "not proven %d; curves %d (closed %d, ending at a suspected tangency %d), segments %zu; "
                "lowest sin(angle) %.4f; worst certified gap bound %.1e; mean %.1f ms, worst %.1f ms.\n\n",
                trials, h_complete, h_incomplete, h_curves, h_closed, h_singular_ends, h_segs, h_lowest_sin,
                h_worst_gap, h_total_ms / trials, h_worst_ms);
    std::printf("Independent cross-check on %d scan lines: kernel crossing count equals the sampled root "
                "count on %d; disagreements in cases the kernel called complete: %d.\n\n",
                lines_checked, lines_agree, lines_disagree_complete);
    hash_double(static_cast<double>(h_complete));
    hash_double(static_cast<double>(h_curves));
}

// ----------------------------------------------------------------------- fill

struct Wave {
    static double dz(int k, int l, double x, double y) {
        const double hp = 1.57079632679489661923;
        double s = 0.2 * std::pow(2.5, k) * std::sin(2.5 * x + 0.3 + k * hp) * std::pow(2.0, l) *
                   std::cos(2.0 * y + l * hp);
        if (k == 0 && l == 0) s += 0.1 * x * y;
        if (k == 1 && l == 0) s += 0.1 * y;
        if (k == 0 && l == 1) s += 0.1 * x;
        if (k == 1 && l == 1) s += 0.1;
        return s;
    }
    static geom::Jet3 jet(double x, double y) {
        geom::Jet3 J{};
        for (int k = 0; k <= 3; ++k)
            for (int l = 0; l <= 3 - k; ++l) J[geom::jet_index(k, l)] = {0, 0, dz(k, l, x, y)};
        J[0].x = x; J[0].y = y; J[1].x = 1.0; J[2].y = 1.0;
        return J;
    }
    static Vec3d side_jet(int side, int j, double s) {
        const bool u_side = side == var::U0 || side == var::U1;
        const double x = u_side ? (side == var::U1 ? 1.0 : 0.0) : s;
        const double y = u_side ? s : (side == var::V1 ? 1.0 : 0.0);
        return u_side ? jet(x, y)[geom::jet_index(j, 0)] : jet(x, y)[geom::jet_index(0, j)];
    }
};

void bench_fill() {
    std::printf("## 3. Variational G3 fill (degree 7, energy order 4)\n\n");
    std::printf("| Control points | free | G1 jump (rad) | G2 jump b(nu,nu) | G3 jump C(nu,nu,nu) | largest other component | solve ms |\n");
    std::printf("|---|---|---|---|---|---|---|\n");
    std::vector<double> hs, g1s, g2s, g3s;
    for (int n : {12, 16, 20, 24, 32}) {
        var::FillSpec spec;
        spec.ncp = n;
        const double t0 = now_ms();
        const var::FillResult r = var::solve_fill(spec, Wave::side_jet);
        const double ms = now_ms() - t0;
        double g1 = 0, g2 = 0, g3 = 0, other = 0;
        for (int side = 0; side < 4; ++side)
            for (int i = 1; i < 40; ++i) {
                const double t = i / 40.0;
                const bool u_side = side == var::U0 || side == var::U1;
                const double x = u_side ? (side == var::U1 ? 1.0 : 0.0) : t;
                const double y = u_side ? t : (side == var::V1 ? 1.0 : 0.0);
                const geom::Jet3 JQ = Wave::jet(x, y), JS = r.surface.jet(x, y);
                const geom::SeamJump j = geom::seam_jump(JQ, JS, u_side ? JQ[2] : JQ[1]);
                g1 = std::max(g1, j.normal_angle);
                g2 = std::max(g2, std::fabs(j.db_nn));
                g3 = std::max(g3, std::fabs(j.dC_nnn));
                for (double v : {j.db_tt, j.db_tn, j.dC_ttt, j.dC_ttn, j.dC_tnn}) other = std::max(other, std::fabs(v));
            }
        std::printf("| %d x %d | %d | %.2e | %.2e | %.2e | %.2e | %.1f |\n", n, n, r.free_dofs, g1, g2, g3, other, ms);
        hs.push_back(1.0 / (n - spec.degree));
        g1s.push_back(g1); g2s.push_back(g2); g3s.push_back(g3);
        for (const auto& c : r.surface.cp) hash_double(c.z);
    }
    std::printf("\nObserved convergence orders in the knot spacing: G1 %.1f, G2 %.1f, G3 %.1f.\n\n",
                diag::fit_rate(hs, g1s).order, diag::fit_rate(hs, g2s).order, diag::fit_rate(hs, g3s).order);

    std::printf("## 4. Energy-order rule (G2 jets prescribed)\n\n");
    std::printf("| Control points | E3 of thin-plate (m=2) solution | E3 of matched (m=3) solution |\n|---|---|---|\n");
    for (int n : {12, 18, 27, 40}) {
        var::FillSpec plate;
        plate.ncp = n;
        plate.continuity = 2;
        plate.energy_weights = {0, 0, 1, 0, 0};
        var::FillSpec matched = plate;
        matched.energy_weights = {0, 0, 0, 1, 0};
        const var::FillResult a = var::solve_fill(plate, Wave::side_jet), b = var::solve_fill(matched, Wave::side_jet);
        std::printf("| %d x %d | %.4f | %.4f |\n", n, n, a.energy[3], b.energy[3]);
        hash_double(a.energy[3]);
        hash_double(b.energy[3]);
    }
    std::printf("\n");
}

// ------------------------------------------------------------------ Booleans

/// An analytic solid: membership test and bounding box, independent of the kernel.
struct Analytic {
    std::function<bool(const Vec3d&)> inside;
    Vec3d lo, hi;
};

struct RandomBox {
    Vec3d o, ex, ey, ez;
    op::PatchSolid solid() const { return op::make_box(o, ex, ey, ez); }
    bool inside(const Vec3d& p) const {
        const Vec3d d = p - o;
        const double x = dot(d, ex) / dot(ex, ex), y = dot(d, ey) / dot(ey, ey), z = dot(d, ez) / dot(ez, ez);
        return x >= 0 && x <= 1 && y >= 0 && y <= 1 && z >= 0 && z <= 1;
    }
};

/// A box with random orientation and size, centred at a random point of [clo, chi].
RandomBox random_box(akt::Rng& rng, const Vec3d& clo, const Vec3d& chi) {
    Vec3d axis;
    double len = 0.0;
    do {
        axis = {rng.range(-1, 1), rng.range(-1, 1), rng.range(-1, 1)};
        len = geom::norm(axis);
    } while (len < 0.2 || len > 1.0);
    axis = axis / len;
    const double ang = rng.range(0.2, 2.9), c = std::cos(ang), s = std::sin(ang);
    const double ax = axis.x, ay = axis.y, az = axis.z;
    const double R[3][3] = {{c + ax * ax * (1 - c), ax * ay * (1 - c) - az * s, ax * az * (1 - c) + ay * s},
                            {ay * ax * (1 - c) + az * s, c + ay * ay * (1 - c), ay * az * (1 - c) - ax * s},
                            {az * ax * (1 - c) - ay * s, az * ay * (1 - c) + ax * s, c + az * az * (1 - c)}};
    const double sx = rng.range(0.4, 1.1), sy = rng.range(0.4, 1.1), sz = rng.range(0.4, 1.1);
    RandomBox b;
    b.ex = {R[0][0] * sx, R[1][0] * sx, R[2][0] * sx};
    b.ey = {R[0][1] * sy, R[1][1] * sy, R[2][1] * sy};
    b.ez = {R[0][2] * sz, R[1][2] * sz, R[2][2] * sz};
    const Vec3d centre{rng.range(clo.x, chi.x), rng.range(clo.y, chi.y), rng.range(clo.z, chi.z)};
    b.o = centre - (b.ex + b.ey + b.ez) * 0.5;
    return b;
}

/// Monte Carlo volume of A and B over A's bounding box, with its standard error.
double monte_carlo(const Analytic& a, const RandomBox& b, int samples, std::uint64_t seed, double& sigma) {
    akt::Rng rng(seed);
    int hits = 0;
    for (int i = 0; i < samples; ++i) {
        const Vec3d p{rng.range(a.lo.x, a.hi.x), rng.range(a.lo.y, a.hi.y), rng.range(a.lo.z, a.hi.z)};
        if (a.inside(p) && b.inside(p)) ++hits;
    }
    const double box = (a.hi.x - a.lo.x) * (a.hi.y - a.lo.y) * (a.hi.z - a.lo.z);
    const double p = static_cast<double>(hits) / samples;
    sigma = box * std::sqrt(std::max(p * (1.0 - p), 1.0 / samples) / samples);
    return box * p;
}

void bench_boolean_family(const char* name, const op::PatchSolid& patch_a, const Analytic& a, const Vec3d& clo,
                          const Vec3d& chi, std::uint64_t seed, int trials) {
    akt::Rng rng(seed);
    const op::Brep A = op::to_brep(patch_a);
    const double va = op::volume(A);
    int ok = 0, audited = 0, empty_inter = 0, beyond3 = 0;
    double worst_ie = 0.0, worst_z = 0.0, worst_gap = 0.0, worst_width = 0.0, total_ms = 0.0, worst_ms = 0.0;
    long curves = 0, vertices = 0;
    std::map<std::string, int> reasons;
    for (int trial = 0; trial < trials; ++trial) {
        const RandomBox rb = random_box(rng, clo, chi);
        const op::Brep B = op::to_brep(rb.solid());
        const double t0 = now_ms();
        op::Imprint im = op::imprint(A, B, with_exact());
        const op::BooleanResult ri = op::select(im, op::Op::Intersection);
        const op::BooleanResult ru = op::select(im, op::Op::Union);
        const op::BooleanResult rd = op::select(im, op::Op::Difference);
        const double ms = now_ms() - t0;
        total_ms += ms;
        worst_ms = std::max(worst_ms, ms);
        if (!im.ok) {
            // Group by the part of the message after the entity names.
            const std::size_t colon = im.error.find(": ");
            ++reasons[colon == std::string::npos ? im.error : im.error.substr(colon + 2)];
            hash_double(-1.0);
            continue;
        }
        ++ok;
        if (ri.ok && ru.ok && rd.ok) ++audited;
        if (ri.empty) ++empty_inter;
        const double vi = op::volume(ri.solid), vu = op::volume(ru.solid), vd = op::volume(rd.solid);
        const double vb = op::volume(B);
        worst_ie = std::max(worst_ie, std::max(std::fabs(vi + vu - va - vb), std::fabs(vd - (va - vi))));
        double sigma = 0.0;
        const double mc = monte_carlo(a, rb, 100000, seed * 1000 + static_cast<std::uint64_t>(trial), sigma);
        const double z = std::fabs(mc - vi) / sigma;
        worst_z = std::max(worst_z, z);
        if (z > 3.0) ++beyond3;
        worst_gap = std::max(worst_gap, op::edge_gap_bound(im));
        worst_width = std::max(worst_width, im.stats.max_vertex_width);
        curves += im.stats.curves;
        vertices += im.stats.new_vertices;
        hash_double(vi);
        hash_double(vu);
        hash_double(vd);
        hash_double(static_cast<double>(im.stats.curves));
        hash_double(static_cast<double>(ri.solid.faces.size() + ru.solid.faces.size() + rd.solid.faces.size()));
        for (const op::BooleanResult* r : {&ri, &ru, &rd})
            if (r->ok && !r->empty)
                tally_tessellation(tally_named("5. random boxes against curved solids"), r->solid,
                                   std::string(name) + ", trial " + std::to_string(trial));
    }
    std::printf("| %s | %d | %d | %d | %d | %ld | %ld | %.1e | %.2f | %d | %.1e | %.1e | %.0f | %.0f |\n", name,
                trials, ok, audited, empty_inter, curves, vertices, worst_ie, worst_z, beyond3, worst_gap,
                worst_width, total_ms / trials, worst_ms);
    for (const auto& kv : reasons) std::printf("|   not solved: %s | %d | | | | | | | | | | | | |\n", kv.first.c_str(), kv.second);
}

void bench_boolean() {
    std::printf("## 5. Boolean operations on patch solids\n\n");

    std::printf("| Case | Curves | New vertices | Regions | Volume of A and B | Reference | Difference | "
                "Inclusion-exclusion residual | All three results closed, genus 0 | Certified edge gap bound | "
                "Widest vertex enclosure | Time, ms |\n|---|---|---|---|---|---|---|---|---|---|---|---|\n");
    struct Case {
        const char* name;
        op::PatchSolid a, b;
        double reference;
        const char* reference_kind;
    };
    auto prim = [](double x) { return 0.5 * (x * std::sqrt(1.0 - x * x) + std::asin(x)); };
    // Cap of the hill above the cutting plane, by midpoint quadrature.
    double cap = 0.0;
    {
        const int N = 2000;
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j) {
                const double x = (i + 0.5) / N, y = (j + 0.5) / N;
                cap += std::max(0.0, 1.0 + 4.8 * x * (1 - x) * y * (1 - y) - 1.15);
            }
        cap /= static_cast<double>(N) * N;
    }
    const std::vector<Case> cases = {
        {"Two offset boxes", op::make_box({0, 0, 0}, {1, 1, 1}), op::make_box({0.5, 0.4, 0.3}, {1.5, 1.4, 1.3}),
         0.21, "closed form"},
        {"Hill cut by a plane (closed loop)", op::make_bump_box({0, 0, 0}, 1, 1, 1, 0.3),
         op::make_box({-0.5, -0.5, 1.15}, {1.5, 1.5, 2.0}), cap, "2D quadrature"},
        {"Cylindrical shell and box (rational faces)", op::make_annular_sector(1.0, 2.0, 0.0, 1.0),
         op::make_box({0.5, 0.6, 0.25}, {1.3, 1.4, 0.75}), (0.64 - (prim(0.8) - prim(0.5) - 0.18)) * 0.5,
         "closed form"},
        {"Box inside a box (no contact)", op::make_box({0, 0, 0}, {1, 1, 1}),
         op::make_box({0.3, 0.3, 0.3}, {0.6, 0.6, 0.6}), 0.027, "closed form"},
    };
    for (const Case& c : cases) {
        const op::Brep A = op::to_brep(c.a), B = op::to_brep(c.b);
        const double t0 = now_ms();
        op::Imprint im = op::imprint(A, B, with_exact());
        const op::BooleanResult ri = op::select(im, op::Op::Intersection);
        const op::BooleanResult ru = op::select(im, op::Op::Union);
        const op::BooleanResult rd = op::select(im, op::Op::Difference);
        const double ms = now_ms() - t0;
        if (!im.ok) { std::printf("| %s | not solved: %s |\n", c.name, im.error.c_str()); continue; }
        const double vi = op::volume(ri.solid), vu = op::volume(ru.solid), vd = op::volume(rd.solid);
        const double va = op::volume(c.a), vb = op::volume(c.b);
        auto closed = [](const op::BooleanResult& r) {
            if (!r.ok || r.empty) return false;
            for (const auto& sh : r.audit.shells)
                if (sh.genus != 0) return false;
            return true;
        };
        std::printf("| %s | %d | %d | %d | %.12f | %.12f (%s) | %.1e | %.1e | %s | %.1e | %.1e | %.0f |\n", c.name,
                    im.stats.curves, im.stats.new_vertices, im.stats.regions, vi, c.reference, c.reference_kind,
                    std::fabs(vi - c.reference),
                    std::max(std::fabs(vi + vu - va - vb), std::fabs(vd - (va - vi))),
                    closed(ri) && closed(ru) && closed(rd) ? "yes" : "no", op::edge_gap_bound(im),
                    im.stats.max_vertex_width, ms);
        hash_double(vi);
        hash_double(vu);
        hash_double(vd);
        hash_double(static_cast<double>(ri.solid.edges.size() + ru.solid.edges.size() + rd.solid.edges.size()));
    }

    std::printf("\nRandom pairs: a fixed solid A against a box of random orientation, size and position "
                "(fixed seeds). Volumes are checked against a Monte Carlo count of 100,000 points per pair "
                "that uses only the analytic definitions of the two solids.\n\n");
    std::printf("| A | Pairs | Solved | All three results pass the audit | Empty intersection | Curves | "
                "New vertices | Worst inclusion-exclusion residual | Worst Monte Carlo deviation, sigma | "
                "Pairs beyond 3 sigma | Worst certified edge gap bound | Widest vertex enclosure | Mean time, ms | "
                "Worst time, ms |\n|---|---|---|---|---|---|---|---|---|---|---|---|---|---|\n");
    const Analytic box{[](const Vec3d& p) {
                           return p.x >= 0 && p.x <= 1 && p.y >= 0 && p.y <= 1 && p.z >= 0 && p.z <= 1;
                       },
                       {0, 0, 0}, {1, 1, 1}};
    const Analytic hill{[](const Vec3d& p) {
                            return p.x >= 0 && p.x <= 1 && p.y >= 0 && p.y <= 1 && p.z >= 0 &&
                                   p.z <= 1.0 + 4.8 * p.x * (1 - p.x) * p.y * (1 - p.y);
                        },
                        {0, 0, 0}, {1, 1, 1.3}};
    const Analytic shell{[](const Vec3d& p) {
                             const double r2 = p.x * p.x + p.y * p.y;
                             return p.x >= 0 && p.y >= 0 && r2 >= 1.0 && r2 <= 4.0 && p.z >= 0 && p.z <= 1;
                         },
                         {0, 0, 0}, {2, 2, 1}};
    bench_boolean_family("Unit box", op::make_box({0, 0, 0}, {1, 1, 1}), box, {0.1, 0.1, 0.1}, {0.9, 0.9, 0.9}, 11, 100);
    bench_boolean_family("Box with a hill on top", op::make_bump_box({0, 0, 0}, 1, 1, 1, 0.3), hill,
                         {0.1, 0.1, 0.5}, {0.9, 0.9, 1.4}, 12, 100);
    bench_boolean_family("Quarter cylindrical shell", op::make_annular_sector(1.0, 2.0, 0.0, 1.0), shell,
                         {0.4, 0.4, 0.1}, {1.6, 1.6, 0.9}, 13, 100);
    std::printf("\n");
}

// --------------------------------------------------------- chained Booleans

/// Chains of operations: each result is the first operand of the next
/// operation, against a fresh box of random orientation, size and position.
/// With `reuse` > 0, a step uses, with that probability, a box already used
/// earlier in the same chain: the operands then share surfaces by construction.
void bench_boolean_chains(int chains, int depth, int first_chain = 0, bool trace = false, double reuse = 0.0,
                          bool grid = false) {
    if (grid) {
        std::printf("## 8. Chained operations on grid-aligned boxes\n\n");
        std::printf("%d chains as in section 6, but every box is axis-aligned with its corners on a grid of "
                    "spacing 0.25 (sides 0.25 to 1, corners from -0.5 to 1.5). Faces of different boxes then lie "
                    "in the same plane, boxes touch face to face, and edges overlap or cross inside shared planes: "
                    "contact that is not shared by construction and has to be decided exactly.\n\n",
                    chains);
    } else if (reuse > 0.0) {
        std::printf("## 7. Chained operations that reuse earlier operands\n\n");
        std::printf("%d chains as in section 6, but at each step after the first, with probability %.0f %%, the "
                    "box is one already used earlier in the same chain, so the two operands share surfaces, "
                    "edges and vertices by construction. Seeds differ from section 6.\n\n",
                    chains, 100.0 * reuse);
    } else {
        std::printf("## 6. Chained Boolean operations\n\n");
        std::printf("%d chains from the unit box, up to %d operations each. The operation is drawn at random "
                    "(difference 45 %%, union 35 %%, intersection 20 %%). After every step the volume is compared "
                    "with a Monte Carlo count of 200,000 points that uses only the analytic definition of the "
                    "chain. A chain stops when a step is not solved or its result is empty.\n\n",
                    chains, depth);
    }
    int reused_steps = 0, shared_steps = 0, nonmanifold_steps = 0;
    struct Row {
        int attempted = 0, solved = 0, audited = 0, empty = 0, beyond3 = 0;
        double worst_z = 0.0, worst_gap = 0.0, worst_width = 0.0, total_ms = 0.0, worst_ms = 0.0;
        long faces = 0, edges = 0, arcs = 0, triples = 0;
        double stage_ms[7] = {0, 0, 0, 0, 0, 0, 0};
        double pairs = 0.0, pairs_solved = 0.0;
        double ssi_ms_empty = 0, ssi_ms_curves = 0, pairs_empty = 0, nodes_empty = 0, nodes_curves = 0;
        double edge_straight = 0, edge_arc = 0, edge_ms_straight = 0, edge_ms_arc = 0, csi_nodes = 0, tsi_nodes = 0;
    };
    std::vector<Row> rows(static_cast<std::size_t>(depth));
    TessTally& tess =
        tally_named(grid ? "8. grid-aligned boxes" : (reuse > 0.0 ? "7. reused operands" : "6. random boxes"));
    std::map<std::string, int> reasons;
    std::vector<std::string> failures;
    const Vec3d lo{-1.2, -1.2, -1.2}, hi{2.2, 2.2, 2.2};
    long total_entities[3] = {0, 0, 0}; // faces, edges and vertices of every committed result
    op::BooleanOptions options = with_exact();
    options.trace = trace;
    options.plane_lines = g_plane_lines;
    for (int chain = first_chain; chain < first_chain + chains; ++chain) {
        // Each chain can be replayed on its own.
        akt::Rng rng((grid ? 5026 : (reuse > 0.0 ? 9026 : 2026)) + static_cast<std::uint64_t>(chain));
        op::Brep current = op::to_brep(op::make_box({0, 0, 0}, {1, 1, 1}));
        std::function<bool(const Vec3d&)> member = [](const Vec3d& p) {
            return p.x >= 0 && p.x <= 1 && p.y >= 0 && p.y <= 1 && p.z >= 0 && p.z <= 1;
        };
        std::vector<RandomBox> used;
        for (int step = 0; step < depth; ++step) {
            Row& row = rows[static_cast<std::size_t>(step)];
            RandomBox rb;
            bool reused = false;
            if (reuse > 0.0 && !used.empty() && rng.unit() < reuse) {
                rb = used[static_cast<std::size_t>(rng.next() % used.size())];
                reused = true;
                ++reused_steps;
            } else if (grid) {
                auto cell = [&](int lo, int hi) { return lo + static_cast<int>(rng.next() % static_cast<std::uint64_t>(hi - lo + 1)); };
                const Vec3d corner{0.25 * cell(-2, 5), 0.25 * cell(-2, 5), 0.25 * cell(-2, 5)};
                const Vec3d size{0.25 * cell(1, 4), 0.25 * cell(1, 4), 0.25 * cell(1, 4)};
                rb.o = corner;
                rb.ex = {size.x, 0, 0};
                rb.ey = {0, size.y, 0};
                rb.ez = {0, 0, size.z};
                used.push_back(rb);
            } else {
                rb = random_box(rng, {0.1, 0.1, 0.1}, {0.9, 0.9, 0.9});
                used.push_back(rb);
            }
            const double pick = rng.unit();
            const op::Op o = pick < 0.45 ? op::Op::Difference : (pick < 0.8 ? op::Op::Union : op::Op::Intersection);
            const op::Brep B = op::to_brep(rb.solid());
            ++row.attempted;
            const double t0 = now_ms();
            const op::Imprint im = op::imprint(current, B, options);
            const op::BooleanResult r = op::select(im, o);
            const double ms = now_ms() - t0;
            if (trace && !im.ok) {
                std::printf("chain %d step %d: not solved: %s\n", chain, step + 1, im.error.c_str());
                for (const std::string& line : im.trace) std::printf("%s\n", line.c_str());
            }
            row.total_ms += ms;
            row.worst_ms = std::max(row.worst_ms, ms);
            for (int i = 0; i < 7; ++i) row.stage_ms[i] += im.stats.ms[i];
            row.pairs += im.stats.carrier_pairs;
            row.pairs_solved += im.stats.carrier_pairs_solved;
            row.ssi_ms_empty += im.stats.ssi_ms_empty;
            row.ssi_ms_curves += im.stats.ssi_ms_curves;
            row.pairs_empty += im.stats.ssi_pairs_empty;
            row.nodes_empty += static_cast<double>(im.stats.ssi_nodes_empty);
            row.nodes_curves += static_cast<double>(im.stats.ssi_nodes_curves);
            row.edge_straight += im.stats.edge_searches_straight;
            row.edge_arc += im.stats.edge_searches_arc;
            row.edge_ms_straight += im.stats.edge_ms_straight;
            row.edge_ms_arc += im.stats.edge_ms_arc;
            row.csi_nodes += static_cast<double>(im.stats.csi_nodes);
            row.tsi_nodes += static_cast<double>(im.stats.tsi_nodes);
            if (!im.ok) {
                const std::size_t colon = im.error.find(": ");
                ++reasons[colon == std::string::npos ? im.error : im.error.substr(colon + 2)];
                failures.push_back("chain " + std::to_string(chain) + ", step " + std::to_string(step + 1) +
                                   (reused ? " (reused box)" : "") + " (stage " + std::to_string(im.stats.stage) +
                                   "): " + im.error);
                hash_double(-1.0);
                break;
            }
            ++row.solved;
            if (im.stats.shared_carriers > 0) ++shared_steps;
            if (r.ok) ++row.audited;
            if (r.ok && !r.empty && !r.manifold) ++nonmanifold_steps;
            if (trace) {
                std::printf("chain %d step %d: op %s, reused %d, ok %d, shared surfaces %d, on-regions %d, regions %d\n",
                            chain, step + 1, o == op::Op::Difference ? "difference" : (o == op::Op::Union ? "union" : "intersection"),
                            reused ? 1 : 0, r.ok ? 1 : 0, im.stats.shared_carriers, im.stats.on_regions,
                            im.stats.regions);
            }
            if (trace && !r.ok) {
                static const char* cls[] = {"out", "in", "on-same", "on-opposite"};
                for (std::size_t ri = 0; ri < im.regions.size(); ++ri) {
                    const op::Region& reg = im.regions[ri];
                    const op::BFace& f = im.solid[reg.solid]->faces[static_cast<std::size_t>(reg.face)];
                    const int g = im.carrier_of[reg.solid][static_cast<std::size_t>(f.carrier)];
                    int fresh = 0, own = 0, foreign = 0;
                    for (const auto& loop : reg.loops)
                        for (const op::Piece& p : loop) {
                            if (p.fresh) ++fresh;
                            else if (im.pieces[static_cast<std::size_t>(p.src)].solid == reg.solid) ++own;
                            else ++foreign;
                        }
                    std::printf("  region %zu: solid %d face %d surface %d%s flipped %d: %s; loops %zu, pieces new %d own %d other %d\n",
                                ri, reg.solid, reg.face, g, " ", f.flipped ? 1 : 0,
                                reg.inside >= 0 ? cls[reg.inside] : "unknown", reg.loops.size(), fresh, own, foreign);
                }
                // Which result edges are used more than twice, and by which regions.
                std::map<int, std::vector<int>> uses;
                for (std::size_t fi = 0; fi < r.solid.faces.size(); ++fi)
                    for (const auto& loop : r.solid.faces[fi].loops)
                        for (const auto& ce : loop) uses[ce.edge].push_back(r.regions[fi]);
                for (const auto& kv : uses)
                    if (kv.second.size() > 2) {
                        std::printf("  result edge %d used by regions", kv.first);
                        for (int x : kv.second) std::printf(" %d", x);
                        const auto& e = r.solid.edges[static_cast<std::size_t>(kv.first)];
                        std::printf(" (%s)\n", e.arc ? "arc" : "side");
                    }
            }
            if (!r.ok) {
                const auto& a = r.audit;
                failures.push_back("chain " + std::to_string(chain) + ", step " + std::to_string(step + 1) +
                                   (reused ? " (reused box)" : "") + ": result rejected by the audit: boundary edges " +
                                   std::to_string(a.boundary_edges) + ", non-manifold " + std::to_string(a.nonmanifold_edges) +
                                   ", pinched vertices " + std::to_string(a.nonmanifold_vertices) +
                                   ", misoriented " + std::to_string(a.misoriented_edges) + ", loops closed " +
                                   (a.loops_closed ? "yes" : "no") + ", shared surfaces " +
                                   std::to_string(im.stats.shared_carriers));
                break;
            }
            const auto prev = member;
            if (o == op::Op::Difference) member = [=](const Vec3d& p) { return prev(p) && !rb.inside(p); };
            else if (o == op::Op::Union) member = [=](const Vec3d& p) { return prev(p) || rb.inside(p); };
            else member = [=](const Vec3d& p) { return prev(p) && rb.inside(p); };
            if (r.empty) { ++row.empty; hash_double(0.0); break; }

            const double vol = op::volume(r.solid);
            akt::Rng mc_rng(static_cast<std::uint64_t>(chain * 100 + step + 1));
            const int samples = 200000;
            int hits = 0;
            for (int i = 0; i < samples; ++i)
                if (member({mc_rng.range(lo.x, hi.x), mc_rng.range(lo.y, hi.y), mc_rng.range(lo.z, hi.z)})) ++hits;
            const double box = (hi.x - lo.x) * (hi.y - lo.y) * (hi.z - lo.z);
            const double p = static_cast<double>(hits) / samples;
            const double sigma = box * std::sqrt(std::max(p * (1.0 - p), 1.0 / samples) / samples);
            const double z = std::fabs(box * p - vol) / sigma;
            row.worst_z = std::max(row.worst_z, z);
            if (z > 3.0) ++row.beyond3;
            row.worst_gap = std::max(row.worst_gap, op::edge_gap_bound(r.solid));
            row.worst_width = std::max(row.worst_width, op::vertex_width(r.solid));
            row.faces += static_cast<long>(r.solid.faces.size());
            row.edges += static_cast<long>(r.solid.edges.size());
            for (const auto& e : r.solid.edges) row.arcs += e.arc ? 1 : 0;
            row.triples += im.stats.triple_vertices;
            hash_double(vol);
            hash_double(static_cast<double>(r.solid.faces.size() * 1000 + r.solid.edges.size()));
            total_entities[0] += static_cast<long>(r.solid.faces.size());
            total_entities[1] += static_cast<long>(r.solid.edges.size());
            total_entities[2] += static_cast<long>(r.solid.vertices.size());
            tally_tessellation(tess, r.solid, "chain " + std::to_string(chain) + ", step " + std::to_string(step + 1));
            current = r.solid;
        }
    }
    std::printf("Over all committed results: %ld faces, %ld edges, %ld vertices.\n\n", total_entities[0], total_entities[1],
                total_entities[2]);
    std::printf("| Step | Attempted | Solved | Result passes the audit | Empty result | Mean faces | Mean edges | "
                "Of those, intersection edges | Three-surface vertices created | Worst Monte Carlo deviation, sigma | "
                "Beyond 3 sigma | Worst certified edge gap bound | Widest vertex enclosure | Mean time, ms | "
                "Worst time, ms |\n|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|\n");
    for (int step = 0; step < depth; ++step) {
        const Row& row = rows[static_cast<std::size_t>(step)];
        const int kept = std::max(1, row.audited - row.empty);
        std::printf("| %d | %d | %d | %d | %d | %.1f | %.1f | %.1f | %ld | %.2f | %d | %.1e | %.1e | %.0f | %.0f |\n",
                    step + 1, row.attempted, row.solved, row.audited, row.empty,
                    static_cast<double>(row.faces) / kept, static_cast<double>(row.edges) / kept,
                    static_cast<double>(row.arcs) / kept, row.triples, row.worst_z, row.beyond3, row.worst_gap,
                    row.worst_width, row.total_ms / std::max(1, row.attempted), row.worst_ms);
    }
    std::printf("\nMean time per stage, ms, and how many surface pairs had to be intersected:\n\n"
                "| Step | Operands | Surface/surface | Edge/surface | Clipping | Edge splitting | Regions | "
                "Classification | Surface pairs | Of those, intersected |\n|---|---|---|---|---|---|---|---|---|---|\n");
    for (int step = 0; step < depth; ++step) {
        const Row& row = rows[static_cast<std::size_t>(step)];
        const double n = std::max(1, row.attempted);
        std::printf("| %d | %.0f | %.0f | %.0f | %.0f | %.0f | %.0f | %.0f | %.0f | %.0f |\n", step + 1,
                    row.stage_ms[0] / n, row.stage_ms[1] / n, row.stage_ms[2] / n, row.stage_ms[3] / n,
                    row.stage_ms[4] / n, row.stage_ms[5] / n, row.stage_ms[6] / n, row.pairs / n, row.pairs_solved / n);
    }
    if (g_profile) {
        std::printf("\nSurface pairs, per step (means): intersected pairs that gave no curve, their time and "
                    "completeness nodes; pairs with curves, likewise.\n\n| Step | Pairs without curves | ms | nodes | "
                    "Pairs with curves | ms | nodes |\n|---|---|---|---|---|---|---|\n");
        for (int step = 0; step < depth; ++step) {
            const Row& row = rows[static_cast<std::size_t>(step)];
            const double n = std::max(1, row.attempted);
            std::printf("| %d | %.1f | %.0f | %.0f | %.1f | %.0f | %.0f |\n", step + 1, row.pairs_empty / n,
                        row.ssi_ms_empty / n, row.nodes_empty / n, row.pairs_solved / n - row.pairs_empty / n,
                        row.ssi_ms_curves / n, row.nodes_curves / n);
        }
        std::printf("\nEdge/surface searches, per step (means): straight edges (curve/surface) and arcs "
                    "(three-surface points along a tube).\n\n| Step | Straight | ms | nodes | Arcs | ms | nodes |\n"
                    "|---|---|---|---|---|---|---|\n");
        for (int step = 0; step < depth; ++step) {
            const Row& row = rows[static_cast<std::size_t>(step)];
            const double n = std::max(1, row.attempted);
            std::printf("| %d | %.1f | %.1f | %.0f | %.1f | %.1f | %.0f |\n", step + 1, row.edge_straight / n,
                        row.edge_ms_straight / n, row.csi_nodes / n, row.edge_arc / n, row.edge_ms_arc / n,
                        row.tsi_nodes / n);
        }
    }
    if (reuse > 0.0)
        std::printf("\nSteps that reused an earlier box: %d. Solved steps whose operands shared at least one "
                    "surface by construction: %d.\n",
                    reused_steps, shared_steps);
    if (grid)
        std::printf("\nSolved steps in which the operands had faces in a common plane: %d.\n", shared_steps);
    std::printf("\nCommitted results that are not manifold (parts touching along an edge or at a point): %d.\n", nonmanifold_steps);
    for (const auto& kv : reasons) std::printf("\nNot solved (%d): %s\n", kv.second, kv.first.c_str());
    std::printf("\n");
    for (const std::string& f : failures) std::printf("- %s\n", f.c_str());
    std::printf("\n");
}

// ------------------------------------------------------------- identities

/// Algebraic identities (blueprint §5.4): two routes to one solid must give the
/// same solid with the same topology (op::compare_solids), not only the same
/// volume. Each route recomputes every operation from scratch.
/// A result the audit refused: the audit's counts, every region of the imprint
/// with its classification, and the edges the selection leaves unbalanced
/// (used a different number of times in each direction). For finding the
/// decision that was wrong.
void dump_refused(const op::Imprint& im, const op::BooleanResult& r) {
    const topo::AuditReport& a = r.audit;
    std::printf("  audit: loops closed %d, chain condition %d, boundary edges %d, misoriented %d, unused %d, more than "
                "two faces %d, open fans %d, vertices with several fans %d, shells %d\n",
                a.loops_closed ? 1 : 0, a.chain_condition ? 1 : 0, a.boundary_edges, a.misoriented_edges, a.unused_edges,
                a.nonmanifold_edges, a.open_fans, a.nonmanifold_vertices, static_cast<int>(a.shells.size()));
    static const char* const kClass[] = {"unknown", "out", "in", "on, same normal", "on, opposite normal"};
    std::vector<char> kept(im.regions.size(), 0);
    for (int ri : r.regions) kept[static_cast<std::size_t>(ri)] = 1;
    std::map<std::pair<bool, int>, std::pair<int, int>> uses; // (fresh, id) -> forward, backward
    auto ends = [&](const op::Piece& p, int& v0, int& v1) {
        if (p.fresh) { v0 = im.arcs[static_cast<std::size_t>(p.id)].v0; v1 = im.arcs[static_cast<std::size_t>(p.id)].v1; }
        else { v0 = im.pieces[static_cast<std::size_t>(p.id)].v0; v1 = im.pieces[static_cast<std::size_t>(p.id)].v1; }
    };
    for (std::size_t ri = 0; ri < im.regions.size(); ++ri) {
        const op::Region& reg = im.regions[ri];
        const op::BFace& face = im.solid[reg.solid]->faces[static_cast<std::size_t>(reg.face)];
        std::printf("  region %d: solid %d face %d surface %d%s, %s%s\n", static_cast<int>(ri), reg.solid, reg.face,
                    im.carrier_of[reg.solid][static_cast<std::size_t>(face.carrier)], face.flipped ? " (flipped)" : "",
                    kClass[reg.inside + 1], kept[ri] ? ", kept" : "");
        for (const auto& loop : reg.loops) {
            std::printf("    loop:");
            for (const op::Piece& p : loop) {
                int v0 = -1, v1 = -1;
                ends(p, v0, v1);
                std::printf(" %s%d%s(v%d>v%d)", p.fresh ? "n" : "p", p.id, p.forward ? "" : "'", p.forward ? v0 : v1,
                            p.forward ? v1 : v0);
                if (kept[ri]) {
                    auto& u = uses[{p.fresh, p.id}];
                    ++(p.forward ? u.first : u.second);
                }
            }
            std::printf("\n");
        }
    }
    for (const auto& [key, u] : uses)
        if (u.first != u.second) {
            int v0 = -1, v1 = -1;
            op::Piece p;
            p.fresh = key.first;
            p.id = key.second;
            ends(p, v0, v1);
            const Vec3d a0 = im.vertices[static_cast<std::size_t>(v0)].point, a1 = im.vertices[static_cast<std::size_t>(v1)].point;
            std::printf("  unbalanced: %s%d used %d forward, %d backward; v%d (%.6f %.6f %.6f) to v%d (%.6f %.6f %.6f)\n",
                        key.first ? "n" : "p", key.second, u.first, u.second, v0, a0.x, a0.y, a0.z, v1, a1.x, a1.y, a1.z);
        }
}

// ------------------------------------------------------------ curved chains

/// Chains with curved operands: cylinder sectors on shared cylinders (exact
/// circle parameterisation, dyadic ranges, so faces of different operands lie
/// on one cylinder or plane), boxes with a hill on top, and rotated boxes.
/// Starts from the quarter cylindrical shell.
void bench_curved_chains(int chains, int depth, int first_chain = 0, bool trace = false) {
    std::printf("## 11. Chained operations with curved operands\n\n");
    std::printf("%d chains of up to %d operations from the quarter cylindrical shell (radii 1 to 2, height 1). Each "
                "step uses, with equal probability, a cylinder sector on the same family of cylinders (radii, angle "
                "parameter and height ranges on dyadic grids, so faces of different operands often lie on one "
                "cylinder or one plane without being the same patch), a box with a hill on top at a random grid "
                "position, or a box of random orientation. Volumes are checked against a Monte Carlo count of "
                "200,000 points that uses only the analytic definition of the chain.\n\n",
                chains, depth);
    struct Operand {
        int kind = 0; // 0 sector, 1 hill box, 2 rotated box
        double r0 = 0, r1 = 0, t0 = 0, t1 = 0, z0 = 0, z1 = 0; // sector
        Vec3d o{}; double sx = 0, sy = 0, sz = 0, bump = 0;    // hill box
        RandomBox box;                                         // rotated box
        op::PatchSolid solid() const {
            if (kind == 0) return op::make_cylinder_sector(r0, r1, t0, t1, z0, z1);
            if (kind == 1) return op::make_bump_box(o, sx, sy, sz, bump);
            return box.solid();
        }
        bool inside(const Vec3d& p) const {
            if (kind == 0) {
                const double rho = std::sqrt(p.x * p.x + p.y * p.y), th = std::atan2(p.y, p.x);
                return rho >= r0 && rho <= r1 && p.z >= z0 && p.z <= z1 && th >= 2.0 * std::atan(t0) &&
                       th <= 2.0 * std::atan(t1);
            }
            if (kind == 1) {
                const double x = (p.x - o.x) / sx, y = (p.y - o.y) / sy;
                if (x < 0 || x > 1 || y < 0 || y > 1 || p.z < o.z) return false;
                return p.z <= o.z + sz + 16.0 * bump * x * (1 - x) * y * (1 - y);
            }
            return box.inside(p);
        }
    };
    auto pick = [](akt::Rng& rng, const std::vector<double>& v) { return v[static_cast<std::size_t>(rng.next() % v.size())]; };
    struct Row {
        int attempted = 0, solved = 0, audited = 0, empty = 0, beyond3 = 0, nonmanifold = 0, coincident = 0;
        double worst_z = 0.0, total_ms = 0.0, worst_ms = 0.0;
        long faces = 0;
    };
    std::vector<Row> rows(static_cast<std::size_t>(depth));
    std::map<std::string, int> reasons;
    std::vector<std::string> failures;
    int kinds[3] = {0, 0, 0};
    long new_vertices = 0, exact_roots = 0, exact_ends = 0; // over the solved steps
    int steps_with_exact_ends = 0;
    TessTally tess; // kept apart from section 9 so its hash stays comparable
    tess.name = "11. curved chains";
    // Covers every operand: sectors reach radius 2.5 and angles down to
    // 2 atan(-0.25) (negative y); hills reach x, y = 1.75 and z = 1.3.
    const Vec3d lo{-0.5, -1.5, -0.5}, hi{2.75, 2.75, 1.75};
    for (int chain = first_chain; chain < first_chain + chains; ++chain) {
        akt::Rng rng(11000 + static_cast<std::uint64_t>(chain));
        Operand start;
        start.kind = 0;
        start.r0 = 1; start.r1 = 2; start.t0 = 0; start.t1 = 1; start.z0 = 0; start.z1 = 1;
        op::Brep current = op::to_brep(start.solid());
        std::function<bool(const Vec3d&)> member = [start](const Vec3d& p) { return start.inside(p); };
        for (int step = 0; step < depth; ++step) {
            Row& row = rows[static_cast<std::size_t>(step)];
            Operand B;
            const double which = rng.unit();
            if (which < 1.0 / 3.0) {
                B.kind = 0;
                B.r0 = pick(rng, {0.75, 1.0, 1.25, 1.5});
                B.r1 = B.r0 + pick(rng, {0.25, 0.5, 1.0});
                B.t0 = pick(rng, {-0.25, 0.0, 0.25, 0.5});
                B.t1 = B.t0 + pick(rng, {0.25, 0.5});
                B.z0 = pick(rng, {-0.25, 0.0, 0.25, 0.5});
                B.z1 = B.z0 + pick(rng, {0.25, 0.5, 1.0});
            } else if (which < 2.0 / 3.0) {
                B.kind = 1;
                B.o = {pick(rng, {-0.25, 0.0, 0.25, 0.5, 0.75}), pick(rng, {-0.25, 0.0, 0.25, 0.5, 0.75}),
                       pick(rng, {-0.25, 0.0, 0.25})};
                B.sx = pick(rng, {0.5, 0.75, 1.0});
                B.sy = pick(rng, {0.5, 0.75, 1.0});
                B.sz = pick(rng, {0.5, 0.75});
                B.bump = pick(rng, {0.1, 0.2, 0.3});
            } else {
                B.kind = 2;
                B.box = random_box(rng, {0.3, 0.3, 0.2}, {1.6, 1.6, 0.8});
            }
            ++kinds[B.kind];
            const double pk = rng.unit();
            const op::Op o = pk < 0.45 ? op::Op::Difference : (pk < 0.8 ? op::Op::Union : op::Op::Intersection);
            const op::Brep Bb = op::to_brep(B.solid());
            ++row.attempted;
            const double t0 = now_ms();
            op::BooleanOptions options = with_exact();
            options.trace = trace;
            options.plane_lines = g_plane_lines;
            const op::Imprint im = op::imprint(current, Bb, options);
            const op::BooleanResult r = op::select(im, o);
            const double ms = now_ms() - t0;
            row.total_ms += ms;
            row.worst_ms = std::max(row.worst_ms, ms);
            const char* kind_name = B.kind == 0 ? "sector" : (B.kind == 1 ? "hill box" : "rotated box");
            if (trace) {
                std::printf("chain %d step %d: %s %s", chain, step + 1,
                            o == op::Op::Difference ? "minus" : (o == op::Op::Union ? "or" : "and"), kind_name);
                if (B.kind == 0) std::printf(" r %g-%g t %g-%g z %g-%g", B.r0, B.r1, B.t0, B.t1, B.z0, B.z1);
                if (B.kind == 1) std::printf(" o (%g %g %g) size %g %g %g bump %g", B.o.x, B.o.y, B.o.z, B.sx, B.sy, B.sz, B.bump);
                std::printf(": imprint %s, result %s\n", im.ok ? "ok" : "failed", r.ok ? (r.empty ? "empty" : "ok") : "refused");
            }
            if (!im.ok) {
                const std::size_t colon = im.error.find(": ");
                ++reasons[colon == std::string::npos ? im.error : im.error.substr(colon + 2)];
                failures.push_back("chain " + std::to_string(chain) + ", step " + std::to_string(step + 1) + " (" +
                                   kind_name + "): " + im.error);
                if (trace)
                    for (const std::string& line : im.trace) std::printf("%s\n", line.c_str());
                hash_double(-1.0);
                break;
            }
            ++row.solved;
            if (im.stats.curved_coincident_pairs > 0) ++row.coincident;
            new_vertices += im.stats.new_vertices;
            exact_roots += im.stats.exact_roots;
            exact_ends += im.stats.exact_ends;
            if (im.stats.exact_ends > 0) ++steps_with_exact_ends;
            if (!r.ok) {
                failures.push_back("chain " + std::to_string(chain) + ", step " + std::to_string(step + 1) + " (" +
                                   kind_name + "): result refused: " + r.error);
                if (trace) {
                    for (const std::string& line : im.trace) std::printf("%s\n", line.c_str());
                    dump_refused(im, r);
                }
                break;
            }
            ++row.audited;
            const auto prev = member;
            const Operand ob = B;
            if (o == op::Op::Difference) member = [=](const Vec3d& p) { return prev(p) && !ob.inside(p); };
            else if (o == op::Op::Union) member = [=](const Vec3d& p) { return prev(p) || ob.inside(p); };
            else member = [=](const Vec3d& p) { return prev(p) && ob.inside(p); };
            if (r.empty) { ++row.empty; hash_double(0.0); break; }
            if (!r.manifold) ++row.nonmanifold;
            const double vol = op::volume(r.solid);
            akt::Rng mc_rng(static_cast<std::uint64_t>(chain * 100 + step + 7));
            const int samples = 200000;
            int hits = 0;
            for (int i = 0; i < samples; ++i)
                if (member({mc_rng.range(lo.x, hi.x), mc_rng.range(lo.y, hi.y), mc_rng.range(lo.z, hi.z)})) ++hits;
            const double boxv = (hi.x - lo.x) * (hi.y - lo.y) * (hi.z - lo.z);
            const double p = static_cast<double>(hits) / samples;
            const double sigma = boxv * std::sqrt(std::max(p * (1.0 - p), 1.0 / samples) / samples);
            const double z = std::fabs(boxv * p - vol) / sigma;
            row.worst_z = std::max(row.worst_z, z);
            if (z > 3.0) {
                ++row.beyond3;
                failures.push_back("chain " + std::to_string(chain) + ", step " + std::to_string(step + 1) + " (" +
                                   kind_name + "): Monte Carlo " + std::to_string(boxv * p) + ", kernel " +
                                   std::to_string(vol) + " (" + std::to_string(z) + " sigma)");
            }
            row.faces += static_cast<long>(r.solid.faces.size());
            hash_double(vol);
            hash_double(static_cast<double>(r.solid.faces.size() * 1000 + r.solid.edges.size()));
            tally_tessellation(tess, r.solid, "chain " + std::to_string(chain) + ", step " + std::to_string(step + 1));
            current = r.solid;
        }
    }
    std::printf("| Step | Attempted | Solved | Result passes the audit | Empty result | Not manifold | Operands with "
                "coincident curved surfaces | Mean faces | Worst Monte Carlo deviation, sigma | Beyond 3 sigma | Mean "
                "time, ms | Worst time, ms |\n|---|---|---|---|---|---|---|---|---|---|---|---|\n");
    for (int step = 0; step < depth; ++step) {
        const Row& row = rows[static_cast<std::size_t>(step)];
        const int kept = std::max(1, row.audited - row.empty);
        std::printf("| %d | %d | %d | %d | %d | %d | %d | %.1f | %.2f | %d | %.0f | %.0f |\n", step + 1, row.attempted,
                    row.solved, row.audited, row.empty, row.nonmanifold, row.coincident,
                    static_cast<double>(row.faces) / kept, row.worst_z, row.beyond3,
                    row.total_ms / std::max(1, row.attempted), row.worst_ms);
    }
    std::printf("\nOperands drawn: %d sectors, %d hill boxes, %d rotated boxes.\n", kinds[0], kinds[1], kinds[2]);
    std::printf("\nExact points, over the solved steps: %ld of %ld new vertices proven exact; %ld roots proven to be an end "
                "vertex of their edge by its exact point, in %d steps.\n",
                exact_roots, new_vertices, exact_ends, steps_with_exact_ends);
    {
        const int ok = std::max(1, tess.meshes - tess.failed);
        std::printf("\nTessellation of every committed result (chord %.0e): meshes %d, not tessellated %d, closed %d "
                    "(not manifold %d), Euler characteristic matches %d, components match %d, mean triangles %.0f, "
                    "worst deviation / chord %.2f, worst volume error / chord %.2f.\n",
                    kTessChord, tess.meshes, tess.failed, tess.watertight, tess.nonmanifold, tess.euler_ok,
                    tess.components_ok, static_cast<double>(tess.triangles) / ok, tess.worst_dev / kTessChord,
                    tess.worst_vol / kTessChord);
        for (const std::string& f : tess.failures) std::printf("- tessellation, %s\n", f.c_str());
    }
    for (const auto& kv : reasons) std::printf("\nNot solved (%d): %s\n", kv.second, kv.first.c_str());
    std::printf("\n");
    for (const std::string& f : failures) std::printf("- %s\n", f.c_str());
    std::printf("\n");
}

// ------------------------------------------------------------------ scaling

/// A unit cube whose six sides are tiled into k x k flat patches (dyadic
/// restrictions, so neighbouring tiles share their sides exactly).
op::PatchSolid tiled(const op::PatchSolid& cube, int k) {
    op::PatchSolid s;
    // Tile corners by bilinear interpolation of the face corners. With
    // dyadic corners and k a power of two no rounding happens, so the tiles
    // of a face are exactly coplanar and neighbours share corners exactly.
    for (const auto& F : cube.faces)
        for (int i = 0; i < k; ++i)
            for (int j = 0; j < k; ++j) {
                geom::PolyPatchD T(1, 1);
                for (int a = 0; a < 2; ++a)
                    for (int b = 0; b < 2; ++b) {
                        const double u = double(i + a) / k, v = double(j + b) / k;
                        const auto& c00 = F.at(0, 0); const auto& c10 = F.at(1, 0);
                        const auto& c01 = F.at(0, 1); const auto& c11 = F.at(1, 1);
                        auto mix = [&](double x00, double x10, double x01, double x11) {
                            return (1 - u) * (1 - v) * x00 + u * (1 - v) * x10 + (1 - u) * v * x01 + u * v * x11;
                        };
                        T.at(a, b) = {mix(c00.x, c10.x, c01.x, c11.x), mix(c00.y, c10.y, c01.y, c11.y),
                                      mix(c00.z, c10.z, c01.z, c11.z), 1.0};
                    }
                s.faces.push_back(T);
            }
    std::string err;
    if (!s.build(err)) { std::printf("tiling %d: %s\n", k, err.c_str()); s.faces.clear(); }
    return s;
}

/// Time of one Boolean against the number of faces of the operand, with the
/// cut kept local: a small rotated box at a corner of a tiled cube, so that
/// the work that must be done stays the same and any growth is overhead.
void bench_scaling(int max_k) {
    std::printf("## 12. Scaling with the number of faces\n\n");
    std::printf("A unit cube whose sides are tiled into k x k flat patches, minus a small parallelepiped at one corner "
                "(dyadic edge vectors, so that tiles stay exactly flat). First with the small solid's six faces whole: "
                "the cut is the same at every k, so growth in time is overhead or the cut crossing more tiles. Then "
                "with the small solid tiled k x k too: both operands large.\n\n");
    // Edge vectors in units of 1/256: a skewed, right-handed parallelepiped.
    // The origin's coordinates are odd multiples of 1/1024, so no tile vertex
    // of it (k <= 64) lies on a plane x, y or z = 0 or 1: general position.
    const double q = 1.0 / 256.0, o = 1.0 / 1024.0;
    const op::PatchSolid small = op::make_box({861 * o, 821 * o, 877 * o}, {73 * q, 37 * q, 0}, {-37 * q, 73 * q, 9 * q},
                                              {5 * q, -7 * q, 77 * q});
    for (int series = 0; series < 2; ++series) {
        std::printf("%s\n\n", series == 0 ? "Small solid whole (6 faces):" : "Small solid tiled k x k as well:");
        std::printf("| k | Faces | Faces of the small solid | Solved | Volume | Surface pairs considered | Of those, intersected | "
                    "Operands, ms | Surface/surface, ms | Edge/surface, ms | Clipping, ms | Edge splitting, ms | "
                    "Regions, ms | Classification, ms | Total, ms |\n|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|\n");
        for (int k = 1; k <= max_k; k *= 2) {
            const op::Brep A = op::to_brep(tiled(op::make_box({0, 0, 0}, {1, 1, 1}), k));
            const op::Brep B = op::to_brep(series == 0 ? small : tiled(small, k));
            const double t0 = now_ms();
            op::BooleanOptions options = with_exact();
            options.plane_lines = g_plane_lines;
            const op::Imprint im = op::imprint(A, B, options);
            const op::BooleanResult r = op::select(im, op::Op::Difference);
            const double ms = now_ms() - t0;
            const double vol = r.ok && !r.empty ? op::volume(r.solid) : 0.0;
            std::printf("| %d | %zu | %zu | %s | %.12f | %d | %d |", k, A.faces.size(), B.faces.size(),
                        im.ok && r.ok ? "yes" : "no", vol, im.stats.carrier_pairs, im.stats.carrier_pairs_solved);
            for (int i = 0; i < 7; ++i) std::printf(" %.0f |", im.stats.ms[i]);
            std::printf(" %.0f |\n", ms);
            if (!im.ok) std::printf("|   not solved: %s |\n", im.error.c_str());
            if (!r.ok) std::printf("|   not selected: %s |\n", r.error.c_str());
            hash_double(vol);
        }
        std::printf("\n");
    }
}
bool g_identity_trace = false;
void bench_identities() {
    std::printf("## 10. Algebraic identities\n\n");
    std::printf("For random operands, seven identities on pairs and three on triples, each computed along both "
                "routes with fresh imprints. Two results count as the same solid when op::compare_solids finds "
                "the same vertices (certified enclosures), edges and faces, after removing seams and vertices "
                "where an edge merely continues; their volumes must also agree. An identity is not evaluated "
                "when a step of either route is not solved.\n\n");
    using Route = std::function<bool(op::Brep&)>;
    std::string last_error, last_full;
    auto apply = [&last_error, &last_full](const op::Brep& a, const op::Brep& b, op::Op o, op::Brep& out) {
        const op::BooleanResult r = op::boolean_op(a, b, o, with_exact());
        if (!r.ok) {
            const std::size_t colon = r.error.find(": ");
            last_error = colon == std::string::npos ? r.error : r.error.substr(colon + 2);
            last_full = r.error;
            return false;
        }
        out = r.solid;
        return true;
    };
    struct Family {
        const char* name;
        int pairs = 0, triples = 0;
        int evaluated = 0, same = 0, different = 0, not_solved = 0;
        double worst_dvol = 0.0, ms = 0.0;
        std::vector<std::string> failures;
        std::map<std::string, int> reasons;
    };
    std::uint64_t hash = 0xcbf29ce484222325ull;
    const op::Brep unit = op::to_brep(op::make_box({0, 0, 0}, {1, 1, 1}));
    const op::Brep hill = op::to_brep(op::make_bump_box({0, 0, 0}, 1, 1, 1, 0.3));
    const op::Brep shell = op::to_brep(op::make_annular_sector(1.0, 2.0, 0.0, 1.0));
    std::vector<Family> fams;
    const int kPairs = 16, kTriples = 6;
    for (int fi = 0; fi < 4; ++fi) {
        Family F;
        F.name = fi == 0 ? "Unit box, random boxes" : fi == 1 ? "Box with a hill, random boxes"
                 : fi == 2 ? "Cylindrical shell, random boxes" : "Grid-aligned boxes (contact)";
        const op::Brep* A0 = fi == 0 ? &unit : fi == 1 ? &hill : fi == 2 ? &shell : nullptr;
        akt::Rng rng(7100 + static_cast<std::uint64_t>(fi));
        const Vec3d clo = fi == 2 ? Vec3d{0.5, 0.5, 0.2} : (fi == 1 ? Vec3d{0.1, 0.1, 0.2} : Vec3d{0.1, 0.1, 0.1});
        const Vec3d chi = fi == 2 ? Vec3d{1.5, 1.5, 0.8} : (fi == 1 ? Vec3d{0.9, 0.9, 1.2} : Vec3d{0.9, 0.9, 0.9});
        auto grid_box = [&]() {
            auto cell = [&](int lo, int hi) { return lo + static_cast<int>(rng.next() % static_cast<std::uint64_t>(hi - lo + 1)); };
            const Vec3d c{0.25 * cell(-2, 3), 0.25 * cell(-2, 3), 0.25 * cell(-2, 3)};
            const Vec3d s{0.25 * cell(2, 5), 0.25 * cell(2, 5), 0.25 * cell(2, 5)};
            return op::to_brep(op::make_box(c, {c.x + s.x, c.y + s.y, c.z + s.z}));
        };
        auto next_box = [&]() { return fi == 3 ? grid_box() : op::to_brep(random_box(rng, clo, chi).solid()); };
        const double t0 = now_ms();
        for (int trial = 0; trial < kPairs; ++trial) {
            const op::Brep A = A0 ? *A0 : grid_box();
            const op::Brep B = next_box();
            const op::Brep C = next_box();
            const bool triple = trial < kTriples;
            ++F.pairs;
            if (triple) ++F.triples;
            auto run = [&](const char* name, const Route& left, const Route& right) {
                op::Brep x, y;
                if (!left(x) || !right(y)) { ++F.not_solved; ++F.reasons[last_error]; if (g_identity_trace) std::printf("  not evaluated: %s, %s, trial %d: %s\n", F.name, name, trial, last_full.c_str()); hash_into(hash, -1.0); return; }
                ++F.evaluated;
                const op::SolidComparison c = op::compare_solids(x, y);
                const double dv = std::fabs((x.faces.empty() ? 0.0 : op::volume(x)) - (y.faces.empty() ? 0.0 : op::volume(y)));
                F.worst_dvol = std::max(F.worst_dvol, dv);
                if (c.same && dv < 1e-9) ++F.same;
                else {
                    ++F.different;
                    F.failures.push_back(std::string(name) + ", trial " + std::to_string(trial) + ": " +
                                         (c.same ? "volumes differ" : c.reason));
                }
                hash_into(hash, c.same ? 1.0 : 0.0);
                hash_into(hash, static_cast<double>(c.vertices[0] * 1000 + c.face_groups[0]));
            };
            using O = op::Op;
            auto r1 = [&](const op::Brep& a, const op::Brep& b, O o) {
                return Route([&, o](op::Brep& out) { return apply(a, b, o, out); });
            };
            auto r2 = [&](const op::Brep& a, O o1, const op::Brep& b, O o2, const op::Brep& c) { // (a o1 b) o2 c
                return Route([&, o1, o2](op::Brep& out) { op::Brep t; return apply(a, b, o1, t) && apply(t, c, o2, out); });
            };
            auto r2r = [&](const op::Brep& a, O o1, const op::Brep& b, O o2, const op::Brep& c) { // a o1 (b o2 c)
                return Route([&, o1, o2](op::Brep& out) { op::Brep t; return apply(b, c, o2, t) && apply(a, t, o1, out); });
            };
            const Route just_a = [&](op::Brep& out) { out = A; return true; };
            const Route split_join = [&](op::Brep& out) {
                op::Brep s, t;
                return apply(A, B, O::Difference, s) && apply(A, B, O::Intersection, t) && apply(s, t, O::Union, out);
            };
            run("A and B = B and A", r1(A, B, O::Intersection), r1(B, A, O::Intersection));
            run("A or B = B or A", r1(A, B, O::Union), r1(B, A, O::Union));
            run("A - (A - B) = A and B", r2r(A, O::Difference, A, O::Difference, B), r1(A, B, O::Intersection));
            run("(A - B) or (A and B) = A", split_join, just_a);
            run("A or (A and B) = A", r2r(A, O::Union, A, O::Intersection, B), just_a);
            run("A and (A or B) = A", r2r(A, O::Intersection, A, O::Union, B), just_a);
            run("(A or B) - B = A - B", r2(A, O::Union, B, O::Difference, B), r1(A, B, O::Difference));
            if (triple) {
                run("(A or B) or C = A or (B or C)", r2(A, O::Union, B, O::Union, C), r2r(A, O::Union, B, O::Union, C));
                run("(A and B) and C = A and (B and C)", r2(A, O::Intersection, B, O::Intersection, C),
                    r2r(A, O::Intersection, B, O::Intersection, C));
                run("A - (B or C) = (A - B) - C", r2r(A, O::Difference, B, O::Union, C),
                    r2(A, O::Difference, B, O::Difference, C));
            }
        }
        F.ms = now_ms() - t0;
        fams.push_back(std::move(F));
    }
    std::printf("| Operands | Pairs | Triples | Identities evaluated | Same solid | Different | Not evaluated | "
                "Worst volume difference | Time, s |\n|---|---|---|---|---|---|---|---|---|\n");
    for (const Family& F : fams)
        std::printf("| %s | %d | %d | %d | %d | %d | %d | %.1e | %.0f |\n", F.name, F.pairs, F.triples, F.evaluated,
                    F.same, F.different, F.not_solved, F.worst_dvol, F.ms / 1000.0);
    std::printf("\n");
    for (const Family& F : fams)
        for (const std::string& f : F.failures) std::printf("- %s, %s\n", F.name, f.c_str());
    for (const Family& F : fams)
        for (const auto& kv : F.reasons) std::printf("- %s, not evaluated (%d): %s\n", F.name, kv.second, kv.first.c_str());
    std::printf("Identity hash: %016llx\n", static_cast<unsigned long long>(hash));
}

} // namespace

int main(int argc, char** argv) {
    // "ak_bench chains [count] [depth] [first] [trace]" runs only the chained-Boolean
    // section; "ak_bench reuse ..." the same with reused operands; "ak_bench scale [k]"
    // the scaling section up to k x k tiles. A last argument "nolines" switches off
    // the direct plane/plane certification (BooleanOptions::plane_lines) in the
    // chain and scaling sections, for comparison.
    if (argc > 1 && std::string(argv[argc - 1]) == "nolines") {
        g_plane_lines = false;
        --argc;
    }
    // Likewise "noexact" (before "nolines" if both): BooleanOptions::exact_points
    // off in every section.
    if (argc > 1 && std::string(argv[argc - 1]) == "noexact") {
        g_exact_points = false;
        --argc;
    }
    if (argc > 1 && std::string(argv[1]) == "scale") {
        bench_scaling(argc > 2 ? std::atoi(argv[2]) : 16);
        std::printf("Scaling result hash: %016llx\n", static_cast<unsigned long long>(g_hash));
        return 0;
    }
    if (argc > 1 && std::string(argv[1]) == "curved") {
        const int count = argc > 2 ? std::atoi(argv[2]) : 40, depth = argc > 3 ? std::atoi(argv[3]) : 6;
        const int first = argc > 4 ? std::atoi(argv[4]) : 0;
        g_dump_faces = argc > 5 && std::string(argv[5]) == "trace";
        bench_curved_chains(count, depth, first, g_dump_faces);
        std::printf("Curved chain result hash: %016llx\n", static_cast<unsigned long long>(g_hash));
        return 0;
    }
    if (argc > 1 && std::string(argv[1]) == "identities") {
        g_identity_trace = argc > 2 && std::string(argv[2]) == "trace";
        bench_identities();
        return 0;
    }
    if (argc > 1 && (std::string(argv[1]) == "chains" || std::string(argv[1]) == "reuse" || std::string(argv[1]) == "grid")) {
        // Optional: first chain index, and "trace" to print the stages of any failing step.
        const int count = argc > 2 ? std::atoi(argv[2]) : 40, depth = argc > 3 ? std::atoi(argv[3]) : 6;
        const int first = argc > 4 ? std::atoi(argv[4]) : 0;
        const double reuse = std::string(argv[1]) == "reuse" ? 0.4 : 0.0;
        g_profile = argc > 5 && std::string(argv[5]) == "profile";
        bench_boolean_chains(count, depth, first, argc > 5 && std::string(argv[5]) == "trace", reuse,
                             std::string(argv[1]) == "grid");
        std::printf("Chain result hash: %016llx\n\n", static_cast<unsigned long long>(g_hash));
        print_tessellation();
        return 0;
    }
    std::printf("# akernel benchmark output\n\n");
    bench_predicates();
    bench_ssi();
    bench_fill();
    std::printf("Result hash: %016llx\n\n", static_cast<unsigned long long>(g_hash));
    // The Boolean section has its own hash so the one above stays comparable
    // with earlier runs.
    g_hash = 0xcbf29ce484222325ull;
    bench_boolean();
    std::printf("Boolean result hash: %016llx\n\n", static_cast<unsigned long long>(g_hash));
    g_hash = 0xcbf29ce484222325ull;
    bench_boolean_chains(40, 6);
    std::printf("Chain result hash: %016llx\n\n", static_cast<unsigned long long>(g_hash));
    g_hash = 0xcbf29ce484222325ull;
    bench_boolean_chains(40, 6, 0, false, 0.4);
    std::printf("Reuse chain result hash: %016llx\n\n", static_cast<unsigned long long>(g_hash));
    g_hash = 0xcbf29ce484222325ull;
    bench_boolean_chains(40, 6, 0, false, 0.0, true);
    std::printf("Grid chain result hash: %016llx\n\n", static_cast<unsigned long long>(g_hash));
    print_tessellation();
    std::printf("\n");
    bench_identities();
    std::printf("\n");
    g_hash = 0xcbf29ce484222325ull;
    bench_curved_chains(40, 6);
    std::printf("Curved chain result hash: %016llx\n\n", static_cast<unsigned long long>(g_hash));
    g_hash = 0xcbf29ce484222325ull;
    bench_scaling(16);
    std::printf("Scaling result hash: %016llx\n", static_cast<unsigned long long>(g_hash));
    return 0;
}
