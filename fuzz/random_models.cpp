// akernel — random model generator with oracles (W10: test corpus and oracles).
//
// Without a foreign-format reader the public CAD corpora are out of reach, so
// the kernel makes its own: every case is a short feature chain grown from one
// 64-bit seed. A case starts from a random primitive and applies 1 to
// `--depth` Boolean operations, each with another random primitive placed so
// that it overlaps the current result. Primitives are boxes (on a dyadic grid,
// so faces of different operands often lie in one plane, or at a random
// orientation), boxes with a hill on top, exact cylinder sectors, and the
// kernel's blends: a box with one edge rounded or chamfered, a box with every
// edge rounded, and a ring with one circular edge rounded or chamfered.
//
// Every step is checked by oracles that do not trust the operation:
//   audit        the topology of every built solid and every Boolean result
//                is a closed, oriented boundary (op::audit, and the edge fans
//                where parts touch);
//   validity     the certified self-intersection check (op/validity.hpp):
//                Invalid on a solid the kernel built is a failure, Undecided
//                is counted;
//   closed form  each primitive's volume against its formula, through the
//                quadrature (op::volume) and, where covered, the certified
//                enclosure (op::certified_volume), which must contain it;
//   identities   from one imprint, vol(A or B) + vol(A and B) = vol A + vol B
//                and vol(A - B) = vol A - vol(A and B), in doubles to a
//                tolerance and, when every volume is certified, exactly (the
//                interval of the defect must contain 0);
//   membership   a Monte Carlo count of points inside the analytic definition
//                of the whole chain (no kernel code) against the volume of the
//                result.
// A refusal (an operation that answers "not solved") is not a failure: the
// kernel never guesses, and refusals are what this program counts, per
// operation and per reason, so that they can be driven to zero. An oracle
// violation or an exception is a failure, and makes the exit code 1.
//
// Every case is reproducible from its seed on every platform (SplitMix64 and
// the deterministic transcendentals of num/detmath.hpp, no <random>):
//
//   ak_random_models --seed 1234 --count 1 --verbose
//
// replays case 1234 alone and prints each step. See fuzz/README.md.
#include <akernel/num/detmath.hpp>
#include <akernel/op/boolean.hpp>
#include <akernel/op/certified_volume.hpp>
#include <akernel/op/fillet.hpp>
#include <akernel/op/surfaces.hpp>
#include <akernel/op/validity.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdarg>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <fstream>
#include <functional>
#include <initializer_list>
#include <map>
#include <numbers>
#include <optional>
#include <sstream>
#include <string>
#include <vector>

#include "check.hpp"

using namespace ak;
using geom::Vec3d;

namespace {

constexpr double kPi = std::numbers::pi;

// ------------------------------------------------------------------ options

struct Options {
    std::uint64_t seed = 1;     ///< first case
    int count = 50;             ///< number of cases
    int depth = 3;              ///< most Boolean steps per case
    double time_limit = 0.0;    ///< seconds; 0: no limit. No new case starts after it.
    int mc_samples = 20000;     ///< Monte Carlo points per checked result
    bool validity_results = true; ///< also run the self-intersection check on Boolean results
    bool verbose = false;
    std::string json;           ///< write the summary as JSON here
    std::string markdown;       ///< write the summary as Markdown here
    std::string progress;       ///< write the seed of the running case here (crash triage)
    std::vector<std::uint64_t> seeds; ///< explicit seeds (--seeds-file), instead of seed..seed+count
    bool seeds_file = false;          ///< --seeds-file given (an empty file runs no case)
};

// ------------------------------------------------------------------ primitives

enum class Kind { GridBox, RotatedBox, HillBox, Sector, EdgeFillet, EdgeChamfer, RoundedBox, RingBlend, Count };

const char* kind_name(Kind k) {
    switch (k) {
    case Kind::GridBox: return "grid box";
    case Kind::RotatedBox: return "rotated box";
    case Kind::HillBox: return "hill box";
    case Kind::Sector: return "cylinder sector";
    case Kind::EdgeFillet: return "edge fillet";
    case Kind::EdgeChamfer: return "edge chamfer";
    case Kind::RoundedBox: return "rounded box";
    case Kind::RingBlend: return "ring blend";
    default: return "?";
    }
}

/// Is (x, y) of the section [0, a] x [0, b] kept by the cut of the corner (a, b)?
inline bool section_kept(double x, double y, double a, double b, bool round, double dA, double dB) {
    if (round) {
        const double cx = a - dA, cy = b - dA;
        if (x <= cx || y <= cy) return true;
        return (x - cx) * (x - cx) + (y - cy) * (y - cy) <= dA * dA;
    }
    // Chamfer: dA on the side x = a (from y = b down), dB on y = b (from x = a back).
    return (a - x) / dB + (b - y) / dA >= 1.0;
}

/// Area removed from the corner and its centroid's distance from the corner along each axis.
inline void cut_moments(bool round, double dA, double dB, double& area, double& dx, double& dy) {
    if (round) {
        area = dA * dA * (1.0 - kPi / 4.0);
        dx = dy = dA * (5.0 / 6.0 - kPi / 4.0) / (1.0 - kPi / 4.0);
    } else {
        area = 0.5 * dA * dB;
        dx = dB / 3.0; // along x, the setback on y = b is dB
        dy = dA / 3.0;
    }
}

/// One random primitive: its parameters, how the kernel builds it, and its
/// analytic definition (membership, volume, bounding box) that the oracles
/// use instead of the kernel.
struct Prim {
    Kind kind = Kind::GridBox;
    // Boxes (all kinds but sectors and rings): o + x X + y Y + z Z,
    // x in [0, a], y in [0, b], z in [0, c], X, Y, Z orthonormal right-handed.
    Vec3d o{}, X{1, 0, 0}, Y{0, 1, 0}, Z{0, 0, 1};
    double a = 1, b = 1, c = 1;
    double bump = 0;   // hill box
    double r = 0;      // fillet / rounded box radius; chamfer setback dA
    double r2 = 0;     // chamfer setback dB
    // Sector and ring: r0 <= rho <= r1, z0 <= z <= z1.
    double r0 = 0, r1 = 0, z0 = 0, z1 = 0;
    double t0 = 0, t1 = 0; // sector: circle parameters
    int quarters = 4;      // ring: quarter turns from +x
    op::RingEdge edge = op::RingEdge::OuterTop;
    bool ring_round = true;

    std::string describe() const {
        char buf[640];
        switch (kind) {
        case Kind::Sector:
            std::snprintf(buf, sizeof buf, "%s r %.17g..%.17g, t %.17g..%.17g, z %.17g..%.17g", kind_name(kind), r0, r1, t0,
                          t1, z0, z1);
            break;
        case Kind::RingBlend:
            std::snprintf(buf, sizeof buf, "%s r %.17g..%.17g, z %.17g..%.17g, edge %d, %s %.17g/%.17g, %d quarters",
                          kind_name(kind), r0, r1, z0, z1, static_cast<int>(edge), ring_round ? "round" : "chamfer", r, r2,
                          quarters);
            break;
        default:
            std::snprintf(buf, sizeof buf,
                          "%s o (%.17g, %.17g, %.17g), X (%.6f, %.6f, %.6f), Y (%.6f, %.6f, %.6f), sides %.17g x %.17g x "
                          "%.17g, bump %.17g, r %.17g/%.17g",
                          kind_name(kind), o.x, o.y, o.z, X.x, X.y, X.z, Y.x, Y.y, Y.z, a, b, c, bump, r, r2);
        }
        return buf;
    }

    /// Build with the kernel. Empty with `error` set when it refuses.
    std::optional<op::PatchSolid> build(std::string& error) const {
        auto from = [&](const op::BlendedSolid& s) -> std::optional<op::PatchSolid> {
            if (!s.solid) error = s.error;
            return s.solid;
        };
        switch (kind) {
        case Kind::GridBox:
        case Kind::RotatedBox: return op::make_box(o, X * a, Y * b, Z * c);
        case Kind::HillBox: return op::make_bump_box(o, a, b, c, bump);
        case Kind::Sector: return op::make_cylinder_sector(r0, r1, t0, t1, z0, z1);
        case Kind::EdgeFillet: return from(op::fillet_box_edge(o, X * a, Y * b, Z * c, r));
        case Kind::EdgeChamfer: return from(op::chamfer_box_edge(o, X * a, Y * b, Z * c, r, r2));
        case Kind::RoundedBox: return from(op::rounded_box(o, a, b, c, r));
        case Kind::RingBlend:
            return from(op::blend_ring_edge(r0, r1, z0, z1, edge, ring_round ? op::EdgeCut::round(r) : op::EdgeCut::flat(r, r2),
                                            quarters));
        default: error = "unknown kind"; return std::nullopt;
        }
    }

    bool inside(const Vec3d& p) const {
        switch (kind) {
        case Kind::Sector:
        case Kind::RingBlend: {
            const double rho = std::sqrt(p.x * p.x + p.y * p.y);
            if (rho < r0 || rho > r1 || p.z < z0 || p.z > z1) return false;
            const double th = num::det::atan2(p.y, p.x);
            if (kind == Kind::Sector) return th >= 2.0 * num::det::atan(t0) && th <= 2.0 * num::det::atan(t1);
            if (quarters < 4) {
                const double t = th < 0 ? th + 2.0 * kPi : th;
                if (t > quarters * kPi / 2.0) return false;
            }
            const bool outer = edge == op::RingEdge::OuterTop || edge == op::RingEdge::OuterBottom;
            const bool top = edge == op::RingEdge::OuterTop || edge == op::RingEdge::InnerTop;
            const double x = outer ? rho - r0 : r1 - rho, y = top ? p.z - z0 : z1 - p.z;
            return section_kept(x, y, r1 - r0, z1 - z0, ring_round, r, r2);
        }
        default: break;
        }
        const Vec3d d = p - o;
        const double x = dot(d, X), y = dot(d, Y), z = dot(d, Z);
        if (x < 0 || x > a || y < 0 || y > b || z < 0) return false;
        switch (kind) {
        case Kind::HillBox: {
            const double u = x / a, v = y / b;
            return z <= c + 16.0 * bump * u * (1 - u) * v * (1 - v);
        }
        case Kind::EdgeFillet: return z <= c && section_kept(x, y, a, b, true, r, r);
        case Kind::EdgeChamfer: return z <= c && section_kept(x, y, a, b, false, r, r2);
        case Kind::RoundedBox: {
            if (z > c) return false;
            const double qx = std::clamp(x, r, a - r), qy = std::clamp(y, r, b - r), qz = std::clamp(z, r, c - r);
            const double ex = x - qx, ey = y - qy, ez = z - qz;
            return ex * ex + ey * ey + ez * ez <= r * r;
        }
        default: return z <= c;
        }
    }

    double volume() const {
        switch (kind) {
        case Kind::HillBox: return a * b * (c + 4.0 * bump / 9.0);
        case Kind::Sector: {
            const double th = 2.0 * (num::det::atan(t1) - num::det::atan(t0));
            return 0.5 * (r1 * r1 - r0 * r0) * th * (z1 - z0);
        }
        case Kind::EdgeFillet: return a * b * c - c * r * r * (1.0 - kPi / 4.0);
        case Kind::EdgeChamfer: return a * b * c - 0.5 * c * r * r2;
        case Kind::RoundedBox: {
            const double A = a - 2 * r, B = b - 2 * r, C = c - 2 * r;
            return A * B * C + 2.0 * r * (A * B + B * C + A * C) + kPi * r * r * (A + B + C) + 4.0 / 3.0 * kPi * r * r * r;
        }
        case Kind::RingBlend: {
            // Pappus: angle times the first moment in rho of the section.
            const double th = quarters * kPi / 2.0, h = z1 - z0;
            const bool outer = edge == op::RingEdge::OuterTop || edge == op::RingEdge::OuterBottom;
            double area = 0, dx = 0, dy = 0;
            cut_moments(ring_round, r, r2, area, dx, dy);
            // The cut corner sits at x = w: rho = r1 for an outer edge, r0 for an inner one.
            const double rho_c = outer ? r1 - dx : r0 + dx;
            return th * (0.5 * h * (r1 * r1 - r0 * r0) - area * rho_c);
        }
        default: return a * b * c;
        }
    }

    void bounds(Vec3d& lo, Vec3d& hi) const {
        if (kind == Kind::Sector || kind == Kind::RingBlend) {
            lo = {-r1, -r1, z0};
            hi = {r1, r1, z1};
            return;
        }
        const double top = kind == Kind::HillBox ? c + bump : c;
        lo = {1e300, 1e300, 1e300};
        hi = {-1e300, -1e300, -1e300};
        for (int i = 0; i < 8; ++i) {
            const Vec3d p = o + X * ((i & 1) ? a : 0.0) + Y * ((i & 2) ? b : 0.0) + Z * ((i & 4) ? top : 0.0);
            lo = {std::min(lo.x, p.x), std::min(lo.y, p.y), std::min(lo.z, p.z)};
            hi = {std::max(hi.x, p.x), std::max(hi.y, p.y), std::max(hi.z, p.z)};
        }
    }
};

/// A value from a list, uniformly.
double pick(akt::Rng& rng, std::initializer_list<double> v) {
    return *(v.begin() + static_cast<std::ptrdiff_t>(rng.next() % v.size()));
}

/// A random rotation (axis and angle), as the columns X, Y, Z.
void random_frame(akt::Rng& rng, Vec3d& X, Vec3d& Y, Vec3d& Z) {
    Vec3d axis;
    double len = 0.0;
    do {
        axis = {rng.range(-1, 1), rng.range(-1, 1), rng.range(-1, 1)};
        len = geom::norm(axis);
    } while (len < 0.2 || len > 1.0);
    axis = axis / len;
    const double ang = rng.range(0.2, 2.9), c = num::det::cos(ang), s = num::det::sin(ang);
    const double ax = axis.x, ay = axis.y, az = axis.z;
    X = {c + ax * ax * (1 - c), ay * ax * (1 - c) + az * s, az * ax * (1 - c) - ay * s};
    Y = {ax * ay * (1 - c) - az * s, c + ay * ay * (1 - c), az * ay * (1 - c) + ax * s};
    Z = cross(X, Y); // right-handed to rounding
}

/// A random primitive whose box overlaps the region [lo, hi].
Prim random_prim(akt::Rng& rng, const Vec3d& lo, const Vec3d& hi) {
    Prim p;
    const double which = rng.unit();
    // Weights: boxes and the curved primitives most, blends often enough that
    // every kind shows in a nightly run of a few hundred cases.
    if (which < 0.18) p.kind = Kind::GridBox;
    else if (which < 0.32) p.kind = Kind::RotatedBox;
    else if (which < 0.46) p.kind = Kind::HillBox;
    else if (which < 0.60) p.kind = Kind::Sector;
    else if (which < 0.72) p.kind = Kind::EdgeFillet;
    else if (which < 0.80) p.kind = Kind::EdgeChamfer;
    else if (which < 0.90) p.kind = Kind::RoundedBox;
    else p.kind = Kind::RingBlend;

    // A dyadic point of the grid of step 1/4 near the region, so faces of
    // different operands often lie in one plane.
    auto grid = [&](double l, double h) { return std::floor(rng.range(l - 0.5, h) * 4.0) / 4.0; };
    const Vec3d g{grid(lo.x, hi.x), grid(lo.y, hi.y), grid(lo.z, hi.z)};

    switch (p.kind) {
    case Kind::GridBox:
        p.o = g;
        p.a = pick(rng, {0.25, 0.5, 0.75, 1.0});
        p.b = pick(rng, {0.25, 0.5, 0.75, 1.0});
        p.c = pick(rng, {0.25, 0.5, 0.75, 1.0});
        break;
    case Kind::RotatedBox: {
        random_frame(rng, p.X, p.Y, p.Z);
        p.a = rng.range(0.4, 1.1);
        p.b = rng.range(0.4, 1.1);
        p.c = rng.range(0.4, 1.1);
        const Vec3d centre{rng.range(lo.x, hi.x), rng.range(lo.y, hi.y), rng.range(lo.z, hi.z)};
        p.o = centre - (p.X * p.a + p.Y * p.b + p.Z * p.c) * 0.5;
        break;
    }
    case Kind::HillBox:
        p.o = g;
        p.a = pick(rng, {0.5, 0.75, 1.0});
        p.b = pick(rng, {0.5, 0.75, 1.0});
        p.c = pick(rng, {0.5, 0.75});
        p.bump = pick(rng, {0.1, 0.2, 0.3});
        break;
    case Kind::Sector:
        p.r0 = pick(rng, {0.5, 0.75, 1.0, 1.25, 1.5});
        p.r1 = p.r0 + pick(rng, {0.25, 0.5, 1.0});
        p.t0 = pick(rng, {-0.5, -0.25, 0.0, 0.25, 0.5});
        p.t1 = p.t0 + pick(rng, {0.25, 0.5});
        p.z0 = grid(lo.z, hi.z);
        p.z1 = p.z0 + pick(rng, {0.25, 0.5, 1.0});
        break;
    case Kind::EdgeFillet:
    case Kind::EdgeChamfer:
        p.a = pick(rng, {0.5, 0.75, 1.0});
        p.b = pick(rng, {0.5, 0.75, 1.0});
        p.c = pick(rng, {0.5, 0.75, 1.0});
        if (p.kind == Kind::EdgeFillet) {
            p.r = std::min(p.a, p.b) * pick(rng, {0.25, 0.5, 0.625});
        } else {
            p.r = p.b * pick(rng, {0.25, 0.5});
            p.r2 = p.a * pick(rng, {0.25, 0.5});
        }
        if (rng.unit() < 0.5) {
            p.o = g;
        } else {
            random_frame(rng, p.X, p.Y, p.Z);
            const Vec3d centre{rng.range(lo.x, hi.x), rng.range(lo.y, hi.y), rng.range(lo.z, hi.z)};
            p.o = centre - (p.X * p.a + p.Y * p.b + p.Z * p.c) * 0.5;
        }
        break;
    case Kind::RoundedBox:
        p.o = g;
        p.a = pick(rng, {0.5, 0.75, 1.0});
        p.b = pick(rng, {0.5, 0.75, 1.0});
        p.c = pick(rng, {0.5, 0.75, 1.0});
        p.r = pick(rng, {0.0625, 0.125, 0.1875});
        break;
    case Kind::RingBlend: {
        p.r0 = pick(rng, {0.5, 0.75, 1.0});
        p.r1 = p.r0 + pick(rng, {0.5, 0.75, 1.0});
        p.z0 = grid(lo.z, hi.z);
        p.z1 = p.z0 + pick(rng, {0.5, 0.75});
        p.edge = static_cast<op::RingEdge>(rng.next() % 4);
        p.ring_round = rng.unit() < 0.6;
        const double m = std::min(p.r1 - p.r0, p.z1 - p.z0);
        p.r = m * pick(rng, {0.25, 0.5});
        p.r2 = p.ring_round ? p.r : m * pick(rng, {0.25, 0.5});
        p.quarters = static_cast<int>(1 + rng.next() % 4);
        break;
    }
    default: break;
    }
    return p;
}

// ------------------------------------------------------------------ tallies

/// Counts for one operation (a primitive kind's construction, or one Boolean op).
struct Tally {
    int attempts = 0, ok = 0, refused = 0, empty = 0, failed = 0;
    int valid = 0, invalid = 0, undecided = 0; // self-intersection check
    int certified = 0;                         // certified volume available
    double ms = 0.0, worst_ms = 0.0;
    std::uint64_t worst_seed = 0;
    std::map<std::string, int> reasons;        // refusal reasons
    void time(double t, std::uint64_t seed) {
        ms += t;
        if (t > worst_ms) { worst_ms = t; worst_seed = seed; }
    }
};

struct Failure {
    std::uint64_t seed = 0;
    std::string op, oracle, detail;
};

struct Run {
    Options opt;
    std::map<std::string, Tally> ops; // keyed "build: <kind>", "boolean: union", ...
    std::vector<Failure> failures;
    int cases = 0, steps = 0, chains_finished = 0;
    double worst_identity = 0.0, worst_mc_sigma = 0.0, worst_closed_form = 0.0;
    std::uint64_t seed = 0; // the running case
    std::string log;        // verbose lines of the running case

    void say(const char* fmt, ...) {
        if (!opt.verbose) return;
        char buf[1024];
        va_list ap;
        va_start(ap, fmt);
        std::vsnprintf(buf, sizeof buf, fmt, ap);
        va_end(ap);
        std::fputs(buf, stdout);
    }
    void fail(const std::string& what, const std::string& oracle, const std::string& detail) {
        failures.push_back({seed, what, oracle, detail});
        ++ops[what].failed;
        std::printf("FAIL seed %llu: %s: %s: %s\n", static_cast<unsigned long long>(seed), what.c_str(), oracle.c_str(),
                    detail.c_str());
        std::fflush(stdout);
    }
};

/// A refusal reason with the numbers taken out, so that similar refusals share a bucket.
std::string bucket(const std::string& why) {
    std::string s;
    bool num = false;
    for (char ch : why) {
        if (ch >= '0' && ch <= '9') {
            if (!num) s += '#';
            num = true;
        } else if (num && (ch == '.' || ch == 'e' || ch == '-' || ch == '+')) {
            continue; // still inside the number
        } else {
            s += ch;
            num = false;
        }
    }
    if (s.size() > 140) s.resize(140);
    return s.empty() ? "(no reason given)" : s;
}

double now_ms() {
    using clock = std::chrono::steady_clock;
    return std::chrono::duration<double, std::milli>(clock::now().time_since_epoch()).count();
}

const char* op_name(op::Op o) {
    return o == op::Op::Union ? "union" : o == op::Op::Intersection ? "intersection" : "difference";
}

/// The closed-boundary oracle on a solid the kernel produced.
bool audited(const op::Brep& b, const topo::AuditReport& a, std::string& why) {
    if (!a.closed()) {
        char buf[200];
        std::snprintf(buf, sizeof buf, "not closed: %d boundary, %d misoriented, %d unused edges, %d open fans, loops %s",
                      a.boundary_edges, a.misoriented_edges, a.unused_edges, a.open_fans, a.loops_closed ? "closed" : "open");
        why = buf;
        return false;
    }
    if (!a.manifold() && !op::check_edge_fans(b, why)) return false;
    return true;
}

void validity(Run& run, const std::string& name, const op::Brep& b) {
    Tally& t = run.ops[name];
    const op::Validity v = op::check_self_intersection(b);
    if (v.state == op::Validity::Valid) ++t.valid;
    else if (v.state == op::Validity::Undecided) ++t.undecided;
    else {
        ++t.invalid;
        char buf[300];
        std::snprintf(buf, sizeof buf, "faces %d and %d cross near (%.4f, %.4f, %.4f): %s", v.face_a, v.face_b, v.where.x,
                      v.where.y, v.where.z, v.why.c_str());
        run.fail(name, "self-intersection", buf);
    }
}

/// Build one primitive and check it. Empty when the kernel refused it or an oracle failed.
std::optional<op::Brep> build(Run& run, const Prim& p) {
    const std::string name = std::string("build: ") + kind_name(p.kind);
    Tally& t = run.ops[name];
    ++t.attempts;
    run.say("  build %s\n", p.describe().c_str());
    const double t0 = now_ms();
    std::string error;
    const auto solid = p.build(error);
    if (!solid) {
        t.time(now_ms() - t0, run.seed);
        ++t.refused;
        ++t.reasons[bucket(error)];
        run.say("    refused: %s\n", error.c_str());
        return std::nullopt;
    }
    op::Brep b = op::to_brep(*solid);
    t.time(now_ms() - t0, run.seed);
    std::string why;
    if (!audited(b, op::audit(b), why)) {
        run.fail(name, "audit", why);
        return std::nullopt;
    }
    validity(run, name, b);
    // Closed-form volume: through the quadrature, and through the certified
    // enclosure where it covers the solid.
    const double exact = p.volume(), q = op::volume(b);
    const double rel = std::fabs(q - exact) / std::max(1e-12, exact);
    run.worst_closed_form = std::max(run.worst_closed_form, rel);
    if (rel > 1e-8) {
        char buf[160];
        std::snprintf(buf, sizeof buf, "quadrature %.15g, closed form %.15g (relative %.1e)", q, exact, rel);
        run.fail(name, "closed-form volume", buf);
        return std::nullopt;
    }
    if (const auto cv = op::certified_volume(b)) {
        ++t.certified;
        const double slack = 1e-13 * std::max(1.0, exact); // the closed form is a double
        if (!(cv->lo <= exact + slack && exact - slack <= cv->hi)) {
            char buf[200];
            std::snprintf(buf, sizeof buf, "certified [%.17g, %.17g] misses the closed form %.17g", cv->lo, cv->hi, exact);
            run.fail(name, "certified volume", buf);
            return std::nullopt;
        }
    }
    ++t.ok;
    return b;
}

/// The analytic chain: membership of a point, from the primitives alone.
struct Chain {
    std::vector<Prim> prims;
    std::vector<op::Op> ops; // ops[i] combines prims[i + 1]
    bool inside(const Vec3d& p) const {
        bool in = prims[0].inside(p);
        for (std::size_t i = 0; i < ops.size(); ++i) {
            const bool b = prims[i + 1].inside(p);
            in = ops[i] == op::Op::Union ? (in || b) : ops[i] == op::Op::Intersection ? (in && b) : (in && !b);
        }
        return in;
    }
    void bounds(Vec3d& lo, Vec3d& hi) const {
        prims[0].bounds(lo, hi);
        for (std::size_t i = 1; i < prims.size(); ++i) {
            Vec3d l, h;
            prims[i].bounds(l, h);
            // Union widens; the other ops can only shrink, but widening is safe.
            lo = {std::min(lo.x, l.x), std::min(lo.y, l.y), std::min(lo.z, l.z)};
            hi = {std::max(hi.x, h.x), std::max(hi.y, h.y), std::max(hi.z, h.z)};
        }
    }
};

/// One case: a chain of Boolean operations grown from `seed`.
void run_case(Run& run, std::uint64_t seed) {
    run.seed = seed;
    ++run.cases;
    akt::Rng rng(seed * 0x9E3779B97F4A7C15ull + 0x5DEECE66Dull);
    run.say("case %llu\n", static_cast<unsigned long long>(seed));
    Chain chain;
    chain.prims.push_back(random_prim(rng, {0, 0, 0}, {1, 1, 1}));
    auto first = build(run, chain.prims[0]);
    if (!first) return;
    op::Brep current = *first;
    double v_current = op::volume(current);
    std::optional<num::Interval> cv_current = op::certified_volume(current);
    const int depth = 1 + static_cast<int>(rng.next() % static_cast<std::uint64_t>(std::max(1, run.opt.depth)));

    for (int step = 0; step < depth; ++step) {
        ++run.steps;
        Vec3d lo, hi;
        chain.bounds(lo, hi);
        // Place the next operand around the current result's box, not beyond it.
        const Prim B = random_prim(rng, lo, hi);
        const double which = rng.unit();
        const op::Op chosen = which < 0.4 ? op::Op::Union : which < 0.6 ? op::Op::Intersection : op::Op::Difference;
        auto bb = build(run, B);
        if (!bb) return;

        // One imprint, all three operations: the identities need them all.
        const double t0 = now_ms();
        const op::Imprint im = op::imprint(current, *bb);
        const double t_imprint = now_ms() - t0;
        {
            Tally& t = run.ops["boolean: imprint"];
            ++t.attempts;
            t.time(t_imprint, seed);
            if (!im.ok) {
                ++t.refused;
                ++t.reasons[bucket(im.error)];
                run.say("  step %d: imprint refused: %s\n", step + 1, im.error.c_str());
                return;
            }
            ++t.ok;
        }
        const op::Op all[3] = {op::Op::Union, op::Op::Intersection, op::Op::Difference};
        op::BooleanResult res[3];
        double vol[3] = {0, 0, 0};
        std::optional<num::Interval> cvol[3];
        bool solved = true;
        for (int k = 0; k < 3; ++k) {
            const std::string name = std::string("boolean: ") + op_name(all[k]);
            Tally& t = run.ops[name];
            ++t.attempts;
            const double s0 = now_ms();
            res[k] = op::select(im, all[k]);
            t.time(now_ms() - s0, seed);
            if (!res[k].ok) {
                ++t.refused;
                ++t.reasons[bucket(res[k].error)];
                run.say("  step %d %s: refused: %s\n", step + 1, op_name(all[k]), res[k].error.c_str());
                solved = false;
                continue;
            }
            if (res[k].empty) {
                ++t.empty;
                ++t.ok;
                continue;
            }
            std::string why;
            if (!audited(res[k].solid, res[k].audit, why)) {
                run.fail(name, "audit", why);
                solved = false;
                continue;
            }
            ++t.ok;
            vol[k] = op::volume(res[k].solid);
            cvol[k] = op::certified_volume(res[k].solid);
            if (cvol[k]) ++t.certified;
            if (run.opt.validity_results) validity(run, name, res[k].solid);
        }
        const double vb = op::volume(*bb);
        run.say("  step %d: vol A %.9f, vol B %.9f, union %.9f, intersection %.9f, difference %.9f (%.0f ms)\n", step + 1,
                v_current, vb, vol[0], vol[1], vol[2], t_imprint);

        if (solved) {
            // Identities among the volumes, in doubles.
            const double scale = std::max(1.0, v_current + vb);
            const double d1 = std::fabs(vol[0] + vol[1] - v_current - vb) / scale;
            const double d2 = std::fabs(v_current - vol[1] - vol[2]) / scale;
            run.worst_identity = std::max({run.worst_identity, d1, d2});
            if (d1 > 1e-7 || d2 > 1e-7) {
                char buf[240];
                std::snprintf(buf, sizeof buf,
                              "vol A %.12g, vol B %.12g, union %.12g, intersection %.12g, difference %.12g: defects %.1e, %.1e",
                              v_current, vb, vol[0], vol[1], vol[2], d1, d2);
                run.fail("boolean: identities", "volume identity", buf);
            }
            // And exactly, when every volume is certified.
            const auto cb = op::certified_volume(*bb);
            if (cv_current && cb && cvol[0] && cvol[1] && cvol[2]) {
                const num::Interval e1 = *cvol[0] + *cvol[1] - *cv_current - *cb;
                const num::Interval e2 = *cv_current - *cvol[1] - *cvol[2];
                if (!(e1.lo <= 0.0 && 0.0 <= e1.hi) || !(e2.lo <= 0.0 && 0.0 <= e2.hi)) {
                    char buf[200];
                    std::snprintf(buf, sizeof buf, "certified defects [%.3e, %.3e] and [%.3e, %.3e] exclude 0", e1.lo, e1.hi,
                                  e2.lo, e2.hi);
                    run.fail("boolean: identities", "certified volume identity", buf);
                }
            }
        }

        const int k = chosen == op::Op::Union ? 0 : chosen == op::Op::Intersection ? 1 : 2;
        if (!res[k].ok) return;
        chain.prims.push_back(B);
        chain.ops.push_back(chosen);
        if (res[k].empty) {
            run.say("  step %d: %s is empty; chain ends\n", step + 1, op_name(chosen));
            ++run.chains_finished;
            return;
        }

        // Membership: Monte Carlo over the analytic chain, against the volume.
        // The deviation is measured in standard errors of the count the
        // Brep's volume predicts (hits ~ Binomial(n, vol / box)), with a floor
        // of one expected hit so that tiny results are not judged on zero
        // hits. A count beyond 4 is drawn again with ten times the points from
        // another stream, and only a second count beyond 6 is a failure: a
        // false alarm needs two independent unlikely draws.
        if (run.opt.mc_samples > 0) {
            Vec3d l, h;
            chain.bounds(l, h);
            const double box = (h.x - l.x) * (h.y - l.y) * (h.z - l.z);
            const double p0 = std::clamp(vol[k] / box, 0.0, 1.0);
            auto count = [&](int n, std::uint64_t stream, double& est, double& z) {
                akt::Rng mc(seed * 0xD1B54A32D192ED03ull + stream);
                int hits = 0;
                for (int i = 0; i < n; ++i)
                    hits += chain.inside({mc.range(l.x, h.x), mc.range(l.y, h.y), mc.range(l.z, h.z)}) ? 1 : 0;
                const double expect = n * p0, sd = std::sqrt(std::max(expect * (1.0 - p0), 1.0));
                est = box * hits / n;
                z = std::fabs(hits - expect) / sd;
            };
            double est = 0, z = 0;
            const std::uint64_t stream = 2 * static_cast<std::uint64_t>(step);
            count(run.opt.mc_samples, stream, est, z);
            run.say("  step %d: %s kept, Monte Carlo %.6f (%.1f sigma)\n", step + 1, op_name(chosen), est, z);
            if (z > 4.0) {
                count(10 * run.opt.mc_samples, stream + 1, est, z);
                run.say("  step %d: drawn again: Monte Carlo %.6f (%.1f sigma)\n", step + 1, est, z);
            }
            run.worst_mc_sigma = std::max(run.worst_mc_sigma, z);
            if (z > 6.0) {
                char buf[200];
                std::snprintf(buf, sizeof buf, "Brep volume %.9f, Monte Carlo %.9f (%.1f sigma) after step %d", vol[k], est, z,
                              step + 1);
                run.fail(std::string("boolean: ") + op_name(chosen), "membership", buf);
                return;
            }
        }
        current = res[k].solid;
        v_current = vol[k];
        cv_current = cvol[k];
    }
    ++run.chains_finished;
}

// ------------------------------------------------------------------ reports

std::string json_escape(const std::string& s) {
    std::string o;
    for (char ch : s) {
        if (ch == '"' || ch == '\\') { o += '\\'; o += ch; }
        else if (ch == '\n') o += "\\n";
        else if (static_cast<unsigned char>(ch) < 0x20) o += ' ';
        else o += ch;
    }
    return o;
}

std::string markdown(const Run& run, double seconds) {
    std::ostringstream out;
    char buf[512];
    std::snprintf(buf, sizeof buf,
                  "## Random models\n\n%d cases (seeds %llu and on), %d Boolean steps, %d chains run to the end, %.0f s.\n\n",
                  run.cases, static_cast<unsigned long long>(run.opt.seeds.empty() ? run.opt.seed : run.opt.seeds.front()),
                  run.steps, run.chains_finished, seconds);
    out << buf;
    out << "| operation | attempts | ok | refused | empty | **failed** | valid | invalid | undecided | certified vol. | mean ms | "
           "worst ms (seed) |\n|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|\n";
    for (const auto& [name, t] : run.ops) {
        std::snprintf(buf, sizeof buf, "| %s | %d | %d | %d | %d | %s | %d | %d | %d | %d | %.1f | %.0f (%llu) |\n", name.c_str(),
                      t.attempts, t.ok, t.refused, t.empty, t.failed ? ("**" + std::to_string(t.failed) + "**").c_str() : "0",
                      t.valid, t.invalid, t.undecided, t.certified, t.attempts ? t.ms / t.attempts : 0.0, t.worst_ms,
                      static_cast<unsigned long long>(t.worst_seed));
        out << buf;
    }
    std::snprintf(buf, sizeof buf,
                  "\nWorst relative volume identity defect %.1e; worst Monte Carlo deviation %.1f sigma; worst closed-form "
                  "volume error of a primitive %.1e.\n",
                  run.worst_identity, run.worst_mc_sigma, run.worst_closed_form);
    out << buf;
    bool any = false;
    for (const auto& [name, t] : run.ops)
        for (const auto& [why, n] : t.reasons) {
            if (!any) out << "\n### Refusals by reason\n\n| operation | count | reason |\n|---|--:|---|\n";
            any = true;
            out << "| " << name << " | " << n << " | " << why << " |\n";
        }
    if (run.failures.empty()) out << "\nNo oracle failures.\n";
    else {
        out << "\n### Failures\n\n| seed | operation | oracle | detail |\n|--:|---|---|---|\n";
        for (const Failure& f : run.failures)
            out << "| " << f.seed << " | " << f.op << " | " << f.oracle << " | " << f.detail << " |\n";
        out << "\nReplay one with `ak_random_models --seed <seed> --count 1 --verbose`.\n";
    }
    return out.str();
}

std::string json(const Run& run, double seconds) {
    std::ostringstream o;
    o << "{\n  \"cases\": " << run.cases << ",\n  \"steps\": " << run.steps << ",\n  \"seconds\": " << seconds
      << ",\n  \"worst_identity_defect\": " << run.worst_identity << ",\n  \"worst_mc_sigma\": " << run.worst_mc_sigma
      << ",\n  \"operations\": {";
    bool first = true;
    for (const auto& [name, t] : run.ops) {
        o << (first ? "\n" : ",\n") << "    \"" << json_escape(name) << "\": {\"attempts\": " << t.attempts << ", \"ok\": " << t.ok
          << ", \"refused\": " << t.refused << ", \"empty\": " << t.empty << ", \"failed\": " << t.failed
          << ", \"valid\": " << t.valid << ", \"invalid\": " << t.invalid << ", \"undecided\": " << t.undecided
          << ", \"certified_volume\": " << t.certified << ", \"ms\": " << t.ms << ", \"worst_ms\": " << t.worst_ms
          << ", \"worst_seed\": " << t.worst_seed << ", \"reasons\": {";
        bool fr = true;
        for (const auto& [why, n] : t.reasons) {
            o << (fr ? "" : ", ") << "\"" << json_escape(why) << "\": " << n;
            fr = false;
        }
        o << "}}";
        first = false;
    }
    o << "\n  },\n  \"failures\": [";
    first = true;
    for (const Failure& f : run.failures) {
        o << (first ? "\n" : ",\n") << "    {\"seed\": " << f.seed << ", \"operation\": \"" << json_escape(f.op)
          << "\", \"oracle\": \"" << json_escape(f.oracle) << "\", \"detail\": \"" << json_escape(f.detail) << "\"}";
        first = false;
    }
    o << "\n  ]\n}\n";
    return o.str();
}

void usage() {
    std::printf(
        "ak_random_models: random feature models checked by oracles (fuzz/README.md)\n\n"
        "  --seed N            first seed (default 1)\n"
        "  --count N           number of cases (default 50)\n"
        "  --seeds-file PATH   run the seeds listed in a file (one per line, # comments) instead\n"
        "  --depth N           most Boolean steps per case (default 3)\n"
        "  --time-limit S      start no case after S seconds\n"
        "  --mc-samples N      Monte Carlo points per kept result (default 20000; 0 = off)\n"
        "  --no-validity-results  skip the self-intersection check on Boolean results\n"
        "  --verbose           print every step\n"
        "  --json PATH         write the summary as JSON\n"
        "  --markdown PATH     write the summary as Markdown\n"
        "  --progress PATH     write the running seed here before each case (to find a crash)\n\n"
        "Exit code 1 when any oracle fails or an exception escapes; refusals are counted, not failures.\n");
}

bool parse(int argc, char** argv, Options& o) {
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        auto value = [&]() -> const char* {
            if (i + 1 >= argc) {
                std::printf("%s needs a value\n", a.c_str());
                std::exit(2);
            }
            return argv[++i];
        };
        if (a == "--seed") o.seed = std::strtoull(value(), nullptr, 10);
        else if (a == "--count") o.count = std::atoi(value());
        else if (a == "--depth") o.depth = std::atoi(value());
        else if (a == "--time-limit") o.time_limit = std::atof(value());
        else if (a == "--mc-samples") o.mc_samples = std::atoi(value());
        else if (a == "--no-validity-results") o.validity_results = false;
        else if (a == "--verbose") o.verbose = true;
        else if (a == "--json") o.json = value();
        else if (a == "--markdown") o.markdown = value();
        else if (a == "--progress") o.progress = value();
        else if (a == "--seeds-file") {
            o.seeds_file = true;
            std::ifstream in(value());
            if (!in) {
                std::printf("cannot read the seeds file\n");
                return false;
            }
            std::string line;
            while (std::getline(in, line)) {
                const auto hash = line.find('#');
                if (hash != std::string::npos) line.resize(hash);
                std::istringstream ls(line);
                unsigned long long s = 0;
                if (ls >> s) o.seeds.push_back(s);
            }
        } else if (a == "--help" || a == "-h") {
            usage();
            std::exit(0);
        } else {
            std::printf("unknown option %s\n", a.c_str());
            usage();
            return false;
        }
    }
    return true;
}

} // namespace

int main(int argc, char** argv) {
    Run run;
    if (!parse(argc, argv, run.opt)) return 2;
    std::vector<std::uint64_t> seeds = run.opt.seeds;
    if (!run.opt.seeds_file)
        for (int i = 0; i < run.opt.count; ++i) seeds.push_back(run.opt.seed + static_cast<std::uint64_t>(i));

    const double start = now_ms();
    for (const std::uint64_t s : seeds) {
        if (run.opt.time_limit > 0 && (now_ms() - start) / 1000.0 > run.opt.time_limit) {
            std::printf("time limit reached after %d cases\n", run.cases);
            break;
        }
        if (!run.opt.progress.empty()) {
            std::ofstream p(run.opt.progress, std::ios::trunc);
            p << s << "\n";
        }
        try {
            run_case(run, s);
        } catch (const std::exception& e) {
            run.seed = s;
            run.fail("case", "exception", e.what());
        } catch (...) {
            run.seed = s;
            run.fail("case", "exception", "unknown exception");
        }
        // With --progress (a supervised run that may crash), keep the reports
        // current so a crash loses only the case that crashed.
        if (!run.opt.progress.empty()) {
            const double so_far = (now_ms() - start) / 1000.0;
            if (!run.opt.json.empty()) std::ofstream(run.opt.json) << json(run, so_far);
            if (!run.opt.markdown.empty()) std::ofstream(run.opt.markdown) << markdown(run, so_far);
        }
    }
    if (!run.opt.progress.empty()) std::remove(run.opt.progress.c_str());
    const double seconds = (now_ms() - start) / 1000.0;

    const std::string md = markdown(run, seconds);
    std::fputs(md.c_str(), stdout);
    if (!run.opt.markdown.empty()) std::ofstream(run.opt.markdown) << md;
    if (!run.opt.json.empty()) std::ofstream(run.opt.json) << json(run, seconds);
    return run.failures.empty() ? 0 : 1;
}
