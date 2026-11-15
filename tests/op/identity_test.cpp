// Algebraic identities of Boolean operations (docs/BLUEPRINT.md §5.4): two
// routes to the same solid must give the same solid, with the same topology,
// not just the same volume. op::compare_solids decides "the same" (see
// op/compare.hpp); the volume is checked as well.
#include <akernel/op/boolean.hpp>
#include <akernel/op/compare.hpp>

#include <cmath>
#include <cstdio>
#include <functional>
#include <map>
#include <string>
#include <vector>

#include "check.hpp"

using namespace ak;
using geom::Vec3d;
using op::Brep;
using op::Op;

namespace {

struct Tally {
    int evaluated = 0, same = 0, not_solved = 0;
};
std::map<std::string, Tally> g_tally;
bool g_verbose = false;

/// A result, or nothing if the kernel did not solve the operation.
bool apply(const Brep& a, const Brep& b, Op o, Brep& out, std::string& why) {
    const op::BooleanResult r = op::boolean_op(a, b, o);
    if (!r.ok) { why = r.error; return false; }
    out = r.solid;
    return true;
}

/// Check that two routes give the same solid. Each route is a function that
/// computes a solid or reports that a step was not solved.
using Route = std::function<bool(Brep&, std::string&)>;
void identity(const std::string& name, const std::string& where, const Route& left, const Route& right) {
    Tally& t = g_tally[name];
    Brep x, y;
    std::string why;
    if (!left(x, why) || !right(y, why)) {
        ++t.not_solved;
        if (g_verbose) std::printf("  %s, %s: not evaluated (%s)\n", name.c_str(), where.c_str(), why.c_str());
        return;
    }
    ++t.evaluated;
    const op::SolidComparison c = op::compare_solids(x, y);
    AK_CHECK(c.same);
    if (!c.same) std::printf("FAIL %s, %s: %s\n", name.c_str(), where.c_str(), c.reason.c_str());
    const double vx = x.faces.empty() ? 0.0 : op::volume(x), vy = y.faces.empty() ? 0.0 : op::volume(y);
    AK_CHECK_NEAR(vx, vy, 1e-10);
    if (c.same) ++t.same;
}

Route op1(const Brep& a, const Brep& b, Op o) {
    return [&a, &b, o](Brep& out, std::string& why) { return apply(a, b, o, out, why); };
}
/// (a o1 b) o2 c
Route op2(const Brep& a, Op o1, const Brep& b, Op o2, const Brep& c) {
    return [&a, &b, &c, o1, o2](Brep& out, std::string& why) {
        Brep t;
        return apply(a, b, o1, t, why) && apply(t, c, o2, out, why);
    };
}
/// a o1 (b o2 c)
Route op2r(const Brep& a, Op o1, const Brep& b, Op o2, const Brep& c) {
    return [&a, &b, &c, o1, o2](Brep& out, std::string& why) {
        Brep t;
        return apply(b, c, o2, t, why) && apply(a, t, o1, out, why);
    };
}
/// (a o1 b) o2 (a o3 b)
Route op_pair(const Brep& a, const Brep& b, Op o1, Op o2, Op o3) {
    return [&a, &b, o1, o2, o3](Brep& out, std::string& why) {
        Brep s, t;
        return apply(a, b, o1, s, why) && apply(a, b, o3, t, why) && apply(s, t, o2, out, why);
    };
}
Route just(const Brep& a) {
    return [&a](Brep& out, std::string&) { out = a; return true; };
}

void pair_identities(const std::string& where, const Brep& A, const Brep& B) {
    identity("A and B = B and A", where, op1(A, B, Op::Intersection), op1(B, A, Op::Intersection));
    identity("A or B = B or A", where, op1(A, B, Op::Union), op1(B, A, Op::Union));
    identity("A - (A - B) = A and B", where, op2r(A, Op::Difference, A, Op::Difference, B), op1(A, B, Op::Intersection));
    identity("(A - B) or (A and B) = A", where, op_pair(A, B, Op::Difference, Op::Union, Op::Intersection), just(A));
    identity("A or (A and B) = A", where, op2r(A, Op::Union, A, Op::Intersection, B), just(A));
    identity("A and (A or B) = A", where, op2r(A, Op::Intersection, A, Op::Union, B), just(A));
    identity("(A or B) - B = A - B", where, op2(A, Op::Union, B, Op::Difference, B), op1(A, B, Op::Difference));
}

void triple_identities(const std::string& where, const Brep& A, const Brep& B, const Brep& C) {
    identity("(A or B) or C = A or (B or C)", where, op2(A, Op::Union, B, Op::Union, C),
             op2r(A, Op::Union, B, Op::Union, C));
    identity("(A and B) and C = A and (B and C)", where, op2(A, Op::Intersection, B, Op::Intersection, C),
             op2r(A, Op::Intersection, B, Op::Intersection, C));
    identity("A - (B or C) = (A - B) - C", where, op2r(A, Op::Difference, B, Op::Union, C),
             op2(A, Op::Difference, B, Op::Difference, C));
}

/// A box of random orientation, size and position (as in the benchmark).
Brep random_box(akt::Rng& rng, const Vec3d& lo, const Vec3d& hi) {
    const double a = rng.range(0, 6.283185307179586), b = rng.range(-1.2, 1.2), c = rng.range(0, 6.283185307179586);
    const double ca = std::cos(a), sa = std::sin(a), cb = std::cos(b), sb = std::sin(b), cc = std::cos(c), sc = std::sin(c);
    // Rotation Rz(a) Rx(b) Rz(c), columns.
    const Vec3d e0{ca * cc - sa * cb * sc, sa * cc + ca * cb * sc, sb * sc};
    const Vec3d e1{-ca * sc - sa * cb * cc, -sa * sc + ca * cb * cc, sb * cc};
    const Vec3d e2{sa * sb, -ca * sb, cb};
    const double sx = rng.range(0.4, 1.1), sy = rng.range(0.4, 1.1), sz = rng.range(0.4, 1.1);
    const Vec3d centre{rng.range(lo.x, hi.x), rng.range(lo.y, hi.y), rng.range(lo.z, hi.z)};
    const Vec3d o = centre - e0 * (0.5 * sx) - e1 * (0.5 * sy) - e2 * (0.5 * sz);
    return op::to_brep(op::make_box(o, e0 * sx, e1 * sy, e2 * sz));
}

} // namespace

int main(int argc, char** argv) {
    g_verbose = argc > 1 && std::string(argv[1]) == "verbose";
    const Brep unit = op::to_brep(op::make_box({0, 0, 0}, {1, 1, 1}));

    // --- the comparator itself -------------------------------------------------------
    {
        const Brep same_again = op::to_brep(op::make_box({0, 0, 0}, {1, 1, 1}));
        AK_CHECK(op::compare_solids(unit, same_again).same);
        // The same cube built mirrored: other patches, same planes, same solid.
        const Brep mirrored = op::to_brep(op::make_box({1, 0, 1}, {-1, 0, 0}, {0, 1, 0}, {0, 0, -1}));
        const op::SolidComparison m = op::compare_solids(unit, mirrored);
        AK_CHECK(m.same);
        if (!m.same) std::printf("mirrored cube: %s\n", m.reason.c_str());
        // Different solids are told apart.
        AK_CHECK(!op::compare_solids(unit, op::to_brep(op::make_box({0, 0, 0}, {1, 1, 2}))).same);
        AK_CHECK(!op::compare_solids(unit, op::to_brep(op::make_bump_box({0, 0, 0}, 1, 1, 1, 0.3))).same);
        AK_CHECK(!op::compare_solids(unit, Brep{}).same);
        AK_CHECK(op::compare_solids(Brep{}, Brep{}).same);
        // A solid and the same solid glued from two pieces: the seam faces and
        // the vertices along the cut are removed before comparing.
        const Brep lower = op::to_brep(op::make_box({0, 0, 0}, {1, 1, 0.5}));
        const op::BooleanResult glued = op::boolean_op(lower, op::to_brep(op::make_box({0, 0, 0.5}, {1, 1, 1})), Op::Union);
        AK_CHECK(glued.ok);
        if (glued.ok) {
            const op::SolidComparison g = op::compare_solids(unit, glued.solid);
            AK_CHECK(g.same);
            if (!g.same) std::printf("glued cube: %s\n", g.reason.c_str());
            std::printf("glued cube: %d faces, comparison %s (seams removed %d, vertices removed %d)\n",
                        static_cast<int>(glued.solid.faces.size()), g.same ? "same" : "different", g.seams[1],
                        g.removed_vertices[1]);
        }
    }

    // --- the empty solid as an operand -------------------------------------------------
    {
        const Brep none;
        const op::BooleanResult u = op::boolean_op(unit, none, Op::Union), u2 = op::boolean_op(none, unit, Op::Union);
        const op::BooleanResult i = op::boolean_op(unit, none, Op::Intersection);
        const op::BooleanResult d = op::boolean_op(unit, none, Op::Difference), d2 = op::boolean_op(none, unit, Op::Difference);
        AK_CHECK(u.ok && !u.empty && op::compare_solids(u.solid, unit).same);
        AK_CHECK(u2.ok && !u2.empty && op::compare_solids(u2.solid, unit).same);
        AK_CHECK(i.ok && i.empty);
        AK_CHECK(d.ok && !d.empty && op::compare_solids(d.solid, unit).same);
        AK_CHECK(d2.ok && d2.empty);
    }

    // --- curated pairs -----------------------------------------------------------------
    const Brep boxB = op::to_brep(op::make_box({0.5, 0.4, 0.3}, {1.5, 1.4, 1.3}));
    const Brep hill = op::to_brep(op::make_bump_box({0, 0, 0}, 1, 1, 1, 0.3));
    const Brep slab = op::to_brep(op::make_box({-0.5, -0.5, 1.15}, {1.5, 1.5, 2.0}));
    const Brep shell = op::to_brep(op::make_annular_sector(1.0, 2.0, 0.0, 1.0));
    const Brep shell_box = op::to_brep(op::make_box({0.5, 0.6, 0.25}, {1.3, 1.4, 0.75}));
    const Brep hill2 = op::to_brep(op::make_bump_box({0.6, 0.6, -0.35}, 1, 1, 1, 0.3));
    const double c = std::cos(0.6), s = std::sin(0.6);
    const Brep rotated = op::to_brep(
        op::make_box({0.4, 0.3, 0.35}, {0.8 * c, 0.8 * s, 0}, {-0.7 * s, 0.7 * c, 0.1}, {0.05, -0.05, 0.75}));
    const Brep face_to_face = op::to_brep(op::make_box({1, 0.2, 0.3}, {2, 1.3, 0.9}));
    const Brep overlapping = op::to_brep(op::make_box({0.5, 0, 0.3}, {1.5, 1, 0.8}));
    const Brep corner_box = op::to_brep(op::make_box({0.2, 0.2, 0.6}, {0.8, 0.7, 1.5}));

    pair_identities("offset boxes", unit, boxB);
    pair_identities("hill and slab (closed loop)", hill, slab);
    pair_identities("cylindrical shell and box", shell, shell_box);
    pair_identities("hill and cylindrical shell", hill2, shell);
    pair_identities("box and rotated box", unit, rotated);
    pair_identities("boxes face to face", unit, face_to_face);
    pair_identities("boxes overlapping, two planes shared", unit, overlapping);
    triple_identities("three boxes", unit, boxB, corner_box);
    triple_identities("hill, slab, box", hill, slab, corner_box);
    triple_identities("box, rotated box, box", unit, rotated, boxB);

    // --- random pairs and triples --------------------------------------------------------
    akt::Rng rng(4242);
    const Brep* fixed[3] = {&unit, &hill, &shell};
    const char* fixed_name[3] = {"unit box", "hill", "cylindrical shell"};
    const Vec3d lo[3] = {{0.1, 0.1, 0.1}, {0.1, 0.1, 0.2}, {0.5, 0.5, 0.2}};
    const Vec3d hi[3] = {{0.9, 0.9, 0.9}, {0.9, 0.9, 1.2}, {1.5, 1.5, 0.8}};
    for (int k = 0; k < 3; ++k)
        for (int trial = 0; trial < 6; ++trial) {
            const Brep B = random_box(rng, lo[k], hi[k]);
            const Brep C = random_box(rng, lo[k], hi[k]);
            const std::string where = std::string(fixed_name[k]) + ", random trial " + std::to_string(trial);
            pair_identities(where, *fixed[k], B);
            if (trial < 2) triple_identities(where, *fixed[k], B, C);
        }

    std::printf("| Identity | Evaluated | Same solid | Not evaluated (a step not solved) |\n|---|---|---|---|\n");
    for (const auto& [name, t] : g_tally)
        std::printf("| %s | %d | %d | %d |\n", name.c_str(), t.evaluated, t.same, t.not_solved);
    return akt::finish("identity_test");
}
