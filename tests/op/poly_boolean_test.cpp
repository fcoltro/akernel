// Booleans between a solid with procedural faces and a flat-faced solid
// (op/poly_boolean.hpp). Every result must be a closed shell that meshes
// watertight; the volumes must satisfy the inclusion-exclusion identities,
// match closed forms where there are some, and match the same region cut
// out by successive plane cuts (op/plane_cut.hpp), an independent route.
#include <akernel/op/plane_cut.hpp>
#include <akernel/op/poly_boolean.hpp>
#include <akernel/op/proc_fillet.hpp>
#include <akernel/op/tessellate.hpp>

#include <cmath>
#include <cstdio>

#include "check.hpp"

using namespace ak;
using geom::Vec3d;

namespace {

op::PolyPatchD top_patch(double a, double L, double h, double c, double x0 = 0, double y0 = 0) {
    op::PolyPatchD T(1, 3);
    static constexpr double s[4] = {0, 0, 1, 1};
    for (int i = 0; i <= 1; ++i)
        for (int j = 0; j <= 3; ++j) T.at(i, j) = {x0 + a * i, y0 + L * j / 3.0, h + c * a * i * s[j], 1.0};
    return T;
}

const char* op_name(op::PolyOp o) {
    return o == op::PolyOp::Union ? "union" : o == op::PolyOp::Intersection ? "intersection" : "difference";
}

/// Check a result: closed, watertight; the volume returned.
double check_result(const char* name, op::PolyOp o, const op::PolyBooleanResult& r) {
    if (!r.solid) std::printf("%s %s: refused: %s\n", name, op_name(o), r.error.c_str());
    AK_CHECK(r.solid.has_value());
    if (!r.solid) return 0.0;
    const auto rep = op::audit(*r.solid);
    AK_CHECK(rep.closed());
    const double v = op::volume(*r.solid);
    const op::TessResult t = op::tessellate(*r.solid);
    if (!t.ok) std::printf("  tessellation failed: %s\n", t.error.c_str());
    AK_CHECK(t.ok);
    double dev = -1.0;
    if (t.ok) {
        AK_CHECK(op::audit(t.mesh).watertight());
        dev = op::mesh_deviation(t.mesh, *r.solid);
        AK_CHECK(dev < 2e-3);
    }
    std::printf("%s %s: %zu faces, %d curves, closed %d, euler %d, volume %.15f, mesh deviation %.2e, %.2f s\n", name,
                op_name(o), r.solid->faces.size(), r.curves, static_cast<int>(rep.closed()), rep.euler, v, dev, r.seconds);
    return v;
}

/// The three operations, and their volume identities.
struct Vols {
    double uni, inter, diff;
};
Vols check_all(const char* name, const op::Brep& A, const op::Brep& P, double tol) {
    const double vA = op::volume(A), vP = op::volume(P);
    Vols w{};
    w.uni = check_result(name, op::PolyOp::Union, op::poly_boolean(A, P, op::PolyOp::Union));
    w.inter = check_result(name, op::PolyOp::Intersection, op::poly_boolean(A, P, op::PolyOp::Intersection));
    w.diff = check_result(name, op::PolyOp::Difference, op::poly_boolean(A, P, op::PolyOp::Difference));
    std::printf("  identities: union + intersection - A - P = %.2e, difference + intersection - A = %.2e\n",
                w.uni + w.inter - vA - vP, w.diff + w.inter - vA);
    AK_CHECK_NEAR(w.uni + w.inter, vA + vP, tol);
    AK_CHECK_NEAR(w.diff + w.inter, vA, tol);
    return w;
}

/// A cut to the box lo..hi by its six planes.
double box_by_cuts(const op::Brep& A, const Vec3d& lo, const Vec3d& hi) {
    op::Brep S = A;
    const Vec3d n[6] = {{1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1}};
    const double d[6] = {hi.x, -lo.x, hi.y, -lo.y, hi.z, -lo.z};
    for (int k = 0; k < 6; ++k) {
        // Planes that miss the solid are skipped by checking the vertices.
        bool cuts = false;
        for (const op::BVertex& v : S.vertices) cuts = cuts || dot(n[k], v.point) > d[k];
        if (!cuts) continue;
        const op::PlaneCutResult r = op::cut_by_plane(S, n[k], d[k]);
        if (!r.solid) {
            std::printf("  plane cut %d refused: %s\n", k, r.error.c_str());
            AK_CHECK(false);
            return 0.0;
        }
        S = *r.solid;
    }
    return op::volume(S);
}

/// Removed volume of the rounded block over y in [0, y_end] by slices (as in
/// proc_fillet_test).
double removed_by_slices(const solve::BlendCarrier& blend, double a, double L, double h, double c, double y_end) {
    auto cut_area = [&](double y0, int n) {
        const double t = y0 / L, sy = t * t * (3 - 2 * t);
        std::vector<Vec3d> poly;
        for (int i = 0; i <= n; ++i) {
            const double th = static_cast<double>(i) / n;
            Vec3d P, Ps, Pt;
            auto f = [&](double sg) {
                blend.eval(sg, th, P, Ps, Pt);
                return P.y - y0;
            };
            double lo = 0.0, hi = 1.0, flo = f(lo), fhi = f(hi);
            int side = 0;
            for (int it = 0; it < 100 && hi - lo > 1e-15; ++it) {
                const double m = (lo * fhi - hi * flo) / (fhi - flo), fm = f(m);
                if (fm == 0.0) { lo = hi = m; break; }
                if ((fm > 0.0) == (flo > 0.0)) { lo = m; flo = fm; if (side == -1) fhi /= 2; side = -1; }
                else { hi = m; fhi = fm; if (side == 1) flo /= 2; side = 1; }
                if (std::fabs(fm) < 1e-14) { lo = hi = m; break; }
            }
            f(0.5 * (lo + hi));
            poly.push_back(P);
        }
        poly.push_back({a, y0, h + c * a * sy});
        double A = 0.0;
        for (std::size_t i = 0; i < poly.size(); ++i) {
            const Vec3d& p = poly[i];
            const Vec3d& q = poly[(i + 1) % poly.size()];
            A += p.x * q.z - q.x * p.z;
        }
        return std::fabs(A) / 2;
    };
    static constexpr double gx[4] = {0.1834346424956498, 0.5255324099163290, 0.7966664774136267, 0.9602898564975363};
    static constexpr double gw[4] = {0.3626837833783620, 0.3137066458778873, 0.2223810344533745, 0.1012285362903763};
    const int cells = 8;
    double total = 0.0;
    for (int k = 0; k < cells; ++k)
        for (int i = 0; i < 4; ++i)
            for (int sgn = -1; sgn <= 1; sgn += 2) {
                const double y = y_end * (k + 0.5 + 0.5 * sgn * gx[i]) / cells;
                const double a1 = cut_area(y, 300), a2 = cut_area(y, 600);
                total += 0.5 * gw[i] * (y_end / cells) * (a2 + (a2 - a1) / 3);
            }
    return total;
}

} // namespace

int main() {
    // A cube and a box through one of its corners' edges: closed forms.
    {
        const op::Brep A = op::to_brep(op::make_box(Vec3d{0, 0, 0}, Vec3d{1, 1, 1}));
        const op::Brep P = op::to_brep(op::make_box(Vec3d{0.5, 0.25, 0.3}, Vec3d{1.5, 0.75, 1.7}));
        const Vols w = check_all("cube and box", A, P, 1e-12);
        AK_CHECK_NEAR(w.inter, 0.5 * 0.5 * 0.7, 1e-12);
        AK_CHECK_NEAR(w.diff, 1.0 - 0.175, 1e-12);
    }

    // The rounded freeform block (a = 1, L = 2, h = 3/4, r = 1/4, c = 1/2)
    // and a box notch through the blend, the top and the side: the
    // intersection is compared with six plane cuts of the block.
    const double a = 1.0, L = 2.0, h = 0.75, r = 0.25, c = 0.5;
    const op::ProcFilletResult rb = op::round_block_edge(top_patch(a, L, h, c), r);
    AK_CHECK(rb.solid.has_value());
    if (rb.solid) {
        const Vec3d lo{0.5, 0.5, 0.3}, hi{1.5, 1.0, 2.0};
        const op::Brep P = op::to_brep(op::make_box(lo, hi));
        const Vols w = check_all("rounded block and notch", *rb.solid, P, 1e-11);
        const double cuts = box_by_cuts(*rb.solid, lo, hi);
        std::printf("  intersection by plane cuts %.15f, difference %.2e\n", cuts, w.inter - cuts);
        // A third route: the box part of the block minus the blend's
        // removed volume over y in [1/2, 1], by slices.
        auto S = [&](double y) { const double t = y / L; return L * (t * t * t - t * t * t * t / 2); };
        const double boxpart = (hi.y - lo.y) * (a - lo.x) * (h - lo.z) + c * (a * a - lo.x * lo.x) / 2 * (S(hi.y) - S(lo.y));
        const double removed = removed_by_slices(*rb.blend, a, L, h, c, hi.y) - removed_by_slices(*rb.blend, a, L, h, c, lo.y);
        std::printf("  intersection by slices %.15f, difference %.2e\n", boxpart - removed, w.inter - (boxpart - removed));
        AK_CHECK_NEAR(w.inter, boxpart - removed, 2e-11);
        AK_CHECK_NEAR(w.inter, cuts, 1e-11);
    }
    return akt::finish("poly_boolean_test");
}
