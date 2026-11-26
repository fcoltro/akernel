// Certified angles across the edges of a Brep (op/edge_continuity.hpp).
//
// Fillets must be proven tangent (G1) to the faces they blend, every other
// edge proven sharp, and every enclosure must hold the angle sampled in
// doubles along the edge. Edges outside the covered class are refused with a
// reason.
#include <akernel/geom/bernstein.hpp>
#include <akernel/op/boolean.hpp>
#include <akernel/op/edge_continuity.hpp>
#include <akernel/op/fillet.hpp>
#include <akernel/op/surfaces.hpp>

#include <cmath>
#include <cstdio>
#include <numbers>
#include <string>

#include "check.hpp"

using namespace ak;
using geom::Vec3d;

namespace {

constexpr double kHalfPi = std::numbers::pi / 2;
/// Below this an edge counts as proven tangent. Bounds on tangent edges are
/// rounding: up to 2e-12 rad on the sphere octants of a rounded box.
constexpr double kTangent = 1e-10;

/// Flat, judged in doubles from normals at 25 points (test bookkeeping only).
bool is_flat(const geom::PolyPatchD& P) {
    const geom::BezierSurface S(P);
    auto n = [&](double u, double v) {
        const geom::Jet3 J = S.jet(u, v);
        const Vec3d N = cross(J[geom::jet_index(1, 0)], J[geom::jet_index(0, 1)]);
        return N * (1.0 / norm(N));
    };
    const Vec3d n0 = n(0.5, 0.5);
    for (int i = 0; i <= 4; ++i)
        for (int j = 0; j <= 4; ++j)
            if (norm(cross(n0, n(0.05 + 0.225 * i, 0.05 + 0.225 * j))) > 1e-12) return false;
    return true;
}

/// Edges proven tangent between a curved face and a flat one: the contacts
/// of a fillet with the plates it blends.
int curved_to_flat_tangent(const op::Brep& b) {
    int n = 0;
    for (int i = 0; i < static_cast<int>(b.edges.size()); ++i) {
        const op::BEdge& e = b.edges[static_cast<std::size_t>(i)];
        if (e.arc || e.proc) continue;
        if (is_flat(b.carriers[static_cast<std::size_t>(e.on[0].carrier)]) ==
            is_flat(b.carriers[static_cast<std::size_t>(e.on[1].carrier)]))
            continue;
        const op::EdgeAngle a = op::certify_edge_angle(b, i);
        if (a.ok && a.hi < kTangent) ++n;
    }
    return n;
}

/// The angle between the two faces' outward normals at edge parameter t, in doubles.
double sampled_angle(const op::Brep& b, int edge, double t) {
    const op::BEdge& e = b.edges[static_cast<std::size_t>(edge)];
    Vec3d n[2];
    int k = 0;
    for (const op::BFace& f : b.faces) {
        bool uses = false;
        for (const auto& loop : f.loops)
            for (const op::BCoedge& ce : loop) uses = uses || ce.edge == edge;
        if (!uses || k == 2) continue;
        const op::BEdge::SideRef* ref = op::side_on(e, f.carrier);
        double u = 0, v = 0;
        op::side_uv(ref->side, ref->reversed ? 1.0 - t : t, u, v);
        const geom::Jet3 J = geom::BezierSurface(b.carriers[static_cast<std::size_t>(f.carrier)]).jet(u, v);
        Vec3d N = cross(J[geom::jet_index(1, 0)], J[geom::jet_index(0, 1)]);
        if (f.flipped) N = N * -1.0;
        n[k++] = N;
    }
    return std::atan2(norm(cross(n[0], n[1])), dot(n[0], n[1]));
}

struct Tally {
    int tangent = 0, right = 0, other = 0, refused = 0;
    double worst_tangent = 0.0, widest = 0.0;
};

/// Certifies every edge; checks soundness against samples and narrowness.
Tally certify_all(const char* name, const op::Brep& b) {
    Tally t;
    int unsound = 0;
    for (int i = 0; i < static_cast<int>(b.edges.size()); ++i) {
        const op::EdgeAngle a = op::certify_edge_angle(b, i);
        if (!a.ok) {
            ++t.refused;
            continue;
        }
        for (int s = 0; s <= 8; ++s) {
            const double x = sampled_angle(b, i, s / 8.0);
            // The sample has rounding of its own, of the order of 1e-15
            // relative to the curvature scale.
            if (x < a.lo - 1e-12 || x > a.hi + 1e-12) ++unsound;
        }
        if (a.hi < kTangent) {
            ++t.tangent;
            t.worst_tangent = std::max(t.worst_tangent, a.hi);
        } else {
            t.widest = std::max(t.widest, a.hi - a.lo);
            if (a.lo <= kHalfPi && kHalfPi <= a.hi && a.hi - a.lo < 1e-12) ++t.right;
            else ++t.other;
        }
    }
    std::printf("  %s: %zu edges, %d proven tangent (worst %.1e rad), %d right angles, %d other, %d refused, "
                "widest sharp enclosure %.1e\n",
                name, b.edges.size(), t.tangent, t.worst_tangent, t.right, t.other, t.refused, t.widest);
    AK_CHECK(unsound == 0);
    if (unsound) std::printf("FAIL %s: %d samples outside their enclosure\n", name, unsound);
    return t;
}

} // namespace

int main() {
    // A box: twelve right angles, nothing tangent.
    {
        const Tally t = certify_all("box", op::to_brep(op::make_box({0, 0, 0}, {2, 1, 0.5})));
        AK_CHECK(t.right == 12 && t.tangent == 0 && t.other == 0 && t.refused == 0);
    }
    // A parallelepiped: sharp edges at other angles, enclosed narrowly.
    {
        const Tally t = certify_all("parallelepiped",
                                    op::to_brep(op::make_box({0, 0, 0}, {1, 0, 0}, {0.5, 1, 0}, {0, 0.25, 1})));
        AK_CHECK(t.tangent == 0 && t.refused == 0 && t.right + t.other == 12);
        AK_CHECK(t.widest < 1e-12);
    }
    // One box edge rounded. The solid is split into cells, so many edges
    // join coplanar or co-cylindrical pieces (tangent too). The fillet's own
    // edges: tangent to both plates, at right angles to the two end caps.
    {
        const op::FilletedSolid f = op::fillet_box_edge({0, 0, 0}, {2, 0, 0}, {0, 1.5, 0}, {0, 0, 1}, 0.5);
        AK_CHECK(f.solid.has_value());
        if (f.solid) {
            const op::Brep B = op::to_brep(*f.solid);
            const Tally t = certify_all("box, one edge rounded", B);
            AK_CHECK(t.refused == 0 && t.other == 0);
            const int cf = curved_to_flat_tangent(B);
            std::printf("    proven tangent between a curved and a flat face: %d\n", cf);
            AK_CHECK(cf == 2);
            AK_CHECK(t.tangent == 17 && t.right == 19);
        }
        // Chamfered instead: no edge tangent, the two chamfer edges at 45 degrees.
        const op::BlendedSolid c = op::chamfer_box_edge({0, 0, 0}, {2, 0, 0}, {0, 1.5, 0}, {0, 0, 1}, 0.5, 0.5);
        AK_CHECK(c.solid.has_value());
        if (c.solid) {
            const op::Brep B = op::to_brep(*c.solid);
            const Tally t = certify_all("box, one edge chamfered", B);
            AK_CHECK(t.other == 2 && t.refused == 0);
            int at45 = 0;
            for (int i = 0; i < static_cast<int>(B.edges.size()); ++i) {
                const op::EdgeAngle a = op::certify_edge_angle(B, i);
                if (a.ok && a.lo <= std::numbers::pi / 4 && std::numbers::pi / 4 <= a.hi) ++at45;
            }
            AK_CHECK(at45 == 2);
        }
    }
    // Rings: the torus-zone fillet is tangent to the cylinder and the cap,
    // and consecutive quarters of one surface meet tangentially.
    for (int quarters : {1, 4}) {
        const auto r = op::blend_ring_edge(1.0, 2.0, 0.0, 1.0, op::RingEdge::OuterTop, op::EdgeCut::round(0.25), quarters);
        AK_CHECK(r.solid.has_value());
        if (!r.solid) continue;
        const Tally t = certify_all(quarters == 4 ? "ring, outer top edge rounded, full turn"
                                                  : "ring, outer top edge rounded, one quarter",
                                    op::to_brep(*r.solid));
        AK_CHECK(t.refused == 0 && t.other == 0);
        // Per quarter: the fillet against the cylinder and against the cap.
        // The torus zone against the flat cap: one contact per quarter.
        const int cf = curved_to_flat_tangent(op::to_brep(*r.solid));
        std::printf("    proven tangent between a curved and a flat face: %d\n", cf);
        AK_CHECK(cf == quarters);
        // The three sharp circles (inner top, inner bottom, outer bottom).
        if (quarters == 4) AK_CHECK(t.right == 3 * 4);
    }
    // Every edge of a rounded box: cylinders and sphere octants, all tangent.
    {
        const auto rb = op::rounded_box({0, 0, 0}, 2.0, 1.5, 1.0, 0.25);
        AK_CHECK(rb.solid.has_value());
        if (rb.solid) {
            const Tally t = certify_all("box, every edge rounded", op::to_brep(*rb.solid));
            AK_CHECK(t.refused == 0 && t.right == 0 && t.other == 0 && t.tangent == static_cast<int>(op::to_brep(*rb.solid).edges.size()));
        }
    }
    // A chamfer that is almost tangent: setbacks 0.5 and 0.001, so the
    // chamfer leaves one plate at atan(0.001 / 0.5), about 2e-3 rad. It must
    // be proven not tangent, with the closed-form angle enclosed.
    {
        const op::BlendedSolid c = op::chamfer_box_edge({0, 0, 0}, {2, 0, 0}, {0, 1.5, 0}, {0, 0, 1}, 0.5, 0.001);
        AK_CHECK(c.solid.has_value());
        if (c.solid) {
            const op::Brep B = op::to_brep(*c.solid);
            const Tally t = certify_all("box, one edge chamfered almost tangentially", B);
            AK_CHECK(t.refused == 0);
            const double shallow = std::atan2(0.001, 0.5), large = kHalfPi - shallow;
            int at_shallow = 0, at_large = 0;
            for (int i = 0; i < static_cast<int>(B.edges.size()); ++i) {
                const op::EdgeAngle a = op::certify_edge_angle(B, i);
                if (!a.ok) continue;
                if (a.lo > 0.0 && a.lo <= shallow + 1e-15 && shallow - 1e-15 <= a.hi && a.hi - a.lo < 1e-12) {
                    ++at_shallow;
                    std::printf("    proven not tangent: [%.17g, %.17g] rad, closed form %.17g\n", a.lo, a.hi, shallow);
                }
                if (a.lo <= large + 1e-15 && large - 1e-15 <= a.hi) ++at_large;
            }
            AK_CHECK(at_shallow == 1 && at_large == 1);
        }
    }
    // Outside the covered class: intersection edges are refused.
    {
        const op::BooleanResult r = op::boolean_op(op::to_brep(op::make_box({0, 0, 0}, {1, 1, 1})),
                                                   op::to_brep(op::make_cylinder_sector(0.25, 0.5, 0.0, 1.0, -1.0, 2.0)),
                                                   op::Op::Difference);
        AK_CHECK(r.ok);
        if (r.ok) {
            int arc = 0, refused_arc = 0;
            for (int i = 0; i < static_cast<int>(r.solid.edges.size()); ++i) {
                if (!r.solid.edges[static_cast<std::size_t>(i)].arc) continue;
                ++arc;
                const op::EdgeAngle a = op::certify_edge_angle(r.solid, i);
                if (!a.ok && a.why == "an edge along a traced intersection curve") ++refused_arc;
            }
            std::printf("  box minus cylinder sector: %d intersection edges, %d refused with the reason\n", arc,
                        refused_arc);
            AK_CHECK(arc > 0 && refused_arc == arc);
        }
    }
    return akt::finish("edge_continuity_test");
}
