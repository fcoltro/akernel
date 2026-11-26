// Certified self-intersection check (op/validity.hpp).
//
// Valid solids must be proven valid or left undecided, never called invalid;
// solids whose faces cross must be proven invalid, with the crossing faces,
// or left undecided, never called valid.
#include <akernel/op/fillet.hpp>
#include <akernel/op/surfaces.hpp>
#include <akernel/op/validity.hpp>

#include <cstdio>
#include <string>

#include "check.hpp"

using namespace ak;
using geom::Vec3d;

namespace {

const char* name_of(op::Validity::State s) {
    return s == op::Validity::Valid ? "valid" : s == op::Validity::Invalid ? "invalid" : "undecided";
}

op::Validity report(const char* name, const op::Brep& b) {
    const op::Validity v = op::check_self_intersection(b);
    std::printf("  %s: %s", name, name_of(v.state));
    if (v.state != op::Validity::Valid) std::printf(" (%s; faces %d and %d)", v.why.c_str(), v.face_a, v.face_b);
    if (v.state == op::Validity::Invalid) std::printf(" near (%.3f, %.3f, %.3f)", v.where.x, v.where.y, v.where.z);
    std::printf("; %d pairs: %d by boxes, %d by planes, %d by intersection, %d undecided\n", v.pairs, v.by_boxes, v.by_planes,
                v.by_ssi, v.undecided);
    return v;
}

/// Two Breps as one, without joining anything: two shells.
op::Brep combine(const op::Brep& a, const op::Brep& b) {
    op::Brep r = a;
    const int nc = static_cast<int>(a.carriers.size()), nv = static_cast<int>(a.vertices.size()),
              ne = static_cast<int>(a.edges.size());
    r.carriers.insert(r.carriers.end(), b.carriers.begin(), b.carriers.end());
    for (op::BVertex v : b.vertices) {
        for (op::Chart& c : v.charts) c.carrier += nc;
        r.vertices.push_back(v);
    }
    for (op::BEdge e : b.edges) {
        e.v0 += nv;
        e.v1 += nv;
        for (auto& s : e.on)
            if (s.carrier >= 0) s.carrier += nc;
        r.edges.push_back(e);
    }
    for (op::BFace f : b.faces) {
        f.carrier += nc;
        for (auto& loop : f.loops)
            for (op::BCoedge& ce : loop) ce.edge += ne;
        r.faces.push_back(f);
    }
    return r;
}

op::Brep box(double x0, double y0, double z0, double x1, double y1, double z1) {
    return op::to_brep(op::make_box({x0, y0, z0}, {x1, y1, z1}));
}

} // namespace

int main() {
    // Valid solids: proven valid.
    AK_CHECK(report("unit box", box(0, 0, 0, 1, 1, 1)).state == op::Validity::Valid);
    AK_CHECK(report("parallelepiped", op::to_brep(op::make_box({0, 0, 0}, {1, 0, 0}, {0.5, 1, 0}, {0, 0.25, 1}))).state ==
             op::Validity::Valid);
    AK_CHECK(report("two boxes apart", combine(box(0, 0, 0, 1, 1, 1), box(2, 0, 0, 3, 1, 1))).state ==
             op::Validity::Valid);
    // Two shells sharing a face region: two faces in one plane overlap, so
    // the union of the boundaries is not the boundary of a solid.
    AK_CHECK(report("two boxes touching along a face", combine(box(0, 0, 0, 1, 1, 1), box(1, 0, 0, 2, 1, 1))).state ==
             op::Validity::Invalid);
    AK_CHECK(report("cylinder sector", op::to_brep(op::make_cylinder_sector(1.0, 2.0, 0.0, 1.0, 0.0, 1.0))).state ==
             op::Validity::Valid);
    AK_CHECK(report("box with a hill on top", op::to_brep(op::make_bump_box({0, 0, 0}, 1, 1, 1, 0.3))).state ==
             op::Validity::Valid);
    {
        const op::FilletedSolid f = op::fillet_box_edge({0, 0, 0}, {2, 0, 0}, {0, 1.5, 0}, {0, 0, 1}, 0.5);
        AK_CHECK(f.solid.has_value());
        if (f.solid) AK_CHECK(report("box with one edge rounded", op::to_brep(*f.solid)).state == op::Validity::Valid);
        const op::BlendedSolid c = op::chamfer_box_edge({0, 0, 0}, {2, 0, 0}, {0, 1.5, 0}, {0, 0, 1}, 0.5, 0.25);
        AK_CHECK(c.solid.has_value());
        if (c.solid) AK_CHECK(report("box with one edge chamfered", op::to_brep(*c.solid)).state == op::Validity::Valid);
    }
    for (int quarters : {1, 4}) {
        const auto r = op::blend_ring_edge(1.0, 2.0, 0.0, 1.0, op::RingEdge::OuterTop, op::EdgeCut::round(0.25), quarters);
        AK_CHECK(r.solid.has_value());
        if (r.solid)
            AK_CHECK(report(quarters == 4 ? "ring with a rounded edge, full turn" : "ring with a rounded edge, one quarter",
                            op::to_brep(*r.solid))
                         .state != op::Validity::Invalid);
    }
    {
        const auto rb = op::rounded_box({0, 0, 0}, 2.0, 1.5, 1.0, 0.25);
        AK_CHECK(rb.solid.has_value());
        if (rb.solid) AK_CHECK(report("box with every edge rounded", op::to_brep(*rb.solid)).state != op::Validity::Invalid);
    }

    // Invalid solids: faces that cross.
    {
        const op::Validity v = report("two overlapping boxes", combine(box(0, 0, 0, 1, 1, 1), box(0.5, 0.25, 0.25, 1.5, 0.75, 0.75)));
        AK_CHECK(v.state == op::Validity::Invalid);
    }
    {
        // Same bottom plane, overlapping there: caught exactly as coplanar.
        const op::Validity v = report("two boxes overlapping in a shared plane", combine(box(0, 0, 0, 1, 1, 1), box(0.5, 0.5, 0, 1.5, 1.5, 1)));
        AK_CHECK(v.state == op::Validity::Invalid);
    }
    {
        const op::Validity v = report("box through a cylinder sector",
                                      combine(op::to_brep(op::make_cylinder_sector(1.0, 2.0, 0.0, 1.0, 0.0, 1.0)),
                                              box(0.2, 1.2, 0.25, 0.6, 3.0, 0.75)));
        AK_CHECK(v.state == op::Validity::Invalid);
    }
    {
        // A hill pushed up through a slab floating above it.
        const op::Validity v = report("hill through a slab", combine(op::to_brep(op::make_bump_box({0, 0, 0}, 1, 1, 1, 0.5)),
                                                                     box(0.25, 0.25, 1.2, 0.75, 0.75, 1.4)));
        AK_CHECK(v.state == op::Validity::Invalid);
    }
    {
        // One shell: a box whose top sags through its floor.
        const op::Validity v = report("box whose top sags through its floor",
                                      op::to_brep(op::make_bump_box({0, 0, 0}, 1, 1, 1, -1.5)));
        AK_CHECK(v.state == op::Validity::Invalid);
    }
    return akt::finish("validity_test");
}
