# akernel

A research geometric modelling kernel built around one rule: **geometric decisions are proven, not guessed**. Every predicate returns a proven answer or says it could not decide, and every intersection curve comes with a certificate.

- Design and theory: [docs/BLUEPRINT.md](docs/BLUEPRINT.md)
- Measured results: [docs/RESULTS.md](docs/RESULTS.md)

## Status

This is the certified core of the kernel described in the blueprint (roadmap Phases 0 and 1, the variational solver, the topology audit) and a first Boolean operation built on it. It is **not** a complete CAD kernel. The order of further work follows a survey of recent research, [docs/SURVEY.md](docs/SURVEY.md) (summarised in the blueprint, §5.3).

The Boolean (intersection, union, difference) works on solids with trimmed faces, and its result is a solid of the same type, so operations can be chained. Operands may share surfaces, edges and vertices by construction (a result combined with one of its own operands, a solid with itself); that sharing is recognised from identical control points and certified parameters, never from distance. Flat faces that lie in the same plane without sharing a construction (boxes placed face to face, stacked, or overlapping with aligned faces) are recognised by exact predicates on their control points. Results may be non-manifold: parts of a solid that touch only along an edge (four faces there) or at a point. Such a result is accepted when its boundary is closed and oriented and the faces around every shared edge alternate between material and empty space, and it can be an operand of further operations. Curved faces that lie on the same surface without being the same patch are recognised exactly when one patch is an affine reparameterisation of the other with dyadic shifts and power-of-two scales, with the two parameters in the same order or exchanged (an inner and an outer face of one cylinder): for example, sectors of one cylinder built with `op::make_cylinder_sector`, which uses an exact rational parameterisation of the circle. `op::on_quadric` decides exactly whether a patch lies on a cylinder, sphere, cone or other quadric. A curved face that touches a flat face's plane along sides of its patch (a box face against a cylinder along one of its lines, a hill on its base), or along a line inside a patch of degree 2 (a plane tangent to a cylinder along a ruling), is decided exactly, from the plane's equation restricted to the patch (`op/contact.hpp`); so are the edges that touch the plane at their ends, and the order of the curves and faces that meet tangentially there. A flat face on a rational carrier (a sector's annular cap) that lies in one plane with a face it is not affinely related to (a box face), or that is crossed in its plane by straight edges of the other solid (the base of a hill), is moved onto an exact parallelogram of its plane before the Boolean, so the exact coplanar overlay applies. A vertex of one solid that lies in a plane of the other at the end of one of its edges is decided exactly, and a circle that is a side of a cylinder patch in one solid and an intersection curve in the other is recognised as one edge. Vertices with round coordinates are known exactly: a new vertex such as (0.75, 1, 0.5), where a cylinder of radius 1.25 meets two planes, is proven to be that point and carries it into later operations, where "does this vertex lie on that surface" is then decided by exact arithmetic, for planes and for curved surfaces alike. Two hill tops (graphs over the same axes) that share a boundary line and cross along a curve leaving it are decided by exact factorisation of their height difference, and the crossing curve is certified against a shifted surface on which it is transversal; a crossing curve that touches the shared line is reported. Still reported as not solved: tangency at a point or between two curved surfaces other than such graphs, coincident curved surfaces related in other ways, flat rational faces in planes that are not axis-aligned, and an edge of one solid lying in a curved surface of the other other than along one of its straight lines (that case is decided exactly, and where the solids only touch along such a line the result reports a line contact), an edge through an edge in space at a point whose coordinates are not round (at a point with round coordinates it is decided exactly), and vertices whose coordinates are not round (their coincidences cannot be proven by the exact points above).

Any solid, operand or result, can be tessellated into a closed triangle mesh whose deviation from the surfaces follows a chosen chord tolerance.

Algebraic identities hold with identical topology, not just equal volume: A − (A − B) and A ∩ B, A ∪ B and B ∪ A, (A − B) ∪ (A ∩ B) and A, associativity and others, computed along both routes and compared by `op::compare_solids`.

| Blueprint area | Built and tested | Not built |
|---|---|---|
| Numerics (§3.2) | Rigorous interval arithmetic; exact expansion arithmetic; filtered exact `orient2d` / `orient3d`; escalation policy; exact points for vertices with round coordinates (exact evaluation of rational patches at dyadic and rational parameters, proven to be the vertex through its root's uniqueness box) | Multi-precision tier; algebraic numbers; lazy expression DAG |
| Carriers (§3.1) | Rational Bezier patches evaluated in doubles or intervals; NURBS surfaces with Bezier extraction; derivatives to order 3; trimmed faces with certified point classification; exact quadric patches (cylinder sectors from the exact rational circle) and an exact patch-on-quadric test | Procedural offsets and blends; subdivision surfaces; implicit carriers as their own type; full spheres and discs (degenerate patches); provenance DAG |
| Intersection (§1.4) | Certified surface/surface intersection: tubes with Krawczyk certificates, loop detection, proof of completeness, explicit report of what could not be proven. Certified curve/surface and three-surface intersection points with the same guarantees. Contact of a patch with a plane along sides of the patch, and of a curve with a plane at its ends, decided exactly; a curve that meets both parameter squares only on their boundaries traced through over enlarged squares | Multi-patch stitching; tangency at a point or between curved surfaces (provenance or exact algebra beyond planes) |
| Booleans (§3.3) | Imprint and intersection / union / difference of solids with trimmed faces, with no tessellation and no tolerance; results are operands of the same type (chained operations); operands that share surfaces, edges and vertices by construction (on-boundary regions classified with the same or opposite orientation); flat faces in contact decided exactly; non-manifold results and operands (contact along an edge or at a point); every result checked by the integer audit and, where not manifold, by the fan check around shared edges; a broad phase over a bounding-volume hierarchy of certified boxes, so the overhead grows linearly with the number of faces; flat faces in different planes intersected along their one line with a direct completeness proof, no subdivision (two 1,536-face solids in about 0.25 s); a curved face touching a flat face along a line, including solids that meet only along it; flat faces on rational carriers (sector caps) moved onto exact parallelograms where they share a plane with an unrelated flat face or are crossed there by edges of the other solid; a vertex in a plane of the other solid at the end of an edge decided exactly; vertices with exact points, identified across operands and decided exactly against planes and curved surfaces they lie on; an undecided point location is always reported, never given a meaning | Coincident curved surfaces related other than by an affine map; tangency at a point or between curved faces; a surface through a vertex inside an edge; healing |
| Continuity (§1.7) | G1 / G2 / G3 jump measurement by parameterisation-independent invariants; certified bounds on the jumps over a whole seam between polynomial or rational patches in a common parameterisation (Bernstein bounds on the jet residuals, carried through the invariants by interval forward differentiation) | Certified bounds for reparameterised joins (gauge functions); tight bounds for rational patches |
| Class-A (§4) | Variational four-sided fill with prescribed jets up to order 3 and energies up to order 4 | N-sided patches; gauge functions (true G^k with reparameterisation); fillets; crown constraints; nonlinear geometric energies |
| Topology (§1.5) | Cell complex with integer audit: loop closure, chain condition, closed and oriented boundary (allowing edges with four or more faces and pinched vertices), manifoldness, Euler-Poincare genus, tolerance ledger; a boundary representation with trimmed faces whose edges and vertices carry certificates instead of tolerances | Euler operators; persistent store; self-intersection check |
| Diagnostics (§2.2) | Conditioning indicator, step and tier telemetry, convergence-rate estimators, tangency early warning | Learned scheduler |
| Tessellation (§5.1, ak::viz) | Watertight meshes of trimmed faces: each edge sampled once and shared by its two faces, faces triangulated by a constrained Delaunay triangulation with exact predicates; mesh audit (closure, Euler characteristic against the Brep's genus) and measured deviation; OBJ output | Certified deviation bounds; curvature-adapted anisotropic meshes; GPU tessellation |
| Other | | GPU paths; STEP / IGES; healing |

## Build and test

Requires the Visual Studio Build Tools (MSVC, CMake and Ninja are bundled).

```powershell
.\build.ps1                 # Release build, run all tests
.\build.ps1 -Config Debug   # Debug build, run all tests
.\build.ps1 -Bench          # also print the benchmark tables
.\build.ps1 -Filter ssi     # run only matching tests
```

The library is header-only. It must be compiled with strict IEEE semantics (`/fp:strict`; `-ffp-contract=off -fno-fast-math` elsewhere), because the rigour of the interval and exact arithmetic depends on it. Only MSVC 19.51 on Windows x64 has been tested.

## Try it

```powershell
.\build\Release\ak_demo.exe demo.obj
```

This intersects two freeform patches and writes both surfaces and the certified intersection curve to an OBJ file that any mesh viewer can open.

```powershell
.\build\Release\ak_boolean_demo.exe build
```

This cuts a box that has a hill on top with a tilted box, then drills a square hole through the result, and writes `boolean_intersection.obj`, `boolean_union.obj`, `boolean_difference.obj` and `boolean_chained.obj` into the `build` folder. Each file is a closed triangle mesh (one group per face) with the model's edges as polylines; the program checks that every mesh is closed before writing it. If Windows refuses to run the Release program, use `.\build\Debug\ak_boolean_demo.exe build` after `.\build.ps1 -Config Debug`.

## Boolean operations

```cpp
#include <akernel/op/boolean.hpp>

using namespace ak;

const op::Brep a = op::to_brep(op::make_bump_box({0, 0, 0}, 1, 1, 1, 0.3));   // box with a hill on top
const op::Brep b = op::to_brep(op::make_box({0.5, 0.4, 0.3}, {1.5, 1.4, 1.5}));
const op::Brep c = op::to_brep(op::make_box({0.2, 0.2, 0.6}, {0.8, 0.7, 1.6}));

const op::BooleanResult r = op::boolean_op(a, b, op::Op::Difference);
if (!r.ok) {
    // Not solved, and r.error says which entities and why
    // (for example two faces that touch). Nothing is guessed.
}
// r.solid     the result: a Brep with trimmed faces, independent of a and b
// r.audit     the integer audit: closed, oriented, genus per shell
// r.manifold  false if parts of the result touch along an edge or at a point

// A result is an operand:
const op::BooleanResult r2 = op::boolean_op(r.solid, c, op::Op::Union);

const double v = op::volume(r2.solid);               // numerical, for checking
const double gap = op::edge_gap_bound(r2.solid);     // certified, along every intersection edge

// Certified test of a parameter point against a trimmed face:
op::FaceGeom fg;
std::string err;
op::face_geometry(r2.solid, 0, fg, err);
const op::Where w = op::locate(fg, num::Interval(0.3), num::Interval(0.7));
// Inside, Outside, Boundary, or Undecided - never a guess.
```

`op::imprint(a, b)` followed by `op::select(im, op)` gives all three operations from one intersection computation.

## Tessellation

```cpp
#include <akernel/op/tessellate.hpp>

op::TessOptions opt;
opt.chord = 1e-3;                                    // wanted deviation from the surfaces
const op::TessResult t = op::tessellate(r2.solid, opt);
if (t.ok) {
    const op::MeshAudit a = op::audit(t.mesh);       // a.watertight(), a.euler, a.components
    const double dev = op::mesh_deviation(t.mesh, r2.solid);   // measured, not bounded
    std::ofstream out("result.obj");
    op::write_obj(out, t.mesh);
}
```

The mesh is closed by construction, not by welding nearby points: every edge of the model is sampled once (side edges on their surface, intersection edges along their certified tube), and both faces along it use exactly those points as the boundary of their triangulations. Its Euler characteristic matches the genus from the model's integer audit.

## Comparing solids

```cpp
#include <akernel/op/compare.hpp>

const op::BooleanResult x = op::boolean_op(a, op::boolean_op(a, b, op::Op::Difference).solid, op::Op::Difference);
const op::BooleanResult y = op::boolean_op(a, b, op::Op::Intersection);
const op::SolidComparison c = op::compare_solids(x.solid, y.solid);   // c.same, c.reason
```

Two solids are the same when their vertices match one to one by certified enclosures, their surfaces are identical (or exactly coplanar planes with the same outward side), and their edges and faces correspond. Seams (two faces of one surface glued along an edge) and vertices where an edge merely continues are removed first, since two correct routes to the same solid may split it differently.

## Surface intersection

```cpp
#include <akernel/geom/primitives.hpp>
#include <akernel/solve/ssi.hpp>

using namespace ak;

const solve::SsiSystem sys(patch_a, patch_b);          // two geom::PolyPatchD
const solve::SsiResult res = solve::intersect(sys);

if (!res.complete) {
    // The kernel could not prove it found everything.
    // res.unverified_hull says where; a tangency is the usual cause.
}
for (const auto& curve : res.curves)
    for (const auto& seg : curve.segs) {
        solve::Box4 enc;                                // certified (u, v, s, t) enclosure
        solve::curve_point(sys, seg, 0.5, enc);
        auto on_a = sys.point_a(enc), on_b = sys.point_b(enc); // interval 3D points
    }
```

## Layout

```
include/akernel/
  num/    interval.hpp  expansion.hpp  predicates.hpp (orient2d, orient3d, incircle)  escalation_policy.hpp  dense.hpp
  geom/   vec.hpp  bernstein.hpp  bspline.hpp  nurbs_surface.hpp  continuity.hpp  seam_certificate.hpp  primitives.hpp  box_index.hpp
  solve/  ssi.hpp  csi.hpp
  var/    fill.hpp
  topo/   complex.hpp
  op/     patch_solid.hpp  brep.hpp  trim.hpp  planar.hpp  surfaces.hpp  contact.hpp  tsi.hpp  boolean.hpp  tessellate.hpp  compare.hpp
  diag/   rate.hpp
tests/    one executable per module (21)
bench/    bench_main.cpp   (source of docs/RESULTS.md)
examples/ intersect_demo.cpp  boolean_demo.cpp
```
