# akernel — measured results

Measured on 2026-11-11, 2026-11-12 (tessellation, the vertex-manifoldness check, algebraic identities, non-manifold results, exact quadrics, curved chains) 2026-11-13 (spatial index, plane pairs, arc searches, tangential contact, flat faces on rational carriers), 2026-11-14 (contact cases left after the flat-face work) 2026-11-15 (exact points, hill-top crossings) and 2026-11-16 (roadmap phases 1 and 2: certified seam bounds, the faster interval step, tangent points) with `ak_bench.exe` and the test suite.

- Build: every number and timing below is from the Release benchmark (`/O2`) unless a row says otherwise. Rows taken from the test programs are marked as such.
- Setup: Windows 11 x64, MSVC 19.51, `/fp:strict`, single thread, fixed random seeds.
- Reproduce: `.\build.ps1 -Bench`. `.\build\Release\ak_bench.exe chains 40 6` runs only the chained-Boolean section, and `curved 40 6` only the chains with curved operands. A last argument `noexact` runs with exact points switched off, and `nolines` without the direct certification of plane pairs (§8 says which hashes each reproduces). Counts and the result hashes are deterministic; timings vary by a few percent between runs.
- Scope: these numbers describe the certified core and Boolean operations on solids with trimmed faces in general position, including chains of operations. No other kernel was measured, so nothing here is a comparison with OCCT, Parasolid, ACIS or CGM.

## Test suite status

Twenty-two test programs, one per module.

| Build | Result on the final sources |
|---|---|
| Release | 22 of 22 pass |
| Debug | 22 of 22 pass |

On this machine Windows refuses some freshly built executables ("An Application Control policy has blocked this file"), in either build, and which file is refused changes from build to build. That is an operating-system policy decision, not a test failure. Every test program has run and passed on the final sources.

## 1. Orientation predicates

One million `orient2d` calls per row.

| Input | Decided by double filter | by intervals | by exact arithmetic | Plain double gives wrong sign | Plain double inconsistent under rotation of arguments | Cascade, ns/call | Plain double, ns/call |
|---|---|---|---|---|---|---|---|
| Uniform random points | 1,000,000 | 0 | 0 | 0 | 0 | 10 | 6 |
| Third point within 2 ulp of the line | 743,917 | 9 | 256,074 | 37,660 | 139,438 | 214 | 6 |

- On ordinary input the cascade costs about 1.7 times a plain evaluation and never leaves the first tier.
- On near-degenerate input plain doubles give the wrong sign in 3.8 % of calls and contradict themselves in 13.9 %. The cascade is exact on all of them, at about 35 times the cost of a plain evaluation. The exact tier allocates memory and has not been optimised.
- On the 64 × 64 one-ulp grid of Kettner et al., plain doubles are wrong on 1,442 of 4,096 points; the cascade is right on all, including the 64 exactly collinear ones (`tests/num/predicates_test.cpp`).

## 2. Certified surface/surface intersection

"Proven complete" means the kernel proved that every intersection point in the domain lies on the certified arc of a returned segment. "Certified gap bound" is a rigorous upper bound on the distance between the two surfaces, evaluated at the certified curve points.

| Case | Curves | Segments | Proven complete | Unresolved boxes | Lowest sin(angle between normals) | Certified gap bound | Time, ms |
|---|---|---|---|---|---|---|---|
| Plane × quarter cylinder | 1 | 22 | yes | 0 | 0.958 | 6.6e-14 | 6 |
| Plane × hill (closed loop) | 1 | 45 | yes | 0 | 0.908 | 8.8e-14 | 26 |
| Plane above hill (no intersection) | 0 | 0 | yes (proven empty) | 0 | n/a | n/a | 0.1 |
| Cylinder × cylinder, radii 1 and 0.8 | 1 | 18 | yes | 0 | 0.600 | 7.7e-14 | 5 |
| Bicubic × bicubic | 1 | 24 | yes | 0 | 0.586 | 2.8e-13 | 25 |
| Plane tangent to hill (one contact point) | 0 | 0 | **no** | 16 | n/a | n/a | 16 |

The tangent case is the intended behaviour. The kernel does not answer "no intersection"; it reports that it could not decide and returns the unresolved region, u and v in [0.484, 0.516], which contains the true contact point (0.5, 0.5).

### Random samples (200 pairs each, fixed seeds)

| Family | Proven complete | Not proven | Curves | Segments | Lowest sin(angle) | Worst certified gap bound | Mean time, ms | Worst time, ms |
|---|---|---|---|---|---|---|---|---|
| Transversal bicubic pairs | 200 | 0 | 200 | 3,710 | 0.889 | 3.6e-13 | 29 | 40 |
| Two random graphs over the same square (loops, several components, near-tangencies) | 200 | 0 | 403 (4 closed) | 12,044 | 0.022 | 7.9e-12 | 161 | 1,462 |

**Independent cross-check.** For the second family the completeness claim was tested without the kernel's own machinery. Along six fixed lines u = const per pair, the roots of f − g were counted by dense sampling (20,001 samples) and compared with the number of times the returned curves cross that line. The counts agree on all 1,200 lines.

**How the second family got to 200 of 200.** The first version proved only 122 of 200. Two fixes, both found by measuring:

1. Axis-aligned bounding boxes cannot separate two nearly parallel, close surface pieces. Adding a separating-axis test along the patch normals (still rigorous) raised the count to 188 and cut the worst time from 4.5 s to 0.6 s.
2. A very short last step onto a domain boundary, or a short step near a tangency, left a certified region too small for the completeness pass. Giving regions a minimum size where the certificate allows, and letting the pass go deeper, proved the remaining 12.

Two of those 12 are kept as regression tests (`tests/solve/ssi_test.cpp`, cases 8 and 9).

### What a fitted edge would give

For the curated cases, the certified junction points of each curve were also joined by interpolating curves in the two parameter domains, as a kernel storing p-curves does, and the distance between the two surfaces was measured at the mid-parameters.

| Case | Certified gap bound | Gap of a cubic interpolant through the same points | Gap of a linear interpolant |
|---|---|---|---|
| Plane × quarter cylinder | 6.6e-14 | 8.6e-04 | 2.3e-03 |
| Plane × hill | 8.8e-14 | 2.8e-03 | 6.2e-03 |
| Cylinder × cylinder | 7.7e-14 | 1.0e-03 | 1.9e-03 |
| Bicubic × bicubic | 2.8e-13 | 9.7e-05 | 2.2e-04 |

This is an illustration of the mechanism, not a measurement of any legacy kernel. A production kernel adds points until the fit meets its tolerance (typically 1e-6 to 1e-8). The point is that a fitted pair of p-curves always has a nonzero gap that must be covered by a tolerance, whereas the certified edge can be evaluated to rounding level anywhere along it.

## 3. Boolean operations

Operands are closed shells whose faces are trimmed Bezier patches, and the result is a solid of the same kind, so a result can be an operand (chained operations, below). In the cases of this first part both operands are solids whose faces are whole patches. One imprint produces all three results (A and B, A or B, A minus B). No surface is tessellated and no tolerance is used.

Checks applied to every result:

- **Audit.** The integer audit of the cell complex: loops close, every edge is used exactly twice in opposite directions, the faces at every vertex form one fan (no pinch point), and each shell has an integer genus.
- **Inclusion-exclusion.** vol(A and B) + vol(A or B) = vol(A) + vol(B), and vol(A minus B) = vol(A) - vol(A and B). The three volumes come from three separately assembled results, so a wrong face selection or a wrong trim loop breaks the identity.
- **Reference volume.** A value the kernel did not compute: a closed form, a plain 2D quadrature, or a Monte Carlo count using only the analytic definitions of the two solids.

Volumes are numerical quadrature in doubles. They test the result; they are not certified bounds.

### Curated cases

| Case | Curves | New vertices | Volume of A and B | Reference | Difference | Inclusion-exclusion residual | Certified edge gap bound | Widest vertex enclosure |
|---|---|---|---|---|---|---|---|---|
| Two offset boxes | 6 | 6 | 0.210000000000 | 0.21 (closed form) | 6.1e-16 | 1.9e-15 | 7.2e-14 | 7.2e-15 |
| Hill cut by a plane (one closed loop, no vertex) | 1 | 0 | 0.030901244622 | 0.030901245191 (2D midpoint quadrature) | 5.7e-10 | 1.0e-14 | 2.8e-13 | n/a |
| Cylindrical shell cut by a box (rational faces) | 4 | 4 | 0.297329064872 | 0.297329064872 (closed form) | 5.6e-17 | 1.0e-15 | 9.8e-14 | 1.2e-14 |
| Box inside a box (no contact) | 0 | 0 | 0.027000000000 | 0.027 (closed form) | 0 | 4.4e-16 | n/a | n/a |
| Box against a rotated box | 12 | 12 | 0.370206 | 0.369885 ± 0.000341 (Monte Carlo, 2,000,000 points) | 0.9 sigma | below 1e-10 | 1.5e-13 | 7.9e-14 |
| Hill solid against the cylindrical shell (curved against curved) | 10 | 10 | 0.718806 | 0.718570 ± 0.000457 (Monte Carlo, 2,000,000 points) | 0.5 sigma | below 1e-10 | 5.2e-13 | 3.1e-14 |

- All three results of every case pass the audit with genus 0. "Box inside a box" gives a difference with two shells (a solid with a cavity).
- The 5.7e-10 in the second row is the accuracy of the reference quadrature, not of the kernel.
- "Certified edge gap bound" is a rigorous bound on the distance between the two surfaces at one certified point per segment of every new edge. "Widest vertex enclosure" is the width, in parameter units, of the certified box around a point where an edge pierces a face.
- The first four rows are from the benchmark, the last two from `tests/op/boolean_test.cpp`, which only asserts the residual to 1e-10.

### Random pairs (100 per family, fixed seeds)

A fixed solid A against a box of random orientation, size (0.4 to 1.1 per side) and position. Each pair is checked against a Monte Carlo count of 100,000 points.

| A | Solved | All three results pass the audit | Curves | New vertices | Worst inclusion-exclusion residual | Worst Monte Carlo deviation | Pairs beyond 3 sigma | Worst certified edge gap bound | Widest vertex enclosure |
|---|---|---|---|---|---|---|---|---|---|
| Unit box | 100 of 100 | 100 | 1,102 | 1,102 | 1.3e-15 | 2.86 sigma | 0 | 6.8e-13 | 2.1e-13 |
| Box with a hill on top | 100 of 100 | 100 | 896 | 893 | 1.7e-15 | 2.71 sigma | 0 | 1.4e-12 | 8.3e-13 |
| Quarter cylindrical shell | 100 of 100 | 100 | 810 | 810 | 4.9e-15 | 3.19 sigma | 1 | 2.0e-12 | 7.3e-13 |

- All 300 pairs were solved, and all 900 results are closed, oriented shells.
- Among 300 independent Monte Carlo comparisons, about 0.8 are expected beyond 3 sigma by chance. One was observed, at 3.19 sigma.
- Every random pair has a planar-faced box as one operand. The only test with two curved operands is the last curated row.
- Random position does not produce contact cases. Those are below.

Time for one imprint plus all three results:

| A | Mean, ms | Worst, ms |
|---|---|---|
| Unit box | 72 | 711 |
| Box with a hill on top | 96 | 453 |
| Quarter cylindrical shell | 68 | 395 |

These operands have six faces, so almost all of this time is the certified intersection of surface pairs that do meet. How the time grows with the number of faces is measured under "Scaling" below.

### Chained operations

A result is used as an operand of the next operation. A new kind of vertex appears here: where a surface of the new operand cuts an edge made by an earlier operation, the vertex lies on three surfaces at once, and is found as a certified root of six equations in six unknowns.

**Curated chains** (`tests/op/boolean_test.cpp`; A, B, C, E are axis-aligned boxes, so every volume has a closed form):

| Chain | Volume | Closed form | Three-surface vertices |
|---|---|---|---|
| (A − B) − C | 0.706000000000 | 0.706 | 2 |
| (A ∪ B) ∩ C | 0.147000000000 | 0.147 | 2 |
| (A ∩ B) ∪ C | 0.444000000000 | 0.444 | 2 |
| (A − B) ∩ (C − E): both operands are results | 0.084000000000 | 0.084 | 2 |
| (hill − slab) − box: the box crosses the closed edge left by the first cut | 1.056638 | 1.056929 ± 0.000358 (Monte Carlo) | 4 |

A fixed chain of five operations on boxes (difference, union, difference, intersection, union) is checked against Monte Carlo at every step. At step 3 a pocket meets a notch and opens a tunnel; the audit reports genus 1 from then on, which is correct. My first version of the test assumed genus 0 throughout and failed on this.

**Random chains.** 40 chains start from the unit box. Each step applies a random operation (difference 45 %, union 35 %, intersection 20 %) with a box of random orientation, size and position, and the result becomes the next operand. Every step is compared with a Monte Carlo count of 200,000 points that uses only the analytic definition of the chain so far.

| Step | Attempted | Solved | Result passes the audit | Mean faces | Mean edges | Of those, intersection edges | Three-surface vertices | Worst Monte Carlo deviation | Beyond 3 sigma | Worst certified edge gap bound | Widest vertex enclosure | Mean time, ms | Worst time, ms |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 40 | 40 | 40 | 12.8 | 31.6 | 11.4 | 0 | 2.02 sigma | 0 | 3.1e-13 | 2.2e-13 | 78 | 478 |
| 2 | 40 | 40 | 40 | 19.6 | 52.0 | 24.8 | 100 | 2.58 sigma | 0 | 7.7e-13 | 3.6e-13 | 128 | 467 |
| 3 | 40 | 40 | 40 | 25.9 | 71.1 | 39.1 | 200 | 2.21 sigma | 0 | 7.7e-13 | 5.0e-13 | 192 | 570 |
| 4 | 40 | 40 | 40 | 31.1 | 86.5 | 51.3 | 289 | 1.99 sigma | 0 | 7.7e-13 | 7.2e-13 | 255 | 820 |
| 5 | 39 | 39 | 39 | 34.1 | 92.6 | 57.4 | 331 | 3.18 sigma | 1 | 7.7e-13 | 4.2e-13 | 271 | 634 |
| 6 | 39 | 39 | 39 | 37.1 | 101.7 | 65.9 | 327 | 3.05 sigma | 1 | 7.7e-13 | 4.6e-13 | 284 | 760 |

- All 238 steps were solved and all 238 results pass the audit. One chain ended at step 4 because its result was empty (an intersection with a box that missed the solid), so 39 chains continued.
- **No growth along the chains.** The certified gap between the two surfaces along any intersection edge stays at or below 7.7e-13 from step 2 to step 6, and the widest certified vertex enclosure stays below 7.2e-13. Nothing accumulates, because a later operation never re-derives geometry from an earlier result: it reuses the certified tubes and vertex enclosures, which still refer to the original surfaces.
- Of 237 Monte Carlo comparisons (the empty result has none), about 0.6 are expected beyond 3 sigma by chance; two were observed, at 3.18 and 3.05 sigma. The comparisons along one chain are not independent of each other, since each step builds on the previous solid.

Time per stage at step 6 (means, ms): surface/surface intersection 160, edge/surface roots 116, everything else under 10. 64 of 140 surface pairs had to be intersected.

**How the random chains got to 238 of 238.** The first run of this benchmark solved 213 of 238 steps (25 not solved) and took up to 89 s per step in Debug. Three causes, all found by measuring:

1. The search for three-surface vertices exhausted its budget of 200,000 nodes. Once a box had been proven to hold exactly one root, other boxes inside it were still searched again whenever that root was already known. Recording every proven box fixed it and cut the time of that stage at step 4 from 33 s to 0.45 s (Debug).
2. A face that touches nothing in the other solid needs one interior point for a ray test, and a 9 × 9 grid search sometimes found none on a narrow face. Any vertex of such a face serves, so the search was removed.
3. A real bug, which ended in a report, not a wrong solid. Each piece of an intersection edge is certified in a box that overlaps its neighbours and reaches a little past the two ends of the curve. A new surface cut the polynomial continuation of an edge beyond the edge's end, outside the parameter square of one of its surfaces. The test "is this point between the edge's two ends" compared only which piece each point was in, and accepted it. The kernel then found that no intersection curve ended at the new vertex and stopped. Positions are now read inside each piece's own range. A regression test places 42 points on the continuation beyond the two ends of a curve, each inside a certified region, and requires all of them to be reported as off the curve.

### Operands that share geometry by construction

A result combined with one of its own operands shares surfaces, edges and vertices with it, as does a solid combined with itself. Sharing is recognised from construction, never from distance: identical control points for surfaces and shared sides, identical certified parameters for vertices, and the same curve between the same two vertices for edge pieces. In a shared surface the faces of both solids are overlaid, and each region is classified as inside, outside, or on the other solid's boundary with the same or the opposite outward normal.

**Curated cases** (`tests/op/boolean_test.cpp`):

| Case | Shared surfaces | Union | Intersection | Difference | Reference |
|---|---|---|---|---|---|
| (A − B) with B | 3 | 1.790000000000 | empty | 0.790000000000 | A ∪ B = 1.79; A − B = 0.79 (closed forms) |
| A with A | 6 | 1.000000000000 | 1.000000000000 | empty | A = 1 |
| (A − R) with R, R a rotated box | all of R's | 1.149778256066 | empty | 0.724350223664 | A ∪ R and A − R computed directly, equal to 1e-11 |
| (hill − slab) with slab: a curved surface and a closed edge shared | 1 | 4.502432088712 | empty | 1.102432088712 | hill ∪ slab and hill − slab computed directly, equal to 1e-11 |

In (A − B) with B, no new vertex is created and no new edge is traced: every vertex and edge of the result already existed in one of the operands. Re-applying the cut to the union, ((A − B) ∪ B) − B, gives 0.79 again.

**Random chains with reused operands.** 40 chains as above, but at each step after the first, with probability 40 %, the box is one already used earlier in the same chain.

| Step | Attempted | Solved | Result passes the audit | Empty result | Mean faces | Worst Monte Carlo deviation | Beyond 3 sigma | Worst certified edge gap bound | Mean time, ms | Worst time, ms |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 40 | 40 | 40 | 0 | 12.8 | 1.93 sigma | 0 | 2.9e-13 | 94 | 824 |
| 2 | 40 | 40 | 40 | 0 | 16.2 | 2.57 sigma | 0 | 3.1e-13 | 173 | 867 |
| 3 | 40 | 40 | 40 | 4 | 22.1 | 2.28 sigma | 0 | 3.4e-13 | 187 | 650 |
| 4 | 36 | 36 | 36 | 2 | 26.8 | 2.03 sigma | 0 | 6.9e-13 | 204 | 663 |
| 5 | 34 | 34 | 34 | 1 | 28.6 | 2.17 sigma | 0 | 4.2e-13 | 231 | 536 |
| 6 | 33 | 33 | 33 | 1 | 28.3 | 2.48 sigma | 0 | 4.8e-13 | 229 | 531 |

- All 223 steps were solved and all 223 results pass the audit. Reused operands mostly cut away or add back material that is already gone or already there, so results are often empty; each empty result ends its chain, which is why fewer steps are attempted than in the plain chains.
- 80 steps reused an earlier box; in 79 of the solved steps the two operands shared at least one surface.
- None of the 215 Monte Carlo comparisons is beyond 3 sigma.

**How the reuse chains got to 223 of 223.** The first run solved 228 of 229 steps, and one solved step produced a result that the audit rejected (14 edges used by more than two faces). Both were real bugs; both were caught by the kernel's own checks, not returned as wrong solids.

1. *Wrong classification across a seam.* A face can be split along a curve where an earlier operation cut it, leaving two faces on the same surface with an edge between them (a seam), and that edge can lie on a surface of the other operand. Classification spread "outside" from one of the two faces to the other across the seam, but the seam was on the other solid's boundary, so one face was in fact inside. The audit rejected the result. Classification now never crosses an edge piece that lies on the other solid's boundary.
2. *A missing edge in a face's map.* An edge of one solid can lie on a surface that carries no face of that solid (a seam again). The other solid's face on that surface needs the edge in its planar map, but edges were only inserted into faces on surfaces both solids have faces on. The map had a dangling curve, and the kernel reported "the planar map is not closed". Edges are now inserted wherever they lie on the face's surface.

Two smaller issues came up while writing this part. Directions along the sides of a parameter square are exact, but interval multiplication widened their exact zero cross product by one unit, so the order of curves at a vertex was reported as undecidable; those products are now computed exactly. And an edge split at a vertex continues straight through it: its two directions there are opposite by construction, which interval tangents alone cannot prove, so curves on the same pair of surfaces are now treated as continuing each other.

### Coincidence decided exactly: flat faces in contact

Two boxes placed face to face each have a face in the plane where they touch, but the two faces are different patches with different parameterisations and extents: nothing is shared by construction. The kernel decides exactly, from the control points, that a patch is a flat parallelogram and that two such patches lie in one plane (exact `orient3d`), and maps parameters between them with interval-certified affine maps. The overlay of the previous section then works across those maps.

**Curated cases** (`tests/op/boolean_test.cpp`; A is the unit cube):

| B | Exactly coplanar surface pairs | A ∩ B | A ∪ B | A − B |
|---|---|---|---|---|
| A unit cube face to face with A | 5 | empty | 2.000000000000 | 1.000000000000 |
| Face to face, B smaller and offset | 1 | empty | 1.660000000000 | 1.000000000000 |
| Stacked on top, smaller | 1 | empty | 1.125000000000 | 1.000000000000 |
| Overlapping, two face planes shared | 2 | 0.250000000000 | 1.250000000000 | 0.750000000000 |
| Inside A, touching its bottom face | 1 | 0.125000000000 | 1.000000000000 | 0.875000000000 |
| The unit cube again, built mirrored (every face parameterised the other way) | 6 | 1.000000000000 | 1.000000000000 | empty |

Every value is the closed form to 1e-12, and every non-empty result is one closed shell. A row of three unit cubes joined by two unions and then notched across both seams has volume 2.500000000000 (closed form 2.5).

**Random chains of grid-aligned boxes.** 40 chains as above, but every box is axis-aligned with its corners on a grid of spacing 0.25 (sides 0.25 to 1). Contact is the rule: faces in common planes, boxes touching face to face, edges overlapping or crossing inside shared planes.

| Step | Attempted | Solved | Result passes the audit | Empty result | Mean faces | Worst Monte Carlo deviation | Beyond 3 sigma | Mean time, ms | Worst time, ms |
|---|---|---|---|---|---|---|---|---|---|
| 1 | 40 | 40 | 40 | 5 | 9.6 | 2.58 sigma | 0 | 22 | 101 |
| 2 | 35 | 35 | 35 | 7 | 12.2 | 2.56 sigma | 0 | 19 | 86 |
| 3 | 28 | 28 | 28 | 4 | 15.5 | 3.01 sigma | 1 | 29 | 106 |
| 4 | 24 | 24 | 24 | 5 | 17.9 | 1.98 sigma | 0 | 28 | 118 |
| 5 | 19 | 19 | 19 | 2 | 19.1 | 2.15 sigma | 0 | 33 | 144 |
| 6 | 17 | 17 | 17 | 0 | 21.0 | 1.97 sigma | 0 | 27 | 93 |

- All 163 steps were solved and all 163 results pass the audit; in 146 steps the two operands had faces in a common plane.
- 17 committed results are not manifold: parts that touch along an edge (four faces there) or at a point (a pinched vertex). The chains continue with them as operands, and every later step still agrees with Monte Carlo.
- One of 140 Monte Carlo comparisons is beyond 3 sigma (3.01 sigma); about 0.4 are expected by chance.
- The faces are planes, so the steps are fast: 19 to 33 ms on average.

Before non-manifold results were supported, the same chains stopped at 5 refused steps (§8 has the earlier hash): 4 unions of boxes touching along an edge, and 1 touching at a corner.

**A pinched result found by the tessellation.** Until the tessellation milestone, the corner-contact result of chain 30 was committed: the audit checked edges, loops and per-shell Euler characteristics, and each of the two shells was a valid genus-0 sphere on its own. The tessellation (§3, below) counted the shared corner once and reported a mesh Euler characteristic of 3 against the Brep's 4. The integer audit now also checks vertex manifoldness: the face corners at each vertex, joined through the edges at that vertex, must form one fan. At that point the result was refused, so that chain stopped one step earlier (147 steps attempted before, 145 after; §8). No other result in any section was affected. Since the next milestone (below), such a result is accepted as a non-manifold solid, with its pinched vertex counted correctly.

### Non-manifold results: parts touching along an edge or at a point

**Definitions.** A union of two boxes that touch only along an edge has four faces at that edge. A union of boxes touching only at a corner has a vertex where the faces form two separate fans. Both are valid solids; they are just not manifold. The integer audit now separates two properties:
- **Closed and oriented:** every loop closes, every edge is used an even number of times and as often in each direction, and the fans at every vertex close.
- **Manifold:** no edge with more than two faces and no vertex with more than one fan.

A non-manifold result is accepted only if it is closed and passes a geometric check, `op::check_edge_fans`. At a point of each edge with more than two faces, the faces' interior directions are sorted by angle around the edge. Consecutive faces must then bound material and empty space alternately: the uses' effective orientations must alternate. Faces at the same angle (tangential contact) are refused. The integer audit alone is not enough. Turning one part of a valid edge-contact solid inside out keeps every edge balanced, and the fan check refuses it (`tests/op/nonmanifold_test.cpp`). An earlier bug in this kernel (§3, reuse chains) also produced edges with more than two faces.

**Using non-manifold solids.** They work as operands. Classification across an operand's own edge was the one place that assumed two faces per edge. Its rule holds for any number of faces: a piece that is not on the other solid's boundary has all regions around it on the same side of that solid. Where such an edge pierces a face of the other solid, one new curve ends per face of the edge, so four instead of two.

Volume, tessellation and `op::compare_solids` handle them too. A mesh edge may be used by four triangles, two in each direction. The comparison treats an edge's faces as a list. The mesh and the Brep agree on the Euler characteristic of the whole complex.

**Curated cases** (`tests/op/nonmanifold_test.cpp`; A is the unit cube; every volume is the closed form to 1e-12):

| Case | Volume | Manifold | Shells | Edges with four faces | Pinched vertices | Euler characteristic (Brep = mesh) |
|---|---|---|---|---|---|---|
| A ∪ box touching it along a whole edge | 2.0 | no | 1 | 1 | 0 | 3 |
| A ∪ box touching it along part of an edge | 1.5 | no | 1 | 1 | 0 | 3 |
| A ∪ box touching it at a corner | 2.0 | no | 2 | 0 | 1 | 3 |
| (edge pair) − slot across the shared edge | 1.9 | no | 1 | 2 | 0 | 2 |
| (edge pair) ∪ block across the shared edge | 2.1 | no | 1 | 2 | 0 | 2 |
| (edge pair) ∩ block across the shared edge | 0.1 | no | 1 | 1 | 0 | 3 |
| (edge pair) − column around the shared edge | 1.92 | yes | 2 | 0 | 0 | 4 |
| slab − (edge pair) | 3.5 | no | 1 | 1 | 0 | −1 |
| (corner pair) ∪ box over the corner | 2.09375 | yes | 1 | 0 | 0 | 2 |
| (corner pair) − box away from the corner | 1.992 | no | 3 | 0 | 1 | 5 |

Commutativity with a non-manifold operand gives the same solid by `op::compare_solids`.

In the benchmark, the grid-aligned chains no longer stop at contact along an edge or at a point. All 163 steps are solved and committed, 17 of them non-manifold, and all agree with Monte Carlo (above). All 130 grid identities are now evaluated, and all give the same solid (below).

**How the grid chains got to 147 of 147.** The first run solved 111 of 145 steps. Each failure was a report, not a wrong answer:

1. Two planes whose overlap is a single point on the boundary of both parameter squares give tracing nothing to follow, so completeness could not be proven. Two planes meet in exactly one line, which is now certified directly as a tube; this is used only when tracing cannot finish, so no earlier result changed.
2. A vertex could carry parameters on a plane only if that plane held a face. An edge can lie on a plane with no face left (a seam), and then the crossing points on it had no parameters there. All planar surfaces now take part.
3. Two vertices at the same grid point, one with exact corner coordinates and one with interval coordinates from an earlier operation, were not recognised as one. Now: two vertices that lie, by construction, on the same three planes whose normals are provably independent are the same point.
4. A crossing of two edges inside a plane was accepted only if both solids had faces in that plane; the edge's own solid may have none there (a seam again).

The work also exposed a soundness gap in point location that predates this milestone. A straight boundary piece was compared against the chord between the midpoints of its end enclosures; when those enclosures are not points, the chord is not exactly the piece, and a query point inside the piece's enclosure box was being decided by the chord. Such a point is now undecided and a different point is tried. No earlier result changed (§8).

### Watertight tessellation

`op::tessellate` turns any solid into a triangle mesh. Each edge is sampled once: side edges on their surface, intersection edges along their certified tube, adaptively until the chord deviation is below the requested tolerance. Both faces along an edge use exactly those points as the boundary of their triangulations. Each face is triangulated in its parameter plane: its loops are the constraint polygon, interior points come from a grid whose spacing follows the surface's second derivatives, and a constrained Delaunay triangulation joins them. The triangulation's decisions use exact `orient2d` and `incircle` (a new predicate, filtered like the others). No point is welded by distance.

**Checks on every mesh:**

- **Closed.** Every mesh edge is used exactly once in each direction, and no triangle is degenerate.
- **Topology agrees with the Brep.** The mesh's Euler characteristic and number of components equal those from the Brep's integer audit.
- **Deviation.** Largest distance from triangle centroids and side midpoints to the surface (closest point by Newton). This is a measurement, not a bound.
- **Volume.** Divergence-theorem volume of the mesh against the Brep's quadrature volume.

**Every committed result of the benchmark** (§5 to §8 of `ak_bench`), chord 1e-3:

| Results of | Meshes | Not tessellated | Closed | Euler characteristic matches | Components match | Mean triangles | Most triangles | Worst deviation / chord | Worst volume error / (chord × area) | Mean time, ms | Worst time, ms |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Random boxes against the box, the hill and the cylindrical shell | 900 | 0 | 900 | 900 | 900 | 3,829 | 12,558 | 0.60 | 0.05 | 19.8 | 77.2 |
| Random chains | 237 | 0 | 237 | 237 | 237 | 723 | 2,252 | 0.00 | 0.00 | 1.7 | 5.0 |
| Chains with reused operands | 215 | 0 | 215 | 215 | 215 | 590 | 1,958 | 0.00 | 0.00 | 1.4 | 4.2 |
| Grid-aligned chains | 140 | 0 | 140 | 140 | 140 | 113 | 452 | 0.00 | 0.00 | 0.8 | 2.2 |

All 1,492 meshes are closed and agree with the Brep's topology. 15 grid-chain meshes are not manifold (parts touching along an edge, so four triangles at a mesh edge). The other 2 non-manifold results touch only at a point, so their meshes have no such edge. In a non-manifold mesh, "closed" means every edge is used as often in each direction. Box results are flat, so their deviation is rounding only (below 1e-15).

**Convergence** (`tests/op/tessellate_test.cpp`):

| Solid | Chord | Triangles | Deviation | Volume error |
|---|---|---|---|---|
| Box with a hill on top | 1e-2 | 582 | 4.3e-03 | 1.1e-03 |
| | 1e-3 | 4,868 | 4.8e-04 | 1.1e-04 |
| | 1e-4 | 48,460 | 5.0e-05 | 1.1e-05 |
| Quarter cylindrical shell | 1e-2 | 1,304 | 2.5e-03 | 1.7e-03 |
| | 1e-3 | 12,410 | 2.5e-04 | 1.8e-04 |
| | 1e-4 | 120,886 | 2.5e-05 | 1.8e-05 |

Deviation and volume error fall in proportion to the chord, and the triangle count grows tenfold per decade: a second-order method. The test also covers a face bounded by a closed intersection loop, a solid with a cavity (2 components), the drilled demo solid (genus 1, Euler characteristic 0), faces in contact, and a length limit on triangle sides. It checks that two runs give bitwise identical meshes.

**What came up while building it:**

1. *A closed curve's seam.* The last piece of a closed intersection curve overlaps the first, and the curve model closes the loop with a short span travelled backwards over the overlap. For point classification this is harmless, because a back-and-forth piece adds nothing to a winding number. A polyline through those points doubles back and crosses itself, and the triangulator reported that. Samples that overshoot the seam are now dropped.
2. *Long triangles on singly curved faces.* On a cylinder the grid was fine across the curved direction and had one line along the straight one, so near corners the Delaunay triangles spanned a large step of the curved direction (deviation 8.8e-3 at any chord). The grid now uses one step size in both directions of the (power-of-two scaled) parameter plane.
3. *A deviation measure that was wrong.* Comparing a mesh point with the surface point at the same parameters overstates the error on rational faces, even on a flat annulus whose triangles lie exactly in its plane. It was replaced by the distance to the closest surface point.
4. *A pinched result* (above): two shells sharing one corner vertex had been committed. The audit now checks vertex manifoldness.

### Algebraic identities

Blueprint §5.4 asks that identities between Boolean expressions give identical topology, not just equal volumes. Each identity below is computed along both of its routes, every operation with a fresh imprint, and the two results are compared by `op::compare_solids`:

- **Vertices** must match one to one. The match is by certified 3D enclosures, which come from the vertices' certified parameters by interval evaluation.
- **Surfaces** must be identical, or exactly coplanar flat faces with the same outward side.
- **Edges and faces** must correspond.
- **Allowed differences.** Seams (two faces of one surface glued along an edge) and vertices where an edge merely continues are removed first, since two correct routes may split the same solid differently. The glued cube in `tests/op/identity_test.cpp` is an example: 10 faces, 4 seams and 4 such vertices, and the same solid as the unit cube.
- **Volumes** must agree as well.

Identities, for pairs: A ∩ B = B ∩ A; A ∪ B = B ∪ A; A − (A − B) = A ∩ B; (A − B) ∪ (A ∩ B) = A; A ∪ (A ∩ B) = A; A ∩ (A ∪ B) = A; (A ∪ B) − B = A − B. For triples: (A ∪ B) ∪ C = A ∪ (B ∪ C); (A ∩ B) ∩ C = A ∩ (B ∩ C); A − (B ∪ C) = (A − B) − C.

Most of these combine a result with one of its own operands, so they exercise the handling of geometry shared by construction on every step.

| Operands (`ak_bench`, §10) | Pairs | Triples | Identities evaluated | Same solid | Different | Not evaluated | Worst volume difference | Time, s |
|---|---|---|---|---|---|---|---|---|
| Unit box, random boxes | 16 | 6 | 130 | 130 | 0 | 0 | 6.7e-16 | 48 |
| Box with a hill, random boxes | 16 | 6 | 130 | 130 | 0 | 0 | 1.1e-15 | 70 |
| Cylindrical shell, random boxes | 16 | 6 | 130 | 130 | 0 | 0 | 4.4e-15 | 57 |
| Grid-aligned boxes (contact) | 16 | 6 | 130 | 130 | 0 | 0 | 1.1e-15 | 7 |

`tests/op/identity_test.cpp` adds curated pairs and triples: offset boxes, the hill and slab (closed loop), the cylindrical shell (rational faces), the hill against the shell (curved against curved), a rotated box, and boxes face to face and overlapping. It also adds 18 more random pairs and 6 triples, and the empty solid as an operand. All 202 identities it evaluates give the same solid.

- **Different: none.** Every identity that could be evaluated holds with identical topology.
- **Not evaluated: none.** Before non-manifold results were supported, 3 grid identities could not be evaluated: in one trial the two boxes touch along an edge.

**What the identities exposed.** Before the fixes below, 7 more identities could not be evaluated: 3 in the test, 4 in the benchmark. Each was a report, not a wrong solid:

1. *An empty operand was refused.* A ∩ B is empty for boxes that only touch, so (A − B) ∪ (A ∩ B) failed on its last step. `boolean_op` now treats an operand without faces as the empty set (A ∪ ∅ = A, A ∩ ∅ = ∅, A − ∅ = A, ∅ − B = ∅).
2. *Two closed edges on one curve were refused.* A − B and A ∩ B of the hill and the slab both carry the closed intersection loop as an edge, with one vertex. Pieces of the two operands on the same curve are unified, and the direction of one relative to the other was read from their end vertices, which says nothing when both ends are the same vertex. The direction is now decided separately:
   - pieces of one certified curve run the same way;
   - otherwise the curve's tangent is enclosed by intervals at a certified point and its sign against each tube's slicing direction is compared.
   The identity uses separate imprints for its two parts, so the second rule was the one exercised.
3. *The three-surface search ran out of depth beside a known root.* An intersection edge of A ∩ B ends at a vertex on the edge between the hill's top and one of its sides; there the side, the top and the box face meet. Combining that result with the hill again intersects the edge with the hill's top, which the edge leaves at about 3°. The root at the vertex was found and proven. But the boxes right beside it (within 3e-5) could be neither excluded nor certified before the depth limit of 48. Now, once a root is proven unique in a box, boxes of doubling size around it are tested with the same Krawczyk operator, and every one that passes is recorded as covering only that root. This is sound: a covered search box holds no root other than one already found. The edge/surface search received the same change. No earlier result changed (§8).

### Exact quadrics and coincident curved surfaces

**Exact quadric patches.** A quadric is x^T A x = 0 in homogeneous coordinates. `op::on_quadric` decides exactly whether a rational patch lies on one: all Bernstein coefficients of A(P, P) must vanish, computed with expansion arithmetic. The circle has the exact rational parameterisation (1 − t², 2t, 1 + t²), whose Bezier data on any interval with dyadic ends are doubles. `op::make_cylinder_sector(r0, r1, t0, t1, z0, z1)` builds cylindrical sectors from it, with faces lying exactly on their cylinders and planes (`tests/op/quadric_test.cpp`):

| Patch | On its quadric (exact) |
|---|---|
| Quarter cylinder r = 1 and r = 2, weights (1, 1, 2) | yes; and not on the other radius |
| The usual quarter circle with middle weight √2/2, rounded to a double | no: a nearby conic |
| Sector over t ∈ [−0.75, 0.5] | inner face on r = 1, outer face on r = 2 |
| Sphere patch, the product of two exact circles (degree 2 × 2) | yes, on the unit sphere; not on radius 2 |
| Cone patch s·(circle, 1), s ∈ [1, 2] | yes, on x² + y² = z²; not on x² + y² = 4z² |

**Coincident curved surfaces.** Two sectors of one cylinder have different patches, so until now an operation in which their curved faces overlap was reported as not solved. `op::affine_reparam(P, Q)` now proves exactly that Q(u, v) = P(a u + b, c v + d):
- **Candidate map:** Q's corners are located on P in doubles, and the candidate map is snapped to a dyadic grid.
- **Exact verification:** P is restricted to the candidate rectangle by blossoming in expansion arithmetic and compared with Q by cross-multiplication.
- **Requirement:** both scales must be powers of two, so that the map and its inverse are exact and exact parameters stay exact in both directions.
- **In the Boolean:** such surfaces join the classes that coplanar faces already used. The overlay, mapped arcs, foreign edges, on-boundary classification and the comparison work through those classes unchanged.

| Case (every value the closed form ½(θ₁ − θ₀)(r₁² − r₀²)(z₁ − z₀) to 1e-12) | Volume | Same solid as |
|---|---|---|
| A = sector(1, 2; t 0–1; z 0–1) ∪ B = sector(1.5, 2; t ½–1; z ¼–¾), B against A's outer wall | 2.356194490192 | A |
| A ∩ B | 0.281531735097 | B |
| A − B (a pocket in the outer wall), then (A − B) ∪ B | 2.074662755095, then 2.356194490192 | A |
| Left half ∪ right half of the quarter (t 0–½ and ½–1) | 2.356194490192 | A |
| Left half ∩ right half (they touch) | empty | |
| Overlapping sectors (t 0–½, z 0–1 and t ¼–¾, z 0–½): ∩, ∪, − | 0.328003418811, 1.660723076691, 1.062939408192 | the closed-form sector for ∩ |
| Fourth-quadrant sector (t −½–0) ∪ the quarter | 3.747137317195 | |
| Two sectors related by a scale of 2/3 | reported as not solved | |

Every result is closed, tessellated into a closed mesh with the Brep's Euler characteristic, and its mesh follows the surfaces within 2e-3.

**What came up while building it:**
1. *Flat faces with uniform weights.* A sector's radial face at t = ½ has weight 1 + t² = 1.25 on every corner. The plane tests required unit weights; they now take any uniform positive weight. Coplanarity is decided projectively (a 4 × 4 determinant of homogeneous points, exact), and corner identity by exact cross-multiplication, with no division by weights.
2. *Interval maps destroyed exactness.* Interval arithmetic widens every result by a unit, even an exact one, so a dyadic map turned an exact corner chart into an interval straddling a side, and "does this vertex lie on that side" became undecidable. Affine maps now evaluate exactly when every input is a point and the result is a double.
3. *A wrong preimage.* The flat annulus, extended beyond its square, reaches the same points again with negative radius. Locating a corner then returned (−4, −2) instead of (2, 0). Among converged starts, the one nearest the parameter square is now taken.
4. *One point found twice.* Where two sectors overlap, a point on three surfaces is found from an edge of each. When the side test could not decide, the two roots stayed separate vertices. They are now identified by proving, with the three-surface Krawczyk test, that a box containing both holds exactly one common point. The vertex then takes the intersection of both sets of certified parameters.

### Chained operations with curved operands

The random chains so far used boxes only. `ak_bench` §11 runs 40 chains of up to 6 operations from the quarter cylindrical shell. Each step uses, with equal probability, one of three operands:
- a cylinder sector on the same family of cylinders, with dyadic radius, angle and height ranges, so that faces of different operands often lie on one cylinder or one plane without being the same patch;
- a box with a hill on top, at a grid position;
- a box of random orientation.

Each step is checked against a Monte Carlo count of 200,000 points that uses only the analytic definition of the chain. The table and the list of refusals below are from before tangential contact was resolved; the current numbers (192 steps, 173 committed) are in "Contact cases left after the flat-face work" below.

| Step | Attempted | Solved and committed | Operands with coincident curved surfaces | Worst Monte Carlo deviation | Beyond 3 sigma |
|---|---|---|---|---|---|
| 1 | 40 | 29 | 8 | 2.20 sigma | 0 |
| 2 | 29 | 18 | 9 | 2.62 sigma | 0 |
| 3 | 18 | 11 | 5 | 1.38 sigma | 0 |
| 4 | 11 | 9 | 3 | 1.40 sigma | 0 |
| 5 | 9 | 8 | 2 | 1.55 sigma | 0 |
| 6 | 8 | 6 | 3 | 0.99 sigma | 0 |

- **Committed results:** all 81 agree with Monte Carlo, and all 81 tessellate into closed meshes with the Brep's Euler characteristic. No result was refused by the audit. (This table first said 70. That count came from a benchmark run made before the last fixes of the curved-chain work: the run's failures include "could not decide whether a vertex lies on it" and "a vertex of the other solid has no parameters on this face", which the finished code no longer produces, and its messages predate the one-plane and tangency wording. Run now, the chains commit 81, and they still do with every pruning step of the spatial-index work switched back off together (same result hash), so the difference is not from that work.)
- **Not solved: 34 steps, each reported.** The geometry here is degenerate by design, because the dyadic grid puts flat faces exactly on cylinders and planes:

| Reason (count) | What it is |
|---|---|
| Tangent along a side (10) | A flat face touching a cylinder along a line: a hill box side at y = 1 against the cylinder r = 1 |
| One plane, parameterisations not related affinely (10) | A sector's annular cap (rational in the angle) and a parallelogram in the same plane. The overlay maps parameters through affine maps only |
| Intersection not proven complete (9) | Tangential contact, or two surfaces meeting only on the boundary of both domains |
| Vertex and ordering cases near such contacts (5) | A vertex that could not be separated from an end of its edge, curves at a vertex that could not be ordered |

**What the curved chains found:**
1. *A benchmark bug that looked like two wrong answers.* The first run reported two results 11 and 14.6 sigma off Monte Carlo. The kernel was right (the union of chain 22 is 0.75π + 0.98 − 0.3216 = 3.0146, exactly what it computed). The Monte Carlo sampling box was too small: sectors reach radius 2.5 and negative y. The box now covers every operand.
2. *A result refused by the audit* (chain 15: the quarter shell minus a sector over its upper right part, on the same cylinders and radial plane). One end of the intersection curve between the sector's radial plane and the shell's top was judged off the curve. The vertex lay exactly at the end of an open curve, but its slice value, measured with the end segment's rounded direction, came out a few units past the end. It was then rejected as a point on the continuation outside the parameter square, so the cut was missing. The join slack (1e-9 in slice) is now allowed at the two ends of an open curve. A point beyond the end by more is still rejected; the regression test for that still passes.
3. *Exact maps between coincident planes.* The map between two flat patches in one plane was computed by interval division, which widens even when the map is exactly u′ = 2u − 1. A vertex at a corner then lost its exact parameters on the other patch, and "does it lie on this edge" became undecidable. Where all corners are exact points (weights powers of two), the coefficients are now computed exactly and kept as points when they are doubles. This also narrowed the vertex enclosures in the grid chains (widest 4.0e-14 → 2.5e-14), which changed that section's hash (§8). Every count and comparison there is unchanged.
4. *Flat faces of higher degree.* The bottom and sides of a box with a hill on top are flat, but of degree (2, 2) and (2, 1). They are recognised exactly as degree-elevated parallelograms: every control point is where degree elevation of the corners puts it, checked with integer-weighted exact sums. Every flat-face test now works on the patch corners.
5. *Naming the one-plane case.* Two flat faces in one plane whose parameterisations are not affinely related used to be reported as tangency. An exact test (every control point of both in one plane, by 4 × 4 determinants) now names that case.

Regression cases for 2 and 3 are in `tests/op/quadric_test.cpp`, with closed-form volumes.

### Scaling and the spatial index

Every operand above has at most a few dozen faces, so the Boolean's time was almost all certified intersection of surface pairs that meet. `ak_bench` §12 measures what happens as models grow: a unit cube whose six sides are tiled into k × k flat patches (k = 1 to 16, so 6 to 1,536 faces), minus a small skewed parallelepiped at one corner. Its edge vectors are dyadic (in units of 1/256), so its own tiles stay exactly flat, and its origin is at odd multiples of 1/1024, so no tile vertex lies on a face plane of the cube. The volume must not depend on k, and it does not: 0.995140302715 in every row.

Small solid whole (6 faces), times in ms:

| k | Faces | Surface pairs considered | Intersected | Operands | Surface/surface | Edge/surface | Splitting | Regions | Total |
|---|---|---|---|---|---|---|---|---|---|
| 1 | 6 | 36 | 6 | 0 | 13 | 0 | 0 | 1 | 15 |
| 4 | 96 | 576 | 8 | 3 | 21 | 1 | 0 | 2 | 28 |
| 8 | 384 | 2,304 | 16 | 13 | 36 | 2 | 2 | 4 | 59 |
| 16 | 1,536 | 9,216 | 27 | 55 | 71 | 5 | 7 | 10 | 157 |

Small solid tiled k × k as well (both operands large):

| k | Faces of each | Surface pairs considered | Intersected | Operands | Surface/surface | Edge/surface | Splitting | Regions | Total |
|---|---|---|---|---|---|---|---|---|---|
| 1 | 6 | 36 | 6 | 0 | 12 | 0 | 0 | 1 | 14 |
| 4 | 96 | 9,216 | 20 | 8 | 44 | 1 | 1 | 3 | 58 |
| 8 | 384 | 147,456 | 48 | 29 | 111 | 4 | 3 | 7 | 160 |
| 16 | 1,536 | 2,359,296 | 96 | 126 | 330 | 11 | 13 | 21 | 519 |

The first version of this benchmark (with a different small box) took 15.5 s at k = 16, 14.8 s of it in preparing the operands, and it failed for every k ≥ 4. What was found:

1. *A needless refusal on tiled planes.* An edge between two tiles of one plane was asked whether a vertex lay on it before checking where along its line the vertex was. A vertex on the same straight line but on another tile's edge far along it has parameters that straddle the line, so the question was undecidable, and the operation was refused. Points beyond the edge's ends are now set aside first.
2. *Plane classes, quadratic.* Every pair of flat surfaces got an exact coplanarity test (1.2 million at k = 16). Coplanarity is an equivalence, so each surface is now tested against one member of each plane found so far, after an interval key (unit normal and distance, up to sign) has ruled most planes out.
3. *Charts on every tile.* A vertex was given parameters on every surface of its plane: 256 per cube side at k = 16, and every later lookup scanned them. A vertex now gets them only on surfaces whose faces or edges are near it. (A first version looked at faces only. That broke one step of the grid chains: an edge can lie along a side of a surface whose own faces are elsewhere, its face being on another surface of the plane. A vertex at the end of such an edge then lacked parameters there, a second copy of it was created, and the step was refused. The benchmark hashes caught it; edges now count too.)
4. *The inside-out check* evaluated a 16 × 12-point quadrature on every face. On polynomial faces with straight sides it now uses Gauss rules sized to the degree, which are exact: 2 × 3 points on a flat bilinear face.
5. *All-against-all loops.* Surface pairs to intersect, surfaces an edge may pierce, vertices of one operand that may be vertices of the other, vertices that may lie on an edge, and edge pieces that may lie in a face all go through a bounding-volume hierarchy over certified boxes (`geom/box_index.hpp`). Disjoint boxes prove that two entities do not meet, so nothing is lost. Results come back in increasing order, so the order of everything downstream is unchanged.
6. *A flat surface against a face entirely on one side of its plane* is no longer passed to the certified surface/surface solver (interval signs at the face's control points). At k = 16 with one operand whole this cut the pairs intersected from 83 to 37 (first benchmark version).

**What these numbers show.** With one operand large the overhead is now about 85 ms for 1,536 faces, growing linearly. With both large, everything except surface/surface intersection takes 190 ms at 1,536 + 1,536 faces. Surface/surface is 330 ms for the 96 pairs whose faces really meet, about 3.5 ms each: output, not overhead. The pairs considered (2.4 million at k = 16) are counted in closed form, not enumerated. On the six-face operands of the earlier sections nothing changed in speed (a chain step on a 37-face solid still takes about 0.3 s, almost all of it certified intersection), and nothing changed in results: every earlier hash is reproduced, except the curved chains', whose recorded value was stale (§8).

Searching surface pairs only over the regions that hold faces (`BooleanOptions::face_domains`) is implemented but off by default: it gave no measurable speed-up on the benchmark chains, and it changes how curves are segmented.

### Plane pairs certified directly

Two flat faces (exact parallelograms, X(u, v) = X0 + u E1 + v E2 and Y(s, t) = Y0 + s F1 + t F2) in planes that are not parallel meet in one straight line. More precisely, the whole solution set of X(u, v) = Y(s, t) in (u, v, s, t) is one line ℓ = P + σ d, and P and d can be enclosed by intervals from the corners. The general solver traced that line and then proved completeness by subdividing the 4D domain (about 1,200 boxes per pair). `BooleanOptions::plane_lines` (on by default) replaces both steps:

1. A certified tube is built along ℓ, a little beyond the part inside the domain, in segments about as long as traced ones.
2. Each certified segment has exactly one solution on each of its slices. Every solution lies on ℓ, so that solution is ℓ's point on the slice, and the segment covers the σ-range between its end slices. Intervals bound that range from inside.
3. If these ranges overlap and contain an outer bound of the σ-range where ℓ is inside the domain, the intersection is complete. That is a proof, without subdivision. If any step cannot be certified (parallel planes, degenerate parallelograms, a gap between ranges), the pair goes to the general solver as before.

It applies only to exact parallelograms. A box in a random orientation has corners rounded to doubles, so its faces are generally not exact parallelograms; the kernel correctly does not treat them as flat, and they still take the general path.

Measured in one run with the option on and off (same machine state; the earlier tables above were measured before this change):

| Case | Off | On |
|---|---|---|
| Unit cube minus the small solid (k = 1) | 15 ms | 5 ms |
| Tiled, one side, k = 16 (1,536 + 6 faces) | 161 ms | 105 ms |
| Tiled, both, k = 16 (1,536 + 1,536 faces) | 524 ms | 245 ms (surface/surface 330 → 52 ms) |
| Grid-box chains (§3, 40 chains, axis-aligned boxes) | 10 s | 6 s |
| Random-box chains (§3, 40 chains, rotated boxes) | 134 s | 126 s |

Every count, audit, Monte Carlo deviation and identity is the same with the option on and off; only the segmentation of the curves differs, so the result hashes differ (§8). With the option off every earlier hash is reproduced exactly.

Found on the way: on the tiled solids, `newton_slice` sometimes stalled at rounding-noise level just above its absolute step threshold (1e-14), because small patches amplify the noise through the inverse Jacobian, and then reported failure for a well-conditioned point (sin of the angle between the surfaces 0.89). A certified point on part of a curve could then not be placed. Newton now also accepts a point whose last step is below 1e-12 after its 25 iterations. It only supplies the guess for the Krawczyk test that certifies the point, so nothing is accepted without proof. Points that converged before are computed exactly as before.

The random-box chain profile shows where their time goes now: surface/surface about 120 to 210 ms per step and edge/surface (three-surface roots along earlier intersection edges) up to about 210 ms per step, both on surfaces that are not exactly flat.

### Edge/surface searches along arcs

`ak_bench chains 40 6 0 profile` now also splits the edge/surface stage: straight edges (curve/surface roots) cost 1 to 2 ms per step; arcs, which are earlier intersection curves whose crossings with a third surface are three-surface points, cost 50 to 230 ms per step, about 2 ms per search. Most of that went into proving that an arc does *not* meet a surface. Within those searches, 7.4 of 9.5 s (12 chains) were spent on one function, `arc_enclosure`: the tight box around a piece of the arc, obtained by contracting the segment's box with four interval Krawczyk iterations. Three changes, none of which changes any result:

1. **One enclosure per piece and edge.** An enclosure depends only on the tube segment and the slice range, not on the surface searched. The search bisects at fixed midpoints, so an arc searched against several surfaces asked for the same enclosures again. They are now cached per edge (`op::ArcEnclosureCache`).
2. **Cheap exclusions first.** Before a tight enclosure is computed, the segment's own box (at the top) or the parent's box cut down by the slice slab (when bisecting) is tested against the surface, by bounding boxes and by extents along the surface's normal. The tight enclosure lies inside that box, so whatever the cheap box excludes the search would have excluded at its first node.
3. **Stop at a fixed point.** The contraction stops when an iteration leaves the box exactly unchanged; the remaining iterations would repeat it bit for bit.

Two other ideas were measured and dropped: restricting the search to the rectangle of the third surface that holds faces (no effect, because on these models the faces span their surfaces), and starting the contraction from a smaller box (the contraction reaches the same box either way).

Every hash of the full benchmark is unchanged (core, Boolean pairs, chains, reuse, grid chains, tessellation, identities, curved chains and scaling). The search visits fewer nodes (in 12 chains, 93,480 instead of 108,066), only because nodes that would have been excluded at once are no longer entered. Timings from two full runs whose unaffected sections (predicates, fill) agree to within a few percent:

| Edge/surface, mean ms per step | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|---|
| Random-box chains, before | 2 | 81 | 157 | 225 | 234 | 245 |
| Random-box chains, after | 1 | 30 | 50 | 65 | 74 | 75 |
| Reused-operand chains, before | 2 | 97 | 143 | 171 | 173 | 198 |
| Reused-operand chains, after | 1 | 38 | 45 | 67 | 64 | 75 |

The whole benchmark went from 1,010 s to 882 s, and the identity test from 156 s to 130 s. Surface/surface on faces that are not exactly flat is now the largest cost of a chain step (about 150 to 290 ms).

### Tangential contact with planes

Tangency was the largest group of refusals in the curved chains: a box face touching a cylinder along one of its lines (the box's side at y = 1 against the cylinder r = 1), where no Krawczyk tube can be certified because the two surfaces have one tangent plane along the whole line. Intervals cannot prove such a contact (any perturbation destroys it), but for a plane it can be decided exactly from the input data (`op/contact.hpp`):

1. **The surface pair.** The plane's equation restricted to a patch, g(u, v) = det(c00, c10, c01, X(u, v)), is a polynomial whose Bernstein coefficients are computed exactly. If they vanish along a side, g is divided by that side's linear factor (exact integer multiples of the coefficients) as often as the side stays zero: twice for a cylinder against a plane along its line. If the quotient has one strict sign over the square (all Bernstein coefficients of one sign, after exact subdivision at 1/2 where needed), the patch meets the plane in that side and nowhere else. The side is already an edge of the patch's faces or outside them, so the pair contributes no curve and needs no completeness search. The same works for several sides at once: a hill on a flat base touches the base's plane along all four.
2. **Edges.** A circular edge ending on the line touches the plane at its end vertex. The same division in one variable proves that it meets the plane only there; the vertex gets parameters on the plane. An edge that lies in the plane meets the touching surface only on the touching line, where the other solid's edges along that line are searched against it; it is not searched again.
3. **Exact parameters.** A vertex on the plane and on a surface that meets the plane only along one side is on that side, so its parameter across the side is set to the exact value (proven on the square enlarged by 1/64, so that enclosures straddling the side are covered), and carried exactly to every surface that shares the side by construction. Order along edges and "is this vertex on this edge" then stay exact.
4. **Order at a vertex.** In a face's parameter plane the plane's trace leaves the touching vertex tangent to the square's side. Its direction enclosure straddles the side's direction, but the curve lies in the face, so its true direction is the side's direction turned inwards by an angle between 0 and the enclosure's width. If every other curve at the vertex is provably outside that narrow cone, the direction is replaced by the side's direction turned inwards by an infinitesimal angle: a symbolic perturbation, compared exactly (a cross product is decided by the exact directions and, on an exact tie, by the first-order term). The same works for a curve tangent to a straight piece inside a face, where the piece lies in a plane that the curve's surface touches only along a side: the curve never crosses that line, so its side is read off any certified point of it. Two curves tangent at one vertex are still refused.
5. **Around an edge.** Where solids touch along the line (a non-manifold result), the curved face and the flat face leave the edge at the same angle. The fan check orders them by the side of the plane on which the curved face lies, from the same exact proof.

One kind of contact that is not tangential also stopped these chains: a transversal curve that meets both parameter squares only on their boundaries (two faces that touch at a corner of each). It cannot be traced inside the squares, and the completeness proof cannot close around the touching point. Such a pair is now searched again over squares enlarged by 1/32: the curve then runs through the point, its tube covers it, and completeness over the larger box implies completeness over the original. Parts outside the squares lie outside every face and are clipped away. The second search is limited to pairs without shared sides (their continuations would duplicate known curves), to 8,000 completeness nodes and to 2,000 tracing steps per branch; the searches that succeed use under 1,000 nodes and a few dozen segments. Without the step limit, one pair near a true tangency traced 50,146 segments and took 28 s before failing.

`tests/op/contact_test.cpp` checks the exact tests (order 2 along a cylinder's line, order 1 for a plane meeting a plane, no single-side contact for a plane through a cylinder or for a hill on its base, a circle touching a plane at one end) and nine Booleans with closed-form volumes, all equal to 1e-11:

| Case | Intersection | Union | Difference |
|---|---|---|---|
| Quarter shell (1 ≤ r ≤ 2) and a box whose face y = 1 touches the inner cylinder inside that face | 0.001314571430 | 2.604879918762 | 2.354879918762 |
| The same with the face x = 1 | 0.038964887548 | 2.817229602644 | 2.317229602644 |
| A box outside the shell whose face x = 2 touches the outer cylinder | empty | 2.481194490192, not manifold (the solids meet along the line) | 2.356194490192 |

Every result passes the audit and tessellates into a closed mesh.

**On the curved chains (§3):**

| | Before | After |
|---|---|---|
| Steps attempted | 115 | 131 |
| Solved and committed | 81 | 100 |
| Refused: tangent along a side | 10 | 0 |
| Refused: curves at a vertex could not be ordered | 4 | 0 |
| Refused: intersection not proven complete | 9 | 8 |
| Refused: a vertex not separated from the ends of its edge | 1 | 4 |
| Refused: flat faces in one plane, parameterisations not affinely related | 10 | 19 |
| Worst Monte Carlo deviation; beyond 3 sigma | 2.62 sigma; 0 | 2.62 sigma; 0 |
| Committed results tessellated closed, Euler characteristic matching | 81 | 98 (and 2 empty results) |

Committed per step: 34, 23, 15, 13, 8, 7. The refusals that are left are reported, never answered. Most are the one-plane case, which the chains now reach more often because they run further. The four vertex-separation refusals include the one degenerate case this work found and did not resolve: a box edge lying exactly on the line along which a cylinder touches that box's face (chain 29), so two solids have collinear, overlapping edges on the contact line. The four earlier "could not be ordered" refusals (rotated boxes) were not tangency: two consecutive pieces of one operand edge, split at a vertex, whose chord directions straddled zero. Two pieces of one edge now continue each other by construction, and those chains run on.

The section's Boolean time went from 47 s (115 steps) to 105 s (131 steps). The slowest steps (4 to 11 s) all end in a refusal after long completeness searches on flat faces in one plane with unrelated parameterisations, or near true tangencies. Every other hash of the benchmark is unchanged (§8): all of this runs only where the earlier code refused.

### Flat faces on rational carriers

After tangency, the largest group of refusals (19 of 31) was a flat face whose carrier is not an affine parallelogram, the annular cap of a cylinder sector (rational in the angle), in one plane with a box face or with another cap parameterised the other way round. The overlay of coincident flat faces works in one parameter plane through affine maps between the carriers, and there is none between a cap and a parallelogram.

Such a face is now moved, before the Boolean, onto an exact parallelogram of its plane (`flat_detail::recarry` in `op/boolean.hpp`); the face itself, its edges and its vertices are unchanged:

- **The new carrier.** The plane of the cap is found exactly (every control point checked with expansion arithmetic; axis-aligned planes only). The parallelogram has corners on a 1/16 grid and a power-of-two extent, so a vertex at an exact dyadic point, such as a corner of the cap with weight 2, gets exact parameters on it, and the parallelogram is related to box faces in the plane by exact affine maps.
- **Straight edges** (the cap's radial sides) stay side edges, of a thin exact parallelogram built on the same side, so that they lie on a surface of the plane's class and vertices on them get parameters there by construction.
- **Curved edges** (the cap's arcs) become arc edges: the transversal intersection of the cylinder with the plane, certified as a tube along the cylinder's side. An arc that an earlier operation cut into the cap (the cap and another surface) is certified again as the curve of the plane and that surface, through points of the old arc.
- **Only where needed.** A cap moves only if a face of the other solid lies in its plane on a carrier it is not affinely related to; then every flat rational carrier in that plane moves, in both operands, and only copies of the operands are changed (`Imprint::owned`). An operand whose cap cannot be moved is left as it was, and the conflict is reported as before.

Three things this needed elsewhere:
1. Arcs of the two solids on different but equivalent surfaces (one cylinder under an exact reparameterisation, one plane) are recognised as one piece through the exact maps between those surfaces, where the surfaces themselves give no match.
2. A straight piece inside a face takes its direction from its edge's exact end points where those have exact parameters (carried from a surface of the same plane where needed), signed by the piece's own chord.
3. An arc on a surface that touches a plane only along a side meets that plane only at its own vertices; it is not searched against it.

Moving the cap also exposed an accidental coincidence: the re-certified curve was first made of 8 chord pieces, whose joints fell exactly on vertices at 1/4 and 3/4 of the side; a piece starting at a joint began with an empty span, and its tessellation had two boundary points at one place. The curves are now made of 7 pieces.

**A fault found by the Debug build.** The first Debug run of `flat_test` did not finish: it stopped in the sector-below case with almost no CPU time, waiting on a "Debug Assertion Failed" dialog. The matching of arc pieces by surface classes (item 1 above) compared an iterator of one map with the end of another, which is undefined behaviour; the checked iterators of the Debug library caught it, while in Release the comparison happened to take the intended branch, so no Release result was affected (the curved-chain hash is the same before and after the fix). Each lookup now has its own iterator. The test programs and the benchmark now also report a failed library assertion on stderr and stop with an error, instead of waiting on a dialog (`tests/check.hpp`), so a fault of this kind fails a test run instead of hanging it.

`tests/op/flat_test.cpp` checks ten Booleans with closed-form volumes, all to 1e-11, each result tessellated closed and used again as an operand: a box standing on the shell's annular cap inside the shell (intersection 0.125, union 2.356194490192, difference 2.231194490192); a box on that plane that also touches the inner cylinder along its line (both contacts at once; 0.001314571430, 2.604879918762, 2.354879918762); a sector below the shell whose top cap meets the shell's bottom cap with swapped parameters (empty, 2.570550820428, 2.356194490192); and a box whose faces share no plane with a cap, where nothing is moved.

**On the curved chains (§3):**

| | Before tangency | After tangency | After flat faces |
|---|---|---|---|
| Steps attempted | 115 | 131 | 156 |
| Solved and committed | 81 | 100 | 128 |
| Refused: flat faces in one plane, parameterisations not affinely related | 10 | 19 | 0 |
| Refused: tangent along a side; curves at a vertex not ordered | 14 | 0 | 1 |
| Refused: intersection not proven complete | 9 | 8 | 13 |
| Refused: other (vertices not separated or ordered, an edge through a vertex or a face boundary, region decisions) | 1 | 4 | 14 |
| Worst Monte Carlo deviation; beyond 3 sigma | 2.62 sigma; 0 | 2.62 sigma; 0 | 2.62 sigma; 0 |
| Committed non-empty results tessellated closed, Euler characteristic matching | 81 | 98 | 124 (and 4 empty results) |

Committed per step: 37, 32, 23, 16, 11, 9. The chains now run deeper, into configurations that were never reached before, and the refusals that are left are spread over many kinds: the largest are pairs not proven complete (13, mostly true tangencies at a point or between curved surfaces) and the four vertex-separation cases of the tangency section. The section takes about 95 s (it took 128 s after the tangency work), because steps that ended in long searches on the one-plane case now solve.

### Contact cases left after the flat-face work

Of the 28 steps still refused after the flat-face work, this round resolved five kinds and found one fault. What changed:

1. **One cylinder, parameters exchanged.** A sector placed in the shell's hole touches the shell's inner cylinder with its outer face. Both faces lie on the cylinder r = 1, but an outer and an inner face exchange the roles of the two parameters, so no axis-aligned affine map relates the patches and the pair was reported as not proven complete. `op::affine_reparam_any` also tests the transposed patch, by the same exact comparison, and the classes of coincident surfaces accept a map that exchanges the parameters (it reverses orientation, which the overlay already reads from the map's determinant). Two caps in one plane with swapped parameters are related the same way, so they are no longer moved onto parallelograms.
2. **Edges in a cap's plane without a face there.** A box with a hill on top has no face in the plane of the hill's base, but the four edges around the hill lie in it. Where that plane holds a cap, the edges cross the cap's face, where the rational parameters bend them. A cap is now also moved when the other solid has a straight edge in its plane that this solid does not have by construction. Three smaller things this needed:
   - a vertex at an exact point of a plane gets exact parameters there (snapped to a dyadic grid, then verified with expansion arithmetic);
   - a vertex on a side gets parameters on every surface that shares the side by construction;
   - an edge that crosses a face's boundary inside a plane holding a side of that face's surface is a crossing of two edges, not a refusal.
3. **A vertex in a plane, at the end of an edge.** Once a cap lies on a plane carrier, an old vertex of the other solid can lie exactly in that plane where one of its edges ends. The curve/surface search finds a regular root that cannot be separated from the edge's end, which was reported as "the surface passes through a vertex of the edge". For a plane this is decided exactly, with the test already used for tangential ends (`op/contact.hpp`): if the plane function of the edge's side vanishes only at its ends, the root is that end vertex, which gets parameters on the plane, and the edge has no other point there.
4. **One circle, two kinds of edge.** The circle of a moved cap is an intersection curve (cylinder and plane). The other solid may have the same circle as a side of a cylinder patch. Pieces of the two solids were compared side with side and curve with curve, never across, so the shared piece stayed two pieces and its face could not be built. A curved side and an intersection curve between the same two vertices are now one piece when the side's midpoint, carried into the curve's two surfaces (through the exact map of its class, and into the plane the side lies in exactly), lies on the certified curve inside the piece. The curve then stands for both, since faces on either of its surfaces can carry it.
5. **Straight edges of faces of higher degree.** The flat faces of a box with a hill on top have degree 2, so their straight edges are degree-2 sides. Only sides of degree 1 counted as straight, and straight pieces are the ones recognised as one edge from their two end vertices alone. A hill box's edge running along an edge of the other solid between the same two vertices therefore stayed a second edge. A side whose control points lie exactly on one line (equal weights, exact collinearity) is now straight.

**A fault found by the audit** (chain 20, step 4: a sector put back into a pocket whose roof is a hill). The audit refused the result: three edges were used once. The sector's top cap lies in the plane of the hill's base, and two base edges cross the cap and meet at a corner inside it, cutting off a small triangle. That triangle was missing from the result. Two causes, both fixed:

- *Point location was needlessly undecided.* A straight boundary piece whose end points are enclosures, not exact points, made every query inside its bounding box undecided. For a diagonal piece that box is large; here it held the midpoints of both base edges. The true piece and the chord used for the winding number are both segments from a point of one end box to a point of the other, so both lie in the convex hull of the two boxes, and the chord can be moved onto the piece inside that hull. A strict sign of (b − a) × (q − a) over all a and b in the end boxes and q in the query box proves the query is outside the hull. That test replaces the bounding box (`op/trim.hpp`, `clear_of_straight`); it also decides query boxes, not only points.
- *An undecided edge was read as a boundary edge.* When a face's planar map is built, a straight edge of the other solid lying in the face's plane was skipped unless its midpoint was proven inside. "Undecided" was taken to mean "along this face's boundary". That is now accepted only where the answer is "on a boundary piece". An edge does not cross the face's boundary (it was split where it does), so any point of it decides: six more points are tried, and an edge that stays undecided stops the operation with a message. With the sharper test switched off, the new regression test ends in that message instead of a missing face.

Reporting undecided edges at first refused five steps, one of which the earlier version had committed (correctly: the skipped edge did run along the boundary). Four were steps with a hill box, whose straight edges ran along an edge of the other solid; item 5 makes those one edge. The fifth is the one left in the table below.

The kernel's own audit caught this one: no wrong solid was committed. But it shows a rule worth keeping: an undecided answer must never be given a meaning.

**A mistake in my first version of item 2.** It looked at the sides of the other solid's surfaces, not at its edges. A result keeps the surfaces of its operands, so in A ∪ (A ∩ B) the shell's own cylinder, kept by A ∩ B, has a side in the plane of the shell's cap, and the cap was moved. The result was the right solid on another carrier, and the structural comparison rightly reported different edges: 12 identities of `tests/op/identity_test.cpp` failed. The rule now looks at edges and leaves out those along sides that the cap's own solid has too.

`tests/op/flat_test.cpp` now checks 29 Booleans with closed-form volumes (10 before), all to 1e-11, each result tessellated closed and used again as an operand. Two of the new cases were checked the other way round as well: with the sharper point location switched off the pocket case is refused, and with degree-2 sides not counted as straight the face-to-face case is refused, each with "could not decide whether an edge of the other solid lies inside it". New: a sector and a box side by side put under the shell; a hill box whose base is the plane of the shell's top cap, which also touches the inner cylinder along its line (0.016270695941, 3.002423794251, 2.339923794251); a sector in the shell's hole against its inner cylinder (empty, 2.380111406147, 2.356194490192); a slot through the shell so marked, with a vertex in a moved cap's plane and the circle shared between a side and a curve (0.298891834250, 2.356194490192, 2.057302655943); the pocket of the fault above, where the volumes before and after putting the sector back differ by the sector's closed-form volume; and a box with a hill on top face to face with a plain box, whose common edges are degree-2 sides in one solid and degree-1 sides in the other (empty, 2.133333333333, 1.133333333333).

**On the curved chains (§3):**

| | After flat faces | Now |
|---|---|---|
| Steps attempted | 156 | 192 |
| Solved and committed | 128 | 173 |
| Refused: intersection not proven complete | 13 | 6 |
| Refused: a surface through a vertex of an edge; an edge through an edge | 4 | 4 |
| Refused: vertices on a curve or curves at a vertex not ordered | 3 | 3 |
| Refused: three-surface search not resolved | 0 | 2 |
| Refused: an edge of the other solid not decided inside or outside a face | 2 | 1 |
| Refused: other (a vertex not separated from the ends of an edge or not decided on it, an edge touching a surface, region and curve-part decisions, a moved cap with an edge left on its old carrier) | 6 | 3 |
| Results refused by the audit | 0 | 0 |
| Worst Monte Carlo deviation; beyond 3 sigma | 2.62 sigma; 0 | 2.62 sigma; 0 |
| Committed non-empty results tessellated closed, Euler characteristic matching | 124 | 168 (and 5 empty results) |

Committed per step: 40, 37, 33, 26, 20, 17. Every first step is now solved, and 17 of the 40 chains run all six steps (9 before). In 80 of the solved steps the operands had coincident curved surfaces (44 before). One committed result is not manifold. The 19 refusals that are left are each reported with the entities involved. The largest group is still pairs not proven complete (6); the rest are single degenerate configurations of vertices and edges. The one undecided edge left is a side of one solid that runs along an intersection curve of the other between the same two vertices and is not recognised as the same edge. The section takes about 68 s for 192 steps.

### Exact points for vertices on round numbers

Most of the 19 refusals left were coincidences of round numbers. A typical one: a vertex of an earlier result at exactly (0.75, 1, 1), where the cylinder r = 1.25 meets two planes, and a face of the next operand in the plane y = 1. Its parameters were enclosures about 1e-13 wide, and no enclosure can prove that a point lies exactly in a plane. Models with round dimensions are full of such points.

**Exact roots.** A certified root comes with a box in which it is the only solution. A point that lies exactly on every surface of the root's system, at parameters inside that box, is therefore the root. Candidates for such a point are cheap to make and are then proven or dropped:

1. *From a surface.* The midpoint of a parameter enclosure is rounded to a multiple of 2^-30 (enclosures are far narrower than that grid, so there is at most one candidate). The surface's point there is computed by de Casteljau steps in expansion arithmetic (`op::exact_point_at`): an exact homogeneous point.
2. *From space.* The root's position is rounded to the same grid. This finds points such as (1, 0.75, 0), which is round in space but sits at one third of its edge and one third around its circle, so that no surface has dyadic parameters for it.
3. *Proof.* The candidate must lie exactly on each of the other surfaces: for a plane, an exact determinant; for any other surface, that surface's exact point at candidate parameters (dyadic, or a fraction with a small denominator found by continued fractions and evaluated exactly in homogeneous form) must be the same point by exact cross-multiplication. Its parameters there must lie inside the uniqueness box. For a straight edge the candidate must lie exactly on the edge's line, whose side must be linearly parameterised (checked exactly).

A root that passes is known exactly: its dyadic parameters become points, the others tight enclosures, and the vertex carries its exact point (`BVertex::xp`), which results pass on to later operations. Nothing is assumed: a candidate that fails any test leaves the root as it was, an enclosure. Vertices of the two operands with the same exact point are one vertex. `BooleanOptions::exact_points` (on by default) switches all of it off, for comparison.

**What it decides.** A root that cannot be separated from an end of its edge was reported as "the surface passes through a vertex of the edge". It is now tested against the edge's end vertices: a vertex with an exact point that lies exactly on the surface (a plane: an exact determinant; any other surface: its exact point at the candidate parameters), with parameters inside the root's uniqueness box, is the root. The vertex gets parameters on the surface and no new vertex is made.

**What it does not decide.** Chain 18 shows the limit. Its vertex (1, 0.75, 0) is now recognised on the next operand's cylinder r = 1.25, and the step then stops at an edge of the first solid that lies entirely in that cylinder, along one of its straight lines. An edge lying in a curved surface of the other solid is a different problem (it is one for planes too, unless the plane holds a face in a known relation), and it is still reported.

**Straight curves between flat surfaces.** One smaller rule went in with it. An intersection curve of two flat surfaces in different planes is a straight line also when one of them is a cap on its rational carrier (flatness and "different planes" by the exact test that names the one-plane case). Such a piece is now matched with a straight edge between the same two vertices, and it continues a straight edge on the same line at a vertex.

`tests/op/exact_test.cpp` (new) checks the candidates (1/3, 5/7, −1/4; none for π), the exact evaluator against the ordinary one on every surface of a sector and a hill box at dyadic and rational parameters, and a constructed case:
- A sector (1.25 ≤ r ≤ 2) with a notch whose edge along y pierces the inner cylinder at (0.75, 1, 0.5), a 3-4-5 point. The volume is the closed form 1.332120426389 to 1e-11, and the new vertex is proven exact.
- A tilted box whose face y + z = 1.5 passes exactly through that vertex, across the three edges that end there. Intersection 0.700327488934, union 24.631792937456, difference 0.631792937456: all three results pass the audit and tessellate closed, and the three volumes, from three separately assembled results, satisfy inclusion-exclusion to 1e-10. With exact points switched off the operation is refused ("the surface passes through a vertex of the edge").

**On the curved chains (§3):**

| | After the contact cases | With exact points |
|---|---|---|
| Steps attempted | 192 | 196 |
| Solved and committed | 173 | 178 |
| Refused: intersection not proven complete | 6 | 6 |
| Refused: a surface through a vertex of an edge; an edge through an edge | 4 | 2 |
| Refused: vertices on a curve or curves at a vertex not ordered | 3 | 3 |
| Refused: three-surface search not resolved | 2 | 2 |
| Refused: an edge touching a surface or lying in it | 1 | 2 |
| Refused: other (an edge of the other solid not decided against a face, a vertex not decided on an edge or not separated from its ends, a moved cap with an edge left on its old carrier) | 3 | 3 |
| Results refused by the audit | 0 | 0 |
| Worst Monte Carlo deviation; beyond 3 sigma | 2.62 sigma; 0 | 2.62 sigma; 0 |
| Committed non-empty results tessellated closed, Euler characteristic matching | 168 | 173 (and 5 empty results) |

Committed per step: 40, 37, 33, 28, 22, 18. Eighteen chains run all six steps, and three committed results are not manifold. Over the solved steps, 178 of 1,144 new vertices were proven exact, and in one step two roots were proven to be an end vertex of their edge. With exact points switched off (`ak_bench curved 40 6 noexact`) the chains commit 175 of 194 steps: the straight curves between flat surfaces alone account for two of the five steps gained.

The gain in steps is small because a chain that passes one coincidence often meets the next. Two chains got further: chain 8 now runs to the end (three more steps, through the exact vertex), and chain 22 runs two more steps (through the straight curve between flat surfaces) and stops at another degenerate case. Chain 18 passes its vertex and stops in the same step, as described above. What changed in kind is that "this vertex lies exactly on that surface" is now a question the kernel can answer for points with round coordinates.

### Edges lying in a curved surface

The last round ended at chain 18: a straight edge of one solid lying entirely in the other's cylinder r = 1.25, along one of its straight lines. The edge/surface search finds no isolated root there and cannot finish. Planes had this case covered; curved surfaces did not.

**Deciding it exactly.** A cylinder patch has degree 1 along its axis, so its curve at a constant angle parameter a/b is the segment between its exact points A and B there (exact evaluation at a rational parameter, as in the previous section). The candidate parameter comes from locating the edge's two end points on the surface in doubles; it must be the same at both ends and near a fraction. The proof: A and B lie exactly on the edge's line (exact cross products), and A and B have equal weights, so the patch is linear along that line. The edge then lies in the surface, along u = a/b or v = a/b (`find_iso_line` in `op/boolean.hpp`).

**Carrying it through the Boolean.** Such an edge is treated as an edge lying in a plane already was:
- *Surface pairs.* A straight side of one surface lying in the other along such a line is a known line of their intersection, certified directly as a tube, the counterpart of a side lying in a plane. Without it, tracing ends at the side and the pair is reported as a suspected tangency.
- *Edges.* The edge is not searched against that surface; its end vertices get parameters on it (exact where the position along the line is dyadic).
- *Faces.* The edge enters the planar map of the surface's faces as a straight piece.
- *Clipping.* The surface's curve with either of the edge's own surfaces is the edge itself along it. A point of the line, in both surfaces' parameters, must lie on the certified curve, which has one solution per slice of its tube.
- *Classification.* Inside/outside is never spread across such an edge.

**Line contacts.** Where a box only touches the cylinder along the line, from inside or from the hole, the edge lies inside the cylinder's face and bounds nothing there; the face is not split. Whether the edge's solid crosses the surface or stays on one side of it is decided with intervals at one point of the edge: both faces along the edge leave it across the sides of their squares, and their directions against the surface's normal must have one strict sign. One point is enough because the sign can change only where the surface is tangent to one of the edge's surfaces, and their curves along the edge are certified regular. A touching piece that nothing continues is taken out of the face's map and recorded. A result that keeps both the face and the edge has parts touching along that line without sharing an edge: `BooleanResult::line_contacts` counts them, and the result is reported as not manifold.

**A point on the contact line.** The ray test that classifies a region starts from one interior point. For the cylinder face of a sector, that point was the middle of its parameter square, which here lies exactly on the box's edge: every ray starts on the other solid. The ray-test step now tries further interior points (up to six) before giving up. A region whose first point works is classified as before.

`tests/op/exact_test.cpp` checks four Booleans against closed forms (inclusion-exclusion to 1e-10, each volume to 1e-11, every result closed and tessellated): a box inside the sector against its inner wall along the line at angle parameter 1/2 (0.125, 1.914408023281, 1.789408023281), the same at 1/3 (a point of the circle with non-dyadic parameter, same volumes), a box in the hole against the wall (empty, 2.039408023281, 1.914408023281), and a box whose edge is where the wall enters it (0.129229563931, 1.972678459351, 1.785178459351). In the first three, the result that keeps the wall and the box's edge reports one line contact and is not manifold.

**On the curved chains (§3):**

| | With exact points | Edges in curved surfaces |
|---|---|---|
| Steps attempted | 196 | 197 |
| Solved and committed | 178 | 180 |
| Refused | 18 | 17 |
| Results refused by the audit | 0 | 0 |
| Worst Monte Carlo deviation; beyond 3 sigma | 2.62 sigma; 0 | 2.62 sigma; 0 |
| Committed non-empty results tessellated closed, Euler characteristic matching | 173 | 175 (and 5 empty results) |

Committed per step: 40, 37, 33, 28, 23, 19. Chain 22 now runs to the end, and chain 18 passes its edge on the cylinder and stops one step later, at a surface through a vertex of an edge. Chain 28 passes its edge as well and stops in the same step at two vertices on a curve that could not be ordered; chain 11 stops in the same step for a different reason (an edge not decided against a face). The 17 refusals left: 6 pairs not proven complete, 3 orderings of vertices or curves, 2 edges through an edge, 2 three-surface searches, and four single cases.

### Ordering cases

Three of the 17 refusals were orderings: two vertices on one curve, and the curves leaving one vertex.

**One vertex reached twice** (chains 28 and 30). In both, the point (0.75, 1, z), where a hill box's edge meets the cylinder r = 1.25, was found twice: once from one edge, once from another edge on other surfaces. Both roots had exact points, and they were the same, but nothing compared them. The two vertices then sat at one place on a curve and could not be ordered. A new root is now compared with every vertex that has an exact point; if the points are the same, the root is that vertex, which receives the root's parameters where it had none. This is under `BooleanOptions::exact_points`.

**A curve continuing a line** (chain 17). At a vertex, a new intersection curve leaves almost exactly opposite to the straight edge that arrives there: it continues the edge, tangent to it. Its direction is an enclosure that straddles the edge's line, so "is it right or left of the way back" cannot be decided, and the clockwise order of the curves at the vertex was refused. The order does not need that answer. The enclosure points strictly backwards and straddles the line, so it contains the exact opposite direction. Against any other curve at the vertex, the cross product of the two directions, if it has one strict sign over the enclosure, gives the order the exact opposite direction would have; and a strict sign also means no direction in the enclosure coincides with the other curve, so the true direction is ordered the same way. `side_of` reports such a direction as "near opposite", and `cw_before` uses that rule; only cases it refused before take this path.

`tests/op/exact_test.cpp` replays chains 30 and 28 of the benchmark (operands on dyadic grids) up to the step that was refused. The step must solve with a merged vertex, satisfy inclusion-exclusion, pass the audit, tessellate closed, and match a Monte Carlo count of 400,000 points that uses only the analytic definitions of the operands: 0.008179360 against 0.008085 ± 0.000073 (1.28 sigma), and 0.380719722 against 0.379619 ± 0.005398 (0.20 sigma). With the merge switched off, both steps are refused as in the benchmark. Chain 17 depends on random rotated boxes and is covered by the benchmark only.

**On the curved chains (§3):**

| | Edges in curved surfaces | Ordering cases |
|---|---|---|
| Steps attempted | 197 | 201 |
| Solved and committed | 180 | 187 |
| Refused | 17 | 14 |
| Results refused by the audit | 0 | 0 |
| Worst Monte Carlo deviation; beyond 3 sigma | 2.62 sigma; 0 | 2.62 sigma; 0 |
| Committed non-empty results tessellated closed, Euler characteristic matching | 175 | 181 (and 6 empty results) |

Committed per step: 40, 37, 34, 29, 26, 21. Chains 17, 28 and 30 now run to the end (chain 30 ends in an empty result). No ordering refusal is left. The 14 refusals: 6 pairs not proven complete (tangency), 2 edges through an edge, 2 three-surface searches not resolved, and four single cases (an edge not decided against a face, a vertex not decided on an edge, a surface through a vertex of an edge, a moved cap with an arc left on its old carrier).

### Edges through edges

Two refusals (chains 9 and 36) were an edge of one solid passing exactly through an edge of the other: (0.5, 1, 0.75) and (0.75, 1, 0.25), each on a hill box's edge between its side y = 1 and its top. The edge/surface search finds the root, but it lies on the boundary of the face it pierces, and unless the edge lies in a surface of the other solid (which the earlier rules cover), that was reported as "edge through edge".

**Deciding it exactly.** The root's exact point (previous sections) is tested against the edges of the other solid that bound faces on the pierced surface. It must lie exactly on both of that edge's surfaces: on a curved one, at parameters found in doubles and then proven (dyadic, or a fraction with a small denominator); on a plane, by an exact determinant. And it must lie strictly inside the edge: a side by its position along the side, from the parameter across it, which must be exact; an arc by its position on the certified curve. Then the root is a vertex on both edges and all four surfaces, and both edges are split there as at any vertex.

On a plane, a vertex's two parameters were made exact together or not at all. At (0.5, 1, 0.75) on the hill box's side the parameter across the side is exactly 1 but the one along it is 1/3, so neither was exact, and "is this vertex on that side" stayed undecided. Each parameter is now made exact alone where it can be: the point must lie exactly on the parallelogram's line at that parameter, between its exact points at the other parameter 0 and 1.

`tests/op/exact_test.cpp` checks the unit box against a parallelepiped B with dyadic edge vectors, one of whose edges passes through the box's edge x = 1, z = 1 at (1, 0.5, 1); B's faces there are y = 0.5 and x + y − z = 0.5, so neither edge lies in a face of the other. In B's own coordinates the intersection is the simplex a + b + c ≤ 1/2, so its volume is 1/192. The three results match the closed forms to 1e-12 (0.005208333333, 1.244791666667, 0.994791666667), pass the audit, tessellate closed, and the intersection agrees with Monte Carlo (0.10 sigma). Without the rule the operation is refused with the old message.

`tests/op/boolean_test.cpp` had kept an "edge crossing edge" case as an expected refusal: a tilted box whose edge crosses the unit box's edge x = 1, y = 1 at (1, 1, 0.5). It is now solved and checked against closed forms: the intersection is 1/48 (in the tilted box's coordinates, the integral of 0.8 (c − b) over its region, 25/1152, times the determinant 0.96), and the union and difference follow, all to 1e-12. The case next to it, a tilted box whose edge leaves the unit box exactly through the unit box's edge y = 1, z = 1, is not a crossing in space: that edge lies in the plane y = z of one of the tilted box's faces. It is still refused, now at the classification of regions ("two regions that share an uncut edge piece disagree"), not at the edge search.

**On the curved chains (§3):**

| | Ordering cases | Edges through edges |
|---|---|---|
| Steps attempted | 201 | 201 |
| Solved and committed | 187 | 189 |
| Refused | 14 | 12 |
| Results refused by the audit | 0 | 0 |
| Worst Monte Carlo deviation; beyond 3 sigma | 2.62 sigma; 0 | 2.62 sigma; 0 |
| Committed non-empty results tessellated closed, Euler characteristic matching | 181 | 183 (and 6 empty results) |

Both chains now run to the end. No edge-through-edge refusal is left. The 12 refusals: 6 pairs not proven complete (tangency), 2 three-surface searches not resolved, and four single cases (an edge not decided against a face, a vertex not decided on an edge, a surface through a vertex of an edge, an arc left on a moved cap's old carrier).

### Single cases

Four refusals were single configurations, each with its own cause.

1. **An arc left behind on a moved cap's plane** (chain 38). A cap face can be bounded by arcs on two different flat carriers of its plane: one on the face's own carrier, one on another carrier of the same plane that an earlier operation left without a face. When the cap was moved onto a parallelogram, only arcs on the face's own carrier were certified again on the new plane. The other kept its old carrier, and the next step refused the operand ("an arc edge does not lie on the carrier of a face it bounds"). Arcs of a moved face on any flat carrier of the same plane now move with it; their vertices get parameters on the new plane from their old carrier.
2. **A root at an end whose parameter is a third** (chain 18). An edge that runs from t = 0 to t = 1/3 of its side ends at (1, 0.75, 0), where the next operand's cylinder r = 1.25 passes. The exact-end rule of "exact points" required the end parameter to be exact. It only needs to lie in the root's uniqueness box, which an enclosure of 1/3 does. The chain then gets one step further and stops at a classification clash (below).
3. **A straight side and a straight intersection curve** (chain 11). A side of one solid along a ruling of a cylinder, and the intersection curve of a plane and that cylinder in the other solid, are the same segment between the same two vertices. Sides and intersection curves are matched by the rule of "contact cases" (a point of the side, carried into the curve's surfaces, lies on the certified curve), but straight sides had been left out of it. They are now included.
4. **A vertex not decided on a side** (chain 19). A vertex with an exact point lies on a side of a hill box, but its parameter across the side was an enclosure. Its parameters on the side's surface are now proven exact from the exact point where they can be (on a plane, one at a time), and the side test is then exact.

`tests/op/exact_test.cpp` replays chains 19 (all six steps) and 38 (three steps), whose operands are on dyadic grids. The last step must solve, satisfy inclusion-exclusion, pass the audit, tessellate closed, and match Monte Carlo: 0.146828251 against 0.147776 ± 0.003380 (0.28 sigma), and 0.184121557 against 0.182118 ± 0.003751 (0.53 sigma). With items 4 and 1 switched off, exactly these two steps are refused, with their old messages. Chain 11 has a random rotated box and is covered by the benchmark only.

**On the curved chains (§3):**

| | Edges through edges | Single cases |
|---|---|---|
| Steps attempted | 201 | 205 |
| Solved and committed | 189 | 196 |
| Refused | 12 | 9 |
| Results refused by the audit | 0 | 0 |
| Worst Monte Carlo deviation; beyond 3 sigma | 2.62 sigma; 0 | 2.62 sigma; 0 |
| Committed non-empty results tessellated closed, Euler characteristic matching | 183 | 190 (and 6 empty results) |

Committed per step: 40, 37, 35, 30, 28, 26. Chains 11, 19 and 38 run to the end; 26 of the 40 chains run all six steps. The 9 refusals: 6 pairs not proven complete (tangency), 2 three-surface searches not resolved, and chain 18's classification clash, which is also what the "edge leaving through an edge" case of `boolean_test` ends in. In that case an edge of one solid lies in a face plane of the other and is crossed there by an edge of it, and the regions on the two sides of that line are not separated correctly; the clash catches it, so nothing wrong is committed. It is the next degenerate case to resolve.

### The classification clash

The last degenerate refusal of the curved chains (chain 18) and the "edge leaving through an edge" case of `tests/op/boolean_test.cpp` ended in the same message: two regions that share an uncut edge piece disagree about inside and outside. The test case shows the cause. A tilted box B, 0.5 ≤ x ≤ 1.3, 1 ≤ y + z ≤ 3, 0 ≤ z − y ≤ 2, has a face in the plane y = z, and the unit box's edge y = 1, z = 1 lies in that face for 0.5 ≤ x ≤ 1. The unit box's two faces at that edge lie on the two sides of B's face: z = 1 inside B, y = 1 outside. The edge is on B's boundary, so classification must not carry inside/outside across it; but it did, because the rule that stops it there (a straight edge lying exactly in a plane that holds faces of the other solid) had been applied only where a curved surface touched a plane somewhere in the operation. "Inside" spread from z = 1 to y = 1, met "outside" along the other part of the edge, and the operation was refused.

The rule now applies to every such edge. Not crossing an edge is always safe: its regions are then decided by their other neighbours or by a ray test. The case solves with the closed forms, intersection 1/8 (x ∈ [0.5, 1] times the triangle z ≥ y, y + z ≥ 1 of the unit square), union 2.475, difference 0.875, all to 1e-12, and `boolean_test` checks them. Chain 18 runs to the end.

**On the curved chains (§3):** 197 of 205 steps committed (196 before), 27 chains run all six steps, and the 8 refusals left are 6 pairs not proven complete and 2 three-surface searches that did not finish. No refusal from a degenerate configuration of vertices and edges is left. Worst Monte Carlo deviation 2.62 sigma, none beyond 3 sigma; all 191 committed non-empty results tessellate closed with the Brep's Euler characteristic.

### Tangency: a plane along a ruling inside a cylinder

Six refusals of the curved chains were tangency, of two kinds:
- **A plane tangent to a cylinder along a ruling inside the patch** (chains 16 and 35). Contact along a side of the patch had been decided exactly since "Tangential contact with planes"; along an interior line it had not.
- **Two hill tops of equal base height that share a boundary line** (chains 7, 29, 32 and 33). The difference of the two height functions factors (for example X(1 − X)·h(y)): the intersection is the shared line together with a curve that crosses it. Where they cross, the surfaces are tangent, and no interval certificate can work there.

This round handles the first kind.

**Deciding it exactly** (`line_contact`, `op/contact.hpp`). The plane function g(u, v) of a patch of degree 2 in one parameter has, on each control row, the form c0 (1 − u)² + 2 c1 u (1 − u) + c2 u². It has a double root exactly when c1² = c0 c2, and is then K (u − r)² with K = c0 − 2 c1 + c2 and r = (c0 − c1)/K. If every row has a double root at the same r (exact cross-multiplication) and every K has one strict sign, then g(u, v) = (u − r)² h(v) with h strictly signed. The patch meets the plane exactly on the line u = r and lies strictly on one side of it everywhere else. The line's parameter r is enclosed by summing the exact expansions' components in interval arithmetic.

**Carrying it through.**
- The surface pair contributes no curve.
- An edge in the plane meets the surface only on the line, and an edge on the surface meets the plane only there; neither crosses, so neither is searched.
- Where the contact line also lies in another plane of the operation (exactly: its end points, computed exactly at the rational parameter r, are in that plane), a curve of that plane and the surface runs along the line. A point of the line, proven on the certified curve, places it on the operands' edge along that line (a side or an earlier intersection curve), which is then an existing edge.
- Regions of the two solids that touch along the line inside both are recorded, and a result that keeps both reports a line contact and is not manifold.

Two smaller rules came up on the way: an edge whose side lies exactly in a plane that a surface touches only along its sides is recognised as lying in that plane even when its carrier is not of that plane's class (the circle of a cap left on the cap's old rational carrier), and the same is accepted where curves at a vertex are ordered against such an edge.

`tests/op/exact_test.cpp` checks two configurations against closed forms: a box whose face x = 1.25 touches the outer cylinder of a sector (0.75 ≤ r ≤ 1.25, t ∈ [−1/4, 1/4]) along the ruling y = 0 from outside (intersection empty, union the sum, with one line contact), and the same box against the sector 1.25 ≤ r ≤ 2, where the tangent plane runs through the sector's material (intersection 0.002096015459, the part of the box with r ≥ 1.25). Without the rule both are refused as tangency.

**On the curved chains (§3):** 200 of 208 steps committed (197 of 205 before). Chain 16 runs to the end. Chain 35 passes its tangency and two further stages and stops where a new curve leaves a vertex tangent to an existing curved edge (the cylinder's circle at the plane the hill top touches): two curves tangent at one vertex, still refused. The 7 refusals left: 4 hill-top contacts of the second kind, chain 35, and 2 three-surface searches.

### Tangency: two hill tops crossing along a shared line

The second kind of tangency above: two hill tops of equal base height whose bumps both vanish along a boundary line they share (chains 7, 29, 32 and 33). Along that line the two surfaces coincide; away from it they cross along a curve that leaves the line. Where the curve meets the line the surfaces are tangent, so the ordinary search cannot finish there, and no interval certificate can.

**Factoring the difference exactly** (`graph_crossing`, `op/boolean.hpp`). Both tops are graphs over the same axes: x affine in u, y affine in v, z a polynomial, weights 1. B's height is restricted to A's square exactly (rational blossoming in expansion arithmetic), giving the height difference D = S·zA − zB with a positive scale S. D vanishes on the shared sides, and each such side is divided out exactly, leaving k: the intersection is the shared lines together with the curve k = 0.

**Certifying the curve.** The curve k = 0 is traced as the intersection of A with a shifted surface C = B + κ e_z, where κ = (D − k)/S is computed exactly and enclosed. A and C differ in height by a positive multiple of k, so their intersection is exactly the k-curve, and it is transversal along it: the ordinary certified search proves it complete (`SsiSystem` now takes the shifted surface as an interval patch; only its enclosures enter the proofs). The curve's ends on the shared sides are the roots of k there, isolated by Bernstein subdivision with a monotonicity proof. They become vertices with parameters on both tops and on the side faces that share those lines.

**Supporting rules.** The shared line is a side of both tops whose control points lie on a line, so it is recognised as an edge lying in the other top (edges in curved surfaces now include such sides of patches of degree above 1). Vertices on it get parameters on every surface that shares the side. A vertex on a face's own boundary whose carrier is planar gets its plane parameters from its box when no earlier step gave it any.

**Refused exactly: a crossing curve that touches the shared line.** If k along a shared side has a double root, the crossing curve touches the line there instead of leaving through it: both tops have the same slope at that point. For a quadratic side this is decided exactly (zero discriminant, root in the side) and refused with that reason. A discriminant that is not zero but tiny against its terms (relative 2^-40) is refused too: the data's rounding blurs the touch into two crossings, or none, closer together than any enclosure can separate. Chain 7 is such a case. Bumps of 0.2 over a width of 0.5 and of 0.3 over 0.75 have equal slopes in decimal, but 0.2 and 0.3 are not exact in binary.

**Tessellation of a lens.** The difference of the two hills contains a thin lens between the tops, pinching to nothing along the shared lines and the crossing curve. Its two faces have the same outline and both chose the same chord between two boundary points, so the mesh touched itself along a segment where the surfaces are apart (one edge used by four triangles, Euler characteristic 3 instead of 2). The tessellator now splits a chord that two faces share and that is not a piece of a Brep edge, at each face's own parameter midpoint (`TessStats::lens_splits`). No earlier mesh had such a chord: the tessellation hash is unchanged.

`tests/op/exact_test.cpp` checks two configurations:
- **Two hills crossing along a shared line**, with the closed form: intersection 0.066965819874, union 0.278867513459, difference 0.142756402348, each to 1e-11. All three results tessellate closed.
- **Dyadic hills whose crossing curve touches the shared line** (bumps 0.25 over 0.5 and 0.375 over 0.75, equal slopes at x = 0.75): refused with the reason above.

Without `graph_crossing` both checks fail (the first is refused as tangency, the second with the older message).

**On the curved chains (§3):** 209 of 214 steps committed (200 of 208 before), and 29 chains run all six steps. Chains 29 and 32 run to the end. Chain 7 is refused with the reason above. Chain 33 passes its crossing and stops at a different case: an arc of the earlier result (a cylinder against the first top) ends exactly on the shared line, so it meets the second top at its end vertex, and the search along arcs cannot place a root at an arc's end. Worst Monte Carlo deviation 2.62 sigma, none beyond 3. All 202 committed non-empty results tessellate closed (4 non-manifold), with the Brep's Euler characteristic and components. The 5 refusals left: the blurred touch of chain 7, two curves tangent at a vertex (chain 35), and three three-surface searches (chains 23, 31 and 33).

### Tangent points, certified (roadmap phase 2, first step)

Where two surfaces touch at a single point, their intersection is a singular root of A(u,v) = B(s,t): interval methods can prove that no root exists but never that one does, since an arbitrarily small perturbation removes it. Every refusal of the kind "a box whose bottom touches the top of the hill in one point" comes from this. The survey ([SURVEY.md](SURVEY.md)) points to deflation, which turns the question into a regular one. `solve/contact_point.hpp` (new) solves the square system in (u, v, s, t, μ)

- A(u,v) − B(s,t) + μ N_A(u,v) = 0 (three equations, N_A = A_u × A_v),
- N_A · B_s = 0 and N_A · B_t = 0 (parallel tangent planes),

whose root at a contact with definite relative curvature is regular, so a 5-dimensional Krawczyk test proves it the only one in a box. μ is the signed distance between the facing points along A's normal. Then:
- **μ proven nonzero.** With the relative curvature D = II_B − II_A definite (decided with intervals over the box), B's height over A near there is μ|N_A| + ½ xᵀDx + …: the surfaces are **separated** there if μ and D have the same sign, and cross in a **small loop** around the point otherwise.
- **μ's enclosure contains 0.** A dyadic candidate near the root is checked in expansion arithmetic, in homogeneous form, so rational patches are covered: the two points are proportional, and B_s, B_t lie in A's projective tangent plane span(P, P_u, P_v) (two 4×4 determinants vanish). The candidate is then a root with μ = 0, the only one in the box: the surfaces **touch** exactly there, at an isolated point if D is definite, at a **singular crossing** of two intersection curves if D is indefinite.
- Otherwise it is **undecided**, with the reason.

`tests/solve/contact_point_test.cpp` (new, the 22nd test program), a paraboloid cap z = 1 − (x² + y²)/2 against horizontal planes:

| Case | Result |
|---|---|
| Plane through the apex | Isolated touch, exactly at (u, v, s, t) = (1/2, 1/2, 3/8, 3/8) |
| Plane 2^-20 above | Separated; μ = 2.384e-07 = 2^-20 / \|N_A\| |
| Plane 2^-20 below | Small loop |
| Plane 2^-44 (6e-14) above | Separated; μ ∈ [1.28e-14, 1.56e-14] |
| Saddle z = (x² − y²)/2 against its tangent plane | Singular crossing, exactly at the saddle point |
| The cap with every weight 2 (rational form), plane through the apex / 2^-20 below | Isolated touch / small loop |
| A rational cap (weights 1, 2, 1 in each direction, apex at height 5/4) against z = 5/4 / 2^-30 above | Isolated touch / separated |
| The cap shifted by 1/3 (apex at non-dyadic parameters) | Undecided: the root is unique in its box, but no exact candidate exists |

With the exact check disabled, 6 of its 15 checks fail. Its statements are local: "separated" and "loop" describe the neighbourhood of the certified pair, and showing that nothing else intersects in a region remains the search's job.

**In the Booleans: a curved face touching a plane at one point.** Where the surface/surface search of a curved patch and a plane cannot finish, and the boxes it leaves are small, the imprint now tries a point contact (`point_touch`, `op/boolean.hpp`), decided exactly on the plane function g of the patch (its exact Bernstein form, homogeneous, so rational patches are covered):
- a dyadic point p* near the boxes left where g, g_u and g_v vanish, in expansion arithmetic;
- g's Hessian definite over a dyadic box containing every box left (its second derivatives restricted exactly to that box, their coefficient hulls giving a positive determinant): g is strictly convex or concave there, so p* is its only zero.

Since the search proved everything outside those boxes, the pair meets in the curves found and at p*. The faces are not split there. If the point lies strictly inside a face of each solid, the two regions are recorded as touching at a point, and a result that keeps both reports a point contact and is not manifold; a point on a face's boundary is refused.

The curated case of §3 ("not in general position") is now solved: a box whose bottom z = 1.3 touches the top of the hill z = 1 + 1.2·B₁(u)B₁(v) at its apex (u, v) = (1/2, 1/2), where the top's height 1 + (2.2 − 1)/4 is the double 1.3 exactly. The intersection is empty, the union is both solids touching at the apex (3.933333333333 = 2.8 + 1.133333333333, two components, one point contact, not manifold), and the difference is the hill (1.133333333333 = 1 + 1.2/9), all to 1e-12 (`tests/op/boolean_test.cpp`). The same box moved so that the apex lies on an edge of its bottom face is still refused. With the point-contact proof disabled, the case is refused as before and 2 checks fail.

**Two curved faces touching at one point.** Two curved patches have no single polynomial whose zeros are their intersection, so the proof goes through the normal offset from one to the other (`curved_touch`, `op/boolean.hpp`; `unique_touch_in`, `solve/contact_point.hpp`):
- the contact is proven exact by `classify_contact` (a regular root of the deflated system, an exact dyadic candidate in homogeneous expansion arithmetic);
- for every (u, v) in a box U around the region the search left, the normal line A(u,v) + λ N_A(u,v) meets B exactly once in a box Y that contains the B-parameters of that region (a Krawczyk test with (u, v) as interval parameters), so every intersection point there has λ(u, v) = 0;
- λ vanishes at the contact with zero gradient (the tangent planes are parallel), and its Hessian, by implicit differentiation with third derivatives of A and second derivatives of B enclosed over U × Y, is definite: λ is strictly convex or concave on U, and the contact is its only zero.

The case of `tests/op/boolean_test.cpp`: a hill whose top rises to 1.25 at (0.5, 0.5) (bump 1/4) and a block over [−0.5, 1.5]² up to z = 2 whose bottom dips from 1.5 down to 1.25 there. The contact is proven at the parameters (1/2, 1/2) on both faces. The intersection is empty, the union 3.555555555556 = (1 + (1/4)(4/9)) + (2 + 4 (1/4)(4/9)), one point contact, not manifold, and its mesh is closed; the difference is the hill, all to 1e-12. With the convexity proof disabled, the case is refused and 2 checks fail.

**Near misses.** Two surfaces that come close without meeting leave the search the same unresolved region. How close it can still prove the gap empty was measured with the hill (apex 1.25) and a box above it, the gap 2^-k: the search alone proves it for k ≤ 46 (gaps down to 1.4e-14) and gives up from k = 48 on (at k = 53 the gap rounds away and the plane touches the apex exactly, a point contact). The interval classifier resolves gaps no finer than the search, so below that only exact arithmetic helps:
- **A curved face against a plane, exactly.** The plane function g is an exact polynomial in the data. Newton in doubles finds its critical point p (a double, hence an exact dyadic point); g, g_u and g_v at p are computed exactly; with a lower bound m on the eigenvalues of g's Hessian over the region (its exact restriction's coefficient hulls), Taylor's theorem gives g ≥ g(p) − |∇g(p)|²/(2m) there for a convex g (mirrored for a concave one). If that keeps g(p)'s sign, the face stays off the plane over the region. The hill against the box 2^-50 (9e-16) above its apex is now solved: the solids are disjoint, union 4.111111111111 to 1e-12.
- **Two curved faces, with intervals.** `solve::separated_in` bounds the normal offset the same way (λ and its gradient enclosed at the facing pair, the eigenvalue bound over the region). It is exact in its logic but limited by interval resolution, so in the Boolean it helps only where the search gives up for other reasons; its unit tests (`contact_point_test`) prove a plane and an upside-down cap 2^-20 apart separated over a box, and refuse to call the 2^-20 loop separated. The dipping block 2^-50 above the hill is still refused.

Red checks: without the exact plane proof the 2^-50 case is refused (2 failures); without the sign test the loop is wrongly called separated (1 failure).

Not covered yet: a contact on a face's boundary (it would need a vertex there), near misses of two curved faces below interval resolution (they need exact algebra for curved pairs), and tiny loops (a crossing too small for the search to trace).

### A compiler fault

While building the quadric tests, `on_quadric` reported that an exact quarter cylinder was not on its cylinder, but only in the `/O2 /fp:strict` build of the test program. With `/Od`, `/O1`, `/Og`, `/fp:precise`, or with inlining off (`/Ob0`), the result was right, and so it was in smaller programs with the same headers and flags. The investigation:

- **Not the arithmetic:** `on_quadric` and every function it calls compiled to identical instructions in the failing program and in a passing one, and the patch passed in was correct.
- **Not the environment:** the SSE control register held its default (0x1F80: round to nearest).
- **Not undefined behaviour in the kernel:** a checked Debug build (`/RTC1`, checked iterators) was clean.
- **The fault itself:** the quadric matrix built in the large `main` held a stale 1.0 in an entry that the source never writes. The struct had a zero default member initialiser and the factory assigned only the nonzero entries. This is a code-generation fault of MSVC 19.51 (`/O2 /fp:strict` with inlining, in a large function).

The factories now write all sixteen entries from one table, and the test checks the matrices they return. The rest of the kernel's exact arithmetic is cross-checked independently throughout (exact predicates against expansion evaluation, volumes against Monte Carlo and closed forms, identities), and none of those checks has failed. But this shows that the rigour depends on the compiler, which is already stated in §10.

### Point classification against a trimmed face

The trimmed top face of the "hill cut by a plane" case, bounded by one closed intersection loop, was queried with `op::locate` and compared with the analytic loop (`tests/op/boolean_test.cpp`).

| Query points | Wrong | Undecided |
|---|---|---|
| 4,000 uniformly random in the parameter square | 0 | 0 |
| 800 placed 1e-6 from the loop, half on each side | 0 | 0 |

The answer is a winding number computed with exact predicates on a polygon that is proven to be equivalent to the true loop for that query point. A point exactly on the curve stays undecided.

### Configurations that are not in general position

Each of these stops with a message naming the two faces involved. None returns a result.

| Configuration | Outcome |
|---|---|
| A box whose bottom touches the top of the hill in one point | Solved since tangent points (§3): the apex is proven the only contact; union with one point contact, not manifold |
| The same, with the apex on an edge of the box's bottom face | Not solved: the edge touches the surface |
| A tilted box whose edge leaves the unit box exactly through one of its edges | Solved since the classification clash (§3): intersection 1/8 |
| An edge crossing an edge at one point, all faces otherwise transversal | Solved since edges through edges (§3): the crossing point (1, 1, 0.5) is exact; intersection 1/48 |
| Two curved faces meeting tangentially along a shared edge | Refused by the fan check around the edge |

Two boxes touching only along an edge or at a point were in this table until the non-manifold milestone; they are now solved (above). A flat face touching a curved one along a line (a plane against a cylinder's line), including a curved and a flat face that meet tangentially around a shared edge, was in it until the tangential-contact work; it is now solved exactly (above). Tangency at a point, and between two curved surfaces, are not, except two hill tops (graphs over the same axes) crossing along a boundary line they share (§3).

The first three are caught by the surface/surface stage. The tilted box's faces are not exactly planar after rounding, so the exact coplanarity tests rightly do not apply to it. Coincident curved surfaces that are not identical patches are also not resolved.

## 4. Variational G3 fill

A square hole in the surface z = 0.2 sin(2.5x + 0.3) cos(2y) + 0.1xy is filled by a degree-7 B-spline patch with jets up to order 3 prescribed on all four sides and the fourth-order energy minimised. Jumps are the largest values over 156 seam points of the parameterisation-independent invariants (§1.7 of the blueprint).

| Control points | Free | G1: normal angle, rad | G2: jump of b(ν,ν) | G3: jump of C(ν,ν,ν) | Largest other component | Solve, ms |
|---|---|---|---|---|---|---|
| 12 × 12 | 16 | 3.8e-08 | 1.0e-06 | 5.3e-05 | 2.5e-05 | 2 |
| 16 × 16 | 64 | 4.3e-10 | 2.3e-08 | 2.1e-06 | 1.1e-06 | 7 |
| 20 × 20 | 144 | 2.5e-11 | 2.2e-09 | 2.9e-07 | 1.6e-07 | 15 |
| 24 × 24 | 256 | 3.8e-12 | 4.0e-10 | 7.0e-08 | 4.3e-08 | 27 |
| 32 × 32 | 576 | 2.6e-13 | 3.2e-11 | 8.2e-09 | 6.3e-09 | 68 |

- Observed convergence orders in the knot spacing: 7.4 (G1), 6.4 (G2), 5.4 (G3), consistent with degree 7.
- The fill is a different surface from the original: at the centre it differs by 4.0e-05.
- A cubic surface has zero fourth-order energy, so the fill must reproduce it. It does, to 5.6e-15 (`tests/var/fill_test.cpp`).
- These jumps are sampled in doubles, on a surface that is not polynomial. Certified bounds, for fills whose neighbour is a polynomial patch, follow below.
- The joins are C^3 in a shared parameterisation. Joins that need a reparameterisation (gauge functions) are not implemented.

### Certified seam bounds

`geom/seam_certificate.hpp` (new) bounds the jumps over a whole seam, rigorously, for two patches in one common parameterisation (the setting of the fill), polynomial or rational. It is the first step of the roadmap in [SURVEY.md](SURVEY.md), and it turns this section's sampled check into a certificate. On each part of the seam:
- **Residuals in Bernstein form.** Every jet entry D^(k,l), k + l ≤ 3, is a polynomial along the seam: of the surface for polynomial patches, of its homogeneous form for rational ones. The residual Δ = J_B − J_A is formed on the control points (exactly zero where the patches agree), differentiated and restricted by de Casteljau steps in interval arithmetic; its largest Bernstein coefficient bounds it on the part.
- **Mean value.** Each jump is H(J_A, J_B) with H(J_A, J_A) = 0, so it equals ∂H/∂J_B · Δ at an intermediate point. The gradient is enclosed by forward differentiation in interval arithmetic, with A's jet over the hull of its coefficients and B's over the hull of both, and the jump is bounded by Σ |∂H/∂J_i| ε_i, rounded upward. For a rational patch the projection from homogeneous coordinates is part of H (40 variables instead of 30).
- **Parts.** The seam is cut into equal parts; the enclosure of the gradient shrinks with the part, so more parts give tighter bounds. A part on which the orientation, a weight or the metric is not proven is halved, up to ten times.

A bound computed this way shrinks with the residual, however small it is. Enclosing each surface's invariants separately could not resolve a jump smaller than the variation of the invariants along a part (here 1e-9 against variations of order 1).

`tests/geom/seam_certificate_test.cpp` (new, the 21st test program) checks:

| Case | Certified | Sampled maximum (401 points) |
|---|---|---|
| Patches equal to third order (rows 0 to 3 shared) | every jump ≤ 5e-162 | 0 |
| Row 1 moved by about 1e-4: G1 jump | normal angle 8.629e-04 rad | 8.616e-04 |
| Row 2 moved: G2 jump | b(ν,ν) 4.357e-03, normal angle 1.5e-14 | 4.224e-03, 0 |
| Row 3 moved: G3 jump | C(ν,ν,ν) 1.736e-02, b(ν,ν) 6.1e-14 | 1.657e-02, 0 |

(64 parts per seam.) On 24 seams between random patches with differences of 1e-6 to 1e-2, on all four sides, every certified bound covers every sampled jump. A neighbour whose normal points the other way is refused with a reason.

Rational patches (weights 0.7 to 1.5, degree 5), with the residual in homogeneous form:

| Case | Certified | Sampled maximum |
|---|---|---|
| Equal to third order | every jump ≤ 5e-162 | 0 |
| G1 jump | normal angle 7.342e-04 rad | 7.180e-04 |
| G2 jump | b(ν,ν) 4.749e-03 | 2.480e-03 |
| G3 jump | C(ν,ν,ν) 3.586e-02 | 1.414e-02 |
| Two adjacent sectors of one cylinder (exact rational carriers, the second continued backwards over the first's parameter) | gap 1.7e-14, angle 9.4e-14, b(ν,ν) 4.0e-13, C(ν,ν,ν) 5.1e-13 | — |

(256 parts for the jumps.) The rational bounds are looser: the projection to the surface combines terms that cancel, which interval enclosures over a part do not see, so the gradient is overestimated; G1 stays within 2 %, G2 and G3 within a factor 2.5, and secondary components can be far looser. On 16 random rational seams (differences in the weights too, all four sides) every bound covers its samples.

The G3 fill of a hole in a polynomial surface of degree 9 (which the fill, reproducing cubics, does not reproduce), certified on every boundary piece of the fill against that surface, 32 parts per seam:

| Fill | Seams | G1, rad: certified (sampled) | G2: certified (sampled) | G3: certified (sampled) | Gap, certified |
|---|---|---|---|---|---|
| 12 × 12 | 20 | 6.37e-06 (6.00e-06) | 2.63e-04 (1.27e-04) | 2.87e-02 (6.03e-03) | 9.98e-07 |
| 16 × 16 | 36 | 1.21e-07 (1.20e-07) | 5.15e-06 (3.72e-06) | 6.55e-04 (3.52e-04) | 9.50e-09 |
| 24 × 24 | 68 | 1.08e-09 (1.05e-09) | 9.39e-08 (8.04e-08) | 1.97e-05 (1.50e-05) | 6.11e-11 |

The certified G1 bound is within 6 % of the sampled maximum at every resolution; G2 and G3 come within 17 % and 31 % at 24 × 24. They are looser on the coarse fill, where the positional gap is largest and the bound adds the contributions of all 30 jet coordinates without cancellation. Red checks: with the derivative residuals left out of the bound, 999 of the test's checks fail; with the weight residual left out of the rational bound, 10 fail. The test (1,603 checks) takes 17.6 s in Release and 43 s in Debug.

An adaptive refinement was tried and dropped: halving the parts whose bound exceeded what the gradient at their midpoint gave spent its budget on parts that never set the maximum (348 s for the test, some bounds looser than with equal parts).

Limits: a common parameterisation (no gauge functions); rational patches must agree in homogeneous form up to the residual (pieces of one rational parameterisation do), and their bounds are looser; the fill's Bezier pieces are taken as polynomial (knot insertion in doubles leaves their weights within rounding of 1, and the pieces are defined with weight 1).

## 5. Energy-order rule

Blueprint §1.7 claims that prescribing jets up to order k is well posed only for the energy of order k + 1. Test: prescribe G2 jets (k = 2) and compare the thin-plate energy (order 2) with the matched energy (order 3), reading off the third-order energy E3 of each solution.

| Control points | E3 of the thin-plate solution | E3 of the matched solution |
|---|---|---|
| 12 × 12 | 21.86 | 11.1847 |
| 18 × 18 | 38.24 | 11.1847 |
| 27 × 27 | 62.79 | 11.1847 |
| 40 × 40 | 98.15 | 11.1847 |

The thin-plate solution's third-order energy grows without bound under refinement (a boundary layer), while the matched energy returns the same value at every resolution. The claim holds on this example.

## 6. The G3 seam reduction

Blueprint §1.7 derives that across a G2 seam the only component of ∇b that can jump is C(ν,ν,ν). Checked on a seam whose two sides are parameterised differently and nonlinearly (`tests/geom/continuity_test.cpp`):

| Seam | Result |
|---|---|
| Same surface on both sides, reparameterised | Every invariant jump ≤ 2.4e-12, while the parametric third derivative differs by 0.79 |
| G1 but not G2 (extra term 0.4x²) | Only b(ν,ν) jumps; it matches the predicted value at all 21 points to 1e-7 |
| G2 but not G3 (extra term x³(0.5 − 0.8y)) | C(ν,ν,ν) jumps by up to 2.76 and matches the prediction to 4.6e-12; all other components ≤ 7.0e-13 |

The closed formula for ∇b also agrees with its finite-difference definition to 3.3e-08 on a random rational bicubic patch. This is numerical confirmation on examples, not a proof.

## 7. Tangency early warning

Blueprint §2.2 proposes the observed convergence order as a warning that precedes failure. The same Newton iteration was run towards a transversal and a tangential intersection (`tests/diag/rate_test.cpp`):

| Target | Observed order |
|---|---|
| Plane cutting the hill | 1.96 |
| Plane touching the hill | 1.00 |

## 8. Determinism

The benchmark ends with a hash over its numeric results: predicate outcomes, the certified curve points of the curated cases, curve and segment counts of both random families, the fill control points and the energies. Timings are excluded.

| Build | Runs | Result hash |
|---|---|---|
| Release (`/O2`) | 3, including one after a rebuild | `95e6ec5ff73a8784` |
| Debug (`/Od`) | 3, two of them after the solver was refactored for the Boolean work | `95e6ec5ff73a8784` |

The results are bitwise identical across runs and across optimisation levels. The refactoring (splitting the Krawczyk test into an operator and a test) left every earlier result unchanged.

The two Boolean sections have their own hashes. The first covers the three volumes and the entity counts of every curated and random pair (304 imprints); the second covers the volume and entity counts of every step of the 40 random chains.

| Hash | Build | Result |
|---|---|---|
| Boolean pairs | Debug, 2 runs of the first implementation (whole-patch operands only) | `3e545d18b11ee35f` |
| Boolean pairs | Release, current implementation | `3e545d18b11ee35f` |
| Random chains | Debug, current implementation before one change that only skips a redundant check | `472f475546f157b1` |
| Random chains | Release, current implementation | `472f475546f157b1` |
| Random chains with reused operands | Release | `3fee6db2fdaa1305` |
| Random chains of grid-aligned boxes | Release, before the vertex-manifoldness check | `85146b6ca8a3ebc6` |
| Random chains of grid-aligned boxes | Release, current (chain 30 now stops at its pinched step 4) | `46aabe1a9bb9dbfd` |
| Tessellation of every committed result (§3) | Release, before non-manifold results | `45afe508e7116fe6` |
| Algebraic identities (§3) | Release, before non-manifold results | `743be27af78abce8` |
| Random chains of grid-aligned boxes | Release, current (non-manifold results committed) | `cfafddb9df9fd21f` |
| Tessellation of every committed result (§3) | Release, current | `4be3b23e245a7be2` |
| Algebraic identities (§3) | Release, current | `896b9c04bb9ba003` |
| Random chains of grid-aligned boxes | Release, after exact maps between coincident planes (narrower vertex enclosures) | `3756858bd81d5b4d` |
| Tessellation of every committed result (§3) | Release, after exact maps between coincident planes | `a237e0c287738fc5` |
| Chains with curved operands (§3) | Release, from a run made before the last curved-chain fixes | `175866c1d5dfd6b3` |
| Chains with curved operands (§3) | Release, current (81 committed; see §3) | `fc823b4c38a56345` |
| Scaling with the number of faces (§3) | Release | `72b9f222ff9ffcb6` |
| Boolean pairs | Release, plane pairs certified directly | `3647b207b5a5c82b` |
| Random chains | Release, plane pairs certified directly | `3b038e88ae34b0fb` |
| Random chains with reused operands | Release, plane pairs certified directly | `93602dd30bf2f0b0` |
| Random chains of grid-aligned boxes | Release, plane pairs certified directly | `4b6ccecf95b83c1c` |
| Tessellation of every committed result (§3) | Release, plane pairs certified directly | `6a2a53712dc4acab` |
| Chains with curved operands (§3) | Release, plane pairs certified directly | `acfb30a14358fea1` |
| Scaling with the number of faces (§3) | Release, plane pairs certified directly | `3a32233dc2ef1304` |
| Chains with curved operands (§3) | Release, tangential contact resolved (100 committed) | `5f149b0ee8cff971` |
| Chains with curved operands (§3) | Release, flat faces on rational carriers moved (128 committed) | `723491a0cd373b12` |
| Chains with curved operands (§3) | Release, contact cases left after the flat-face work (173 committed) | `5116a8b5b834b82d` |
| Chains with curved operands (§3) | Release, straight curves between flat surfaces, exact points off (175 committed) | `63e18895eb39c6a4` |
| Chains with curved operands (§3) | Release, exact points (178 committed) | `ab03e5ee9b0cf590` |
| Boolean pairs | Release, exact points | `c30173d6104eb6c1` |
| Random chains of grid-aligned boxes | Release, exact points | `d0ef9ac8a08dbc49` |
| Tessellation of every committed result (§3) | Release, exact points | `9956a4033e5aafd7` |
| Chains with curved operands (§3) | Release, edges in curved surfaces (180 committed) | `cadd9f55d7dd34c2` |
| Chains with curved operands (§3) | Release, ordering cases (187 committed) | `506cade8c71b1a1a` |
| Chains with curved operands (§3) | Release, ordering cases, exact points off (177 committed) | `bb6becadcb232d18` |
| Chains with curved operands (§3) | Release, edges through edges (189 committed) | `54bed3f149721d1d` |
| Random chains of grid-aligned boxes | Release, edges through edges (one exact parameter on planes) | `c5f2137f04d933f9` |
| Chains with curved operands (§3) | Release, single cases (196 committed) | `ca48d5915d737087` |
| Chains with curved operands (§3) | Release, single cases, exact points off (181 committed) | `076cc79a80b0285b` |
| Chains with curved operands (§3) | Release, classification clash (197 committed) | `004720a80878c448` |
| Chains with curved operands (§3) | Release, tangency along a ruling (200 committed) | `099834461ccbd100` |
| Chains with curved operands (§3) | Release, hill-top crossings (209 committed) | `be8b56d82aa17b20` |

After the coplanar work, again after the tessellation work and the vertex-manifoldness check, and again after the identity work (grown uniqueness regions in the edge/surface and three-surface searches, closed-edge unification, empty operands), the core, single-operation, chain and reuse hashes were re-run in Release and are unchanged, and so are the grid-chain and tessellation hashes of the previous milestone. The grid-chain hash changed first because the newly refused pinched result ended its chain, and then again when non-manifold results became valid and the grid chains ran on. With non-manifold results the core, single-operation, chain and reuse hashes are again unchanged: none of those sections produces a contact of that kind. After the exact-quadric work every hash in this table is unchanged from the non-manifold milestone (core, Boolean pairs, chains, reuse, grid chains `cfafddb9df9fd21f`, tessellation `4be3b23e245a7be2`, identities `896b9c04bb9ba003`).

After the spatial-index work (§3, scaling) every hash in this table is reproduced in Release: core `95e6ec5ff73a8784`, Boolean pairs `3e545d18b11ee35f`, chains `472f475546f157b1`, reuse `3fee6db2fdaa1305`, grid chains `3756858bd81d5b4d`, tessellation `a237e0c287738fc5`, identities `896b9c04bb9ba003`. The curved-chain value recorded before came from a run made before the last fixes of that work; the current value is the same with or without the pruning of the spatial index (§3). All 17 tests pass in Release and Debug.

Certifying plane pairs directly (§3) changes the Boolean hashes, because the curves between flat faces are now segmented differently; the core and identity hashes are unchanged. Every count and check in those sections is identical, and with `BooleanOptions::plane_lines` off (`ak_bench ... nolines`) the earlier hashes above are reproduced exactly. The identity test (`tests/op/identity_test.cpp`) now evaluates all 25 cases per identity (417 checks) where it evaluated 23 (388 checks); the two pairs that were not solved before solve, and all agree. Which change of the spatial-index work made them solvable was not isolated.

The faster three-surface searches (§3) reproduce every hash of the "plane pairs certified directly" rows exactly.

Resolving tangential contact (§3) changes only the curved-chain hash, because more of those steps now solve. Every other hash of the "plane pairs certified directly" rows is reproduced exactly in Release (core `95e6ec5ff73a8784`, Boolean pairs, chains, reuse, grid chains, tessellation, identities `896b9c04bb9ba003`, scaling): the new code runs only where the earlier code refused. Limiting the tracing of the enlarged second search left the curved-chain hash unchanged.

Moving flat faces off rational carriers (§3) again changes only the curved-chain hash; every other hash is reproduced exactly in Release. The full benchmark takes 925 s (1,014 s after the tangency work, 882 s before it). The curved chains were also run in the Debug build (checked iterators, 578 s): no library assertion, and the same hash `723491a0cd373b12` as in Release.

The contact cases resolved after that (§3) again change only the curved-chain hash. The full benchmark was run on the final sources and reproduces every other hash exactly: core `95e6ec5ff73a8784`, Boolean pairs `3647b207b5a5c82b`, chains `3b038e88ae34b0fb`, reuse `93602dd30bf2f0b0`, grid chains `4b6ccecf95b83c1c`, tessellation `6a2a53712dc4acab`, identities `896b9c04bb9ba003`, scaling `3a32233dc2ef1304`. That includes the two changes that are not confined to steps the earlier code refused: the sharper point location decides more queries, and straight sides of degree 2 are now matched by their end vertices. Neither changed any earlier result. All 19 test programs pass in Release and in Debug on those sources. This time the curved chains were run in Release only; the Debug run covered the new code through `flat_test`, which repeats several of the chain steps this round fixed.

Exact points (§3) change four hashes, because vertices with round coordinates now have exact parameters where they had enclosures, and volumes then differ in their last bits: the Boolean pairs (the curated case of two offset boxes; the 300 random pairs have no such vertex), the grid-aligned chains, the tessellation of their results, and the curved chains. Two full benchmark runs on the final sources, one with `BooleanOptions::exact_points` on and one with it off (`ak_bench noexact`), show that nothing else changed:

- **Off:** every hash outside the curved chains is the recorded one, bit for bit: core `95e6ec5ff73a8784`, Boolean pairs `3647b207b5a5c82b`, chains `3b038e88ae34b0fb`, reuse `93602dd30bf2f0b0`, grid chains `4b6ccecf95b83c1c`, tessellation `6a2a53712dc4acab`, identities `896b9c04bb9ba003`, scaling `3a32233dc2ef1304`.
- **On:** the core, chain, reuse, identity and scaling hashes are again those. In the sections whose hash changed, every count, audit and Monte Carlo deviation printed is the same as with the option off, and so are the totals over all committed results: 6,323 faces, 17,139 edges and 11,426 vertices in the random chains, 4,720, 12,425 and 8,268 in the reuse chains, 2,073, 4,685 and 3,051 in the grid chains. The widest vertex enclosure in the grid chains falls from 4.4e-14 to 8.9e-15.
- **One visible difference:** the meshes of the grid results have a few more triangles with exact points (largest 388 against 386, mean 94 against 93). The solids are the same, as the totals show; the triangulation of a face starts from parameters that are now exact. All 140 meshes are closed, with the Brep's Euler characteristic and components, either way.

All 20 test programs pass in Release and in Debug on the final sources.

Edges in curved surfaces (§3) change only the curved-chain hash. Run again on the final sources, with exact points on, every other hash is the one of the exact-points rows above; with them off (`noexact`), every hash is that of the "exact points off" run, the curved chains' included (`63e18895eb39c6a4`): the new rules come into play only after an exact vertex has let a chain get that far. All 20 test programs pass in Release (108 s) and in Debug (699 s).

The ordering cases (§3) again change only the curved-chain hash, with exact points on and off; every other hash of both runs is the same as above. With exact points off, the curved chains differ from the previous "off" run in exactly one step: chain 17, step 5, which the ordering rule for a curve continuing a line now solves (it does not depend on exact points); 177 steps are committed, one more than before. All 20 test programs pass in Release (116 s) and in Debug (874 s).

Edges through edges (§3): with exact points off, every hash is the same as in the ordering run, the curved chains' included. With them on, the curved-chain hash changes as expected, and so does the grid-chain hash: vertices on planes now get one parameter exact where only one is dyadic, which moves volumes in their last bits. Every count, audit, Monte Carlo deviation and vertex width printed for the grid chains is the same, and so are the totals of faces, edges and vertices; the tessellation hash is unchanged. All other hashes are those of the ordering run. All 20 test programs pass in Debug (857 s) and in Release on the final sources (the updated `boolean_test` was run on its own after the suite: 652 checks).

Single cases (§3) change only the curved-chain hash, with exact points on and off; every other hash of both runs is that of the edges-through-edges run, and the grid chains print the same table. With exact points off, the curved chains differ in one step: chain 38, step 3 (the arc left on another flat carrier, a rule that does not depend on exact points) now solves, and the chain runs on to commit 181 steps. All 20 test programs pass in Release (129 s) and in Debug (837 s).

The classification clash (§3) widened a rule that applies in every operation, so all sections were compared. With exact points on, only the curved-chain hash changes; with them off, every hash is the same as in the single-cases run, the curved chains' included. No earlier result changed: in the benchmark the wider rule acts only in chain 18, which needs exact points to reach that step. All 20 test programs pass in Release (132 s) and in Debug (826 s).

Tangency along a ruling (§3): with exact points on, only the curved-chain hash changes. With them off, every hash is the same as before, the curved chains' included: the tangencies are now decided, but the two steps that reach them (chains 16 and 35) stop at a later stage without exact points, so no committed result differs; only those two refusal messages do. All 20 test programs pass in Release (136 s) and in Debug (963 s).

Hill-top crossings (§3): with exact points on, only the curved-chain hash changes. With them off, every hash is the same as before, the curved chains' included (`076cc79a80b0285b`): no step that the new rules decide is reached without exact points. The tessellation hash is unchanged in both runs, so the new splitting of chords shared by two faces changed no earlier mesh. All 20 test programs pass in Release (151 s) and in Debug (916 s).

Phase 1 of the survey's roadmap (v0.17) changes no result. The outward rounding step of the interval arithmetic (`num/interval.hpp`) is now a step of the bit pattern instead of a call to `std::nextafter`. It returns the same double for every input; `interval_test` compares the two on 100,015 values, including both zeros, subnormals, the largest finite values, both infinities and NaN, with no mismatch. Both benchmark runs, with exact points on and off, reproduce every hash of the hill-top crossing run bit for bit. The certified seam bounds (§4) are new code that no benchmark section calls (and their rational extension, added afterwards, changes no other file). Speed, measured on an otherwise idle machine, best of three runs each:

| Program | Before, ms | After, ms | Change |
|---|---|---|---|
| `ssi_test` | 1,685 | 1,243 | −26 % |
| `csi_test` | 491 | 359 | −27 % |
| `boolean_test` | 4,449 | 3,560 | −20 % |
| `quadric_test` | 4,479 | 3,689 | −18 % |

The Release suite takes 110 s (151 s before; `identity_test` 64 s against 98 s), and the Debug suite 782 s (916 s). The full benchmark took 484 s with exact points and 465 s without; the earlier runs ran beside the Debug suite, so those totals are not compared. All 21 test programs pass in Release and in Debug.

Point contacts in the Booleans (§3, "Tangent points, certified"), a curved face against a plane, then two curved faces, then near misses: both benchmark runs, after each step,, with exact points on and off, reproduce every hash above bit for bit. The new proofs run only where the surface/surface search has already failed, and no Boolean of the benchmark has a point contact; the five refusals left in the curved chains are the same. All 22 test programs pass in Release and in Debug (that run's timings, 189 s and 1,123 s, came from a busier machine: `identity_test` took 114 s in Release against 64 s before on unchanged code).

Timings of the tangency-along-a-ruling round are not comparable with the earlier ones. The full benchmark took 503 s in one run and 679 s in another (925 s before), both with the Debug test suite running beside it; sections whose code and results are unchanged ran in about half their earlier time in the first of those runs (a random-box chain step at depth 6: 204 ms against 362 ms), so the difference is the state of the machine, not the code.

So the Boolean results are bitwise identical between Debug and Release; the rewrite for trimmed operands, and the later rewrite for shared geometry, did not change any earlier result. The reuse chains have been run in Release only. All of this covers one machine, one compiler and single-threaded execution; other platforms and parallel execution are untested.

## 9. Against the blueprint's targets (§5.4)

| Target | State |
|---|---|
| Edge gap bounded by evaluation accuracy, no stored tolerance | Met for single patch pairs (certified bounds from 6.6e-14 to 7.9e-12) and for every edge created by the Booleans of §3 (at most 2.0e-12); no tolerance is stored on any entity |
| No silent wrong answers | Met in everything tested: undecidable and contact cases are reported as such; the independent cross-checks of §2 and §3 found no disagreement. The inconsistencies found during development ended in a report or in a result refused by the audit, not in a wrong result (§3: chained operations, cause 3; reuse chains; the face that was not split, in "Contact cases left after the flat-face work"). That last one came from giving an undecided answer a meaning at one place, which now reports it |
| 100 % of committed results pass the topology audit | Met on what exists: all 900 single-operation results of the random pairs, all 238 chain results, all 223 reuse-chain results and all 163 grid-chain results pass (§3). Non-manifold results (17 in the grid chains) are committed only when closed and oriented and when the faces around every edge with more than two faces alternate between material and empty space; results that fail are refused |
| Algebraic identities give identical topology | Met: all 722 identities (520 in the benchmark, 202 in the tests) on flat, curved, rational and touching operands were evaluated and give the same solid after removing seams and pass-through vertices; volumes equal to 4.4e-15 (§3) |
| Crack-free tessellation, shared edge samples (§5.1) | Met on all 1,492 committed non-empty results of the benchmark, 17 of them non-manifold: every mesh closed, Euler characteristic and components equal to the Brep's (§3). Deviation from the surfaces is measured (at most 0.60 × chord), not certified |
| Tolerance growth over chained Booleans | Met on 40 random chains of up to 6 operations: the certified edge gap stays at or below 7.7e-13 and vertex enclosures below 7.2e-13, with no growth from step 2 to step 6 (§3) |
| Bitwise reproducibility | See §8; one compiler and one machine only |
| Speed and memory relative to a legacy kernel | Not measured |
| Certified G1/G2/G3 bounds | Met for seams between patches in a common parameterisation (§4, certified seam bounds): on the finest fill the bounds are within 6 % (G1) and 31 % (G3) of the sampled maxima; rational patches are covered with looser bounds (G2, G3 within a factor 2.5). Reparameterised joins are not covered |
| Healing success rate | Not built |

## 10. What these results do not show

- Behaviour on real CAD models. Every input solid here has six faces, each one Bezier patch, except the tiled boxes of the scaling test (up to 1,536 flat faces); chained results grow to about 37 trimmed faces.
- Curved chains at scale. §3 has 40 chains with cylinder sectors and hill boxes (214 steps, 209 committed). A crossing curve that touches the line two hill tops share, blurred by rounding (1), two curves tangent at a vertex (1) and three-surface searches that do not finish (3, one of them an arc ending on such a shared line) stop the rest; all are reported, never answered.
- Exact points beyond round numbers. A vertex is known exactly only if its homogeneous coordinates are doubles and it has parameters that are dyadic, or fractions with denominators up to 4096, on the surfaces it must be proven on. Vertices with irrational coordinates (most intersections of curved surfaces) stay enclosures, and their coincidences stay unprovable. Exact membership in a curved surface was exercised by one benchmark step (chain 18) and by the evaluator's unit tests, not by a Boolean test with a closed-form volume.
- An edge lying in a curved surface other than along a line of constant parameter (a circle of a sphere, an edge along a cylinder that is not straight in its parameterisation), and an edge through an edge in space at a point whose coordinates are not round: reported. Line contacts are counted and reported; the fan check, which orders faces around shared edges, does not examine them.
- Flat rational faces in planes that are not axis-aligned. A sector's cap is moved onto a parallelogram of its plane only when that plane is x, y or z = constant (exact points on it are then easy to choose); other flat rational faces in one plane with an unrelated flat face are still reported.
- Speed at scale on curved models. The spatial index makes the overhead grow linearly, and flat faces that meet are intersected without subdivision (§3: 0.25 s for two 1,536-face solids), but that was measured on flat tiles only. A step on a 37-face solid with faces that are not exactly flat (rotated boxes, curved operands) still takes about 0.3 s, almost all of it certified surface/surface and edge/surface work, which has not been tuned.
- Multi-patch surfaces (a face spanning several Bezier patches).
- Resolution of every kind of contact. Resolved: sharing by construction (a result with its own operand, a solid with itself), coincident flat faces decided exactly (boxes in contact, and a sector's rational cap against a box face or against edges of the other solid in its plane, by moving the cap onto a parallelogram of its plane), and coincident curved faces whose patches are exact affine reparameterisations of each other, with the parameters in the same order or exchanged (sectors of one cylinder; two caps in one plane). A vertex of one solid in a plane of the other at the end of one of its edges is decided exactly. Non-manifold results (parts touching along an edge or at a point) are represented and used as operands. Tangential contact of a curved surface with a plane along sides of the curved patch (a cylinder along its line, a hill along its base) is decided exactly. A plane touching a degree-2 patch along a line inside it (a cylinder along a ruling) is decided exactly. Detected and reported, not resolved: tangency at a point, tangency between two curved surfaces (two hill tops crossing where they share a boundary line), other coincidences of curved surfaces, and an edge through an edge in general 3D position.
- Non-manifold contact on curved surfaces. One curved case is tested: a box touching the quarter shell along a line of its outer cylinder, where the cylinder and the box face leave the edge at one angle and are ordered by the exact contact proof. Otherwise the fan check computes the angles around an edge in doubles and refuses any two closer than 1e-9 rad. That is a safeguard, not a certificate like the kernel's other decisions.
- Coincident curved surfaces in general. Recognised: one patch an affine reparameterisation of another with dyadic shifts and power-of-two scales (sectors of one cylinder built from its common exact parameterisation). Not recognised: other relations (a scale of 2/3, a rotated parameterisation of the same cylinder, a sphere built in two ways), which are reported as not solved. No quadric solid with a degenerate patch (a full sphere, a disc cap) has been built.
- Exact coplanarity for faces built from rounded data. A rotated box's faces are not exactly planar, so two rotated boxes placed face to face are not recognised as touching.
- Shared curved surfaces in random chains: the reuse chains use boxes only; the curved shared case is the single curated hill-and-slab row.
- Validity of the input. The kernel does not check that an input solid is free of self-intersection.
- Certified volumes or certified continuity jumps. Both are sampled in doubles.
- Certified mesh deviation. The tessellation is closed by construction, but its distance from the surfaces is measured at sample points, not bounded. Its triangles are near-isotropic in a scaled parameter plane, not adapted to principal curvature directions, so meshes of singly curved faces have more triangles than needed.
- Mesh self-intersection. Two faces that meet at a very small angle could give triangles that cross each other at a coarse chord; this is not checked.
- Rigour of the completeness proof under a compiler that breaks strict IEEE semantics; only MSVC with `/fp:strict` was tested.
- One limit of the surface/surface certificate on its own: where a curve ends on a domain boundary, it does not prove that the curve crosses the boundary instead of touching it and turning back. Inside a Boolean this is now proven, because each such end is matched to a regular edge/face root (§3); a stand-alone call to `solve::intersect` still has the limit.
