# Features on complex geometry: fillets, chamfers, sweeps and lofts

Status: research notes and plan, version 0.19 (2026-11-24). Steps F1, F2a and F3a are implemented and tested, and F3b is well under way: procedural carriers with certified enclosures, rolling-ball spines between freeform faces found without a guess with a completeness proof, and the blend itself as a carrier that can be meshed and intersected; F5 is started with fillets and chamfers applied to swept solids (boxes and rings) (RESULTS §11–12). Everything else is design.

Booleans decide *where* two solids meet. Features decide *what new geometry* to put there, and they are where commercial kernels fail most visibly: a fillet that will not build, a blend that self-intersects, a sweep that folds over itself, a loft that twists. This document collects what the literature and the commercial kernels do, says where akernel's certified approach changes the picture and where it does not, and turns it into an ordered plan that extends the roadmap of [BLUEPRINT.md](BLUEPRINT.md) §5.3.

Citations are from the author's reading and from searches made on 2026-11-23. Entries marked **[checked]** were confirmed against a publisher, library or author page during that search; the rest are cited from memory, with venue and year to verify before publication.

---

## 1. How the features fail

A feature succeeds only if four things hold: its surface exists (the defining equations have a solution everywhere along it), it is regular (no fold, cusp or self-intersection), it meets its neighbours along curves the kernel can represent and trim, and the new topology is consistent. Each common failure breaks one of them.

| Failure | Feature | What breaks | Where it is decided |
|---|---|---|---|
| Radius larger than the concave curvature radius | fillet, offset, shell | the ball does not fit: the spine (offset ∩ offset) has a cusp, the blend folds | local offset validity, $1 - r\kappa_i > 0$ (F2a) |
| Ball reaches a third face | fillet | *overflow*: the contact curve must leave the support face and continue on a neighbour | spine ∩ offset of the third face (F6) |
| Support face ends before the blend does | fillet | *cliff edge*: the contact curve hits a boundary that is not a blendable edge | contact curve ∩ face boundary (F6) |
| Several blended edges meet | fillet, chamfer | *vertex blend*: the end caps of the blends overlap or leave a gap | setback vertex blend, N-sided fill (F7) |
| Mixed convex and concave edges at a vertex | fillet | the blends cross each other | same, plus classification of edge convexity (F7) |
| Radius varies to zero, or the radius law is not smooth | variable fillet | the blend degenerates or loses continuity | radius-law validity (F8) |
| Profile larger than the path's curvature radius | sweep | the swept surface folds (local) | sweep Jacobian, analogue of F2a (F1b) |
| Path comes back near itself | sweep, offset, shell | global self-intersection (the offset passes the medial axis) | reach: bottleneck pairs (F2b) |
| Sections twist or their parameterisations disagree | loft | the loft folds or creases | section correspondence and Jacobian (F1c) |
| Interpolated weights cross zero | rational loft | a pole in the surface | refused at construction (F1, done) |
| Offset of a general NURBS | offset, fillet, shell | not rational, so it cannot be a NURBS exactly | procedural carrier (F3), the central architectural question (§5) |

The point of the table: **every failure is a predicate.** A kernel that decides those predicates with proof, and reports when it cannot, has no silent bad fillets. The predicates themselves are not new; certifying them is.

---

## 2. What the commercial kernels do

### 2.1 Parasolid, ACIS, C3D, OCCT

The commercial kernels converge on one pattern. Their internals are proprietary, so this is what their documentation, APIs and published papers show, not their source.

1. **Procedural surfaces.** A rolling-ball blend, a variable-radius blend, an offset or a sweep is stored as its definition: the support surfaces, the radius law, the spine. ACIS names them in its API (`rb_blend_spl_sur`, `var_blend_spl_sur`, `off_spl_sur`, `sweep_spl_sur`); Parasolid's blend surfaces and offset surfaces are likewise defined by reference to their supports. Points are computed on demand by solving the defining equations.
2. **Approximants for everything else.** Each procedural surface carries a B-spline approximation within a fit tolerance, used for intersection, display and export.
3. **Tolerant modelling.** Where two approximants do not meet, the edge gets a local tolerance (Parasolid's tolerant edges and vertices, ACIS's `tedge`/`tvertex`). This absorbs the gap and lets the operation finish. It is also how tolerances grow through a feature history, the problem that akernel's thesis T3 removes.
4. **A large library of special cases.** Overflow onto neighbouring faces, cliff edges, caps, mitres, setback vertex blends, blend-order resolution on chains of edges, and fallbacks that retry with a different strategy. This catalogue, built over decades, is most of what makes a blend "work".
5. **Local operations to rebuild topology.** A blend is put in by removing the edge's neighbourhood and patching faces back in, the same machinery as direct modelling ("remove face", "replace face").

OCCT's filleting (`BRepFilletAPI_MakeFillet`, package `ChFi3d`) is the one whose code is open: it marches along the spine solving the blend equations (`BRepBlend_Walking`), fits B-splines to the result within a tolerance, and has a limited case library. Its well-known failures on overflow and on vertices where several fillets meet are what a missing case library looks like.

**C3D Toolkit** (C3D Labs, extracted from KOMPAS-3D in 2012, licensed since 2013) lists filleting, chamfering, shelling and sheet-metal among its operations [checked: Wikipedia]. Its fillets follow the same procedural-plus-approximation pattern.

### 2.2 Plasticity

Plasticity (Nick Kallen, first released 2023) is a NURBS modeller for concept artists whose selling point is that fillets and Booleans on freeform shapes just work. It currently runs on Parasolid [checked: Wikipedia]. It started on C3D Toolkit and moved to Parasolid; the public pages consulted on 2026-11-23 do not state the date or the reason, so the commonly given reason, better fillet and Boolean robustness on hard cases, is unverified here.

What Plasticity adds on top of the kernel is instructive, because it is the part that makes the kernel look good to users. As recalled from its manual (not re-checked on 2026-11-23):

- **Fillet shape options exposed directly**: round (circular), chamfer, and curvature-continuous ("G2") profiles; variable radius by dragging points on the edge.
- **Automatic handling of the hard cases** without asking the user: the kernel's overflow, cliff-edge and setback strategies, with a preview that updates while dragging.
- **Fillets that "go around" corners** by chaining tangent edges, so the user selects one edge and gets the whole tangent chain.

The lesson for akernel: robustness that users see comes from the case library and from chaining and live previews above the kernel, as much as from the core numerics. akernel can certify each case; it still has to have the cases.

### 2.3 Where certification changes the picture

| Commercial approach | akernel approach |
|---|---|
| Approximant used in decisions, gap absorbed by a tolerance | Decisions read the definition through certified enclosures; no tolerance is stored (theses T2, T3) |
| "Blend too large" found when marching fails | Proven before building: offset validity on the whole support (F2a) |
| Overflow detected by intersecting approximants | Overflow is a certified intersection of the spine with the third face's offset (F6) |
| Self-intersection of offsets checked by sampling or not at all | Reach decided from local curvature bounds and certified bottleneck pairs (F2b) |
| Failure modes discovered by users | Each failure is a predicate that answers Valid, Invalid with a witness, or Undecided |

What certification does **not** remove: the case library (§1), the vertex-blend constructions, and the cost of procedural evaluation.

---

## 3. The literature, by topic

### 3.1 Rolling-ball blends and their exact forms

- Rossignac, Requicha. Constant-radius blending in solid modelling. *Computers in Mechanical Engineering*, 1984. Blends as offsets of offsets; the spine on the intersection of two offset surfaces.
- Choi, Ju. Constant-radius blending in surface modelling. *Computer-Aided Design* 21(4), 1989. The rolling-ball formulation used by most kernels.
- Vida, Martin, Várady. A survey of blending methods that use parametric surfaces. *Computer-Aided Design* 26(5), 1994. The standard taxonomy: rolling ball, cross-sectional, variable radius, chordal width; and the failure cases.
- Pottmann, Peternell. Computing rational parametrizations of canal surfaces. *Journal of Symbolic Computation* 23, 1997; Peternell, Pottmann. A Laguerre geometric approach to rational offsets. *CAGD* 15, 1998. A rolling-ball blend is a canal surface (envelope of spheres); it is rational whenever its spine is rational, and for natural quadrics the spine often is.
- Rational fixed-radius rolling-ball blends between natural quadrics, and Krasauskas et al., Rational rolling ball blending of natural quadrics, *Mathematical Modelling and Analysis*, 2010s [checked: portal.mardi4nfdi.de; journals.vilniustech.lt]. Exact rational blends for plane, cylinder, cone and sphere pairs in common configurations, often via Dupin cyclides.
- Allen, Dutta. Cyclides in pure blending I and II. *CAGD* 14, 1997. Cyclide blends between quadrics; exact and rational.
- Kós, Martin, Várady. Methods to recover constant radius rolling ball blends in reverse engineering. *CAGD* 17, 2000. Useful because it states the blend geometry in a form the kernel can test (spine, radius, contact curves).
- Rolling ball blends and self-intersection. *SPIE* 1830, 1992 [checked: SPIE digital library]. Self-intersection of blends from the radius and support curvature.

Consequence for akernel: blends between planes and natural quadrics (most mechanical fillets) are **exact rational surfaces**, so they fit the current carrier (rational Bézier patches) with no procedural layer. Plane–plane gives a cylinder, plane–cylinder and cylinder–cylinder with parallel or intersecting axes give tori and cyclides, plane–cone and plane–sphere give tori. That is the first fillet step (F3a).

### 3.2 Variable-radius and curvature-continuous blends

- Pegna, Wolter. Geometrical criteria to guarantee curvature continuity of blend surfaces. *ASME J. Mech. Design*, 1992. The linkage-curve theorem already used in BLUEPRINT §1.7.
- Hartmann. G^n-continuous connections between normal ringed surfaces. *CAGD* 18, 2001. Gn blends built on the same offset/spine structure as rolling balls, with a cross-section that meets the supports to order n.
- Bloor, Wilson. Generating blend surfaces using partial differential equations. *Computer-Aided Design* 21(3), 1989. PDE blends; the ancestor of the variational approach in BLUEPRINT §4.
- BLUEPRINT §4.2 already specifies the G3 fillet: a degree-7 cross-section with 8 Hermite conditions per rail, coupled along the spine by the variational solver. SURVEY.md found no 2015–2026 paper on G3 or certified curvature-continuous rolling-ball fillets; that remains open and distinctive.

### 3.3 Vertex blends

- Várady, Rockwood. Vertex blending: problems and solutions. In *Mathematical Methods for Curves and Surfaces II*, 1998 [checked: SZTAKI eprints, Purdue CGVLab]. The catalogue of vertex configurations (all convex, all concave, mixed) and why the naive three-patch corner fails.
- Várady, Rockwood. Geometric construction for setback vertex blending. *Computer-Aided Design* 29(6), 1997 [checked: SZTAKI eprints]. *Setbacks*: stop each edge blend short of the vertex and fill the hole with an N-sided patch whose boundary has 2N sides (N blend ends, N setback curves). This is what Parasolid and ACIS expose as setback vertex blends.
- Vertex blending based on the setback split [checked: SZTAKI eprints]; Salvi, Várady, multi-sided patches (SURVEY.md). The N-sided fill is exactly what akernel's var/fill.hpp and SURVEY's "Route A" build.

Consequence: a vertex blend is an N-sided fill with prescribed boundary jets, which is the existing surfacing track. Fillets and the surfacing track meet here (F7).

### 3.4 Offsets and their self-intersection

- Farouki. The approximation of non-degenerate offset surfaces. *CAGD* 3, 1986.
- Maekawa. An overview of offset curves and surfaces. *Computer-Aided Design* 31, 1999. Local (curvature) and global (medial axis) self-intersection, and why the offset of a NURBS is not rational.
- Pottmann. Rational curves and surfaces with rational offsets. *CAGD* 12, 1995. Pythagorean-normal (PN) surfaces; the exact cases.
- Seong, Elber, Kim. Trimming local and global self-intersections in offset curves/surfaces using distance maps. *Computer-Aided Design* 38, 2006.
- Federer. Curvature measures. *Trans. AMS* 93, 1959. The *reach* of a smooth closed surface is the smaller of $1/\max|\kappa|$ and half the length of its shortest *bottleneck* (a pair $p \ne q$ with $p - q$ parallel to both normals). Below the reach every offset is a regular, embedded surface.

Consequence: global offset validity splits into a local test (F2a, done) and a search for bottleneck pairs, which is a square polynomial system in four parameters, the same kind the kernel already certifies with Krawczyk for surface/surface intersection (F2b). The same test checks an input solid for self-intersection, which RESULTS §10 lists as missing.

### 3.5 Sweeps

- Wang, Jüttler, Zheng, Liu. Computation of rotation minimizing frames. *ACM Transactions on Graphics* 27(1), 2008 [checked: HKU, Microsoft Research]. The double-reflection method: fourth-order accurate, cheap, the standard way to sweep a profile along a space curve without twist.
- Farouki, Giannelli, Manni, Sestini. Quintic space curves with rational rotation-minimizing frames. *CAGD* 26, 2009; Farouki, *Pythagorean-Hodograph Curves*, 2008. Curves whose RMF is rational (RRMF): sweeps along them are exact rational surfaces.
- Jüttler, Mäurer. Rational approximation of rotation minimizing frames using Pythagorean-hodograph cubics. *J. Geometry and Graphics* 3, 1999.
- Abdel-Malek, Yeh. Geometric representation of the swept volume using Jacobian rank-deficiency conditions. *Computer-Aided Design* 29(6), 1997; Blackmore, Leu, Wang. The sweep-envelope differential equation algorithm. *Computer-Aided Design* 29(9), 1997. Self-intersection of sweeps from the rank of the sweep map's Jacobian.

Consequence: extrusion, revolution and translational sweeps are exact (F1, done). A sweep with a moving frame is exact only on RRMF paths; elsewhere it is an approximation whose deviation from the true RMF sweep must be bounded, or the sweep must be a procedural carrier (F1b).

### 3.6 Lofts

- Woodward. Skinning techniques for interactive B-spline surface interpolation. *Computer-Aided Design* 20(8), 1988.
- Piegl, Tiller. *The NURBS Book*, 2nd ed., 1997, §10.3. Compatibility (common degree and knots), section parameters by averaged chord length, knots by averaging. Implemented here (F1).
- Hohmeyer, Barsky. Skinning rational B-spline curves to construct an interpolatory surface. *CVGIP: Graphical Models* 53, 1991. Rational sections and the weight-crossing problem; akernel refuses a loft whose interpolated weights are not positive instead of returning a surface with a pole.

The robust-loft problems are not in the interpolation but in **correspondence** (which point of section k matches which point of section k+1: seam alignment, twist) and in **fairness between sections** (unwanted oscillation with many sections). The first is a user-facing choice plus a Jacobian check for folds (F1c); the second is the variational solver of BLUEPRINT §4.4 with sections as constraints, the same as crown lofting (§4.3).

### 3.7 Certified computation for all of the above

Every predicate in §1 reduces to the signs of polynomials in the parameters of a few patches, sometimes after a square root (the unit normal). The kernel's existing tools cover them: interval Bernstein enclosures and subdivision for "positive everywhere" (F2a uses them with a mean-value form to beat dependency), and Krawczyk certification for "exactly one root here" (spine curves, bottleneck pairs, overflow points). The literature on certified offsets is thin; the closest are the self-intersection papers in SURVEY.md (Botbol–Busé–Chardin on multiple points; the 2025 TOG paper on fast NURBS self-intersection).

---

## 4. The xNURBS energy

xNURBS markets "an energy-minimization method to generate the smoothest NURBS surface among all possible solutions" for given constraints, solving "in milliseconds", with no published functional [checked: McNeel blog, 2022]. SURVEY.md already notes that it documents G2 only and quotes accuracy in thousandths of a millimetre and tenths of a degree, which implies constraints met approximately.

Its lineage is plausibly the physics-based NURBS line of Terzopoulos and Qin (Dynamic NURBS with geometric constraints for interactive sculpting, *ACM TOG* 13(2), 1994 [checked: author's page]; D-NURBS, *IEEE TVCG*, 1996). That is an inference from the field, not a claim of the vendor. D-NURBS uses the *thin plate under tension* energy

$$
E = \int\!\!\int \alpha_{11}\lVert S_u\rVert^2 + \alpha_{22}\lVert S_v\rVert^2 + \beta_{11}\lVert S_{uu}\rVert^2 + 2\beta_{12}\lVert S_{uv}\rVert^2 + \beta_{22}\lVert S_{vv}\rVert^2 \,du\,dv ,
$$

with spatially varying tension $\alpha$ and rigidity $\beta$, quadratic in the control points ($\tfrac12 p^\top K p$), Lagrangian dynamics $M\ddot p + D\dot p + Kp = f$ for interaction, linear constraints eliminated by a null-space basis $p = Gq + q_0$, and nonlinear constraints by Lagrange multipliers with Baumgarte stabilisation.

Against akernel:

| | Thin plate under tension (D-NURBS, plausibly xNURBS) | akernel (BLUEPRINT §1.7, §4.4) |
|---|---|---|
| Energy order | 2 (bending) plus 1 (tension) | matched to continuity: order 3 for G2, order 4 for G3 |
| G2 boundary data | over-determined in the continuum; boundary layers under refinement | well posed; RESULTS §5 measures the thin-plate solution's E3 growing 21.9 → 98.2 under refinement while the matched energy stays at 11.18 |
| Constraints | null-space elimination (same as akernel) | null-space elimination, Lagrange multipliers reported |
| Parameterisation dependence | yes (parametric energy) | the linearised energies too; geometric $E_2, E_3, E_4$ for the final polish |
| Continuity guarantee | measured | certified bounds on the jumps (RESULTS §4) |

What to take from it: the **tension term**. A small $\alpha$ term keeps the parameterisation even (it penalises stretched control nets), which the higher-order energies alone do not. It should be added to akernel's energy as a regulariser with a small weight, and the time-stepping form is a good model for interactive dragging (the prefactored KKT of SURVEY's surfacing track). The energy-order rule stays: for a G3 fillet cross-section the energy must be of order 4.

---

## 5. The architectural question: procedural carriers

The kernel's carriers are rational Bézier patches. That covers extrusions, revolutions, translational sweeps, lofts, and rolling-ball blends between planes and natural quadrics, all exactly. It does not cover the offset of a freeform face, a rolling-ball blend on a freeform face, a variable-radius blend, or a sweep with a moving frame on a general path. Those surfaces are not rational.

Two designs are possible:

- **A. Procedural carriers as the definition.** A carrier type defined by an equation (offset: $S + d\,n$; blend: the canal surface of a spine; sweep: profile in a frame), evaluated with certified enclosures (intervals through the square root of the normal). Intersection, classification and tessellation work through a carrier interface instead of `PolyPatch`. This is thesis T2 ("one definition, many certified evaluators") applied to features, and it is what keeps features free of stored tolerances. Cost: the certified solvers (solve/ssi.hpp, csi.hpp, contact_point.hpp) must be generalised from Bézier patches to the interface, and enclosures of procedural surfaces are looser than Bernstein hulls, so they are slower.
- **B. Rational approximants with a certified deviation.** Fit a rational patch to the procedural surface and attach a certified bound on the distance between them. Decisions read the approximant. Cost: the fitted surface is the definition from then on, so the blend's contact with its supports is no longer exact; the edge between a support and its fillet must carry the deviation, which is a recorded tolerance (thesis T3 allows this, budgeted and local, but it is the commercial pattern).

Decision (2026-11-23, project owner): **A, procedural carriers as the definition**, introduced gradually. Step F3a needs neither (its blends are rational). F2b and F3b build on A.

First layer (2026-11-24). A carrier provides three evaluators: point and first partials in doubles (for Newton steps and seeds), a rigorous enclosure of the point at exact parameters, and rigorous enclosures of the point and first partials over a parameter box. That is exactly what the certified Newton–Krawczyk machinery of `solve/ssi.hpp` reads, so its routines are now templates over the system and work unchanged on any pair of carriers. The offset carrier's enclosures come from the base patch's homogeneous partials up to order 3 through the normal's derivative, intersected with a mean-value form so they converge as the box shrinks. The rolling-ball blend is then the certified spine (offset ∩ offset) plus the circular section at each spine point, with nothing fitted. Second layer (2026-11-24): the completeness verifier now runs on carriers (exclusion by disjoint enclosures and by separation along the normals), so spines are found without a guess; and the blend is a carrier whose enclosures come from the Krawczyk box of the spine, the implicit-function derivative of the spine, and the section arc in forward-mode interval arithmetic, tightened by mean-value forms.

---

### 5.1 Procedural faces in the B-rep (decided 2026-11-24)

Decision (2026-11-24, project owner): a solid may have faces whose surface is a procedural carrier (`geom::Carrier`: offsets, rolling-ball blends, later sweeps). This is how a freeform fillet becomes part of a solid. The plan is to make checks, volume and meshing work on such solids first, and Booleans on those faces after.

Data model, added to `op::Brep` without changing the rational path:

- `Brep::procs`: the procedural carriers. A face with `proc >= 0` is that carrier over its whole parameter square. A blend is parametrised by sigma along the spine and theta across the section, so its square sides are the two contact curves and the two end sections.
- Procedural edges. An edge that is a side of a procedural face is defined once, as that side. Each other face using the edge reads it through the carrier, in its own parameters:
  - on a support it is a contact curve, whose parameters are the spine's own coordinates, with the same certified tubes;
  - on a planar end face it is the end section, mapped by the plane's exact inverse.
  
  So the edge has one definition, as everywhere in the kernel, and no fitted trimming curve.
- What reads which part:
  - the audit reads only incidence, so it is unchanged;
  - volume integrates S · n over procedural faces with the carrier's evaluator, and over trimmed rational faces by Green's theorem along the procedural edges' parameter curves;
  - meshing samples each procedural edge once and shares the samples between its two faces, so the mesh stays watertight.
- Booleans refuse procedural faces until their intersection code reads carriers. The certified carrier intersection already exists (`intersect_carriers`).

First target: a block whose top is freeform, with the edge between the top and a side rounded. The top is chosen with zero slope across the two end planes, so each end section of the blend lies in its end plane, and the end faces are trimmed by a planar arc rather than by a cap.

State (2026-11-24): done for checks, volume and meshing (`op/proc_fillet.hpp`, `round_block_edge`; RESULTS §12, "A freeform fillet inside a solid"). The rounded freeform block audits closed, its volume matches the closed form on a flat top and an independent slice integration on a freeform top to 5e-12, and it meshes watertight. Booleans still refuse procedural faces. Trimmed procedural faces are done too: a blend cut by an end face that is not flat across it is trimmed by a traced intersection edge (`op/proc_edge.hpp`), so the zero-slope restriction at the ends is gone (sloped ends, RESULTS §12). Any convex side edge of a Brep is done too (`round_edge`, RESULTS §12, "Fillets on general solids"). It works on any patch faces and end faces, with flipped faces, and on solids that are already rounded. Next: chains of edges, which meet at rounded vertices.

Fillet corners (2026-11-24). Two fillets that meet at a corner of a shared face G touch each other at one point A. Each blend is tangent to G along its contact curve, and the two contact curves cross at A, so the blends' tangent planes agree there. The intersection curve of the two blends therefore has a singular point at A. Two branches cross there, and only one of them lies inside both arcs. `intersect_carriers` needs transversal crossings, so it traces this curve but cannot prove it complete (checked on a box: fillets of radius 0.25 and 0.15 at a corner). Two ways were put to the project owner: a mitre whose touching point is proven (recommended), or a corner patch that keeps the fillets apart. Plan for the mitre:

1. A, proven: a 2-D Krawczyk on the contact curves in G's parameters, w1(sigma1) = w2(sigma2). This is a transversal crossing.
2. The local picture, proven by a polar blow-up in the second blend's parameters around (sigma2*, 0). Write h = (distance to the first spine) − r1 on the second blend. Then h and its gradient vanish exactly at A, so h / rho^2 is enclosed by its Hessian (Taylor's integral form, with the Hessian of the distance to a curve). On a half-disc rho ≤ rho_max, exclusion over angle intervals and Krawczyk near the two roots prove exactly two branches through A. Each branch is a graph over rho. The sign of (P − c1) · T1 / rho tells which branch lies inside the first blend's arc.
3. Away from A, the ordinary certified tracer on the second blend restricted to theta2 ≥ rho_max (a reparameterised carrier), with exclusion of the rest of the strip near G.
4. Topology in `round_edge`: an end face that is a procedural face, edges at the end vertex that are procedural (split with `ProcEdge::part`), and the new vertex where the second blend's far contact meets the first blend's end curve.

Equal radii are a further special case: there the second blend's far contact runs through the first blend's end vertex. State (2026-11-25): steps 1 to 4 are done for a second radius smaller than the first, on a box and on a freeform top (RESULTS §12, "Fillet corners"). Then (2026-11-25) equal radii, where the mitre ends at the first fillet's end vertex, and the larger radius second, where it ends on the first fillet's far contact and the new fillet goes on to meet the top (RESULTS §12). Both orders and both sides of the shared face are tested. Then a clearance check: the ball is proven clear of every other face of the solid (RESULTS §12, "Clearance"). A third edge at the same corner already works when its radius is small enough that its far contact meets what is left of the first fillet's end curve. Then three fillets at one corner (RESULTS §12, "three fillets at one corner"): the third fillet's mitre crosses both earlier fillets through the point where all three meet, and ends at a second touching point. Box volumes agree with inclusion and exclusion to 2e-12. Open: three equal radii (a vertex of four faces), and orders where the blow-up at the touching point stalls (two close branches). A smooth vertex blend (a setback patch instead of mitres) is the other option for a corner of three fillets.

Decision (2026-11-24, project owner): Booleans on procedural faces start with cuts by solids whose faces are flat. Procedural faces against curved faces come later. The order: trimmed procedural faces (done), a plane cut of any solid with procedural faces (done, `op/plane_cut.hpp`), then Booleans whose other operand has flat faces (done: union, intersection and difference, `op/poly_boolean.hpp`; RESULTS §12). Since 2026-11-24 the same code also takes an other operand with curved rational faces (a cylinder wall through the blend is tested). Then (2026-11-24) both operands with procedural faces: two rounded blocks with crossing fillets (`op/brep_transform.hpp` places solids exactly and re-indexes procedural edges). Next: inner loops in the operands, speed of blend enclosures.

## 6. Plan

The steps extend BLUEPRINT §5.3 as a seventh track, "Features". Cross-platform work (Linux, macOS, CI) is deferred by decision of 2026-11-23.

| Step | Goal | Depends on | State |
|---|---|---|---|
| F1 | Exact generators: NURBS curves (Bézier extraction, degree elevation, compatibility), exact rational circle, extrusion, revolution by quarter turns (cylinder, cone, sphere zone, torus), translational sweep, ruled loft (also the chamfer carrier), skinned loft with refusal of poles | — | **Done** (`geom/nurbs_curve.hpp`, `geom/sweep.hpp`; RESULTS §11) |
| F1b | Sweeps with a moving frame: exact on RRMF paths; double-reflection RMF elsewhere, with certified deviation or as a procedural carrier; sweep Jacobian check (local self-intersection) | F2a, §5 | Open |
| F1c | Loft correspondence (seam alignment, twist control) and fold check; lofts with fairness as a variational problem with sections as constraints | F1, var | Open |
| F2a | Certified local offset validity ($1 - d\kappa_i > 0$ on a whole patch: Valid, Invalid with witness, or Undecided) and certified bounds on $\max\lvert\kappa\rvert$, the largest safe fillet radius | — | **Done** (`geom/offset_validity.hpp`; RESULTS §11) |
| F2b | Global offset validity: reach from F2a and certified bottleneck pairs; the same test for input-solid self-intersection | F2a, §5 for freeform | Open |
| F3a | Constant-radius rolling-ball fillets between planes and natural quadrics, as exact rational surfaces | F1, F2a | **Done** for plane/plane (cylinder sectors) and for planes, cylinders, cones and spheres about one axis (torus zones), `geom/blend.hpp`. Open: quadrics in general position (cyclides), chamfers between quadrics |
| F3b | Rolling-ball fillets on freeform faces: spine as the certified intersection of two offset surfaces; contact curves as varieties in the supports' parameters | F2b, procedural carriers (§5, decided) | **Mostly done**: carrier interface with patch and offset carriers (`geom/carrier.hpp`); the Newton, Krawczyk and tracing routines of `solve/ssi.hpp` made generic over the system; certified spine tracing between freeform faces and the blend section at each spine point (`solve/carrier_ssi.hpp`). Then (2026-11-24): `intersect_carriers` finds carrier intersections without a guess and proves completeness, `find_spines` applies it to two offsets, and `BlendCarrier` makes the blend of a spine a carrier itself (certified enclosures, intersected with a plane in 0.65 s). Trimming and topology: first case done (F5, `op/proc_fillet.hpp`) |
| F4 | Chamfers on freeform faces: distance–distance and distance–angle contact curves, ruled carrier | F3b | **Started**: straight chamfers (equal and unequal setbacks) on edges of swept solids, the same construction as F5's fillets. Open: freeform supports, distance–angle |
| F5 | Blend topology: replace the edge by the blend face, trim the supports along contact curves, chain tangent edges, end caps | F3a, topo | **Started**: `op/fillet.hpp` fillets and chamfers edges of swept solids. The section's corner is cut, then extruded (box edges) or revolved by exact quarter turns (circular edges of rings, closed chains of four edges, torus and cone blends). Concave edges, where the blend adds material, use a seven-cell layout of the step section: a step block, a boss on a plate and a counterbore. Vertex blends: `rounded_box` rounds all twelve edges, the corners exact sphere octants of three regular patches built by inverse stereographic projection, split symmetrically along great circles through the corner's centre (three near-congruent pieces). The supports are trimmed at the contact curves, and end faces are split into five regular cells. Results are closed, exact and usable by the Booleans. Then (2026-11-24) the first freeform fillet in a solid: `op/proc_fillet.hpp` rounds the edge between a freeform top and a side, the blend a procedural face (§5.1), closed, volume checked two ways, meshed watertight. Then plane cuts, and Booleans against solids with flat or curved rational faces (`op/plane_cut.hpp`, `op/poly_boolean.hpp`). Then faces with holes as Boolean operands, so results chain. Open: other freeform support pairs and end conditions (caps, chains), non-swept solids with trimmed faces, vertex blends of unequal radii (setback blends) |
| F6 | The case library: overflow onto neighbouring faces, cliff edges, mitres, each detected by a certified predicate and reported when undecided | F5 | Open |
| F7 | Vertex blends: convexity classification, setbacks, N-sided fill of the vertex hole with the surfacing track's fill | F5, surfacing track | Open |
| F8 | Variable radius, chordal width, conic and G2/G3 cross-sections (BLUEPRINT §4.2) with the tension-regularised, order-matched energy of §4 | F3b, var | Open |
| F9 | Offsets and shells of solids | F2b, F3b | Open |

The order puts first what is exact with the current carriers and what every later step needs (the validity predicates). F3a gives the mechanical fillets most models use, exactly, before the procedural layer exists. The case library (F6) and vertex blends (F7) are the long tail where commercial kernels have spent decades; they come after the geometry is right.
