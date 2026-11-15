# akernel — Theoretical Blueprint for a Unified B-Rep / Freeform Geometric Kernel

Status: design document, version 0.16 (2026-11-15). The certified core is implemented and tested: numerics, carriers, certified intersection, continuity measurement, variational fill and the topology audit. Boolean operations (§3.3) are implemented for solids with trimmed faces; results can be operands, so operations chain; operands that share geometry by construction (step 2 of §3.3, by identity rather than full provenance) are handled, and so is coincidence of flat faces decided by exact predicates. Coincident curved surfaces that are not identical are detected and reported, not resolved. Tangential contact of a curved patch with a plane along sides of the patch is decided by exact algebra (step 3 of §3.3, "singular" case, for planes); other tangency is reported. Non-manifold results (parts touching along an edge or at a single vertex) are represented, checked by the geometric order of the faces around each shared edge, and usable as operands. Quadric patches with exact data (cylinder sectors from the exact rational circle) are built and tested on their quadrics exactly, and curved faces of one surface that are not the same patch are recognised when one is an exact affine reparameterisation of the other (step 2 of §3.3 for that case); other coincidences of curved surfaces remain open (§5.5, item 4). Random chains with curved operands (cylinder sectors, hill boxes, rotated boxes) are benchmarked: every committed result agrees with Monte Carlo and tessellates closed, and tangential and other degenerate contact is reported. The crack-free tessellation of §5.1 (`ak::viz`) is implemented: each edge is sampled once and shared by both faces, and every mesh is checked closed and topologically equal to its Brep; its deviation from the surfaces is measured, not yet bounded. The algebraic-identity target of §5.4 is tested: identities such as A − (A − B) = A ∩ B and associativity give the same solid with identical topology, decided by a structural comparison after removing seams and pass-through vertices. The broad phase (step 1 of §3.3) is a bounding-volume hierarchy over certified boxes, and every all-against-all loop of the Boolean goes through it. Two flat faces in different planes are intersected along their one line with a direct completeness proof instead of subdivision; a Boolean of two 1,536-face solids takes about a quarter of a second. [README.md](../README.md) lists what is and is not built, and [RESULTS.md](RESULTS.md) holds the measurements. Everything else here is still design, and the figures in §5.4 remain targets, not results.

Math is written in LaTeX (`$…$`, `$$…$$`); view it in a renderer that supports math (GitHub, VS Code preview).

---

## 0. Position and scope

### 0.1 Four theses

The whole architecture follows from four commitments.

| # | Thesis | What it replaces |
|---|--------|------------------|
| T1 | **Identity by construction, not by proximity.** Two entities are the same because the kernel can show they have the same definition, never because they are closer than some ε. | Global tolerance comparisons (`dist < resabs`). |
| T2 | **One definition, many certified evaluators.** Every geometric entity has exactly one defining expression. Parametric, implicit and discrete forms are evaluators of it, and each returns an enclosure (a value plus a guaranteed error bound). | Three independent representations of an edge (3D curve + two p-curves) that disagree. |
| T3 | **Tolerance is a local, recorded, budgeted quantity.** A tolerance exists only where a decision could not be proven, is stored on that entity, and is never inherited by entities derived from it. | Tolerances that grow monotonically through a feature history. |
| T4 | **Smoothness is a constraint in an optimisation problem.** G1/G2/G3 are equality constraints on jets; shape is the minimiser of an energy. | Manual control-point manipulation. |

### 0.2 Corrections to the brief

Some goals in the brief cannot be met as literally stated. The blueprint is built on the corrected versions.

| Brief says | What is actually achievable |
|---|---|
| "Permanently eliminate gaps at trimming boundaries" | An edge shared by two faces is a single object (an algebraic variety, §1.4), so there is no gap *by definition*. Evaluating a point on it is numerical, with a certified bound that can be driven below any requested ε. |
| "Deterministic Booleans on freeform surfaces" | Achievable for transversal configurations. Tangential and coincident configurations are decided by provenance (T1) or exact algebra where possible; the remainder return an explicit `Indeterminate` instead of a silent guess (§3.3). No kernel can decide every degenerate real-algebraic configuration in bounded time. |
| "Without resorting to tessellation" | Legacy B-Rep kernels do not tessellate for Booleans either. Their destructive step is refitting intersection curves to splines and creating tolerant edges. That is the step removed here. |
| "Exact trigonometric rational splines" | Conics and surfaces of revolution are exact as rational quadratics. Uniform angular speed (arc-length parameterisation) is not rational for any curve except a line, so "exact" refers to the point set, not the parameterisation (§1.6). |
| "Mirroring and surpassing XNurbs" | XNurbs is publicly described as an energy-minimisation NURBS solver; its algorithm is unpublished. This document specifies measurable targets (§5.4) instead of claiming superiority. |
| "Micchelli et al. error bounds" | Interpreted as the optimal-recovery framework (Micchelli–Rivlin) plus concrete a-posteriori bounds (Filip et al.; Nairn–Peters–Lutterkort). No single named theorem is assumed (§2.2). |
| "G3 across multi-patch boundaries" | Exact G3 against arbitrary given neighbours generally forces degree elevation. The kernel offers exact G3 where degrees allow and otherwise ε-G3 with certified sup-norm bounds on the jumps (§1.7). |

---

## 1. Multi-domain mathematical correlation and unification

### 1.1 One object, four views

A face carrier is a smooth 2-manifold $M \subset \mathbb{R}^3$. The four mathematical domains are four descriptions of the same $M$, each used for the job it does best.

| View | Description of $M$ | Used for |
|---|---|---|
| Differential geometry | Immersion $S:\Omega\subset\mathbb{R}^2\to\mathbb{R}^3$ with fundamental forms $(g,b)$ | Continuity classes, fairness energies, offsets |
| Algebraic geometry | Zero set $V(F)$ of a polynomial, or image of a rational map | Identity, intersection, exact predicates |
| Algebraic topology | A cell in a chain complex $C_2\to C_1\to C_0$ | Validity, watertightness, Euler invariants |
| Field theory | Level set $\{\phi = 0\}$ of a scalar field | Point classification, offsets, healing energies |

The bridges between them are explicit maps, given in §1.2–§1.6.

### 1.2 Differential-geometric core

With $S_\alpha = \partial_\alpha S$ and unit normal $n = (S_1\times S_2)/\lVert S_1\times S_2\rVert$:

$$
g_{\alpha\beta} = S_\alpha\cdot S_\beta,\qquad
b_{\alpha\beta} = n\cdot S_{\alpha\beta},\qquad
W^\alpha{}_\beta = g^{\alpha\gamma}b_{\gamma\beta},\qquad
\partial_\alpha n = -W^\beta{}_\alpha S_\beta .
$$

Principal curvatures $\kappa_1,\kappa_2$ are the eigenvalues of the shape operator $W$; $H=\tfrac12\operatorname{tr}W$, $K=\det W$.

By Bonnet's theorem, a pair $(g,b)$ satisfying the Gauss–Codazzi equations determines the surface up to rigid motion. So $(g,b)$ are the complete invariants, and geometric continuity can be stated without reference to any parameterisation:

| Class | Invariant statement on the glued surface $M = M_1\cup_\Gamma M_2$ |
|---|---|
| G1 | The tangent plane (equivalently $n$) is continuous across $\Gamma$. |
| G2 | G1 and the second fundamental form $b$ is continuous. |
| G3 | G2 and the covariant derivative $\nabla b$ is continuous. |

The Codazzi equation says $C_{\alpha\beta\gamma} := \nabla_\gamma b_{\alpha\beta}$ is totally symmetric. This gives the key reduction used in §1.7.

### 1.3 Bridge to field theory

Let $\phi$ be the signed distance to $M$, positive on the $+n$ side. In the tubular neighbourhood where the closest-point projection is unique:

$$
\lVert\nabla\phi\rVert = 1,\qquad \nabla\phi\big|_M = n,\qquad
\operatorname{Hess}\phi\big|_{TM} = -W,\qquad \Delta\phi\big|_M = -2H .
$$

The offset $S_d = S + d\,n$ is the level set $\{\phi = d\}$. In principal directions $\partial_i S_d = (1-d\kappa_i)\,\partial_i S$, so

$$
\kappa_i(d) = \frac{\kappa_i}{1-d\,\kappa_i}.
$$

Two predicates fall out directly:

- **Local offset validity:** $1 - d\,\kappa_i > 0$ for $i=1,2$ everywhere on the face.
- **Global offset validity:** $\lvert d\rvert < \operatorname{reach}(M)$, the distance from $M$ to its medial axis.

A rolling-ball blend of radius $r$ has its spine on the intersection of two offsets at distance $r$. So blend validity is offset validity plus transversality of the two offsets. Offsets, shells and blends therefore share one failure predicate (§2.2).

### 1.4 Bridge to algebraic geometry: the edge as a variety

A rational patch of bidegree $(m,n)$ is the projection of a polynomial map into homogeneous coordinates, $S = \pi(S^w)$, $S^w:\Omega\to\mathbb{R}^4$. Generically its implicit equation has degree $2mn$, which is 18 for a bicubic. The intersection of two bicubic patches is therefore an algebraic space curve of degree up to $18\times18 = 324$. Such a curve is generally not rational, so **no NURBS curve can represent it exactly**. This is the root of the trimming-gap problem: legacy kernels store a fitted 3D spline and two fitted p-curves, three objects that agree only to a tolerance.

akernel defines the edge as the variety itself:

$$
E = \{(u,v,s,t)\in\Omega_1\times\Omega_2 \;:\; F(u,v,s,t) := S_1(u,v) - S_2(s,t) = 0\}\subset\mathbb{R}^4 .
$$

The Jacobian is $J = [\,S_{1,u}\;\; S_{1,v}\;\; -S_{2,s}\;\; -S_{2,t}\,]\in\mathbb{R}^{3\times4}$.

- $E$ is a regular curve at a point exactly when $\operatorname{rank}J = 3$, which holds exactly when $n_1\times n_2\neq0$ (transversal intersection).
- The tangent of $E$ is the null vector of $J$.
- The condition number of a point on $E$ scales as $1/\sin\angle(n_1,n_2)$.

Both faces reference this single object. Its projection to $(u,v)$ is the trim curve in face 1; its projection to $(s,t)$ is the trim curve in face 2. There is nothing to disagree.

**Certified representation.** $E$ is stored as its defining pair $(S_1,S_2)$ plus a *tube*: a chain of boxes $X_k\subset\mathbb{R}^4$, each certified to contain exactly one arc of $E$. Certification uses the Krawczyk operator on the square system $G=(F,\ell)$, where $\ell$ is a slicing hyperplane transversal to the tangent:

$$
K(X) = \hat x - C\,G(\hat x) + \bigl(I - C\,G'(X)\bigr)(X-\hat x),\qquad C\approx G'(\hat x)^{-1}.
$$

If $K(X)\subset\operatorname{int}X$, then $G$ has exactly one root in $X$. Any point of $E$ can then be refined by Newton iteration to any requested accuracy. Spline approximants of $E$ are kept as caches for display and export only; no topological decision reads them.

Legacy kernels already have procedural intersection curves (for example ACIS `int_cur`). What they lack is the certificate and the rule that decisions use only the definition.

**Implicit side.** For classification the kernel needs the implicit form of a carrier:

- Planes, natural quadrics, tori: closed-form polynomials, evaluated exactly.
- Rational patches: matrix representations (Busé), where a matrix $\mathbb{M}(x)$ drops rank exactly on the surface. This avoids expanding a degree-18 polynomial.
- Everything else: the oriented closest-point field $\phi(x) = (x-\pi_M(x))\cdot n(\pi_M(x))$ with interval evaluation.

The trim region of face 1 against surface 2 is then the sign set of the composite field $g(u,v) = \phi_2(S_1(u,v))$ on $\Omega_1$. This is the coupling the brief asks for: a parametric patch trimmed by the level set of an algebraic or distance field.

### 1.5 Bridge to algebraic topology

A B-Rep is a cell complex with oriented incidence matrices $\partial_2\in\{-1,0,1\}^{E\times F}$ and $\partial_1\in\{-1,0,1\}^{V\times E}$. Validity conditions:

- **Chain condition:** $\partial_1\partial_2 = 0$.
- **Watertightness of a solid shell:** $\partial_2\,\mathbf{1}_F = 0$ (every edge is used exactly twice, with opposite orientation).
- **Euler–Poincaré:** $V - E + F - (L-F) - 2(S-G) = 0$, with $L$ loops, $S$ shells, $G$ total genus.
- **Betti numbers** of the boundary of a solid: $\beta_0=\beta_2=S$, $\beta_1 = 2G$.

These are integer computations and are checked after every operation. Under T2 they are also *sufficient* for watertightness: because both uses of an edge reference one variety, combinatorial closure implies geometric closure. In a legacy kernel the same check passes while the geometry still has a gap of up to the edge tolerance.

For imported data, persistent homology supplies the scale analysis: thicken the model by ε, track $\beta_0,\beta_1$ as ε grows, and read intended topology from the features that persist over a wide ε range (§3.4).

### 1.6 The trigonometric layer

Trigonometry enters in three places, and each has an exact rational form.

1. **Conics.** With $t=\tan(\theta/2)$, $(\cos\theta,\sin\theta) = \bigl(\tfrac{1-t^2}{1+t^2},\tfrac{2t}{1+t^2}\bigr)$. A circular arc of half-angle $\alpha$ is a rational quadratic Bézier with middle weight $w=\cos\alpha$. With rational control data, evaluation at rational parameters is exact in $\mathbb{Q}$.
2. **Rotations.** Rigid motions are the main source of irrational coordinates. Rotations are stored as unit quaternions from the rational Euler–Rodrigues parameterisation. An arbitrary user angle is either kept symbolic (a pair $(c,s)$ with the relation $c^2+s^2=1$) or replaced by a rational rotation within a recorded bound.
3. **Curvature in a direction.** Euler's theorem, $\kappa_n(\theta) = \kappa_1\cos^2\theta + \kappa_2\sin^2\theta$, and its third-order analogue (§1.7) are trigonometric polynomials. They give the compatibility conditions for curve networks at a vertex.

A known limit: polynomial curves with rational offsets are the Pythagorean-hodograph curves (Farouki). General NURBS offsets are not rational and are kept procedural.

### 1.7 G2 and G3: exact conditions

**Jet form.** Let $S_1,S_2$ share the boundary $\Gamma(t) = S_1(0,t) = S_2(0,t)$. $G^k$ means there is a reparameterisation $\varphi$ with $\varphi(0,t)=(0,t)$ such that $S_2$ and $S_1\circ\varphi$ have equal $k$-jets along $\Gamma$. Writing the transversal Taylor coefficients of $\varphi$ as gauge functions $a,b,c,d,e,f$ of $t$ (with $a>0$), the chain rule gives, at $u=0$:

$$
\begin{aligned}
\text{G1:}\quad \partial_u S_2 &= a\,S_{1,u} + b\,S_{1,t}\\[2pt]
\text{G2:}\quad \partial_u^2 S_2 &= a^2 S_{1,uu} + 2ab\,S_{1,ut} + b^2 S_{1,tt} + c\,S_{1,u} + d\,S_{1,t}\\[2pt]
\text{G3:}\quad \partial_u^3 S_2 &= a^3 S_{1,uuu} + 3a^2b\,S_{1,uut} + 3ab^2 S_{1,utt} + b^3 S_{1,ttt}\\
&\quad + 3ac\,S_{1,uu} + 3(ad+bc)\,S_{1,ut} + 3bd\,S_{1,tt} + e\,S_{1,u} + f\,S_{1,t}
\end{aligned}
$$

For curves the same construction gives the β-constraints: $c_2' = \beta_1c_1'$, $c_2''=\beta_1^2c_1''+\beta_2c_1'$, $c_2''' = \beta_1^3c_1''' + 3\beta_1\beta_2c_1'' + \beta_3c_1'$.

**Invariant form.** Let $\tau$ be the unit tangent of $\Gamma$ and $\nu$ the unit conormal (tangent to $M$, normal to $\Gamma$). Write $[\![\cdot]\!]$ for the jump across $\Gamma$.

- **G2 reduces to one scalar.** Given G1, tangential data already agree, so $[\![b(\tau,\tau)]\!] = [\![b(\tau,\nu)]\!] = 0$ follows from differentiating $n$ along $\Gamma$. G2 holds if and only if $[\![b(\nu,\nu)]\!]=0$. This is the linkage-curve theorem of Pegna and Wolter.
- **G3 reduces to one scalar.** Given G2, $b$ is continuous on $\Gamma$, so its tangential derivative is continuous: $[\![C(\cdot,\cdot,\tau)]\!]=0$. Total symmetry of $C$ (Codazzi) then forces every component with at least one $\tau$ index to be continuous. The only free component is

$$
[\![\,C(\nu,\nu,\nu)\,]\!] = [\![\,(\nabla_\nu b)(\nu,\nu)\,]\!] = 0 .
$$

  So G3 across a G2 seam is the continuity of the rate of change of normal curvature in the cross-boundary direction. This second reduction is derived here from Codazzi symmetry and is not quoted from a reference. It has since been checked numerically on seams that are reparameterised nonlinearly (`tests/geom/continuity_test.cpp`, RESULTS §6): the predicted jump of $C(\nu,\nu,\nu)$ is matched to about $10^{-11}$ and every other component stays below $10^{-12}$. A numerical check on examples is not a proof, so a literature check is still advisable before publication.

These two scalars, together with $[\![n]\!]$, are the quantities the kernel bounds and reports. Reflection and highlight lines are $C^{k-1}$ on a $G^k$ surface, which is why G3 is what removes the visible kink in the *curvature* of a highlight.

**Vertex compatibility.** Where $N$ boundary curves meet at a vertex, a $G^k$ fill exists only if the curves' $k$-jets lie in the $k$-jet of a single surface. In Monge form over the tangent plane,

$$
z = \tfrac12\,b_{\alpha\beta}x^\alpha x^\beta + \tfrac16\,c_{\alpha\beta\gamma}x^\alpha x^\beta x^\gamma + O(4),
$$

a curve leaving in direction $\theta_i$ with geodesic curvature $\kappa_{g,i}$ must satisfy

$$
\begin{aligned}
z_i'' &= b_{11}\cos^2\theta_i + 2b_{12}\cos\theta_i\sin\theta_i + b_{22}\sin^2\theta_i,\\
z_i''' - 3\kappa_{g,i}\,b(T_i,T_i^\perp) &= c_{111}\cos^3\theta_i + 3c_{112}\cos^2\theta_i\sin\theta_i + 3c_{122}\cos\theta_i\sin^2\theta_i + c_{222}\sin^3\theta_i .
\end{aligned}
$$

The quadratic form has 3 coefficients and the cubic form 4. A network with $N$ curves at a vertex therefore carries $N-2$ conditions for G1 (coplanar tangents), $N-3$ for G2 and $N-4$ for G3. The solver checks these before attempting a fill and reports which curve violates them and by how much. When full cross-boundary ribbons are prescribed, twist (mixed-partial) compatibility at corners adds further conditions.

**Governing differential equations.** Linearising about a parameter plane, the energy $E_m(S) = \int_\Omega\lVert D^mS\rVert^2$ has Euler–Lagrange equation $\Delta^mS=0$, whose well-posed Dirichlet data are $\partial_\nu^jS$ for $j=0,\dots,m-1$. A $G^k$ boundary prescribes $j=0,\dots,k$. So:

| Target | Order $m=k+1$ | PDE | Energy | Basis regularity needed |
|---|---|---|---|---|
| G1 | 2 | $\Delta^2S=0$ (biharmonic) | thin plate | $C^1$ |
| G2 | 3 | $\Delta^3S=0$ (triharmonic) | variation of curvature | $C^2$ |
| G3 | 4 | $\Delta^4S=0$ | variation of curvature rate | $C^3$ |

This rule has a practical consequence. Imposing G2 or G3 boundary conditions on a pure thin-plate energy is over-determined in the continuum. It appears to work in a fixed spline space but produces boundary layers that worsen under refinement. The energy order must match the continuity order. The implementation confirms this (RESULTS §5): with G2 jets prescribed, the third-order energy of the thin-plate solution grows roughly linearly with the number of control points, while the matched energy gives the same value at every resolution.

The parameterisation-independent versions used for final convergence are

$$
E_2 = \int_M \lVert b\rVert^2\,dA = \int_M(\kappa_1^2+\kappa_2^2)\,dA,\qquad
E_3 = \int_M\lVert\nabla b\rVert^2dA,\qquad
E_4 = \int_M\lVert\nabla^2 b\rVert^2dA .
$$

**Discrete structure.** For a spline patch, $G^k$ across a boundary constrains the $k+1$ control rows next to it. For fixed gauge functions the constraints are linear in control points, $C(\gamma)\,x = d(\gamma)$; jointly in $(\gamma,x)$ they are bilinear. Independence of opposite boundaries needs at least $2(k+1)$ control rows per direction: degree 5 for G2 and degree 7 for G3 in a single Bézier span, or lower degree with interior knots. If $S_1$ has degree $p$ along $\Gamma$ and the gauge has degree $q$, exact satisfaction needs degree $p+q$ in $S_2$. Otherwise the kernel solves the ε-$G^k$ problem and certifies the residual jumps with interval arithmetic.

---

## 2. Historical error learning and fault analysis

### 2.1 Root causes in legacy kernels

The published tolerance constants give the scale: Parasolid uses a linear resolution of $10^{-8}$ in a size box of 1000 units; ACIS uses `resabs` $=10^{-6}$ and `resnor` $=10^{-10}$; OCCT uses `Precision::Confusion()` $=10^{-7}$ with per-sub-shape tolerances. These are quoted from memory of public documentation and should be checked against current releases.

| Failure | Mechanism | Countermeasure here |
|---|---|---|
| Floating-point inconsistency | The damage is not the size of rounding error. Independent numerical decisions contradict each other, and tolerant equality is not transitive ($a\approx b$, $b\approx c$, $a\not\approx c$). Combinatorial topology built on contradictory answers is invalid. | Filtered exact predicates. Every predicate returns a proven sign or `Undecided`; it never returns a guess. |
| Tolerance inflation ("epsilon bubbles") | A tolerant edge or vertex gets a tolerance covering its fit error. Each later operation takes the maximum of its inputs plus new fit error, so tolerance only grows. Vertex spheres eventually swallow short edges and overlap each other. | T3: derived entities are recomputed from carrier definitions, so tolerance is not inherited. Merge is refused when $\tau(e)\ge\tfrac12\operatorname{lfs}(e)$ (local feature size). |
| Self-intersecting offsets | Local: $1-d\kappa_i\le0$. Global: $\lvert d\rvert\ge$ reach. Offsets of NURBS are not NURBS, so they are approximated, and the approximation hides the singularity. | Offsets are procedural carriers. Both predicates of §1.3 are evaluated with intervals before any geometry is built. |
| Boolean topology corruption | (a) Near-tangent intersection: error scales as $1/\sin\theta$. (b) Missed small loops when marching. (c) Coincident faces detected by sampling. (d) Sliver faces from tolerant merging. (e) Inconsistent point classification. | (a) Stability indicator and escalation. (b) Loop-free certificates by Gauss-map separation (Sederberg–Meyers). (c) Coincidence by provenance or exact algebra. (d) No merging without a ledger entry. (e) Classification from the same certified predicates. |
| Blend and sweep failures | Ball radius exceeds local curvature radius; sweep profile crosses the path's centre of curvature ($1-\kappa\,x\le0$ for a profile point at offset $x$ along the path normal); degenerate poles. | Same offset predicate; sweep-Jacobian predicate; explicit pole handling in the evaluator. |
| Exchange damage | Each system writes STEP at its own tolerance; the reader re-interprets edges at a different one. | Healing by re-intersection and variational snapping (§3.4), with deviation reported. |

### 2.2 Error-learning feedback loop

The loop has a certified layer that guarantees correctness and a learned layer that only improves speed and warns early. The learned layer never decides topology.

**(a) Numerical provenance.** Every geometric evaluation returns a triple: value, certified enclosure, condition estimate. These propagate through the operation graph the way a derivative propagates through automatic differentiation.

**(b) Stability indicators.** Each operation has a dimensionless indicator that tends to 0 at collapse.

| Operation | Indicator $\sigma$ | Collapse |
|---|---|---|
| Surface–surface intersection | $\sigma_{\min}(J)$, proportional to $\sin\angle(n_1,n_2)$ | tangency |
| Newton refinement | Kantorovich $h=\beta L\eta$, or Smale $\alpha=\beta\gamma$ | $h>\tfrac12$; $\alpha>\alpha_0=\tfrac{13-3\sqrt{17}}{4}\approx0.1577$ |
| Offset, shell, blend | $\min_i(1-d\kappa_i)$ and $1-\lvert d\rvert/\operatorname{reach}$ | cusp, global overlap |
| Sweep | $\min(1-\kappa_{\text{path}}\,x)$ over the profile | swallowtail |
| Arrangement in a face domain | separation ratio $\rho$ = (distance between distinct features) / (enclosure width) | $\rho\le1$ |
| Any entity | tolerance ledger $\tau(e)/\operatorname{lfs}(e)$ | $\ge\tfrac12$ |

**(c) Residual estimator.** For an approximate intersection point $\tilde x$ with residual $r=F(\tilde x)$, the Newton–Kantorovich theorem turns backward error into certified forward error. With $\lVert G'(\tilde x)^{-1}\rVert\le\beta$, $\lVert G'(\tilde x)^{-1}G(\tilde x)\rVert\le\eta$ and Lipschitz constant $L$ of $G'$ (from interval second derivatives):

$$
h=\beta L\eta\le\tfrac12\;\Longrightarrow\;\text{a root exists within } \frac{1-\sqrt{1-2h}}{\beta L}\ \text{of }\tilde x .
$$

For approximants (tessellation, cached splines) the bounds are a-posteriori. Filip et al. give, for a linear interpolant on a triangle with longest edge $\ell$, $\lVert S-\tilde S\rVert_\infty\le\tfrac18\ell^2(M_1+2M_2+M_3)$ with $M_1,M_2,M_3$ bounding $\lVert S_{uu}\rVert,\lVert S_{uv}\rVert,\lVert S_{vv}\rVert$. Nairn–Peters–Lutterkort bound the distance from a Bézier piece to its control polygon. In optimal-recovery terms (Micchelli–Rivlin), these are bounds on the radius of information: the worst-case error over all surfaces consistent with the samples and the derivative bound.

**(d) Convergence-rate regression.** At refinement levels $h_\ell$ the kernel observes estimators $\eta_\ell$ and fits $\log\eta_\ell = \log C + p\log h_\ell$. Smooth problems converge at the theoretical order (2 for chordal error, 4 for cubic approximants, quadratic for Newton). An observed $\hat p$ below the theoretical order means a singularity is close. This is a cheap, model-free detector that fires *before* the indicator in (b) reaches its threshold, and it also predicts the refinement level needed to meet a budget.

**(e) Learned scheduling.** A regression model predicts, from features of an operation, the precision tier it will need, its cost, and its failure probability.

- Features: the indicators in (b), degree and knot statistics, bounding-volume overlap, dihedral-angle histogram, ledger state, provenance depth.
- Models: gradient-boosted trees per predicate class; a graph network over the face-adjacency graph per operation.
- Training data: telemetry from the benchmark corpus (§5.4) plus fuzzing by controlled perturbation (rotate-and-union, near-tangent sweeps).
- Use: choose the starting tier, order the work, and raise a pre-flight warning ("this blend is predicted to fail at radius 4.2; largest safe radius ≈ 3.7").
- Guard: a wrong prediction costs time, never correctness, because the certified layer still decides.

**(f) Transactions.** Every operation runs against a persistent (structurally shared) topology and either commits or rolls back with a diagnostic naming the entity, the indicator, and the tier reached.

---

## 3. Unifying solids and NURBS in one heterogeneous engine

### 3.1 Hybrid data structure

Three layers, with strict direction of reference: topology → carriers → provenance.

```
Topology (non-manifold cell complex, persistent)
  Body ─ Region ─ Shell ─ Face ─ Loop ─ Coedge ─ Edge ─ Vertex
                           │               │        │       │
Carriers (one definition each)             │        │       │
  SurfaceDef ◄─────────────┘               │        │       │
  EdgeVariety(S1,S2) + tube ◄──────────────┴────────┘       │
  PointVariety(system) + isolating box ◄────────────────────┘
        │
Evaluators of a carrier (each returns value + enclosure)
  parametric : NURBS | analytic | SubD (exact evaluation) | procedural (offset, blend, sweep)
  implicit   : polynomial | matrix representation | certified distance field
  discrete   : tessellation cache with error bound
        │
Provenance DAG
  node = constructive definition with exact inputs (rationals, symbols)
  edges = "defined from"; equality of nodes is structural (hash-consed)
```

Key properties:

- **Co-parameterisation.** A coedge does not own a p-curve. It owns the projection of the edge variety onto its face's two coordinates of $\mathbb{R}^4$. The 3D curve and both p-curves are three views of one tube.
- **Vertices** are zero-dimensional varieties (three surfaces, or curve and surface) with an isolating box. Two vertices are equal if their provenance nodes are equal, or if exact algebra proves it.
- **Fields are carriers.** An implicit body (lattice, level set, distance field) is a `SurfaceDef` whose primary evaluator is implicit. It participates in Booleans through the same edge-variety mechanism, with $F = (\phi_1(x),\phi_2(x))$ on $\mathbb{R}^3$ or $\phi_2(S_1(u,v))$ on $\Omega_1$.
- **Subdivision surfaces** are parametric carriers: regular regions are bicubic B-splines, and extraordinary regions are replaced by $G^2$ caps (§4.1) so that the carrier meets the kernel's continuity contract.
- **Tolerance ledger.** $\tau(e)\ge0$ is stored per entity and is 0 unless a `MergeRecorded` decision was taken. Invariants: $\tau(e)<\tfrac12\operatorname{lfs}(e)$, and $\tau$ is never propagated to derived entities.

### 3.2 Arithmetic model

Predicates use a filter cascade. Each tier either proves the sign or passes the question up.

| Tier | Arithmetic | Typical share of calls (hypothesis) |
|---|---|---|
| 0 | `double` with static error bound | > 95 % |
| 1 | Interval arithmetic, outward rounding | most of the rest |
| 2 | Double-double / quad-double expansions (error-free transforms) | rare |
| 3 | Multi-precision ball arithmetic, precision doubling | very rare |
| 4 | Symbolic: provenance identity, exact rational or algebraic computation | degenerate cases |

What happens when tier 3 reaches its precision cap and tier 4 cannot prove identity is a policy decision, isolated in one function (`on_undecided`, see `include/akernel/num/escalation_policy.hpp`).

### 3.3 Boolean operations

Result faces keep their **original carriers**. No surface is refitted, so Class-A quality is preserved exactly.

1. **Broad phase.** Bounding-volume hierarchy over certified hulls of control nets; conservative, so no true pair is missed.
2. **Coincidence by provenance.** Pairs of faces whose carriers are the same provenance node, or provably the same after exact knot insertion and degree elevation, are routed to a 2D overlap computation in the shared parameter domain.
3. **Certified intersection topology.** Subdivide patch pairs until each pair is one of:
   - *separated*: hulls disjoint;
   - *simple transversal*: Gauss maps separated, which excludes closed loops, and boundary crossings counted exactly;
   - *singular*: tangency suspected. Resolved by provenance (a fillet is tangent to its supports by construction), then exact algebra for analytic pairs, then the escalation policy.
4. **Tube construction.** Trace each branch with the Krawczyk certificate of §1.4. Branch points are roots of the augmented system $F=0,\ n_1\times n_2=0$.
5. **Arrangement per face.** Order intersection points along each curve using isolating intervals. Build the planar arrangement in $\Omega$. A pair of points that cannot be separated and cannot be proven identical goes to `on_undecided`.
6. **Classification.** Each arrangement cell is classified in/out/on against the other body using the composite field sign of §1.4, with local neighbourhood analysis at edges for regularisation.
7. **Assembly and audit.** Build the result complex. Check $\partial_1\partial_2=0$, shell closure and Euler–Poincaré. Commit or roll back.

**As implemented (v0.4, `op/boolean.hpp`).** The implementation covers steps 3 to 7 for solids with trimmed faces in general position. Operands and results are the same type (`op/brep.hpp`), so a result can be an operand. It differs from the plan above in these places:

- *Vertices.* A new vertex is a certified root of "edge of one solid against surface of the other", kept if it is proven to lie inside a face on that surface. For an edge that is a patch boundary this is three equations in three unknowns. For an edge that is itself an intersection curve of two surfaces (every edge made by an earlier Boolean) it is a three-surface point: six equations in six unknowns, searched along the edge's certified tube (`op/tsi.hpp`). Both use a Krawczyk test. An intersection curve ends at a vertex only if the vertex is proven to lie on the curve's certified arc. A regular root also proves that the curve crosses the face boundary there, which the surface/surface certificate alone does not.
- *Clipping.* Surfaces are intersected whole and the curves are then cut at the vertices on them. A curve can enter or leave a face only at such a vertex, so one certified point per part decides whether that part lies in a face of each solid.
- *What an entity stores.* A vertex stores certified parameter enclosures on every surface it lies on; an intersection edge stores its tube segments and its two certified end positions. No entity stores a tolerance, and a chained operation reuses these certificates instead of re-deriving geometry from fitted data.
- *Completeness.* Loops are found by the subdivision pass of §1.4 (exclude or cover every sub-box), not by Gauss-map separation.
- *Classification.* At one certified point of each intersection curve, an interval evaluation of $(\nu\cdot n_{\text{other}})$ decides which side of the curve lies inside the other solid; the sign is constant along a regular curve. Cells without curves inherit across uncut edge pieces, and a shell that meets nothing is classified by counting certified ray crossings. The composite field of §1.4 is not used.

- *Coincidence (step 2), v0.5.* Implemented by identity of construction, without a provenance DAG: two surfaces are one if their control points are identical, a side is shared if its boundary control points are, a vertex is shared if its certified parameters on a shared surface are identical, and two edge pieces are one if they lie on the same curve between the same two vertices. A side shared by two surfaces lies on the boundary of both parameter squares, where subdivision cannot separate them; it is certified directly as a known tube and handed to the completeness proof. In a shared surface the faces of both solids are overlaid in its one parameter square (the other solid's edges that lie inside a face are inserted into that face's planar map, ordered at each vertex by interval cross products), and each region is classified as in, out, on with the same outward normal, or on with the opposite one. Selection keeps an on-region once, from the first operand: with the same normal for union and intersection, with the opposite normal for difference. - *Coincidence decided by exact algebra, v0.6.* For flat faces the kernel decides coincidence that is not by construction, by exact predicates on the input data (`op/planar.hpp`): a patch is an affine parallelogram if its degree is (1,1), its weights are 1 and c00 + c11 = c10 + c01 exactly; two such patches are one plane if the corners of one are exactly coplanar with three corners of the other (exact `orient3d`). Planes are grouped into classes, with interval-certified affine maps between their parameter planes, and the overlay above works across those maps. Further exact facts used: a side of one patch lying exactly in another's plane is a known line of their intersection; two vertices that lie, by construction, on the same three planes with independent normals are one point; two straight edge pieces between the same two vertices are one piece. A patch whose corners were rounded off a plane is, correctly, not planar and is never merged. Not covered: coincident curved surfaces that are not identical patches (they would need exact algebra for quadrics and beyond), and non-manifold results.

*Broad phase (step 1), v0.12.* A static bounding-volume hierarchy (`geom/box_index.hpp`) over certified boxes: of each face (its carrier over a parameter rectangle that holds the face), each edge and each vertex. Boxes enclose their entities for certain, so two entities whose boxes are disjoint cannot meet, and pruning by the index changes no answer. Queries report ids in increasing order, so code that scanned ids in order keeps its order and its results. It is used wherever one solid's entities were tested against all of the other's: surface pairs to intersect, surfaces an edge may pierce, vertices of one operand that may be vertices of the other, vertices that may lie on an edge, and edge pieces that may lie in a face. Two cheap certified filters come before the certified narrow phase: a flat surface is skipped against a face whose control points all lie strictly on one side of its plane (interval signs; the face lies in the convex hull of its control points), and planes are grouped into classes by testing each against one member per class, after an interval key (unit normal and distance, up to sign) has ruled most classes out. A vertex carries parameters only on surfaces of its plane whose faces or edges are near it, not on all of them (a plane tiled into n faces would otherwise give each vertex n charts). The index is built over faces, not over the patches of the control net as planned above; with one Bezier patch per carrier the two coincide. Measurements are in [RESULTS.md](RESULTS.md) §3 (scaling).

*Plane pairs (steps 3 and 4 for two flat faces), v0.12.* For two exact parallelograms in non-parallel planes the solution set in $(u,v,s,t)$ is one line $\ell = P + \sigma d$, with $P$ and $d$ enclosed by intervals from the corners. A certified tube segment has exactly one solution on each slice, and every solution lies on $\ell$, so each segment covers a $\sigma$-range bounded from inside by intervals; if those ranges overlap and contain an outer bound of the $\sigma$-range inside the domain, the intersection is complete. This replaces the completeness subdivision of §1.4 for that case (`BooleanOptions::plane_lines`); anything that cannot be certified falls back to it.

*Three-surface searches, v0.12.* The search for points where an earlier intersection edge meets a third surface (`op/tsi.hpp`) keeps the enclosures of pieces of the edge's tube per edge, since they do not depend on the surface. It tries cheap exclusions (bounding boxes and extents along the surface normal, on boxes known to contain the piece) before computing a tight enclosure, and stops the Krawczyk contraction at a fixed point. None of this changes a result; it removes about two thirds of that stage's time ([RESULTS.md](RESULTS.md) §3).

*Tangential contact with a plane (step 3, singular case), v0.13.* For a plane through the corners c00, c10, c01 of a flat face and a patch X, the plane function $g(u,v)=\det(c_{00},c_{10},c_{01},X(u,v))$ is a polynomial with exact Bernstein coefficients (`op/contact.hpp`). Where it vanishes on sides of the square it is divided by their linear factors, exactly, as often as they stay zero (the order of contact: 2 for a cylinder along its line, 1 for a transversal crossing). If the quotient has one strict sign over the square, proven from its coefficients after exact subdivision, the patch meets the plane in those sides and nowhere else: the pair contributes no curve and needs no completeness search, and vertices on those sides get parameters on the plane. The same division in one variable proves that an edge meets a plane only at its end vertices. Consequences used downstream: a vertex on the plane and on such a patch has an exact parameter across the side; in a face's parameter plane, a curve leaving a vertex tangent to a side of the square is ordered by its known side (it lies in the face) once no other curve can lie in its direction cone; around an edge where the solids touch along the line, the curved and the flat face are ordered by the side of the plane on which the curved one lies. Not covered: contact along a line inside a patch, tangency at a point, and tangency between two curved surfaces. A transversal curve that meets both parameter squares only on their boundaries is handled separately: the pair is searched again over squares enlarged by 1/32, where the curve can be traced through the touching point; completeness over the enlarged box implies it over the original. Measurements: [RESULTS.md](RESULTS.md) §3, "Tangential contact with planes".

*Flat faces on rational carriers (step 2), v0.14.* The coincidence overlay needs an affine map between the parameter planes of two flat carriers, which a rational flat carrier (a sector's annular cap) does not have with a parallelogram. Since a face may lie on any carrier of its surface, such a face is moved before the Boolean onto an exact parallelogram of its plane (axis-aligned planes; corners on a dyadic grid, power-of-two extent). Its straight edges become sides of thin exact parallelograms on the same lines; its curved edges become arc edges, the transversal intersections of the neighbouring surfaces with the plane, certified as tubes. This happens only where the conflict occurs and only on copies of the operands, so every result that could be computed before is unchanged. Measurements: [RESULTS.md](RESULTS.md) §3, "Flat faces on rational carriers".

*Contact left by coincidence handling (steps 2, 4 and 5), v0.15.* Four rules, each a proof (the first three by exact arithmetic, the last by intervals):

- *Coincidence with the parameters exchanged.* An inner and an outer face of one cylinder are related by an affine map that exchanges the two parameters. The exact comparison of step 2 is applied to the transposed patch as well, and the classes of coincident surfaces carry such maps (they reverse orientation).
- *A vertex in a plane at the end of an edge.* A regular root of the edge/surface system that cannot be separated from an end of the edge is, for a plane, decided by the plane function of the edge's side: if it vanishes only at the side's ends, the root is the end vertex and the edge has no other point in the plane.
- *One curve, two representations.* A curved side of a patch and an intersection curve between the same two vertices are one edge when a point of the side, carried exactly into the curve's two surfaces, lies on the certified curve. The same holds for straight sides of any degree: control points exactly on one line make a side a segment, fixed by its end vertices.
- *Point location against enclosed ends.* The winding number of a loop about a point is computed on a polygon of chords. It is valid when each chord can be moved onto its true piece without crossing the point. For a straight piece whose ends are enclosures $A$ and $B$, chord and piece both lie in $\operatorname{conv}(A\cup B)$, and a strict sign of $(b-a)\times(q-a)$ over $a\in A$, $b\in B$ proves that $q$ is outside it.

One rule of method came out of a fault the audit caught (a face that was not split by edges lying in it): *an undecided answer is never given a meaning.* "Undecided" had been read as "on the boundary" at one place; it is now reported. Measurements: [RESULTS.md](RESULTS.md) §3, "Contact cases left after the flat-face work".

*Exact points (a first step of the algebraic tier of §3.2), v0.16.* Interval enclosures cannot prove a coincidence, and models with round dimensions consist of them. For rational data many are provable without algebraic numbers, from one lemma: a certified root comes with a box in which it is the only solution, so *a point that lies exactly on every surface of the system at parameters inside that box is the root*. Candidates are cheap (an enclosure's midpoint rounded to a dyadic grid, on a surface or in space; a fraction with a small denominator by continued fractions), and the proof is exact: de Casteljau steps in expansion arithmetic, in homogeneous form for rational parameters, compared by cross-multiplication; an exact determinant for a plane. A vertex proven this way carries its exact homogeneous point, and results pass it on. With it the kernel decides, exactly and for curved surfaces as well as planes, whether a vertex lies on a surface (which resolves "the surface passes through a vertex of the edge" where the vertex is an end of the edge), and whether two vertices are the same point. A candidate that fails leaves an ordinary enclosure, so nothing depends on guessing well. The limits are those of the representation: coordinates that are doubles, and parameters that are dyadic or fractions with denominators up to 4096. Points with irrational coordinates need the algebraic numbers of §3.2. Measurements: [RESULTS.md](RESULTS.md) §3, "Exact points for vertices on round numbers".

*Edges in curved surfaces (steps 2 to 6), v0.16.* The same exact evaluation proves that a straight edge lies in a curved surface along a line of constant parameter (a ruling of a cylinder): the patch is linear in the other parameter, and its exact points at the line's two ends lie on the edge's line. The edge is then handled as an edge in a plane: a known line of the surface pair, no edge/surface search, a straight piece in the surface's faces, no classification across it. Where the edge's solid stays on one side of the surface (decided with intervals at one point of the edge, which suffices because the curves along the edge are certified regular), the edge bounds nothing in the face, and a result that keeps both has a *line contact*, reported as a non-manifold feature. Measurements: [RESULTS.md](RESULTS.md) §3, "Edges lying in a curved surface".

*Orderings (steps 3 and 5), v0.16.* Two orderings that intervals alone cannot settle are decided without new arithmetic. Two roots with the same exact point are one vertex, whichever edges they were found from. And at a vertex, a curve that continues the way back (its direction enclosure points strictly backwards and straddles that line) is ordered against each other curve by the strict sign of their cross product over the enclosure: the enclosure contains the exact opposite direction, whose order is fixed by the other curve's side, and a strict sign means the true direction cannot pass the other curve. Measurements: [RESULTS.md](RESULTS.md) §3, "Ordering cases".

*Edges through edges (step 2), v0.16.* An edge that pierces a surface exactly on the boundary of a face there passes through an edge of the other solid. Where the root has an exact point, that point is tested against the other solid's edges: it must lie exactly on both of an edge's surfaces and strictly inside it. The root is then a vertex on both edges and four surfaces. On planes, each parameter of an exact point is made exact on its own where it is dyadic, so that a vertex at a third of the way along a side is still exactly on it. Measurements: [RESULTS.md](RESULTS.md) §3, "Edges through edges".

*Single cases (steps 2 to 4), v0.16.* Four more uses of the same tools: a moved cap takes along arcs that lie on other flat carriers of its plane; an end vertex whose parameter is an enclosure (1/3) is the root when that enclosure lies in the root's uniqueness box; straight sides are matched with straight intersection curves (a plane and a cylinder along a ruling); and a vertex with an exact point is decided on a side through parameters proven exact from that point. Measurements: [RESULTS.md](RESULTS.md) §3, "Single cases".

*Classification across edges in the other solid's planes (step 6), v0.16.* Inside/outside is carried from region to region across an operand's own edges, which is valid only where the edge is not on the other solid's boundary. A straight edge lying exactly in a plane that holds faces of the other solid may be on it, and is never crossed; its regions are decided by other neighbours or by a ray test. (This had been applied only in operations with a tangential contact somewhere.) With it, no refusal from a degenerate configuration of vertices and edges is left in the curved benchmark chains. Measurements: [RESULTS.md](RESULTS.md) §3, "The classification clash".

*Tangency along an interior line (step 3, singular case), v0.16.* The plane function of a patch of degree 2 in one parameter has a double root on every control row at the same r exactly when the plane is tangent to the patch along the line u = r; if the rows' leading coefficients share one strict sign, g = (u − r)² h with h strictly signed, and the contact is exactly that line (a plane tangent to a cylinder along a ruling). The pair contributes no curve, edges in the plane or on the surface are not searched against the other, curves that run along the line are matched to the operands' edges there, and touching regions are reported as a line contact. Tangency between two curved surfaces (two hill tops crossing where they share a boundary line, where the difference of their height functions factors) remains open. Measurements: [RESULTS.md](RESULTS.md) §3, "Tangency: a plane along a ruling inside a cylinder".

Why tangential cases cannot be handled by precision alone: a tangency is a root of an over-determined system. Interval methods can prove that a root does *not* exist, but cannot prove that one does, because existence is destroyed by arbitrarily small perturbation. Proof needs exact algebra or knowledge of how the geometry was constructed. That is the mathematical reason for T1.

**Determinism.** Results are bitwise reproducible across thread counts and platforms: canonical ordering of entities, deterministic parallel reductions, no transcendental library calls inside predicates (rational trigonometry instead), fixed contraction rules for fused multiply-add. $A\cup B$ and $B\cup A$ produce identical output.

**Alternative mode.** For downstream consumers that need untrimmed, analysis-suitable geometry, a separate *watertight rebuild* operation reparameterises faces so that trim curves become knot lines (Urick et al.). This changes surfaces within a stated bound and is never applied implicitly.

### 3.4 Healing imported data

Imported STEP/IGES models arrive with gaps, because the writer stored fitted curves. Healing is posed as inference followed by optimisation.

1. **Topology inference by persistence.** Build the candidate edge-pairing graph at scale ε. Track $\beta_0$ and $\beta_1$ of the sewn complex as ε increases. Choose pairings from the range where the invariants are stable. This replaces a single global stitch tolerance.
2. **Re-intersection.** Where two paired faces meet transversally, discard the imported edge curve and recompute the edge as the variety of the two (slightly extended) carriers. The gap becomes zero by construction and both surfaces are untouched.
3. **Variational snapping.** Where faces meet tangentially (blends), re-intersection is ill-conditioned. Solve for control-point displacements $\delta x$:

$$
\min_{\delta x}\;\sum_{\text{edges}}\int\lVert S_i(c_i(t)) - S_j(c_j(t))\rVert^2dt
\;+\;\lambda\,\lVert\delta x\rVert_M^2
\;+\;\mu\,\bigl(E_3(S+\delta S)-E_3(S)\bigr)
\quad\text{s.t.}\quad C_{G^k}(x+\delta x)=0 .
$$

   The first term closes gaps, the second keeps the change minimal in a mass-matrix norm, the third forbids loss of fairness, and the constraints restore the continuity class inferred from the import (detected from the statistics of normal and curvature jumps). Analytic carriers move only through their defining parameters, so a plane stays a plane.
4. **Report.** Per face: maximum deviation from the imported geometry, continuity class achieved, and remaining ledger entries. Entries that could not be healed stay as recorded tolerances.

---

## 4. Class-A surface engine

### 4.1 Multi-sided blend patches

Problem: given $N$ boundary curves with $G^3$ jet ribbons

$$
R_i(s,\delta) = \sum_{j=0}^{3}\frac{\delta^j}{j!}\,\partial_\nu^jS_i(s),
$$

find a fill surface. Algorithm:

1. **Compatibility check** at each corner using the conditions of §1.7. Report or repair violations (minimal-norm correction of the offending curve's end jets).
2. **Layout.** Either a single trimmed four-sided patch when $N\le4$ or one side is short, or a split into $N$ quads around a central point with caps of bidegree 7. Published $G^2$ constructions use bidegree 7 (Loop–Schaefer) or lower with multiple rings (Karčiauskas–Peters). A $G^3$ cap of that kind is a research item (§5.5); the fallback is ε-$G^3$ at the central point with a certified bound.
3. **Linear solve.** Minimise $w_3E_3^{\text{lin}} + w_4E_4^{\text{lin}}$ subject to $C(\gamma)x=d(\gamma)$ with gauges initialised from boundary parameter speeds.
4. **Gauge update and geometric polish.** Alternate between updating $\gamma$ and Gauss–Newton steps on the geometric $E_3+wE_4$ until the certified jumps meet the specification.

### 4.2 Fillet transitions

A rolling-ball fillet is only G1: curvature jumps from the support's normal curvature to $1/r$. For G3, each cross-section must match position, tangent, curvature and curvature derivative at both rails. That is 8 Hermite conditions per coordinate, so the cross-section is degree 7 (or lower degree with interior knots), constrained by the β-conditions of §1.7. The remaining freedom is used to shape the curvature profile: monotone rise to a single interior extremum, no interior inflection. Cross-sections are then coupled along the spine by the same variational solver, so the fillet is fair along its length and not just in each section. Validity is checked first with the offset predicate of §1.3.

### 4.3 Crown lofting

Crown is deliberate slight convexity that keeps a panel from reading as hollow. It is expressed as an inequality on directional curvature, not as hand-moved control points:

$$
\min_x\;E_3+wE_4\quad\text{s.t.}\quad S(\text{sections})\ \text{interpolated},\qquad
\kappa_n(u,v;\,\mathbf{d})\ge\kappa_{\min}(u,v),\qquad \text{no inflection of highlight lines along }\mathbf d .
$$

Linearised, the curvature bound is linear in control points, so each step is a convex quadratic programme.

### 4.4 Variational spline solver

**Energy.** The mechanical model is the Kirchhoff–Love thin shell relative to a reference surface $(\bar g,\bar b)$:

$$
W=\frac12\int_M\Bigl(h\,A^{\alpha\beta\gamma\delta}\varepsilon_{\alpha\beta}\varepsilon_{\gamma\delta}
+\frac{h^3}{12}A^{\alpha\beta\gamma\delta}\rho_{\alpha\beta}\rho_{\gamma\delta}\Bigr)dA,
\qquad
\varepsilon_{\alpha\beta}=\tfrac12(g_{\alpha\beta}-\bar g_{\alpha\beta}),\quad
\rho_{\alpha\beta}=b_{\alpha\beta}-\bar b_{\alpha\beta},
$$

$$
A^{\alpha\beta\gamma\delta}=\frac{E}{1-\nu^2}\Bigl[\nu\,g^{\alpha\beta}g^{\gamma\delta}+\tfrac{1-\nu}{2}\bigl(g^{\alpha\gamma}g^{\beta\delta}+g^{\alpha\delta}g^{\beta\gamma}\bigr)\Bigr].
$$

The membrane term controls stretching of the parameterisation; the bending term is the thin-plate strain energy. For $G^2$ and $G^3$ work the higher-order terms $E_3,E_4$ are added, as required by the order rule of §1.7. Using the design basis for the mechanics is the isogeometric approach (Hughes et al.; Kiendl et al.), so the same assembly code serves design and shell analysis.

**Discrete problem.** With control vector $x$:

$$
\min_x\ \tfrac12x^\top Kx-f^\top x\quad\text{s.t.}\quad Cx=d,\qquad
\begin{bmatrix}K&C^\top\\C&0\end{bmatrix}\begin{bmatrix}x\\\lambda\end{bmatrix}=\begin{bmatrix}f\\d\end{bmatrix}.
$$

- Assembly by Gauss quadrature on Bézier elements (Bézier extraction); elements are independent, so assembly parallelises on CPU or GPU.
- Constraint elimination by null-space basis, $x=x_p+Zy$, giving the positive-definite system $(Z^\top KZ)\,y=Z^\top(f-Kx_p)$.
- Linear solve by sparse Cholesky; for interactive dragging, the factor is reused with low-rank updates as constraint values change.
- Nonlinear polish by Gauss–Newton with Levenberg–Marquardt damping on the geometric functionals.
- Inequalities (crown, curvature bounds) by an active-set quadratic programme around the same factorisation.
- Lagrange multipliers $\lambda$ are reported to the user: a large multiplier identifies the constraint that is fighting the shape.

---

## 5. System blueprint and roadmap

### 5.1 Modules

```
┌──────────────────────────────────────────────────────────────────────────┐
│ ak::api      C++20 API, C ABI, journal/replay, STEP AP242 / IGES I/O     │
├───────────────┬──────────────────┬──────────────────┬────────────────────┤
│ ak::ops       │ ak::var          │ ak::heal         │ ak::viz            │
│ boolean,      │ energies,        │ persistence,     │ watertight         │
│ offset, blend,│ G^k constraints, │ re-intersection, │ tessellation, GPU  │
│ sweep, loft   │ QP / GN solvers  │ snapping         │ eval and analysis  │
├───────────────┴──────────────────┴──────────────────┴────────────────────┤
│ ak::topo     cell complex, Euler operators, persistent store, audit      │
├──────────────────────────────────────────────────────────────────────────┤
│ ak::solve    subdivision / Bézier clipping, Krawczyk, Kantorovich, SSI   │
├──────────────────────────────────────────────────────────────────────────┤
│ ak::geom     carriers, evaluators, provenance DAG, tubes                 │
├──────────────────────────────────────────────────────────────────────────┤
│ ak::num      double+bound, interval, expansions, ball MP, rational,      │
│              algebraic, lazy expression DAG, escalation policy           │
├──────────────────────────────────────────────────────────────────────────┤
│ ak::diag     telemetry, indicators, rate regression, learned scheduler   │
└──────────────────────────────────────────────────────────────────────────┘
```

| Subsystem | Responsibility | Key design points |
|---|---|---|
| **Topology Manager** (`ak::topo`) | Cell complex, incidence, transactions | Non-manifold; persistent with structural sharing (lock-free readers, cheap undo); Euler operators as the only mutators; integer audit after each commit |
| **Exact/Interval Evaluator** (`ak::num`, `ak::solve`) | Number types, predicates, root certification | Evaluators are templates over the number type, so one NURBS routine serves double, interval, jet and ball arithmetic; predicates return `Sign` or `Undecided` |
| **NURBS/Variational Solver** (`ak::geom`, `ak::var`) | Carriers and optimisation | Carriers immutable and hash-consed; energies and constraints assembled per Bézier element; solvers as in §4.4 |
| **Visualisation Tessellator** (`ak::viz`) | Display meshes, analysis shading | Each edge is sampled once from its tube and shared by both faces, so meshes are crack-free by construction; sampling density from the Filip bound |
| **Diagnostics** (`ak::diag`) | §2.2 loop | Advisory only; can be disabled without changing any result |

### 5.2 GPU compute

| Workload | GPU role | Constraint |
|---|---|---|
| Batch evaluation of points, normals, curvature | compute shader over Bézier elements | results are filters |
| Broad-phase overlap | BVH traversal | conservative boxes, widened by a fixed ulp margin |
| Energy and constraint assembly | per-element quadrature | double precision where the API offers it |
| Linear solve | preconditioned conjugate gradient on the null-space system | CPU Cholesky remains the reference |
| Display | hardware tessellation or ray casting of trimmed patches; zebra, isophote, curvature-comb shaders | uses the shared edge samples |

Rule: **the GPU proposes, the CPU certifies.** GPUs give no directed rounding, drivers may fuse or reorder operations, and some shading languages have no 64-bit float (WGSL, for example). No topological decision depends on an uncertified GPU value. The kernel builds and passes its tests with the GPU disabled.

### 5.3 Roadmap

Commercial kernels embody decades of work by large teams. The plan narrows scope first and widens it only at gates. Durations assume a team of 6–10 and are rough estimates.

| Phase | Months | Deliverable | Exit gate |
|---|---|---|---|
| 0 Numerics | 0–6 | `ak::num`, templated NURBS evaluation, filtered predicates | Predicates agree with a multi-precision oracle on $10^9$ random and adversarial cases |
| 1 Certified intersection | 6–15 | `ak::solve`, edge varieties, tubes | All branches found on a curated SSI suite, each with a certificate |
| 2 Analytic solids | 12–24 | `ak::topo`, Booleans on planes/quadrics/tori, tessellation, STEP read | Chained-Boolean stress test with zero ledger entries |
| 3 Freeform Booleans | 20–34 | NURBS and SubD carriers in Booleans, provenance coincidence | §5.4 watertightness and identity targets on the freeform corpus |
| 4 Variational and Class-A | 24–40 | `ak::var`, G2 then G3 fills, fillets, crown | Certified jump bounds meet the Class-A specification |
| 5 Offsets, blends, healing | 34–50 | procedural offsets, shelling, `ak::heal` | Healing success rate on the import corpus |
| 6 Hardening | 46–60 | learned scheduler, GPU paths, API freeze | Determinism across three platforms |

### 5.4 Validation benchmark suite

**Corpora.** ABC dataset (about one million B-Rep models); Fusion 360 Gallery; NIST and CAx-IF STEP test cases for exchange; a synthetic adversarial set (near-tangent cylinders, rotate-and-union by small angles, nested thin shells, coincident offsets); a Class-A set of curve networks with known analytic solutions (sphere, torus and Monge-patch fills where the exact surface is known).

**Baselines.** OCCT (open source, directly comparable). CGAL exact Nef polyhedra as an exactness oracle for planar models. Commercial kernels where licence terms allow publishing results; this must be checked per licence.

| Axis | Metric | Definition | Target (hypothesis) |
|---|---|---|---|
| Watertightness | Edge gap | $\sup_t\lVert S_1(c_1(t))-S_2(c_2(t))\rVert$ over all edges | bounded by evaluation ε, user-selectable; no stored tolerance |
| | Closure rate | share of outputs with $\partial_2\mathbf 1=0$ and Euler–Poincaré satisfied | 100 % of committed results |
| | Tolerance growth | $\max\tau$ after 1000 chained random Booleans | 0 on transversal inputs |
| Correctness | Inclusion–exclusion | $\lvert\mathrm{vol}(A\cup B)+\mathrm{vol}(A\cap B)-\mathrm{vol}A-\mathrm{vol}B\rvert$ | within certified quadrature bound |
| | Algebraic identities | $A\setminus(A\setminus B)=A\cap B$; commutativity; associativity | identical topology |
| | Honest failure | share of failures that are explicit `Indeterminate` | 100 % (no silent invalid bodies) |
| Determinism | Reproducibility | bitwise equality across thread counts and platforms | 100 % |
| Speed | Boolean time | median and 99th percentile against face count | median within 2–3× of OCCT; better 99th percentile |
| | Escalation rate | share of predicate calls leaving tier 0 and tier 1 | < 5 % and < 0.1 % |
| Memory | Footprint | bytes per face including tubes and provenance | within 2× of a p-curve B-Rep |
| Smoothness | G1 | $\sup\angle(n_1,n_2)$ on seams | certified ≤ specification (order of 0.01°–0.05°) |
| | G2 | $\sup\lvert[\![b(\nu,\nu)]\!]\rvert/\max\lvert\kappa\rvert$ | certified ≤ specification |
| | G3 | $\sup\lvert[\![C(\nu,\nu,\nu)]\!]\rvert$, normalised | certified ≤ specification |
| | Fairness | $E_3$, $E_4$; number of curvature-profile extrema; highlight-line inflections | no extrema beyond those forced by the data |
| | Accuracy | deviation from the analytic solution on the Class-A set | converges at the theoretical order |
| Healing | Success | share of import-corpus models closed with zero ledger entries; maximum deviation | reported, no fixed target yet |

The speed and memory targets express an expected trade: the certified kernel is likely slower in the median case and more predictable in the tail. Both expectations need measurement.

### 5.5 Open problems and risks

1. **Tangential intersections without provenance** (imported data, user-drawn tangencies). No complete, fast decision procedure exists. The design contains the problem, it does not solve it.
2. **A $G^3$ cap at extraordinary vertices** with good shape and moderate degree is not an established construction.
3. **Cost of tubes and provenance** on models with $10^5$–$10^6$ faces is unknown until Phase 3.
4. **Coincidence of differently parameterised freeform surfaces** is decidable exactly only in special cases (common refinement, sub-patch).
5. **The Codazzi reduction for G3** in §1.7 is confirmed numerically but has not been checked against the literature.
6. **Interoperability.** Other systems expect p-curves and tolerances. Export must generate them from tubes at a stated accuracy, and round-trip tests must confirm that re-import heals back to zero ledger entries.

---

## References

Cited from memory; verify bibliographic details before publication.

- Busé, L. Implicit matrix representations of rational Bézier curves and surfaces. *Computer-Aided Design*, 2014.
- Celniker, G., Gossard, D. Deformable curve and surface finite-elements for free-form shape design. *SIGGRAPH*, 1991.
- Farouki, R. *Pythagorean-Hodograph Curves*. Springer, 2008.
- Filip, D., Magedson, R., Markot, R. Surface algorithms using bounds on derivatives. *CAGD*, 1986.
- Hoffmann, C. *Geometric and Solid Modeling*. Morgan Kaufmann, 1989.
- Hughes, T., Cottrell, J., Bazilevs, Y. Isogeometric analysis. *CMAME*, 2005.
- Kiendl, J., Bletzinger, K.-U., Linhard, J., Wüchner, R. Isogeometric shell analysis with Kirchhoff–Love elements. *CMAME*, 2009.
- Koch, S. et al. ABC: a big CAD model dataset for geometric deep learning. *CVPR*, 2019.
- Loop, C., Schaefer, S. G2 tensor product splines over extraordinary vertices. *Computer Graphics Forum*, 2008.
- Micchelli, C., Rivlin, T. A survey of optimal recovery. In *Optimal Estimation in Approximation Theory*, 1977.
- Moreton, H., Séquin, C. Functional optimization for fair surface design. *SIGGRAPH*, 1992.
- Nairn, D., Peters, J., Lutterkort, D. Sharp, quantitative bounds on the distance between a polynomial piece and its Bézier control polygon. *CAGD*, 1999.
- Patrikalakis, N., Maekawa, T. *Shape Interrogation for Computer Aided Design and Manufacturing*. Springer, 2002.
- Pegna, J., Wolter, F.-E. Geometrical criteria to guarantee curvature continuity of blend surfaces. *ASME J. Mech. Design*, 1992.
- Peters, J. Geometric continuity. In *Handbook of Computer Aided Geometric Design*, 2002.
- Requicha, A. Representations for rigid solids. *ACM Computing Surveys*, 1980.
- Sederberg, T., Meyers, R. Loop detection in surface patch intersections. *CAGD*, 1988.
- Shewchuk, J. Adaptive precision floating-point arithmetic and fast robust geometric predicates. *Discrete & Computational Geometry*, 1997.
- Smale, S. Newton's method estimates from data at one point. 1986.
- Urick, B., Marussig, B., Cohen, E., Crawford, R., Hughes, T., Riesenfeld, R. Watertight Boolean operations. *Computer-Aided Design*, 2019.
- Várady, T., Salvi, P., Karikó, G. A multi-sided Bézier patch with a simple control structure. *Computer Graphics Forum*, 2016.
- Welch, W., Witkin, A. Variational surface modeling. *SIGGRAPH*, 1992.
- Yap, C. Towards exact geometric computation. *Computational Geometry*, 1997.
