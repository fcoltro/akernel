# Random models

`ak_random_models` (fuzz/random_models.cpp) is the kernel's own test corpus.
aKernel reads no foreign formats, so public CAD collections are out of reach;
instead every test model is a short feature chain grown from one 64-bit seed,
and every step is checked by oracles that do not trust the operation.

## What a case is

A random primitive, then 1 to `--depth` (default 3) Boolean operations, each
with another random primitive placed to overlap the current result. The
operation kept for the next step is a union (40%), an intersection (20%) or a
difference (40%). The primitives:

| kind | built by | why it is in the mix |
|---|---|---|
| grid box | `make_box`, corners on the 1/4 grid | faces of different operands in one plane, shared edges and corners |
| rotated box | `make_box`, random orientation | generic position |
| hill box | `make_bump_box`, on the grid | a curved polynomial face |
| cylinder sector | `make_cylinder_sector`, dyadic parameters | rational faces, exact cylinders shared between operands |
| edge fillet | `fillet_box_edge`, on the grid or rotated | a rolling-ball blend |
| edge chamfer | `chamfer_box_edge`, on the grid or rotated | a flat blend at an angle |
| rounded box | `rounded_box`, all twelve edges | cylinders and sphere corners meeting |
| ring blend | `blend_ring_edge`, 1 to 4 quarter turns | torus zones and cones |

## Oracles

| oracle | applies to | a failure means |
|---|---|---|
| audit | every built primitive and every Boolean result | the boundary is not closed and oriented (`op::audit`, and `op::check_edge_fans` where parts touch) |
| self-intersection | every primitive and every Boolean result | `op::check_self_intersection` proved two faces cross. *Undecided* is counted, not failed |
| closed form | every primitive | the quadrature volume is off its formula by more than 1e-8 (relative), or the certified enclosure misses it |
| identities | every Boolean step | vol(A∪B) + vol(A∩B) ≠ vol A + vol B, or vol(A−B) ≠ vol A − vol(A∩B), by more than 1e-7 (relative); when all five volumes are certified, the interval of each defect must contain 0 |
| membership | the result kept at each step | a Monte Carlo count of points inside the analytic chain (no kernel code) is more than 6 standard errors from the Brep's volume on a second, ten times larger draw |

A *refusal* (the kernel answers "not solved" for an imprint, or a blend
constructor rejects its parameters) is not a failure: the kernel does not
guess. Refusals are counted per operation and grouped by reason, numbers
removed, so the reasons can be driven to zero one by one. The exit code is 1
when an oracle fails or an exception escapes, 0 otherwise.

## Running it

```
ak_random_models                         # seeds 1 to 50
ak_random_models --seed 5000 --count 400 --markdown out.md --json out.json
ak_random_models --seed 1234 --count 1 --verbose   # replay one case, step by step
ak_random_models --seeds-file fuzz/regression_seeds.txt
```

Every case is the same on every platform: the random source is SplitMix64
(tests/check.hpp) and the angles go through num/detmath.hpp, as the benchmark
hashes do. A seed that fails on the nightly Linux run fails the same way on
Windows.

`fuzz/nightly.sh EXE OUTDIR FIRST COUNT SECONDS` runs a range of seeds and
survives crashes: the program writes the seed it is on to a progress file, and
after a crash the script records that seed in `OUTDIR/crashes.txt` and goes on
from the next one.

## In CI

- `ctest` runs `random_models_smoke` (seeds 1 to 12) and
  `random_models_regressions` (fuzz/regression_seeds.txt) with every other test.
- `.github/workflows/nightly.yml` runs every night at 03:23 UTC on Windows
  (MSVC) and Linux (GCC): the regression seeds, seeds 1 to 300 (the same every
  night, so the counts compare from one night to the next), and up to 600 new
  seeds starting at `YYYYMMDD0000`. The tables are on the run's summary page and
  the JSON and Markdown reports are kept as artifacts. It can also be started by
  hand with another first seed, count or time budget.

When a failing seed's bug is fixed, add the seed to `fuzz/regression_seeds.txt`
with a comment naming the fix.
