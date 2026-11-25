# Helmer

**Whitney stratifications of algebraic varieties and maps, in C++ with Python bindings.**

Helmer is a reimplementation of Martin Helmer's Macaulay2 package
[`WhitneyStratifications`](http://martin-helmer.com/Software/WhitStrat/WhitneyStratifications.m2)
(version 2.23, June 2025), named in his memory. It needs no Macaulay2: Gröbner bases
come from [msolve](https://msolve.lip6.fr)'s F4 library, polynomial arithmetic and
factorisation from [FLINT](https://flintlib.org), and everything Macaulay2's
top-level packages used to provide (dimension, elimination, saturation, minimal
primes, conormal varieties, …) is implemented here. On the benchmarks of Helmer and
Mohr it is 38× to 5000× faster than the Macaulay2 package, with the same strata.

```python
import helmer

R  = helmer.Ring("x, y, z")
X  = helmer.Ideal(R, ["y^2*z - x^2"])        # Whitney umbrella
WS = helmer.whitney_stratify(X)
# {0: [ideal(z, y, x)], 1: [ideal(y, x)], 2: [ideal(y^2*z - x^2)]}
```

`WS[i]` lists the closures of the `i`-dimensional strata as prime ideals over ℚ; the
open strata are the connected components of `V(P)` minus the lower-dimensional
closures, exactly as in the Macaulay2 package. Over ℝ the same polynomials describe a
valid stratification (some strata may be empty or drop dimension), see [HLN].

The same from the command line, and from C++:

```
$ helmer "y^2*z - x^2"
Whitney stratification in C^3 with coordinates x, y, z (0.00 s)
  dim 2:  y^2*z - x^2
  dim 1:  y, x
  dim 0:  z, y, x
$ helmer --projective "x*y*z + x*y*w + x*z*w + y*z*w" --vars x,y,z,w     # the Cayley cubic: four nodes
```

```cpp
auto R = helmer::Ring::make({"x", "y", "z"});
helmer::Ideal X(R, {helmer::Poly(R, "y^2*z - x^2")});
helmer::Stratification S = helmer::whitney_stratify(X);   // S.strata[d]: prime ideals
```

## Install

```
pip install helmer-strat
```

Binary wheels for Linux (x86-64) and macOS 15 or later (Apple silicon and Intel), CPython 3.9–3.13, bundle
everything Helmer needs (GMP, MPFR, FLINT, msolve). The package is imported as `helmer` and installs the
`helmer` command.

Elsewhere (Windows, Linux on ARM), or to build from a clone of this repository, you need FLINT ≥ 3.0 with
GMP, msolve ≥ 0.6 (its `libneogb` library), a C++20 compiler and CMake ≥ 3.18. On Ubuntu 24.04:

```
sudo apt install libflint-dev libmsolve-dev cmake g++
pip install .                         # in a clone; or pip install helmer-strat --no-binary helmer-strat
```

On other systems, install FLINT from your package manager and build msolve from
[its source](https://github.com/algebraic-solving/msolve); `tools/build_wheel_deps.sh` does both.

## Getting started

[`examples/`](examples) has six short Python scripts, each runnable as it is, and a C++ program:

| | shows |
|---|---|
| [01_first_stratification.py](examples/01_first_stratification.py) | rings, ideals, reading a stratification; curves, surfaces, several equations, projective varieties |
| [02_maps.py](examples/02_maps.py) | `map_stratify`: the roots of a x² + b x + c as (a, b, c) varies; where a map is not proper |
| [03_algebra.py](examples/03_algebra.py) | Gröbner bases, dimension, minimal primes, elimination, saturation, conormal varieties |
| [04_certified.py](examples/04_certified.py) | certified elimination: polynomial identities that prove where the critical values of a function lie |
| [05_linkages.py](examples/05_linkages.py) | singular configurations of four-bar linkages; free translation directions |
| [06_conormal_or_polar.py](examples/06_conormal_or_polar.py) | the two algorithms, side by side |
| [cpp/stratify.cpp](examples/cpp/stratify.cpp) | the C++ API (`cmake -DHELMER_EXAMPLES=ON`) |

The functions are documented in [docs/API.md](docs/API.md).

## What is implemented

| Macaulay2 (`WhitneyStratifications`)   | Helmer                                   | Notes |
|---|---|---|
| `whitneyStratify(I, AssocPrimes=>false)` | `whitney_stratify(I, projective=False)` | Conormal algorithm of [HN] in the associated-prime-free form of [HM]. |
| `mapStratify(F, X, Y, AssocPrimes=>false, isProper=>…)` | `map_stratify(F, X, Y, proper=True)` | Stratifies source and target so Thom's isotopy lemma applies. |
| `nonProperSet(F, X, Y)` | `non_proper_set(F, X, Y)` | Jelonek set, via the closure of the graph in ℙⁿ × kᵐ. |
| `whitneyStratifyPol(I)` | `whitney_stratify(I, method="polar")` | Polar-variety route: randomised, non-minimal output. |
| `polarVars(I)` | `polar_varieties(I)` | |
| | `whitney_filtration(I)` | The closed dimension filtration, deferring prime decomposition; see below. |
| `minCoarsenWS(V)` | `min_coarsen(S)`, or `whitney_stratify(I, minimal=True)` | The minimal Whitney stratification; see below. |
| `polarSequence(X, Z)` | `polar_sequence(X, Z)`, `polar_multiplicity(X, Z, w)` | |
| `multiplicity(X, Y)` (`SegreClasses`) | `multiplicity(X, Y)` | Computed locally, without Segre classes. |
| `containedInSingularLocus(X, Y)` (`SegreClasses`) | `in_singular_locus(X, Y)` | |
| | `contained(X, Y)` | V(X) ⊆ V(Y); for ideals, `X.contains(Y)`. |
| `conormal(I)` | `conormal(I)` | In kⁿ × ℙⁿ⁻¹, variables `v_0..v_{n-1}`. |
| `decompose(I)` | `minimal_primes(I)` | New implementation; see below. |
| `eliminate`, `saturate`, `intersect`, `gb`, `dim`, `codim`, `independentSets` | same names / `Ideal` methods | |

Not yet ported: `mapStratifyPol`, the Euler-obstruction / characteristic-cycle functions
(`eulerObsMatrix`, `multCC`) and Segre classes; these come next, together with
Chern–Schwartz–MacPherson classes of strata. `AssocPrimes=>true` (the original [HN] route
via associated primes) is deliberately not ported: [HM] shows it is slower, and minimality
is better obtained by minimisation.

## How it works

* **Gröbner bases over ℚ.** msolve's F4 modulo 31-bit primes, Chinese remaindering
  and rational reconstruction; a reconstruction is accepted only after a
  fresh prime confirms it and every input generator reduces to zero over ℚ.
  Each call samples distinct random primes and groups results by leading monomials. Only two monomial orders are ever
  needed (DRL, and a two-block elimination order), both native to msolve.
* **Minimal primes.** Factorisation splitting, then a GTZ-style reduction: for a
  maximal independent set *u*, a block basis gives `h` with
  `V(I) = V(I : h^∞) ∪ V(I + h)`; the first part is the contraction of a
  zero-dimensional ideal over k(u), whose primes come from factoring the minimal
  polynomial of a linear form (FLINT multivariate factorisation over ℚ). Each candidate
  is **certified prime** by checking `dim_{k(u)} k(u)[x∖u]/P = deg μᵢ`; if a candidate
  fails, the ideal is made radical (Seidenberg) and/or the linear form is re-drawn.
* **Conormal varieties: hyperplane charts.** For a prime `I`, the conormal
  computation first looks for an affine linear form ℓ vanishing on every component of
  the singular locus but not on `V(I)`. On the dense open set ℓ ≠ 0, `V(I)` is smooth
  and the conormal equations `I + K` (K the (c+1)-minors of the Jacobian bordered by the
  dual coordinates) already define the conormal bundle, so `I(Con X) = (I + K) : ℓ^∞`.
  Helmer computes that saturation in one Gröbner basis: in coordinates where ℓ is a
  variable, homogenise and use graded reverse lexicographic order with ℓ last
  (Bayer's trick). It prefers a coordinate (in reverse input order), then a sparse
  linear form; `saturation="minor"` skips the chart.
* **Conormal varieties: one minor.** When no chart exists (the singular locus spans
  the affine hull), Helmer saturates by a *single* sparse *c*-minor `g ∉ I` instead of
  the whole minor ideal. Over the dense open set `g ≠ 0` the conormal is smooth and
  reduced, so `(I + K) : g^∞` is exactly `I(Con X)`. Saturations use the squarefree part
  of `g`. Charts and single minors are deterministic, and are the main reason for the
  speed-up on larger examples.
* **Reduced varieties.** `whitney_stratify` stratifies the reduced variety even when
  its input ideal is not radical: it uses the minimal primes for the component strata
  and their intersection for the singular locus of an equidimensional union.
* **Free directions.** On affine inputs Helmer first looks for constant directions v
  with ∂_v f = 0 for every generator f, so that `V(I)` is a product with an affine
  space (equations in the differences of joint coordinates of a mechanism, say). It
  stratifies a complementary linear slice, pulls the strata back and adds the number of
  free directions to their dimensions.

### What is probabilistic

Stratification has one Monte-Carlo step, the modular lift of Gröbner bases over ℚ (as in msolve
itself and Macaulay2's modular methods): a wrong answer requires a reconstruction that
is confirmed by a fresh 31-bit prime and passes the generator-reduction check.
The latter establishes that the reconstructed ideal contains the input ideal. The prime certificate in
`minimal_primes` and the conormal saturation are exact. Ideal membership and
elimination can be certified (next section); verifying every lifted Gröbner basis over
ℚ, which would certify a whole stratification, is not implemented yet.
Multiplicities, polar sequences and minimisation (next section) are Monte Carlo by
nature, being values at a *general* point; each is repeated with independent primes and
random choices until two runs agree.

## Filtrations without prime decomposition

`whitney_filtration(I)` returns the Whitney stratification as its closed dimension
filtration, `{d: J_d}` with V(J_d) the union of the stratum closures of dimension ≤ d,
and avoids prime decomposition wherever it can. Mechanisms are where this matters:
their configuration spaces are often reduced complete intersections whose prime
decomposition is the expensive step.

* **Reduced complete intersections stay whole.** If the number of equations equals the
  height, V(I) is a complete intersection, hence unmixed (S1); if moreover the Jacobian
  locus has smaller dimension, it is generically smooth on every component (R0), so I is
  reduced by Serre's criterion. Such an I is processed as one piece.
* **Finite loci are made reduced without decomposing.** Adding the squarefree part of
  each coordinate's eliminant gives the radical of a zero-dimensional ideal
  (Seidenberg's lemma).
* **Symmetries bound Whitney failures along curves.** For a lower curve Z, Helmer finds
  the affine-linear vector fields D with D(I_X) ⊆ I_X and D(I_Z) ⊆ I_Z. Their flows
  preserve the pair, so Whitney regularity at a general point of Z holds along every orbit
  on which some D is nonzero: if the common zeros of the fields on Z are finite, they are
  the only candidate failures, and no conormal variety is needed. (A family whose singular
  axis carries a ℂ*-action, like a weighted-homogeneous family with a shifted parameter,
  is the typical case.)

Elsewhere it falls back to the prime-by-prime route; `defer=False` uses it throughout.
The filtration is a valid Whitney filtration but need not be minimal, and it carries no
prime decomposition, so `min_coarsen` (which needs prime closures) does not apply to it
directly. On the test inputs every stratum of the minimal Whitney stratification lies in
the filtration (`tests/test_filtration.py`). The gain is largest where prime decomposition
was the bottleneck: the pinned five-bar and the rotating four-bar go from over two minutes
to milliseconds, and so does a Briançon–Speder-type quartic family with a shifted parameter
(eager route over two minutes; polar route then `min_coarsen`, 14 s; filtration, 4 ms, equal
to the minimal stratification).

## Multiplicities and minimal stratifications

`multiplicity(X, Y)` is the Samuel multiplicity of V(Y) along the irreducible V(X), which
equals the multiplicity of Y at a general point x of X. Helmer computes it at one such
point, over F_p for a random 31-bit prime p:

* **A general point.** Slice X by dim X random affine hyperplanes, take an F_p-root of
  the eliminant in one slice coordinate (FLINT), and solve for the others. Strata of
  conjugate points (an irreducible zero-dimensional prime of degree > 1) have F_p-points
  only for some primes; those primes are skipped.
* **A curve germ.** Cut Y by w − 1 general hyperplanes through x (w = dim_x Y), which keeps
  the multiplicity at x, and saturate by one more linear form through x. That removes any
  embedded component at x, so what is left near x is a Cohen–Macaulay curve, whose
  multiplicity is the length at x of its section by a further general hyperplane through
  x. Non-Cohen–Macaulay Y are therefore fine: for two planes meeting in a point of ℂ⁴ the
  answer is 2, where the length of a general linear section through the point is 3.

Every system solved has length at most deg Y. Macaulay2's `SegreClasses` instead counts
d^(dim Y − dim X) · deg Y points globally (d the degree of the equations of X) and
subtracts.

`polar_sequence(X, Z)` lists the multiplicities at a general point of X of the generic
polar varieties of Z, by dimension, ending with e_X(Z). Each is computed the same way, on
the slice: the polar locus is cut by the maximal minors of a random matrix bordered by
combinations of the Jacobian rows, saturated by one combination of c-minors (which turns the minor
equations into Lagrange-type equations, much cheaper for F4 than any low-degree equation
of the singular locus). Entries of dimension at most dim X vanish and are not computed: the
polar variety is too small to contain X, or equals the closure of a subset of Z_reg; so do
all entries below the top when e_X(Z) = 1.

`min_coarsen(S)` (and `whitney_stratify(I, minimal=True)`) returns the coarsest Whitney
stratification refined by S, using Teissier's criterion. Strata are processed from the
top down. A stratum X merges into the open stratum of Y when Y is the unique smallest
remaining stratum whose closure contains X, every remaining stratum whose closure contains
X also contains Y, Y-bar is smooth along X (e_X(Y-bar) = 1), and every larger such stratum
Z has the same polar multiplicities along X as along Y. Polar sequences are compared entry
by entry and the comparison stops at the first difference. On the benchmarks coarsening
costs about as much as stratifying: 0.07 s for [HM] Example 1.1 (10 strata to 7), 0.2 s
for X₂ (4 to 3), 0.9 s for `briancon_speder_small` and 19 s for Briançon–Speder, which
keeps its origin: there the polar curve of the hypersurface has multiplicity 2 against 0
along the t-axis, and the polar surface 32 against 30.

## Certified membership and elimination

`certify_membership(I, h, block)` returns cofactors with h = Σ c_j g_j over the
generators g_j of I, and `eliminate_certified(I, vars)` returns `eliminate(I, vars)`
together with such a certificate in I for each generator. A certificate is *found*
heuristically: the cofactors' coefficients are the unknowns of a linear system over a
bounded support (degree ≤ D_A in the `block` variables, ≤ D_B in the others, for a
small range of bounds ordered by system size). One elimination modulo a ~26-bit prime
picks a nonsingular square subsystem (pivot columns and independent rows), which is
solved by p-adic lifting from a single LU factorisation, with rational reconstruction
attempted as the precision grows. It is *returned* only after the
identity has been checked exactly in FLINT's `fmpq_mpoly`, so it proves membership
whatever the reliability of the modular steps or of the Gröbner basis; `None` proves
nothing. `verify_membership(I, h, cofactors)` repeats the exact check against the
generators of the ideal I (a list of generators and a ring also works).
`max_unknowns` (default 4,000) caps the linear systems tried: the first elimination is
dense, so memory grows like rows × unknowns and time like unknowns³.

For an elimination generator P this proves P ∈ I, hence P vanishes on V(I) over ℂ,
which is the direction needed to bound critical values: with
I = (∂f/∂x_1, …, ∂f/∂x_n, f − ε), every critical value of f is a root of the
certified P(ε), including when the critical points are not isolated. Only a nonzero
certified generator gives a bound; if the elimination ideal comes back zero, both
returned lists are empty.

```python
R = helmer.Ring("x, y, epsilon")
f = "(x^2+y^2-1)*(x^2+y^2-4)"            # critical points form a circle
I = helmer.Ideal(R, ["2*x*(2*x^2+2*y^2-5)", "2*y*(2*x^2+2*y^2-5)", f + "-epsilon"])
J, (cert,) = helmer.eliminate_certified(I, ["x", "y"])
# J = (epsilon^2 - 7/4*epsilon - 9): critical values 4 (origin) and -9/4 (a circle)
```

Measured on critical-value ideals (single equation f, or F = Σ f_i² for systems):
the certificates found have block degree deg f − 1 to deg f + 1 and degree
deg P − 1 in ε, so the linear systems stay small.

Timings of `eliminate_certified` (beyond `eliminate`), from `bench/certify_bench.py`;
"before" is the CRT version (a full dense elimination for every prime).

| input | deg P | unknowns | before (s) | after (s) |
|---|---|---|---|---|
| torus in ℝ³ (deg 4) | 3 | 240 | 0.07 | 0.02 |
| two spheres meeting in a circle (deg 4) | 3 | 240 | 0.08 | 0.02 |
| Clifford torus (F, deg 4, 4 variables) | 3 | 1,890 | 17.2 | 3.6 |
| 5 concentric circles (deg 10, critical circles) | 5 | 990 | 1.5 | 0.28 |
| 5 separated circles (deg 10) | 10 | 1,980 | 22.7 | 3.5 |
| 6 separated circles (deg 12) | 12 | 2,808 | 122.8 | 16.8 |
| 8 concentric circles (deg 16, critical circles) | 8 | 3,264 | 144.0 | 15.1 |

Almost all of the remaining time is the one dense elimination per bound tried
(including the bounds that turn out too small); the lifting itself takes 1–2 s on the
largest cases. A sparse first elimination, or skipping bounds that never succeed, is
the next step.

## Conormal or polar?

In Macaulay2 the polar route is far faster than the conormal route (X₁ below: 2.7 s vs
210 s), because M2 saturates the conormal ideal by the entire ideal of minors. Helmer
uses a hyperplane chart or one sparse minor (exact for prime input, see above), which removes
that bottleneck. As a result the conormal route is the default: on every benchmark so far it
is faster than the polar route, and it is deterministic.

| Example | M2 conormal | M2 polar | Helmer conormal | Helmer polar |
|---|---|---|---|---|
| [HM] X₁ | 210 s | 2.7 s | 0.04 s | 1.2 s |
| [HM] X₂ | 7.2 s | 208 s | 0.19 s | 16.0 s |
| [HM] Example 1.1 | 14.1 s | 22.0 s | 0.18 s | 8.6 s |

Helmer times measured 25 Nov 2026 (best of three, `set_seed(1)` for the polar route); Macaulay2
times as in the table below.

The polar route draws random linear projections, so it returns valid but generally
non-minimal stratifications with random extra strata (exactly as in M2). Neither
route is minimal in general: on Example 1.1 the conormal route returns three
non-minimal curves, and on X₂ both conormal and the associated-primes route return a
non-minimal point. The tests check every polar output against the *minimal*
stratification computed by M2's `minCoarsenWS`. `polar_bound` (default 1000) sets the
size of the random coefficients: smaller is faster but raises the small probability of
an unlucky projection.

The Briançon–Speder family `x⁵ + y⁷z + z¹⁵ + t·xy⁶` (μ-constant, but not Whitney
equisingular along the t-axis) uses the chart `z ≠ 0`. Helmer finds the point at the
origin where Whitney's condition fails, in 35 s; its degree-nine analogue
`briancon_speder_small` takes 0.3 s (`python bench/run_helmer.py briancon_speder`).

## Validation and speed

`tests/test_helmer.py` checks Helmer against Macaulay2 1.22 + WhitneyStratifications
2.23 (with `AssocPrimes=>false`) on the package's documentation examples, the
benchmark varieties of [HM], and assorted singular surfaces, curves and projective
hypersurfaces; plus both root-classification examples for `mapStratify` and a battery
of `decompose` comparisons. All agree: same strata, same dimensions. Polar outputs
are checked to refine the minimal stratification. Further tests cover nonreduced input,
the hyperplane charts against the single-minor saturation, free directions (against the
direct computation), and certified membership and elimination. `tests/test_local.py`
checks multiplicities and polar sequences against closed forms (cones over plane curves of
degree d: polar sequence (0, d(d−1), d), Euler obstruction 2d − d²), `min_coarsen` against
M2's `minCoarsenWS` where recorded, the polar route against the conormal route (the polar
route adds random extra strata; after coarsening the two agree on every affine example),
and projective varieties against their affine cones. `tests/test_filtration.py` checks that `whitney_filtration` contains the minimal stratification on every example (156 tests, about 40 s).
Independently of the fixtures, on every example of `bench/examples.py` the chart
saturation agrees with the exact saturation by all minors, and adding a free variable
(with a random rational shear) shifts every stratum up by one dimension.

Wall-clock times on the same 2-core machine (M2 `cpuTime`; Helmer via Python, best of
three, 25 Nov 2026):

| Example | Macaulay2 | Helmer |
|---|---|---|
| Whitney cusp | 0.09 s | < 0.01 s |
| Projective quartic (`quartic_P`) | 0.39 s | 0.01 s |
| [HM] Example 1.1 (hypersurface in ℂ⁴) | 14.1 s | 0.18 s |
| [HM] X₁ (complete intersection in ℂ⁴) | 210 s | 0.04 s |
| [HM] X₂ | 7.2 s | 0.19 s |
| `mapStratify` Feynman example (package docs) | 1.38 s | 0.07 s |

## Differences from the Macaulay2 package

* Components of a non-equidimensional `X` are filed under their own dimension (the M2
  package files every component under `dim X`), and all pairs of components are
  intersected up front.
* In `mapStratify`, the M2 code pulls Whitney points of `Y` back to `X` with the
  *forward* elimination (`elimIm`), which is only correct for coordinate projections;
  Helmer uses the preimage `f⁻¹(K) ∩ X` throughout.
* `minCoarsenWS` (2.23) differs from Teissier's criterion as implemented here in three
  places, which can make it return a stratification that is not Whitney:
  (a) for affine input the test for Y-bar being singular along X computes the
  multiplicity of X along *itself*, which is always 1, so the test never fires;
  (b) the polar test loops over the strata Z containing X but reads them as `Zlist_b`
  instead of `Zlist_(contY_b)`, i.e. by position in the whole level, so it can compare
  strata that do not contain X and skip those that do (umbrella ∪ the plane x = 5 loses
  the umbrella's origin if the plane comes first);
  (c) a stratum contained in no stratum closure of the next nonempty level up is dropped
  unchecked (cone ∪ the plane z = 5 loses the cone's vertex). Both examples are tests.
* The Gröbner prime stream is seeded from system randomness. For reproducible
  runs, `helmer.set_seed` seeds both the prime stream and the separate stream
  for geometric choices.

## Building and testing

With the dependencies of [Install](#install) (for a build from source):

```
pip install .                      # the Python module and the helmer command
pytest tests/                      # 85 tests, about 20 s
sh examples/run_all.sh             # every example
cmake -S . -B build -DHELMER_PYTHON=OFF -DHELMER_EXAMPLES=ON && cmake --build build   # C++ library and example
```

## Citing

If Helmer helps your research, please cite the software (GitHub's "Cite this repository"
button gives the entry from `CITATION.cff`; version 1.0.0 is
[doi:10.5281/zenodo.24531435](https://doi.org/10.5281/zenodo.24531435) on Zenodo) and the papers behind the algorithms: [HN] for the conormal algorithm, [HM] for the form
implemented here and for the polar route, and [HLN] for stratifications of real
varieties and maps.

## Licence

GPL-3.0-or-later. msolve is GPL-2.0-or-later and FLINT LGPL-3.0-or-later, so any
distributed binary of Helmer is GPL regardless.

## References

* [HN] M. Helmer, V. Nanda. *Conormal spaces and Whitney stratifications.* Found. Comput. Math. (2022), with correction (2023); arXiv:2106.14555.
* [HLN] M. Helmer, A. Leykin, V. Nanda. *Effective Whitney stratification of real algebraic varieties.* Math. Comp. (2025); arXiv:2307.05427.
* [HM] M. Helmer, R. Mohr. *Faster computation of Whitney stratifications and their minimization.* arXiv:2406.17122.
* [HT] M. Helmer, F. Tellander. *Spectral decomposition of Euler–Mellin integrals.* arXiv:2505.12458.
* msolve: J. Berthomieu, C. Eder, M. Safey El Din. *msolve: a library for solving polynomial systems.* ISSAC 2021.
* FLINT: The FLINT team. *FLINT: Fast Library for Number Theory.*
