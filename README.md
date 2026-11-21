# Helmer

**Whitney stratifications of algebraic varieties and maps, in C++ with Python bindings.**

Helmer is a reimplementation of Martin Helmer's Macaulay2 package
[`WhitneyStratifications`](http://martin-helmer.com/Software/WhitStrat/WhitneyStratifications.m2)
(version 2.23, June 2025), named in his memory. It needs no Macaulay2: Gröbner bases
come from [msolve](https://msolve.lip6.fr)'s F4 library, polynomial arithmetic and
factorisation from [FLINT](https://flintlib.org), and everything Macaulay2's
top-level packages used to provide (dimension, elimination, saturation, minimal
primes, conormal varieties, …) is implemented here.

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

## What is implemented

| Macaulay2 (`WhitneyStratifications`)   | Helmer                                   | Notes |
|---|---|---|
| `whitneyStratify(I, AssocPrimes=>false)` | `whitney_stratify(I, projective=False)` | Conormal algorithm of [HN] in the associated-prime-free form of [HM]. |
| `mapStratify(F, X, Y, AssocPrimes=>false, isProper=>…)` | `map_stratify(F, X, Y, proper=True)` | Stratifies source and target so Thom's isotopy lemma applies. |
| `nonProperSet(F, X, Y)` | `non_proper_set(F, X, Y)` | Jelonek set, via the closure of the graph in ℙⁿ × kᵐ. |
| `whitneyStratifyPol(I)` | `whitney_stratify(I, method="polar")` | Polar-variety route: randomised, non-minimal output. |
| `polarVars(I)` | `polar_varieties(I)` | |
| `conormal(I)` | `conormal(I)` | In kⁿ × ℙⁿ⁻¹, variables `v_0..v_{n-1}`. |
| `decompose(I)` | `minimal_primes(I)` | New implementation; see below. |
| `eliminate`, `saturate`, `intersect`, `gb`, `dim`, `codim`, `independentSets` | same names / `Ideal` methods | |

Not yet ported: `mapStratifyPol`, minimisation (`minCoarsenWS`), and the Euler-obstruction /
characteristic-cycle functions (`eulerObsMatrix`, `multCC`, `polarSequence`). These all
rest on Samuel multiplicities (Macaulay2's `SegreClasses`), which is the next piece of
algebra to build. `AssocPrimes=>true` (the original [HN] route via associated primes)
is deliberately not ported: [HM] shows it is slower, and minimality is better obtained
by minimisation.

## How it works

* **Gröbner bases over ℚ.** msolve's F4 modulo 31-bit primes, Chinese remaindering
  and rational reconstruction; a reconstruction is accepted only after an
  independent prime confirms it, and primes whose leading-monomial pattern differs
  from the majority are discarded as unlucky. Only two monomial orders are ever
  needed (DRL, and a two-block elimination order), both native to msolve.
* **Minimal primes.** Factorisation splitting, then a GTZ-style reduction: for a
  maximal independent set *u*, a block basis gives `h` with
  `V(I) = V(I : h^∞) ∪ V(I + h)`; the first part is the contraction of a
  zero-dimensional ideal over k(u), whose primes come from factoring the minimal
  polynomial of a linear form (FLINT multivariate factorisation over ℚ). Each candidate
  is **certified prime** by checking `dim_{k(u)} k(u)[x∖u]/P = deg μᵢ`; if a candidate
  fails, the ideal is made radical (Seidenberg) and/or the linear form is re-drawn.
* **Conormal varieties.** For a prime `I` of codimension *c*, Helmer saturates by a
  *single* sparse *c*-minor `g ∉ I` instead of the whole minor ideal. Over the dense
  open set `g ≠ 0` the conormal is smooth and reduced, so `(I + K) : g^∞` is exactly
  `I(Con X)`. This is deterministic and is the main reason for the speed-up on larger
  examples.

### What is probabilistic

The only Monte-Carlo step is the modular lift of Gröbner bases over ℚ (as in msolve
itself and Macaulay2's modular methods): a wrong answer requires a reconstruction that
is confirmed by an independent random-looking 31-bit prime. The prime certificate in
`minimal_primes` and the conormal saturation are exact. A fully certified mode
(verifying the lifted basis over ℚ) is on the roadmap.

## Certified membership and elimination

`certify_membership(I, h, block)` returns cofactors with h = Σ c_j g_j over the
generators g_j of I, and `eliminate_certified(I, vars)` returns `eliminate(I, vars)`
together with such a certificate in I for each generator. A certificate is *found*
heuristically: the cofactors' coefficients are the unknowns of a linear system over a
bounded support (degree ≤ D_A in the `block` variables, ≤ D_B in the others, for a
small range of bounds ordered by system size), solved modulo 62-bit primes with
Chinese remaindering and rational reconstruction. It is *returned* only after the
identity has been checked exactly in FLINT's `fmpq_mpoly`, so it proves membership
whatever the reliability of the modular steps or of the Gröbner basis; `None` proves
nothing. `verify_membership(I, h, cofactors)` repeats the exact check against the
generators of the ideal I (a list of generators and a ring also works).
`max_unknowns` (default 4,000) caps the linear systems tried: the solve is dense, so
memory grows like rows × unknowns and time like unknowns³.

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

| input | deg P | unknowns | seconds |
|---|---|---|---|
| torus in ℝ³ (deg 4) | 3 | 240 | 0.04 |
| two spheres meeting in a circle (deg 4) | 3 | 240 | 0.04 |
| genus-2 surface (deg 8) | 2 | 960 | 1.0 |
| ℝP² sextic in ℝ⁴ | 2 | 1,260 | 6 |
| three great circles (F, deg 6) | 4 | 1,920 | 15 |
| Clifford torus (F, deg 4, 4 variables) | 3 | 1,890 | 15 |
| 5 separated circles (deg 10) | 10 | 1,650 | 20 |
| 6 separated circles (deg 12) | 12 | 2,808 | 103 |
| 8 concentric circles (deg 16, critical circles) | 8 | 3,264 | 126 |

The time is spent in dense modular row reduction, one per prime and per bound tried;
a square-subsystem solve, threads across primes, or p-adic lifting would cut it.

## Conormal or polar?

In Macaulay2 the polar route is far faster than the conormal route (X₁ below: 2.7 s vs
210 s), because M2 saturates the conormal ideal by the entire ideal of minors. Helmer
saturates by one sparse minor instead (exact for prime input, see above), which removes
that bottleneck. As a result the conormal route is the default: on every benchmark so far it
is at least as fast as the polar route, and it is deterministic.

| Example | M2 conormal | M2 polar | Helmer conormal | Helmer polar |
|---|---|---|---|---|
| [HM] X₁ | 210 s | 2.7 s | 1.2 s | 0.7 s |
| [HM] X₂ | 7.2 s | 208 s | 2.7 s | 13.8 s |
| [HM] Example 1.1 | 14.1 s | 22.0 s | 0.8 s | 8.3 s |

The polar route draws random linear projections, so it returns valid but generally
non-minimal stratifications with random extra strata (exactly as in M2). Neither
route is minimal in general: on Example 1.1 the conormal route returns three
non-minimal curves, and on X₂ both conormal and the associated-primes route return a
non-minimal point. The tests check every polar output against the *minimal*
stratification computed by M2's `minCoarsenWS`. `polar_bound` (default 1000) sets the
size of the random coefficients: smaller is faster (X₂: 5.0 s at bound 5) but raises
the small probability of an unlucky projection.

Known hard cases for both systems: the Briançon–Speder family
`x⁵ + y⁷z + z¹⁵ + t·xy⁶` and a quartic hypersurface in ℂ⁵. A single Gröbner basis
(the conormal saturation) dominates; this is where the remaining speed has to come from.

## Validation and speed

`tests/test_helmer.py` checks Helmer against Macaulay2 1.22 + WhitneyStratifications
2.23 (with `AssocPrimes=>false`) on the package's documentation examples, the
benchmark varieties of [HM], and assorted singular surfaces, curves and projective
hypersurfaces; plus both root-classification examples for `mapStratify` and a battery
of `decompose` comparisons. All agree: same strata, same dimensions. Polar outputs
are checked to refine the minimal stratification (46 tests, ~20 s).

Wall-clock times on the same 2-core machine (M2 `cpuTime`; Helmer via Python):

| Example | Macaulay2 | Helmer |
|---|---|---|
| Whitney cusp | 0.09 s | 0.04 s |
| Projective quartic (`quartic_P`) | 0.39 s | 0.21 s |
| [HM] Example 1.1 (hypersurface in ℂ⁴) | 14.1 s | 0.65 s |
| [HM] X₁ (complete intersection in ℂ⁴) | 210 s | 0.93 s |
| [HM] X₂ | 7.2 s | 2.3 s |
| `mapStratify` Feynman example (package docs) | 1.38 s | 0.06 s |

## Differences from the Macaulay2 package

* Components of a non-equidimensional `X` are filed under their own dimension (the M2
  package files every component under `dim X`), and all pairs of components are
  intersected up front.
* In `mapStratify`, the M2 code pulls Whitney points of `Y` back to `X` with the
  *forward* elimination (`elimIm`), which is only correct for coordinate projections;
  Helmer uses the preimage `f⁻¹(K) ∩ X` throughout.
* Random choices (independent sets) are made deterministically; `helmer.set_seed`
  controls the remaining randomness.

## Building

Requires a C++20 compiler, CMake ≥ 3.18, FLINT ≥ 3.0, GMP and msolve ≥ 0.6 (its
`libneogb`). On Ubuntu 24.04: `apt install libflint-dev libmsolve-dev`, then

```
pip install .
pytest tests/
```

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
