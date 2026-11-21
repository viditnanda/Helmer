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
| `conormal(I)` | `conormal(I)` | In kⁿ × ℙⁿ⁻¹, variables `v_0..v_{n-1}`. |
| `decompose(I)` | `minimal_primes(I)` | New implementation; see below. |
| `eliminate`, `saturate`, `intersect`, `gb`, `dim`, `codim`, `independentSets` | same names / `Ideal` methods | |

Not yet ported: the polar-variety route (`whitneyStratifyPol`, `mapStratifyPol`,
`polarVars`), minimisation (`minCoarsenWS`), and the Euler-obstruction /
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

## Validation and speed

`tests/test_helmer.py` checks Helmer against Macaulay2 1.22 + WhitneyStratifications
2.23 (with `AssocPrimes=>false`) on the package's documentation examples, the
benchmark varieties of [HM], and assorted singular surfaces, curves and projective
hypersurfaces; plus both root-classification examples for `mapStratify` and a battery
of `decompose` comparisons. All agree: same strata, same dimensions (32 tests, ~4 s).

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
