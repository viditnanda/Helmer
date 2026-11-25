# Changes

## Unreleased

* Minimal Whitney stratifications: `min_coarsen` and `whitney_stratify(..., minimal=True)`
  (Macaulay2 `minCoarsenWS`), by Teissier's polar criterion. On every affine benchmark the
  polar route's stratification (random extra strata) and the conormal route's coarsen to the
  same answer; the M2 results recorded for [HM] X₂ and Example 1.1 are reproduced. Three
  defects of `minCoarsenWS` 2.23 are avoided (README, "Differences").
* Local invariants at a general point, computed over F_p on curve germs cut out by slices
  through one general F_p-point (Monte Carlo, two agreeing runs): `multiplicity` (Samuel
  multiplicity, any X ⊆ Y; no Segre classes needed), `polar_sequence`, `polar_multiplicity`,
  `in_singular_locus`, and `contained` (V(X) ⊆ V(Y)).
* Euler obstructions, characteristic cycles and Chern classes of strata: `chern_classes`
  (per stratum: χ, CSM class, Chern–Mather class, CC multiplicity), `euler_characteristic`,
  `csm_class`, `euler_obstruction_matrix` and `cc_multiplicities` (M2 `eulerObsMatrix`,
  `multCC`), `euler_obstruction`, `global_polar_degrees`, `global_euler_obstruction`. Local
  Euler obstructions by Lê–Teissier from polar sequences; Euler characteristics of all strata
  from global polar degrees (Seade–Tibăr–Verjovsky) by the unitriangular Euler-obstruction
  matrix; CSM and Chern–Mather classes from general linear sections by Schürmann's
  transversality theorem. No Segre classes; affine and projective.
* `whitney_filtration`: the closed dimension filtration, deferring prime decomposition
  (reduced complete intersections kept whole by Serre's criterion, finite loci made radical
  by Seidenberg's lemma, Whitney failures along curves bounded by affine-linear symmetries).
  Mechanisms such as the pinned five-bar and the rotating four-bar go from minutes to
  milliseconds. Chart saturations also prune their bases to minimal leading monomials.
* Internals: a mod-p layer over msolve's F4 and FLINT's `nmod_mpoly` (`helmer/modp.hpp`).

## 1.0.0 (2026-11-25)

First release.

* `whitney_stratify` (conormal and polar routes), `map_stratify`, `non_proper_set`, `conormal`,
  `polar_varieties`: the Whitney-stratification functions of Martin Helmer's Macaulay2 package
  `WhitneyStratifications` 2.23 (with `AssocPrimes=>false`), checked against it on its documentation
  examples and the benchmarks of Helmer–Mohr.
* The algebra they rest on, without Macaulay2: Gröbner bases over ℚ (msolve's F4, multi-modular), dimension,
  elimination, saturation, intersection, and minimal primes certified prime.
* Speed: hyperplane charts and single-minor saturation for conormal varieties, squarefree saturation, the
  reduced variety for singular loci, free translation directions split off before stratifying. On the
  Helmer–Mohr benchmarks, 38× to 5000× faster than the Macaulay2 package.
* Certified ideal membership and elimination (`certify_membership`, `verify_membership`,
  `eliminate_certified`): cofactor identities checked exactly over ℚ.
* The `helmer` command, a C++ example, six Python examples, and `docs/API.md`.
* Binary wheels for Linux (x86-64) and macOS 15 or later (arm64, x86-64), CPython 3.9–3.13, with GMP, MPFR, FLINT and
  msolve bundled and built without CPU-specific instructions.
