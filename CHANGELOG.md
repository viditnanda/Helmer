# Changes

## 0.1.0 (unreleased)

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
* Binary wheels for Linux (x86-64) and macOS (arm64, x86-64), CPython 3.9–3.13, with GMP, MPFR, FLINT and
  msolve bundled and built without CPU-specific instructions.
