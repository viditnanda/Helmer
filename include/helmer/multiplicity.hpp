// Local invariants of varieties at a general point of a subvariety, and the minimal
// (coarsest) Whitney stratification.
//
// Everything here is computed at one *general* F_p-point of the subvariety, for a
// random 31-bit prime p and random affine slices; it is Monte Carlo, and each call
// repeats the computation with independent primes and slices until `trials` runs
// agree (default 2). An unlucky choice makes a run disagree, never two in the same way
// except with negligible probability.
//
// How: e_X(Y) = e_x(Y) for a general point x of X. Such an x is found over F_p by
// slicing X with dim X random affine hyperplanes and taking an F_p-root of the
// eliminant. Then Y is cut by w-1 general hyperplanes through x (w = dim_x Y), which
// preserves the multiplicity at x; saturating by a further general linear form through
// x removes any embedded component at x, leaving a Cohen–Macaulay curve germ, whose
// multiplicity is the length at x of its section by one more general hyperplane through
// x. Every system solved has length at most deg Y. (Macaulay2's SegreClasses obtains
// e_X(Y) instead from a global Bezout count with d^(dim Y - dim X) deg Y points.)
#pragma once

#include "helmer/whitney.hpp"

#include <vector>

namespace helmer {

struct LocalOptions {
  int trials = 2;     // independent runs that must agree
  int max_runs = 8;   // give up (throw) after this many runs without agreement
  // Primes skipped because the subvariety has no F_p-point there (cheap; a stratum of
  // conjugate points defined over a field K has one at a proportion of primes that is at
  // least 1/[K:QQ] by Chebotarev).
  int max_pointless_primes = 400;
};

// Samuel multiplicity e_X(Y) of Y along X: the multiplicity of Y at a general point of
// X (0 if X is not contained in V(Y)). X should be irreducible (a prime ideal); Y is
// any ideal in the same ring, and only the components of Y of largest dimension through
// X count. Projective input (homogeneous X and Y) is handled as affine cones, which
// have the same multiplicities.
long multiplicity(const Ideal& X, const Ideal& Y, const LocalOptions& o = {});

// V(X) ⊆ V(Y) (radical containment: some power of every generator of Y lies in X),
// tested modulo random primes with a random combination of the generators of Y.
// For containment of ideals, Y ⊆ X, use Ideal::contains.
bool contained(const Ideal& X, const Ideal& Y, const LocalOptions& o = {});

// X ⊆ Sing V(Y), for X irreducible and Y prime: e_X(Y) > 1 (Y is regular at a point of
// multiplicity one). Needs no equations for the singular locus.
bool in_singular_locus(const Ideal& X, const Ideal& Y, const LocalOptions& o = {});

// Multiplicities at a general point of X of the generic polar varieties of Z (Z prime
// of dimension d): entry j is the multiplicity of the j-dimensional polar variety
// (the closure of the critical locus of a general linear projection Z_reg -> k^(j+1)),
// for j = 0..d; entry d is e_X(Z). Zero where X is not contained in the polar variety.
// By Teissier, (Z_reg, Y) satisfies Whitney's conditions at general points of Y iff
// the sequence is the same at general points of all strata involved.
std::vector<long> polar_sequence(const Ideal& X, const Ideal& Z, const LocalOptions& o = {});
// One entry of polar_sequence: the multiplicity at a general point of X of the
// w-dimensional generic polar variety of Z.
long polar_multiplicity(const Ideal& X, const Ideal& Z, int w, const LocalOptions& o = {});

// The coarsest Whitney stratification refined by S (S itself a Whitney stratification,
// such as whitney_stratify's output). A stratum X is merged into the open stratum of
// Y when Y is the unique smallest stratum whose closure contains X, every stratum
// whose closure contains X also contains Y, Y-bar is smooth along X (e_X(Y-bar) = 1),
// and every such larger stratum Z has the same polar sequence along X as along Y.
Stratification min_coarsen(const Stratification& S, const LocalOptions& o = {},
                           const std::function<void(const std::string&)>& log = {});

}  // namespace helmer
