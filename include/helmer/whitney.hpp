// Whitney stratification of complex (hence real) algebraic varieties.
//
// whitney_stratify implements the conormal algorithm of Helmer–Nanda (FoCM 2022)
// in the associated-prime-free form of Helmer–Mohr (arXiv:2406.17122), i.e.
// Macaulay2's  whitneyStratify(I, AssocPrimes=>false).
#pragma once

#include "helmer/primes.hpp"

#include <functional>
#include <map>

namespace helmer {

// Strata closures indexed by (complex) dimension. The open strata of dimension i
// are the connected components of V(P) minus the union of the lower-dimensional
// closures, for P in strata[i].
struct Stratification {
  RingPtr ring;
  bool projective = false;
  std::map<int, std::vector<Ideal>> strata;
};

struct WhitneyOptions {
  bool projective = false;
  int level = -1;                         // -1: all levels (M2 default n - codim)
  std::string saturation = "random";      // "random" | "exact" (conormal saturation)
  // "conormal": Helmer–Mohr conormal algorithm (deterministic up to the modular lift).
  // "polar": polar-variety algorithm (M2 whitneyStratifyPol): usually faster, uses
  //          random projections, and returns a valid but generally non-minimal
  //          stratification with probability ~1 (see docs).
  std::string method = "conormal";
  long polar_bound = 1000;  // random coefficients in [-B, B] for the polar projections
  std::function<void(const std::string&)> log;  // progress messages, may be empty
};

// Polar varieties of a prime I of codimension c in k^n: for l = 1..dim(I)+1, the
// closure over the smooth locus of { rank [Jac I ; A_l] < c + l } with A_l a random
// l x n matrix (i.e. the critical locus of a random linear projection to k^{l-1}...),
// saturated by the singular-locus ideal Jsing (M2: polarVars).
std::vector<Ideal> polar_varieties(const Ideal& I, const Ideal& Jsing, long bound = 1000);

// QQ[x_0..x_{n-1}, v_0..v_{n-1}]: the ring of the conormal variety in k^n x P^{n-1}.
RingPtr conormal_ring(const RingPtr& R);
// Conormal variety Con(V(I)) for I of codimension c: (I + (c+1)-minors[v; Jac I]) : (c-minors)^∞.
// saturation: "minor" (exact when I is prime: saturate by one sparse c-minor not in I),
// "random" (one random combination of all c-minors), or "exact" (intersection).
Ideal conormal(const Ideal& I, const RingPtr& S, const std::string& saturation = "random");
// Helmer–Mohr "WhitPoints": for Z ⊆ X, an ideal (in R) whose zero set contains the
// points of Z at which (X_reg, Z) can fail Whitney (B) (outside Sing Z).
Ideal whitney_points(const Ideal& Z, const Ideal& ConX, const RingPtr& R);

Stratification whitney_stratify(const Ideal& I, const WhitneyOptions& opts = {});

// ---------------------------------------------------------------- maps
// Stratify a polynomial map f = (F_1..F_m): X = V(IX) ⊆ k^n -> Y = V(IY) ⊆ k^m so that
// each stratum of X maps submersively into a stratum of Y (Thom's isotopy lemma
// then applies over each stratum of Y). Macaulay2: mapStratify(F, X, Y,
// AssocPrimes=>false).
struct MapStratification {
  Stratification source, target;
};
struct MapOptions {
  bool proper = true;              // if false, add the Jelonek non-properness flag
  std::string strata = "all";      // "all" | "singular_only"
  std::function<void(const std::string&)> log;
};
MapStratification map_stratify(const std::vector<Poly>& F, const Ideal& IX, const Ideal& IY,
                               const MapOptions& opts = {});
// Jelonek set of points of Y over which f|X is not proper (as an ideal of k[y]).
Ideal non_proper_set(const std::vector<Poly>& F, const Ideal& IX, const Ideal& IY);

}  // namespace helmer
