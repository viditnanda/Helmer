// Euler obstructions, characteristic cycles, Euler characteristics and Chern classes of
// the strata of a Whitney stratification.
//
// Local data: Eu_{T-bar} at a general point of a stratum S, by the Lê–Teissier formula
// from the polar sequence (multiplicity.hpp). Global data: for an affine variety Z
// (pure dimension d), the global Euler obstruction of Seade–Tibăr–Verjovsky,
//     Eu(Z) = sum over strata S of Z of chi(S) Eu_Z(S) = sum_{i=0}^{d} (-1)^{d-i} gamma_i(Z),
// where gamma_i is the degree of the i-dimensional global polar variety (the closure
// of the critical locus of a general linear projection Z_reg -> C^{i+1}), a point count.
// Since the Euler-obstruction matrix of a stratification is unitriangular, the global
// Euler obstructions of the stratum closures determine the Euler characteristics of
// all strata. Applied to general linear sections, this gives chi(S ∩ L_r) for general
// codimension-r linear spaces L_r, which by Schürmann's transversality theorem
// (c_*(phi · 1_L) = c(N_L)^{-1} ∩ c_*(phi) · [L]) determine Chern–Schwartz–MacPherson
// classes in A_*(P^N) by a triangular binomial system; likewise Chern–Mather classes,
// from Eu instead of 1. No Segre classes are involved.
//
// Monte Carlo, like multiplicity.hpp: every invariant is computed modulo random primes
// with random choices until `trials` runs agree.
#pragma once

#include "helmer/multiplicity.hpp"

#include <vector>

namespace helmer {

// Eu_Z at a general point of X (Z prime, X irreducible): 0 if X is not in V(Z).
long euler_obstruction(const Ideal& X, const Ideal& Z, const LocalOptions& o = {});

// Degrees gamma_0, ..., gamma_d of the global polar varieties of the affine variety
// V(Z) (Z radical, pure dimension d), gamma_d = deg Z.
std::vector<long> global_polar_degrees(const Ideal& Z, const LocalOptions& o = {});

// The global Euler obstruction sum_i (-1)^(d-i) gamma_i of the affine variety V(Z); it
// equals chi(V(Z)) when V(Z) is smooth.
long global_euler_obstruction(const Ideal& Z, const LocalOptions& o = {});

// The strata of a stratification in a fixed order, with their Euler obstructions:
// eu[s][t] = Eu of the closure of stratum t at a general point of stratum s.
struct EulerData {
  RingPtr ring;
  bool projective = false;
  std::vector<int> dims;        // as in Stratification (projective dimensions if projective)
  std::vector<Ideal> strata;    // stratum closures (prime)
  std::vector<std::vector<long>> eu;
};
EulerData euler_obstruction_matrix(const Stratification& S, const LocalOptions& o = {},
                                   const std::function<void(const std::string&)>& log = {});

// Characteristic-cycle multiplicities of the constant function on the variety, as
// Macaulay2's multCC: the m_t with sum_t (-1)^(dim X - dim t) m_t Eu_{closure t} = 1.
std::vector<long> cc_multiplicities(const EulerData& E);

struct ChernData {
  EulerData E;
  int ambient = 0;                          // N: classes live in A_*(P^N) (A^n ⊂ P^n if affine)
  std::vector<std::vector<long>> sections;  // [s][r] = chi(open stratum s ∩ general L_r)
  std::vector<long> chi;                    // chi of each open stratum
  std::vector<std::vector<long>> csm;       // [s][k] = degree of the dimension-k part of c_SM(1_s)
  std::vector<std::vector<long>> mather;    // [t][k] = same for c_M(closure t) = c_*(Eu_{closure t})
};
ChernData chern_classes(const Stratification& S, const LocalOptions& o = {},
                        const std::function<void(const std::string&)>& log = {});

}  // namespace helmer
