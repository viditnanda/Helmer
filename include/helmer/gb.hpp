// Gröbner bases over QQ via msolve's F4 (modulo 31-bit primes) + CRT + rational
// reconstruction. Two orders are supported, which is all Helmer needs:
//   * DRL on all variables (ebl = 0), and
//   * a two-block elimination order: DRL on the first `ebl` variables of `perm`,
//     then DRL on the rest.
#pragma once

#include "helmer/poly.hpp"

namespace helmer {

struct GBStats {
  long calls = 0;        // groebner() calls
  long f4_runs = 0;      // msolve F4 runs (one per prime)
  long unlucky = 0;      // primes discarded as unlucky
  double seconds = 0.0;  // wall time inside groebner()
};
GBStats& gb_stats();

struct GB {
  // Reduced Gröbner basis, each element monic w.r.t. the order used.
  std::vector<Poly> polys;
  // Leading exponent vector of each element, in *ring* variable order.
  std::vector<std::vector<int>> lead;
  bool is_unit() const;  // basis == {1}
};

// perm[k] = ring index of the k-th variable in the order; must be a permutation.
GB groebner(const RingPtr& R, const std::vector<Poly>& F, const std::vector<int>& perm, int ebl);
GB groebner_drl(const RingPtr& R, const std::vector<Poly>& F);

// Options for the modular lift.
struct GBOptions {
  int min_primes = 2;        // primes agreeing before accepting a reconstruction
  int max_primes = 4000;     // hard cap
  int msolve_threads = 1;
};
GBOptions& gb_options();

}  // namespace helmer
