// Minimal primes over QQ (Macaulay2's `decompose`).
//
// Algorithm (GTZ-style, factorisation driven):
//  1. If some element of the reduced DRL basis factors, split along the factors.
//  2. Otherwise pick a maximal independent set u, compute a block basis for
//     [x\u] > [u] and h = lcm of its leading coefficients in k[u]; then
//     V(I) = V(I : h^∞) ∪ V(I + h).
//  3. I : h^∞ is the contraction of the zero-dimensional ideal I·k(u)[x\u]. Its
//     minimal primes come from factoring the minimal polynomial mu(T) of a
//     linear form z over k(u); each factor mu_i gives the candidate
//     (I + mu_i(z)) : h_i^∞, and is *certified* prime by checking
//     dim_{k(u)} k(u)[x\u]/(I + mu_i(z)) = deg mu_i (otherwise z is retried).
//  4. Recurse on I + (irreducible factors of h); discard non-minimal primes.
#pragma once

#include "helmer/ideal.hpp"

namespace helmer {

struct PrimeStats {
  long calls = 0;
  long splits = 0;
  long z_retries = 0;
};
PrimeStats& prime_stats();

std::vector<Ideal> minimal_primes(const Ideal& I);

// Keep only inclusion-minimal ideals, removing duplicates.
std::vector<Ideal> minimalise(std::vector<Ideal> P);

}  // namespace helmer
