// Certified ideal membership: h = sum_j c_j g_j with explicit cofactors, checked by
// exact polynomial arithmetic over QQ.
//
// The cofactors are *found* by linear algebra modulo word-size primes with Chinese
// remaindering and rational reconstruction, which is heuristic; the certificate is
// *accepted* only after the identity has been verified exactly in fmpq_mpoly. A
// returned certificate is therefore a proof that h lies in the ideal (g_1, ..., g_m),
// independent of the modular steps and of the Gröbner basis computations.
#pragma once

#include "helmer/ideal.hpp"

#include <optional>

namespace helmer {

struct MembershipCertificate {
  Poly target;                  // h
  std::vector<Poly> gens;       // g_1, ..., g_m
  std::vector<Poly> cofactors;  // c_1, ..., c_m with h == sum_j c_j g_j
  int bound_block = 0;          // degree bound used in the block variables
  int bound_rest = 0;           // degree bound used in the other variables
};

// Exact check of h == sum_j c_j g_j (all in one ring).
bool verify(const MembershipCertificate& C);
// The same, and C.gens are exactly the generators of I (so C proves C.target in I).
bool verify(const Ideal& I, const MembershipCertificate& C);

struct CertifyOptions {
  long max_unknowns = 4000;  // largest linear system tried (columns); the solve is dense,
                             // so memory ~ 8 * rows * columns bytes and time ~ columns^3
  int max_primes = 64;       // primes per support before giving up on it
  int extra_block = 4;       // block bounds tried: d_A - 2, ..., max(d_A, deg_A h) + extra_block
  int extra_rest = 2;        // other bounds: d_B(h) - 1, ..., max(d_B(h), deg h - e) + extra_rest
};

// Search for cofactors c_j whose terms have degree <= D_A in the `block` variables
// and <= D_B in the others, over a small range of (D_A, D_B) ordered by system
// size. d_A is the largest block degree of the g_j, deg_A h and d_B(h) the degrees of
// h in the block and other variables, e the smallest total degree of a g_j. Returns a verified certificate, or nothing if none was found
// within the limits (which proves nothing either way).
std::optional<MembershipCertificate> certify_membership(const std::vector<Poly>& gens,
                                                        const Poly& h,
                                                        const std::vector<int>& block,
                                                        const CertifyOptions& opt = {});

// The same with explicit bounds.
std::optional<MembershipCertificate> certify_membership_bounded(
    const std::vector<Poly>& gens, const Poly& h, const std::vector<int>& block, int DA,
    int DB, int max_primes = 64);

// eliminate(I, elim) together with a membership certificate in I for each generator
// of the elimination ideal (nullopt where none was found). A certified generator P
// vanishes on V(I) whatever the reliability of the elimination itself; that is the
// direction needed to bound critical values. Only a *nonzero* certified generator
// bounds anything: if eliminate() returns the zero ideal there is nothing to certify
// and `certificates` is empty, so callers must check for one.
struct CertifiedElimination {
  Ideal ideal;
  std::vector<std::optional<MembershipCertificate>> certificates;
};
CertifiedElimination eliminate_certified(const Ideal& I, const std::vector<int>& elim,
                                         const CertifyOptions& opt = {});

}  // namespace helmer
