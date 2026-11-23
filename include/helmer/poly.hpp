// Helmer — Whitney stratification in C++.
// Polynomial layer: thin RAII wrappers around FLINT's fmpq_mpoly (exact rational
// multivariate polynomials). Arithmetic, parsing, factorisation and gcds are FLINT's;
// Gröbner bases are msolve's (see gb.hpp).
#pragma once

#include <flint/fmpq_mpoly.h>
#include <flint/fmpq_mpoly_factor.h>

#include <cstdint>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace helmer {

class Ring;
using RingPtr = std::shared_ptr<Ring>;

// A polynomial ring QQ[x_0, ..., x_{n-1}] with named variables.
// FLINT ordering is degrevlex with x_0 > x_1 > ... ; this agrees with msolve's DRL,
// so FLINT division by an msolve DRL basis is a valid normal form (tested in tests/).
class Ring : public std::enable_shared_from_this<Ring> {
 public:
  static RingPtr make(std::vector<std::string> names);
  ~Ring();
  Ring(const Ring&) = delete;
  Ring& operator=(const Ring&) = delete;

  int nvars() const { return static_cast<int>(names_.size()); }
  const std::vector<std::string>& names() const { return names_; }
  int index_of(const std::string& name) const;  // -1 if absent
  fmpq_mpoly_ctx_struct* ctx() const { return ctx_; }
  const char** cnames() const { return const_cast<const char**>(cnames_.data()); }

 private:
  explicit Ring(std::vector<std::string> names);
  std::vector<std::string> names_;
  std::vector<const char*> cnames_;
  mutable fmpq_mpoly_ctx_t ctx_;
};

class Poly {
 public:
  explicit Poly(RingPtr R);                          // zero
  Poly(RingPtr R, const std::string& s);             // parse, e.g. "y^2*z - x^2"
  static Poly constant(RingPtr R, long c);
  static Poly var(RingPtr R, int i);
  // Build from exponent vectors (ring variable order) and rational coefficients.
  static Poly from_terms(RingPtr R, const std::vector<std::vector<long>>& exps,
                         const std::vector<fmpq*>& coeffs);

  Poly(const Poly& o);
  Poly(Poly&& o) noexcept;
  Poly& operator=(const Poly& o);
  Poly& operator=(Poly&& o) noexcept;
  ~Poly();

  const RingPtr& ring() const { return R_; }
  fmpq_mpoly_struct* raw() { return p_; }
  const fmpq_mpoly_struct* raw() const { return p_; }

  bool is_zero() const;
  bool is_constant() const;
  bool is_one() const;
  long nterms() const;
  long total_degree() const;
  long degree_in(int var) const;
  bool involves(int var) const;
  std::vector<long> term_exp(long i) const;  // exponent vector of i-th term (FLINT order)
  std::string str() const;

  Poly operator+(const Poly& o) const;
  Poly operator-(const Poly& o) const;
  Poly operator*(const Poly& o) const;
  Poly operator-() const;
  Poly scale_si(long c) const;
  Poly derivative(int var) const;
  bool equal(const Poly& o) const;

  // Map into ring S: variable i of this ring goes to variable map[i] of S.
  Poly embed(const RingPtr& S, const std::vector<int>& map) const;
  // Substitute: variable i of this ring -> vals[i] (all in ring S).
  Poly substitute(const RingPtr& S, const std::vector<Poly>& vals) const;

  // Normalisations
  Poly primitive() const;  // integer-primitive, positive leading coefficient (FLINT order)
  Poly monic() const;      // leading coefficient 1 (FLINT order)

 private:
  RingPtr R_;
  fmpq_mpoly_struct* p_;
};

Poly gcd(const Poly& a, const Poly& b);
Poly lcm(const Poly& a, const Poly& b);
Poly divexact(const Poly& a, const Poly& b);
// Distinct irreducible non-constant factors with multiplicities.
std::vector<std::pair<Poly, long>> factor(const Poly& a);
// Squarefree "radical" part: product of distinct irreducible factors (primitive).
Poly radical_part(const Poly& a);

// Process-wide RNG for the probabilistic steps (seedable for reproducibility).
std::mt19937_64& rng();
std::mt19937_64& gb_prime_rng();
void set_seed(uint64_t s);
long random_coeff(long bound);  // uniform in [-bound, bound] \ {0}

}  // namespace helmer
