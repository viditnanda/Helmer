// Ideals over QQ and the commutative-algebra operations the Whitney algorithms use.
#pragma once

#include "helmer/gb.hpp"

#include <optional>

namespace helmer {

class Ideal {
 public:
  Ideal(RingPtr R, std::vector<Poly> gens = {});
  static Ideal parse(RingPtr R, const std::vector<std::string>& gens);

  const RingPtr& ring() const { return R_; }
  const std::vector<Poly>& gens() const { return gens_; }

  // Reduced DRL Gröbner basis (cached). After canonicalise(), gens() == gb().polys.
  const GB& gb() const;
  Ideal& canonicalise();

  bool is_unit() const;
  bool is_zero() const;
  // Krull dimension of R/I; -1 for the unit ideal.
  int dim() const;
  // nvars - dim; nvars + 1 for the unit ideal ("infinity").
  int codim() const;
  // A maximal independent set of variables of size dim() (indices into the ring).
  std::vector<int> independent_set() const;

  bool contains(const Poly& f) const;      // ideal membership (exact normal form)
  bool contains(const Ideal& J) const;     // J ⊆ this
  bool equals(const Ideal& J) const;       // same ideal
  Poly normal_form(const Poly& f) const;

  bool is_homogeneous() const;
  std::string str() const;
  std::string key() const;  // canonical text of the reduced GB (for hashing / dedupe)

 private:
  RingPtr R_;
  std::vector<Poly> gens_;
  mutable std::shared_ptr<GB> gb_;
};

// Largest set of variables containing no leading monomial's support (the Krull
// dimension of a monomial ideal is its size). Leading exponents in ring order.
std::vector<int> max_independent_set(int n, const std::vector<std::vector<int>>& leads);

Ideal operator+(const Ideal& I, const Ideal& J);
Ideal operator+(const Ideal& I, const Poly& f);

// I ∩ k[vars not in elim], returned in the same ring (generators avoid `elim`).
Ideal eliminate(const Ideal& I, const std::vector<int>& elim);
// (I : f^∞) by Rabinowitsch's trick.
Ideal saturate(const Ideal& I, const Poly& f);
// (I : J^∞). Strategy "random": saturate by one random combination of J's
// generators (equal to I:J^∞ with probability ~1 - 1/2^30 when it matters; see
// docs). Strategy "exact": intersect the saturations by each generator.
Ideal saturate(const Ideal& I, const Ideal& J, const std::string& strategy = "random");
Ideal intersect(const Ideal& I, const Ideal& J);

// Embed an ideal into a bigger ring (variable i -> map[i]); and pull an ideal whose
// generators only involve the image variables back down.
Ideal embed(const Ideal& I, const RingPtr& S, const std::vector<int>& map);
Ideal restrict_to(const Ideal& I, const RingPtr& R, const std::vector<int>& map);

// Jacobian matrix: rows = generators, columns = variables.
std::vector<std::vector<Poly>> jacobian(const std::vector<Poly>& F);
// All k x k minors (nonzero ones) of a matrix of polynomials.
std::vector<Poly> minors(int k, const std::vector<std::vector<Poly>>& M);
// Random QQ-linear combination of polynomials (integer coefficients in [-B, B]).
Poly random_combination(const std::vector<Poly>& F, long bound = (1L << 30));

// Singular-locus ideal of an equidimensional radical ideal of codimension c:
// I + c-minors of Jac(I).
Ideal jacobian_ideal(const Ideal& I);

// Block-order data used by the Helmer–Mohr / GTZ steps: Gröbner basis w.r.t. the
// order [vars not in u] > [u], and lcm over its elements of their leading
// coefficients viewed in k[u].
struct BlockData {
  GB gb;
  Poly h;  // lcm of leading coefficients in k[u] (primitive)
};
BlockData block_lc_lcm(const Ideal& I, const std::vector<int>& u);

}  // namespace helmer
