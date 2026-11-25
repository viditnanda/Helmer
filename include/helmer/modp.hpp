// Polynomials over F_p and the few operations on them that the local invariants
// (multiplicity.hpp) need: Gröbner bases by msolve's F4 modulo a single prime, block
// elimination, saturation by one element, and the length of a zero-dimensional
// quotient. Nothing here is lifted back to QQ; these computations are Monte Carlo
// by design, with the random choices drawn from F_p for a random ~31-bit prime p.
#pragma once

#include "helmer/poly.hpp"

#include <flint/nmod_mpoly.h>

#include <optional>

namespace helmer::modp {

// A reduction of rational data hit a denominator divisible by p: choose another prime.
struct Unlucky : std::runtime_error {
  Unlucky() : std::runtime_error("helmer: unlucky prime") {}
};
// A random choice turned out non-generic (detected): retry with fresh choices.
struct Retry : std::runtime_error {
  explicit Retry(const std::string& why) : std::runtime_error("helmer: retry: " + why) {}
};

// The subvariety has no F_p-point on the slices tried (e.g. conjugate points whose
// field of definition does not split at p): try another prime.
struct NoPoint : std::runtime_error {
  NoPoint() : std::runtime_error("helmer: no F_p-point found") {}
};

// F_p[z_0, ..., z_{n-1}], degrevlex.
class Ctx {
 public:
  Ctx(ulong p, int n);
  ~Ctx();
  Ctx(const Ctx&) = delete;
  Ctx& operator=(const Ctx&) = delete;
  ulong p;
  int n;
  nmod_mpoly_ctx_t ctx;
};
using CtxPtr = std::shared_ptr<Ctx>;
CtxPtr make_ctx(ulong p, int n);

class MPoly {
 public:
  explicit MPoly(CtxPtr C);
  MPoly(const MPoly& o);
  MPoly(MPoly&& o) noexcept;
  MPoly& operator=(const MPoly& o);
  MPoly& operator=(MPoly&& o) noexcept;
  ~MPoly();

  static MPoly constant(const CtxPtr& C, ulong c);
  static MPoly var(const CtxPtr& C, int i);

  const CtxPtr& ctx() const { return C_; }
  nmod_mpoly_struct* raw() { return p_; }
  const nmod_mpoly_struct* raw() const { return p_; }

  bool is_zero() const;
  long nterms() const;
  long total_degree() const;
  long degree_in(int i) const;
  ulong constant_term() const;
  ulong eval(const std::vector<ulong>& pt) const;

  MPoly operator+(const MPoly& o) const;
  MPoly operator-(const MPoly& o) const;
  MPoly operator*(const MPoly& o) const;
  MPoly scale(ulong c) const;
  MPoly derivative(int i) const;
  // Substitute variable i -> vals[i] (all in the context of vals).
  MPoly compose(const CtxPtr& target, const std::vector<MPoly>& vals) const;

 private:
  CtxPtr C_;
  nmod_mpoly_struct* p_;
};

ulong random_prime();          // uniform-ish 31-bit prime from Helmer's seeded RNG
ulong random_elem(ulong p);    // uniform in F_p
ulong random_unit(ulong p);    // uniform in F_p^*

// Reduce a rational polynomial (same number of variables) modulo p.
MPoly reduce(const Poly& f, const CtxPtr& C);  // throws Unlucky

struct ModGB {
  std::vector<MPoly> polys;              // reduced basis, leading term first
  std::vector<std::vector<int>> lead;    // leading exponents
  bool unit = false;
};
// Reduced Gröbner basis for DRL (ebl = 0) or the block order eliminating the first
// ebl variables (DRL on each block).
ModGB groebner(const CtxPtr& C, const std::vector<MPoly>& F, int ebl = 0);

// Krull dimension of F_p[z]/(lead terms); -1 for the unit ideal.
int dim(const ModGB& G, int n);
// Number of standard monomials (= dim_{F_p} F_p[z]/I) of a zero-dimensional
// basis; -1 if the ideal is not zero-dimensional.
long staircase(const ModGB& G, int n);
long length(const CtxPtr& C, const std::vector<MPoly>& F);
// Length of the part of a zero-dimensional F away from V(h): len F_p[z,T]/(F, 1 - T h).
long length_off(const CtxPtr& C, const std::vector<MPoly>& F, const MPoly& h);
// F : h^infinity (generators of a Gröbner basis), by eliminating T from (F, 1 - T h).
std::vector<MPoly> saturate(const CtxPtr& C, const std::vector<MPoly>& F, const MPoly& h);

// h vanishes on V(F), i.e. h is in the radical of F: (F, 1 - T h) is the unit ideal.
bool in_radical(const CtxPtr& C, const std::vector<MPoly>& F, const MPoly& h);

// Affine map z -> x + B z from F_p^s into F_p^n with B random; images of x_0..x_{n-1}.
std::vector<MPoly> random_affine_map(const CtxPtr& Cz, const std::vector<ulong>& x, int n);

// Maximal minors of a k x n matrix of polynomials (k <= n), all C(n, k) of them.
std::vector<MPoly> maximal_minors(const CtxPtr& C, const std::vector<std::vector<MPoly>>& M);

}  // namespace helmer::modp
