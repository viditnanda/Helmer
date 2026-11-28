#include "helmer/poly.hpp"

#include <flint/flint.h>

namespace helmer {

// ---------------------------------------------------------------- Ring
Ring::Ring(std::vector<std::string> names) : names_(std::move(names)) {
  if (names_.empty()) throw std::invalid_argument("helmer: a ring needs at least one variable");
  for (auto& s : names_) cnames_.push_back(s.c_str());
  fmpq_mpoly_ctx_init(ctx_, static_cast<slong>(names_.size()), ORD_DEGREVLEX);
}
Ring::~Ring() { fmpq_mpoly_ctx_clear(ctx_); }

RingPtr Ring::make(std::vector<std::string> names) {
  return RingPtr(new Ring(std::move(names)));
}

int Ring::index_of(const std::string& name) const {
  for (int i = 0; i < nvars(); ++i)
    if (names_[i] == name) return i;
  return -1;
}

// ---------------------------------------------------------------- Poly
Poly::Poly(RingPtr R) : R_(std::move(R)) {
  p_ = static_cast<fmpq_mpoly_struct*>(flint_malloc(sizeof(fmpq_mpoly_struct)));
  fmpq_mpoly_init(p_, R_->ctx());
}

Poly::Poly(RingPtr R, const std::string& s) : Poly(std::move(R)) {
  if (fmpq_mpoly_set_str_pretty(p_, s.c_str(), R_->cnames(), R_->ctx()) != 0)
    throw std::invalid_argument("helmer: cannot parse polynomial '" + s + "'");
}

Poly Poly::constant(RingPtr R, long c) {
  Poly p(std::move(R));
  fmpq_mpoly_set_si(p.p_, c, p.R_->ctx());
  return p;
}

Poly Poly::var(RingPtr R, int i) {
  Poly p(std::move(R));
  fmpq_mpoly_gen(p.p_, i, p.R_->ctx());
  return p;
}

Poly Poly::from_terms(RingPtr R, const std::vector<std::vector<long>>& exps,
                      const std::vector<fmpq*>& coeffs) {
  Poly p(std::move(R));
  const int n = p.R_->nvars();
  std::vector<ulong> e(n);
  for (size_t t = 0; t < exps.size(); ++t) {
    for (int i = 0; i < n; ++i) e[i] = static_cast<ulong>(exps[t][i]);
    fmpq_mpoly_push_term_fmpq_ui(p.p_, coeffs[t], e.data(), p.R_->ctx());
  }
  fmpq_mpoly_sort_terms(p.p_, p.R_->ctx());
  fmpq_mpoly_combine_like_terms(p.p_, p.R_->ctx());
  return p;
}

Poly::Poly(const Poly& o) : Poly(o.R_) { fmpq_mpoly_set(p_, o.p_, R_->ctx()); }
Poly::Poly(Poly&& o) noexcept : R_(std::move(o.R_)), p_(o.p_) { o.p_ = nullptr; }
Poly& Poly::operator=(const Poly& o) {
  if (this == &o) return *this;
  Poly tmp(o);
  *this = std::move(tmp);
  return *this;
}
Poly& Poly::operator=(Poly&& o) noexcept {
  if (this == &o) return *this;
  if (p_) {
    fmpq_mpoly_clear(p_, R_->ctx());
    flint_free(p_);
  }
  R_ = std::move(o.R_);
  p_ = o.p_;
  o.p_ = nullptr;
  return *this;
}
Poly::~Poly() {
  if (p_) {
    fmpq_mpoly_clear(p_, R_->ctx());
    flint_free(p_);
  }
}

bool Poly::is_zero() const { return fmpq_mpoly_is_zero(p_, R_->ctx()); }
bool Poly::is_constant() const { return fmpq_mpoly_is_fmpq(p_, R_->ctx()); }
bool Poly::is_one() const { return fmpq_mpoly_is_one(p_, R_->ctx()); }
long Poly::nterms() const { return fmpq_mpoly_length(p_, R_->ctx()); }
long Poly::total_degree() const { return fmpq_mpoly_total_degree_si(p_, R_->ctx()); }
long Poly::degree_in(int var) const { return fmpq_mpoly_degree_si(p_, var, R_->ctx()); }
bool Poly::involves(int var) const { return !is_zero() && degree_in(var) > 0; }

std::vector<long> Poly::term_exp(long i) const {
  std::vector<slong> e(R_->nvars());
  fmpq_mpoly_get_term_exp_si(e.data(), p_, i, R_->ctx());
  return std::vector<long>(e.begin(), e.end());
}

std::string Poly::str() const {
  char* c = fmpq_mpoly_get_str_pretty(p_, R_->cnames(), R_->ctx());
  std::string s(c);
  flint_free(c);
  return s;
}

#define HELMER_SAME_RING(o) \
  if (R_.get() != (o).R_.get()) throw std::invalid_argument("helmer: polynomials in different rings")

Poly Poly::operator+(const Poly& o) const {
  HELMER_SAME_RING(o);
  Poly r(R_);
  fmpq_mpoly_add(r.p_, p_, o.p_, R_->ctx());
  return r;
}
Poly Poly::operator-(const Poly& o) const {
  HELMER_SAME_RING(o);
  Poly r(R_);
  fmpq_mpoly_sub(r.p_, p_, o.p_, R_->ctx());
  return r;
}
Poly Poly::operator*(const Poly& o) const {
  HELMER_SAME_RING(o);
  Poly r(R_);
  fmpq_mpoly_mul(r.p_, p_, o.p_, R_->ctx());
  return r;
}
Poly Poly::operator-() const {
  Poly r(R_);
  fmpq_mpoly_neg(r.p_, p_, R_->ctx());
  return r;
}
Poly Poly::scale_si(long c) const {
  Poly r(R_);
  fmpq_mpoly_scalar_mul_si(r.p_, p_, c, R_->ctx());
  return r;
}
Poly Poly::derivative(int var) const {
  Poly r(R_);
  fmpq_mpoly_derivative(r.p_, p_, var, R_->ctx());
  return r;
}
bool Poly::equal(const Poly& o) const {
  HELMER_SAME_RING(o);
  return fmpq_mpoly_equal(p_, o.p_, R_->ctx());
}

Poly Poly::embed(const RingPtr& S, const std::vector<int>& map) const {
  if (static_cast<int>(map.size()) != R_->nvars())
    throw std::invalid_argument("helmer: embed map has wrong length");
  std::vector<slong> c(map.begin(), map.end());
  Poly r(S);
  fmpq_mpoly_compose_fmpq_mpoly_gen(r.p_, p_, c.data(), R_->ctx(), S->ctx());
  return r;
}

Poly Poly::substitute(const RingPtr& S, const std::vector<Poly>& vals) const {
  if (static_cast<int>(vals.size()) != R_->nvars())
    throw std::invalid_argument("helmer: substitute needs one value per variable");
  std::vector<fmpq_mpoly_struct*> C;
  for (auto& v : vals) {
    if (v.ring().get() != S.get()) throw std::invalid_argument("helmer: substitute ring mismatch");
    C.push_back(const_cast<fmpq_mpoly_struct*>(v.raw()));
  }
  Poly r(S);
  if (!fmpq_mpoly_compose_fmpq_mpoly(r.p_, p_, C.data(), R_->ctx(), S->ctx()))
    throw std::runtime_error("helmer: substitution overflow");
  return r;
}

Poly Poly::primitive() const {
  if (is_zero()) return *this;
  Poly r(R_);
  fmpq_t c;
  fmpq_init(c);
  fmpq_mpoly_content(c, p_, R_->ctx());
  fmpq_mpoly_scalar_div_fmpq(r.p_, p_, c, R_->ctx());
  fmpq_mpoly_get_term_coeff_fmpq(c, r.p_, 0, R_->ctx());
  if (fmpq_sgn(c) < 0) fmpq_mpoly_neg(r.p_, r.p_, R_->ctx());
  fmpq_clear(c);
  return r;
}

Poly Poly::monic() const {
  if (is_zero()) return *this;
  Poly r(R_);
  fmpq_mpoly_make_monic(r.p_, p_, R_->ctx());
  return r;
}

// ---------------------------------------------------------------- free functions
Poly gcd(const Poly& a, const Poly& b) {
  Poly g(a.ring());
  if (!fmpq_mpoly_gcd(g.raw(), a.raw(), b.raw(), a.ring()->ctx()))
    throw std::runtime_error("helmer: gcd failed");
  return g;
}

Poly divexact(const Poly& a, const Poly& b) {
  Poly q(a.ring());
  if (!fmpq_mpoly_divides(q.raw(), a.raw(), b.raw(), a.ring()->ctx()))
    throw std::runtime_error("helmer: inexact division");
  return q;
}

Poly lcm(const Poly& a, const Poly& b) {
  if (a.is_zero() || b.is_zero()) return Poly(a.ring());
  Poly g = gcd(a, b);
  return (divexact(a, g) * b).primitive();
}

std::vector<std::pair<Poly, long>> factor(const Poly& a) {
  std::vector<std::pair<Poly, long>> out;
  if (a.is_zero() || a.is_constant()) return out;
  const auto& R = a.ring();
  fmpq_mpoly_factor_t f;
  fmpq_mpoly_factor_init(f, R->ctx());
  if (!fmpq_mpoly_factor(f, a.raw(), R->ctx())) {
    fmpq_mpoly_factor_clear(f, R->ctx());
    throw std::runtime_error("helmer: factorisation failed");
  }
  for (slong i = 0; i < f->num; ++i) {
    Poly p(R);
    fmpq_mpoly_set(p.raw(), f->poly + i, R->ctx());
    out.emplace_back(p.primitive(), static_cast<long>(fmpz_get_si(f->exp + i)));
  }
  fmpq_mpoly_factor_clear(f, R->ctx());
  return out;
}

Poly radical_part(const Poly& a) {
  if (a.is_zero()) return a;
  // In characteristic zero the gcd with all partials removes one copy
  // of each irreducible factor from a.
  Poly repeated = a;
  for (int i = 0; i < a.ring()->nvars() && !repeated.is_constant(); ++i) {
    Poly partial = a.derivative(i);
    if (!partial.is_zero()) repeated = gcd(repeated, partial);
  }
  return divexact(a, repeated).primitive();
}

std::mt19937_64& rng() {
  static std::mt19937_64 g(0x48656c6d6572ULL);  // "Helmer"
  return g;
}
std::mt19937_64& gb_prime_rng() {
  static std::mt19937_64 g(std::random_device{}());
  return g;
}
void set_seed(uint64_t s) {
  rng().seed(s);
  gb_prime_rng().seed(s);
}
long random_coeff(long bound) {
  std::uniform_int_distribution<long> d(1, bound);
  long c = d(rng());
  return (rng()() & 1) ? c : -c;
}

}  // namespace helmer

namespace helmer {
Poly Poly::compose(const std::vector<Poly>& images) const {
  if (static_cast<int>(images.size()) != R_->nvars())
    throw std::invalid_argument("compose: need one image per variable");
  if (images.empty()) throw std::invalid_argument("compose: no variables");
  const RingPtr& T = images[0].ring();
  for (auto& g : images)
    if (g.ring() != T) throw std::invalid_argument("compose: images in different rings");
  Poly r(T);
  std::vector<fmpq_mpoly_struct*> v;
  for (auto& g : images) v.push_back(const_cast<fmpq_mpoly_struct*>(g.raw()));
  if (!fmpq_mpoly_compose_fmpq_mpoly(r.raw(), p_, v.data(), R_->ctx(), T->ctx()))
    throw std::runtime_error("compose: exponent overflow");
  return r;
}
}  // namespace helmer
