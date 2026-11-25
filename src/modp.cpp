#include "helmer/modp.hpp"

#include "helmer/ideal.hpp"

#include <flint/fmpq.h>
#include <flint/fmpz.h>
#include <flint/ulong_extras.h>

#include <cstdlib>
#include <functional>
#include <unordered_map>

extern "C" {
// msolve's neogb (see gb.cpp).
int64_t export_f4(void* (*mallocp)(size_t), int32_t* bld, int32_t** blen, int32_t** bexp,
                  void** bcf, const int32_t* lens, const int32_t* exps, const void* cfs,
                  const uint32_t field_char, const int32_t mon_order, const int32_t elim_block_len,
                  const int32_t nr_vars, const int32_t nr_gens, const int32_t ht_size,
                  const int32_t nr_threads, const int32_t max_nr_pairs,
                  const int32_t reset_hash_table, const int32_t la_option,
                  const int32_t reduce_gb, const int32_t pbm_file, const int32_t info_level);
void free_f4_julia_result_data(void (*freep)(void*), int32_t** blen, int32_t** bexp, void** bcf,
                               const int64_t ngens, const int64_t field_char);
}

namespace helmer::modp {

// ------------------------------------------------------------------ rings and polys

Ctx::Ctx(ulong p_, int n_) : p(p_), n(n_) { nmod_mpoly_ctx_init(ctx, n > 0 ? n : 1, ORD_DEGREVLEX, p); }
Ctx::~Ctx() { nmod_mpoly_ctx_clear(ctx); }
CtxPtr make_ctx(ulong p, int n) { return std::make_shared<Ctx>(p, n); }

MPoly::MPoly(CtxPtr C) : C_(std::move(C)), p_(new nmod_mpoly_struct) { nmod_mpoly_init(p_, C_->ctx); }
MPoly::MPoly(const MPoly& o) : C_(o.C_), p_(new nmod_mpoly_struct) {
  nmod_mpoly_init(p_, C_->ctx);
  nmod_mpoly_set(p_, o.p_, C_->ctx);
}
MPoly::MPoly(MPoly&& o) noexcept : C_(std::move(o.C_)), p_(o.p_) { o.p_ = nullptr; }
MPoly& MPoly::operator=(const MPoly& o) {
  if (this != &o) {
    MPoly t(o);
    std::swap(C_, t.C_);
    std::swap(p_, t.p_);
  }
  return *this;
}
MPoly& MPoly::operator=(MPoly&& o) noexcept {
  std::swap(C_, o.C_);
  std::swap(p_, o.p_);
  return *this;
}
MPoly::~MPoly() {
  if (p_) {
    nmod_mpoly_clear(p_, C_->ctx);
    delete p_;
  }
}

MPoly MPoly::constant(const CtxPtr& C, ulong c) {
  MPoly r(C);
  nmod_mpoly_set_ui(r.p_, c % C->p, C->ctx);
  return r;
}
MPoly MPoly::var(const CtxPtr& C, int i) {
  MPoly r(C);
  nmod_mpoly_gen(r.p_, i, C->ctx);
  return r;
}

bool MPoly::is_zero() const { return nmod_mpoly_is_zero(p_, C_->ctx); }
long MPoly::nterms() const { return nmod_mpoly_length(p_, C_->ctx); }
long MPoly::total_degree() const { return nmod_mpoly_total_degree_si(p_, C_->ctx); }
long MPoly::degree_in(int i) const { return nmod_mpoly_degree_si(p_, i, C_->ctx); }
ulong MPoly::constant_term() const {
  std::vector<ulong> z(C_->n, 0);
  return eval(z);
}
ulong MPoly::eval(const std::vector<ulong>& pt) const {
  std::vector<ulong> v(pt);
  v.resize(std::max(C_->n, 1), 0);
  return nmod_mpoly_evaluate_all_ui(p_, v.data(), C_->ctx);
}

MPoly MPoly::operator+(const MPoly& o) const {
  MPoly r(C_);
  nmod_mpoly_add(r.p_, p_, o.p_, C_->ctx);
  return r;
}
MPoly MPoly::operator-(const MPoly& o) const {
  MPoly r(C_);
  nmod_mpoly_sub(r.p_, p_, o.p_, C_->ctx);
  return r;
}
MPoly MPoly::operator*(const MPoly& o) const {
  MPoly r(C_);
  nmod_mpoly_mul(r.p_, p_, o.p_, C_->ctx);
  return r;
}
MPoly MPoly::scale(ulong c) const {
  MPoly r(C_);
  nmod_mpoly_scalar_mul_ui(r.p_, p_, c % C_->p, C_->ctx);
  return r;
}
MPoly MPoly::derivative(int i) const {
  MPoly r(C_);
  nmod_mpoly_derivative(r.p_, p_, i, C_->ctx);
  return r;
}
MPoly MPoly::compose(const CtxPtr& target, const std::vector<MPoly>& vals) const {
  MPoly r(target);
  if (static_cast<int>(vals.size()) != C_->n) throw std::logic_error("modp::compose: arity");
  if (C_->n == 0) {
    nmod_mpoly_set_ui(r.p_, constant_term(), target->ctx);
    return r;
  }
  std::vector<nmod_mpoly_struct*> v(vals.size());
  for (size_t i = 0; i < vals.size(); ++i) v[i] = const_cast<nmod_mpoly_struct*>(vals[i].raw());
  if (!nmod_mpoly_compose_nmod_mpoly(r.p_, p_, v.data(), C_->ctx, target->ctx))
    throw std::runtime_error("modp::compose: exponent overflow");
  return r;
}

ulong random_prime() {
  // Primes in [2^30, 2^31): msolve's F4 takes characteristics below 2^31.
  std::uniform_int_distribution<ulong> d(1UL << 30, (1UL << 31) - 1);
  ulong p;
  do { p = d(gb_prime_rng()); } while (!n_is_prime(p));
  return p;
}
ulong random_elem(ulong p) { return std::uniform_int_distribution<ulong>(0, p - 1)(rng()); }
ulong random_unit(ulong p) { return std::uniform_int_distribution<ulong>(1, p - 1)(rng()); }

MPoly reduce(const Poly& f, const CtxPtr& C) {
  const auto* ctx = f.ring()->ctx();
  const slong len = fmpq_mpoly_length(f.raw(), ctx);
  const int n = f.ring()->nvars();
  if (n != C->n) throw std::logic_error("modp::reduce: arity");
  MPoly r(C);
  fmpq_t c;
  fmpq_init(c);
  std::vector<ulong> e(std::max(n, 1));
  for (slong i = 0; i < len; ++i) {
    fmpq_mpoly_get_term_coeff_fmpq(c, f.raw(), i, ctx);
    ulong den = fmpz_fdiv_ui(fmpq_denref(c), C->p);
    if (den == 0) {
      fmpq_clear(c);
      throw Unlucky();
    }
    ulong num = fmpz_fdiv_ui(fmpq_numref(c), C->p);
    fmpq_mpoly_get_term_exp_ui(e.data(), f.raw(), i, ctx);
    nmod_mpoly_push_term_ui_ui(r.raw(), n_mulmod2_preinv(num, n_invmod(den, C->p), C->p,
                                                         n_preinvert_limb(C->p)),
                               e.data(), C->ctx);
  }
  fmpq_clear(c);
  nmod_mpoly_sort_terms(r.raw(), C->ctx);
  nmod_mpoly_combine_like_terms(r.raw(), C->ctx);
  return r;
}

// ------------------------------------------------------------------ Gröbner bases

ModGB groebner(const CtxPtr& C, const std::vector<MPoly>& F, int ebl) {
  ModGB out;
  const int n = C->n;
  std::vector<int32_t> lens, exps, cfs;
  std::vector<ulong> e(std::max(n, 1));
  for (auto& f : F) {
    if (f.is_zero()) continue;
    if (f.total_degree() == 0) {  // nonzero constant
      out.unit = true;
      out.polys = {MPoly::constant(C, 1)};
      out.lead = {std::vector<int>(n, 0)};
      return out;
    }
    const slong len = f.nterms();
    lens.push_back(static_cast<int32_t>(len));
    for (slong i = 0; i < len; ++i) {
      nmod_mpoly_get_term_exp_ui(e.data(), f.raw(), i, C->ctx);
      for (int v = 0; v < n; ++v) exps.push_back(static_cast<int32_t>(e[v]));
      cfs.push_back(static_cast<int32_t>(nmod_mpoly_get_term_coeff_ui(f.raw(), i, C->ctx)));
    }
  }
  if (lens.empty()) return out;  // zero ideal
  if (n == 0) throw std::logic_error("modp::groebner: no variables");
  int32_t bld = 0;
  int32_t* blen = nullptr;
  int32_t* bexp = nullptr;
  void* bcf = nullptr;
  int64_t nt = export_f4(std::malloc, &bld, &blen, &bexp, &bcf, lens.data(), exps.data(), cfs.data(),
                         static_cast<uint32_t>(C->p), 0, ebl, n, static_cast<int32_t>(lens.size()),
                         12, 1, 0, 0, 2, 1, 0, 0);
  gb_stats().f4_runs++;
  if (nt <= 0 || bld <= 0) {
    if (blen || bexp || bcf) free_f4_julia_result_data(std::free, &blen, &bexp, &bcf, bld, C->p);
    throw std::runtime_error("modp::groebner: msolve failed");
  }
  const int32_t* c = static_cast<const int32_t*>(bcf);
  int64_t off = 0;
  for (int32_t g = 0; g < bld; ++g) {
    MPoly p(C);
    for (int32_t k = 0; k < blen[g]; ++k) {
      for (int v = 0; v < n; ++v) e[v] = static_cast<ulong>(bexp[(off + k) * n + v]);
      if (k == 0) out.lead.emplace_back(e.begin(), e.begin() + n);
      nmod_mpoly_push_term_ui_ui(p.raw(), static_cast<ulong>(c[off + k]), e.data(), C->ctx);
    }
    nmod_mpoly_sort_terms(p.raw(), C->ctx);
    nmod_mpoly_combine_like_terms(p.raw(), C->ctx);
    off += blen[g];
    out.polys.push_back(std::move(p));
  }
  free_f4_julia_result_data(std::free, &blen, &bexp, &bcf, bld, C->p);
  out.unit = out.polys.size() == 1 && out.polys[0].total_degree() == 0;
  return out;
}

int dim(const ModGB& G, int n) {
  if (G.unit) return -1;
  return static_cast<int>(max_independent_set(n, G.lead).size());
}

long staircase(const ModGB& G, int n) {
  if (G.unit) return 0;
  // zero-dimensional iff every variable has a pure power among the leads
  for (int v = 0; v < n; ++v) {
    bool pure = false;
    for (auto& l : G.lead) {
      bool ok = l[v] > 0;
      for (int u = 0; u < n && ok; ++u)
        if (u != v && l[u] != 0) ok = false;
      if (ok) { pure = true; break; }
    }
    if (!pure) return -1;
  }
  // Each standard monomial m != 1 has the standard parent m / x_{max var of m}:
  // enumerate the staircase as a tree.
  long count = 0;
  std::vector<int> m(n, 0);
  auto standard = [&]() {
    for (auto& l : G.lead) {
      bool div = true;
      for (int u = 0; u < n; ++u)
        if (l[u] > m[u]) { div = false; break; }
      if (div) return false;
    }
    return true;
  };
  std::function<void(int)> rec = [&](int maxvar) {
    ++count;
    if (count > 50000000) throw std::runtime_error("modp::staircase: too large");
    for (int v = maxvar; v < n; ++v) {
      m[v]++;
      if (standard()) rec(v);
      m[v]--;
    }
  };
  rec(0);
  return count;
}

long length(const CtxPtr& C, const std::vector<MPoly>& F) { return staircase(groebner(C, F), C->n); }

namespace {
// Images of z_0..z_{n-1} in F_p[T, z_0..z_{n-1}] (T first, the eliminated block).
std::vector<MPoly> shift_vars(const CtxPtr& CT, int n) {
  std::vector<MPoly> v;
  for (int i = 0; i < n; ++i) v.push_back(MPoly::var(CT, i + 1));
  return v;
}
std::vector<MPoly> rabinowitsch(const CtxPtr& C, const CtxPtr& CT, const std::vector<MPoly>& F,
                                const MPoly& h) {
  auto up = shift_vars(CT, C->n);
  std::vector<MPoly> G;
  for (auto& f : F) G.push_back(f.compose(CT, up));
  G.push_back(MPoly::constant(CT, 1) - MPoly::var(CT, 0) * h.compose(CT, up));
  return G;
}
}  // namespace

long length_off(const CtxPtr& C, const std::vector<MPoly>& F, const MPoly& h) {
  auto CT = make_ctx(C->p, C->n + 1);
  return length(CT, rabinowitsch(C, CT, F, h));
}

std::vector<MPoly> saturate(const CtxPtr& C, const std::vector<MPoly>& F, const MPoly& h) {
  auto CT = make_ctx(C->p, C->n + 1);
  auto G = groebner(CT, rabinowitsch(C, CT, F, h), /*ebl*/ 1);
  std::vector<MPoly> down{MPoly::constant(C, 0)};
  for (int i = 0; i < C->n; ++i) down.push_back(MPoly::var(C, i));
  std::vector<MPoly> out;
  for (auto& g : G.polys)
    if (g.degree_in(0) <= 0) out.push_back(g.compose(C, down));
  return out;
}

bool in_radical(const CtxPtr& C, const std::vector<MPoly>& F, const MPoly& h) {
  auto CT = make_ctx(C->p, C->n + 1);
  return groebner(CT, rabinowitsch(C, CT, F, h)).unit;
}

std::vector<MPoly> random_affine_map(const CtxPtr& Cz, const std::vector<ulong>& x, int n) {
  std::vector<MPoly> img;
  for (int i = 0; i < n; ++i) {
    MPoly f = MPoly::constant(Cz, x[i]);
    for (int j = 0; j < Cz->n; ++j) f = f + MPoly::var(Cz, j).scale(random_elem(Cz->p));
    img.push_back(std::move(f));
  }
  return img;
}

std::vector<MPoly> maximal_minors(const CtxPtr& C, const std::vector<std::vector<MPoly>>& M) {
  const int k = static_cast<int>(M.size());
  if (k == 0) return {MPoly::constant(C, 1)};
  const int n = static_cast<int>(M[0].size());
  if (k > n) return {};
  if (n > 30) throw std::runtime_error("modp::maximal_minors: too many columns");
  // Laplace expansion along the last row, building minors of the first r rows on
  // every r-subset of columns.
  std::unordered_map<uint32_t, MPoly> cur;
  cur.emplace(0u, MPoly::constant(C, 1));
  for (int r = 0; r < k; ++r) {
    std::unordered_map<uint32_t, MPoly> next;
    for (auto& [mask, det] : cur) {
      if (det.is_zero()) continue;
      for (int j = 0; j < n; ++j) {
        if (mask & (1u << j) || M[r][j].is_zero()) continue;
        const uint32_t m2 = mask | (1u << j);
        // sign: position of column j in m2, plus row index r
        const int pos = __builtin_popcount(m2 & ((1u << j) - 1));
        MPoly term = M[r][j] * det;
        if ((pos + r) & 1) term = term.scale(C->p - 1);
        auto it = next.find(m2);
        if (it == next.end()) next.emplace(m2, std::move(term));
        else it->second = it->second + term;
      }
    }
    cur = std::move(next);
  }
  std::vector<MPoly> out;
  for (auto& [mask, det] : cur)
    if (!det.is_zero()) out.push_back(det);
  return out;
}

}  // namespace helmer::modp
