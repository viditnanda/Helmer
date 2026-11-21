#include "helmer/gb.hpp"

#include <flint/fmpq.h>
#include <flint/fmpz.h>
#include <flint/ulong_extras.h>

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <unordered_map>

extern "C" {
// From msolve's neogb (libneogb). Declared here to avoid pulling in its headers.
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

namespace helmer {

GBStats& gb_stats() {
  static GBStats s;
  return s;
}
GBOptions& gb_options() {
  static GBOptions o;
  return o;
}

bool GB::is_unit() const { return polys.size() == 1 && polys[0].is_constant() && !polys[0].is_zero(); }

namespace {

struct Z {  // RAII fmpz
  fmpz v;
  Z() { fmpz_init(&v); }
  Z(const Z& o) { fmpz_init(&v); fmpz_set(&v, &o.v); }
  Z& operator=(const Z& o) { fmpz_set(&v, &o.v); return *this; }
  ~Z() { fmpz_clear(&v); }
};
struct Q {  // RAII fmpq
  fmpq v;
  Q() { fmpq_init(&v); }
  Q(const Q& o) { fmpq_init(&v); fmpq_set(&v, &o.v); }
  Q& operator=(const Q& o) { fmpq_set(&v, &o.v); return *this; }
  ~Q() { fmpq_clear(&v); }
};

using Mono = std::vector<int32_t>;  // msolve variable order
struct MonoHash {
  size_t operator()(const Mono& m) const {
    size_t h = 1469598103934665603ULL;
    for (int32_t e : m) h = (h ^ static_cast<size_t>(e)) * 1099511628211ULL;
    return h;
  }
};

// One modular result: per element, terms (mono, coeff mod p); element 0..k-1, lead first.
struct ModResult {
  std::vector<std::vector<std::pair<Mono, uint32_t>>> elems;
  std::vector<int32_t> signature;  // concatenated lead monomials
};

// Integer input, pre-permuted to msolve order.
struct IntInput {
  int nvars = 0;
  std::vector<int32_t> lens;
  std::vector<int32_t> exps;
  std::vector<Z> coeffs;
};

bool run_f4(const IntInput& in, uint32_t p, int ebl, ModResult& out) {
  const size_t nt = in.coeffs.size();
  std::vector<int32_t> cfs(nt);
  for (size_t i = 0; i < nt; ++i) {
    ulong r = fmpz_fdiv_ui(&in.coeffs[i].v, p);
    cfs[i] = static_cast<int32_t>(r);
  }
  int32_t bld = 0;
  int32_t* blen = nullptr;
  int32_t* bexp = nullptr;
  void* bcf = nullptr;
  const int32_t ngens = static_cast<int32_t>(in.lens.size());
  int64_t nterms = export_f4(std::malloc, &bld, &blen, &bexp, &bcf, in.lens.data(), in.exps.data(),
                             cfs.data(), p, /*mon_order DRL*/ 0, ebl, in.nvars, ngens,
                             /*ht_size*/ 12, gb_options().msolve_threads, /*max_nr_pairs*/ 0,
                             /*reset_ht*/ 0, /*la_option*/ 2, /*reduce_gb*/ 1, /*pbm*/ 0,
                             /*info*/ 0);
  gb_stats().f4_runs++;
  if (nterms <= 0 || bld <= 0) {
    if (blen || bexp || bcf) free_f4_julia_result_data(std::free, &blen, &bexp, &bcf, bld, p);
    return false;
  }
  out.elems.assign(bld, {});
  out.signature.clear();
  const int32_t* c = static_cast<const int32_t*>(bcf);
  int64_t off = 0;
  for (int32_t g = 0; g < bld; ++g) {
    auto& el = out.elems[g];
    el.reserve(blen[g]);
    for (int32_t k = 0; k < blen[g]; ++k) {
      Mono m(bexp + (off + k) * in.nvars, bexp + (off + k + 1) * in.nvars);
      uint32_t cf = static_cast<uint32_t>(c[off + k]);
      if (k == 0) out.signature.insert(out.signature.end(), m.begin(), m.end());
      el.emplace_back(std::move(m), cf);
    }
    out.signature.push_back(-1);
    off += blen[g];
  }
  free_f4_julia_result_data(std::free, &blen, &bexp, &bcf, bld, p);
  return true;
}

// CRT accumulator for one lead-monomial signature.
struct Accum {
  int nprimes = 0;
  Z modulus;
  std::vector<std::vector<Mono>> monos;  // per element
  std::vector<std::unordered_map<Mono, size_t, MonoHash>> index;
  std::vector<std::vector<Z>> res;
  bool has_candidate = false;
  std::vector<std::vector<Q>> cand;  // reconstruction candidate, aligned with monos
  int next_check = 1;

  void add(const ModResult& r, uint32_t p) {
    if (nprimes == 0) {
      fmpz_one(&modulus.v);
      monos.assign(r.elems.size(), {});
      index.assign(r.elems.size(), {});
      res.assign(r.elems.size(), {});
    }
    Z tmp;
    for (size_t g = 0; g < r.elems.size(); ++g) {
      std::unordered_map<Mono, uint32_t, MonoHash> cur;
      for (auto& [m, c] : r.elems[g]) {
        cur[m] = c;
        if (!index[g].count(m)) {  // new monomial: previous residue is 0
          index[g][m] = monos[g].size();
          monos[g].push_back(m);
          res[g].emplace_back();
        }
      }
      for (size_t j = 0; j < monos[g].size(); ++j) {
        auto it = cur.find(monos[g][j]);
        ulong r2 = it == cur.end() ? 0 : it->second;
        if (nprimes == 0) {
          fmpz_set_ui(&res[g][j].v, r2);
        } else {
          fmpz_CRT_ui(&tmp.v, &res[g][j].v, &modulus.v, r2, p, 0);
          fmpz_set(&res[g][j].v, &tmp.v);
        }
      }
    }
    fmpz_mul_ui(&modulus.v, &modulus.v, p);
    nprimes++;
  }

  bool reconstruct() {
    std::vector<std::vector<Q>> c(monos.size());
    for (size_t g = 0; g < monos.size(); ++g) {
      c[g].resize(monos[g].size());
      for (size_t j = 0; j < monos[g].size(); ++j)
        if (!fmpq_reconstruct_fmpz(&c[g][j].v, &res[g][j].v, &modulus.v)) return false;
    }
    cand = std::move(c);
    has_candidate = true;
    return true;
  }

  // Does the candidate reduce to r modulo p?
  bool confirms(const ModResult& r, uint32_t p) const {
    if (r.elems.size() != monos.size()) return false;
    for (size_t g = 0; g < monos.size(); ++g) {
      std::unordered_map<Mono, uint32_t, MonoHash> cur;
      for (auto& [m, cf] : r.elems[g]) cur[m] = cf;
      size_t matched = 0;
      for (size_t j = 0; j < monos[g].size(); ++j) {
        const fmpq* q = &cand[g][j].v;
        ulong den = fmpz_fdiv_ui(fmpq_denref(q), p);
        if (den == 0) return false;
        ulong num = fmpz_fdiv_ui(fmpq_numref(q), p);
        ulong val = n_mulmod2_preinv(num, n_invmod(den, p), p, n_preinvert_limb(p));
        auto it = cur.find(monos[g][j]);
        ulong got = 0;
        if (it != cur.end()) {
          got = it->second;
          matched++;
        }
        if (got != val) return false;
      }
      if (matched != cur.size()) return false;  // r has a monomial the candidate lacks
    }
    return true;
  }
};

}  // namespace

GB groebner(const RingPtr& R, const std::vector<Poly>& F, const std::vector<int>& perm, int ebl) {
  auto t0 = std::chrono::steady_clock::now();
  gb_stats().calls++;
  const int n = R->nvars();
  if (static_cast<int>(perm.size()) != n) throw std::invalid_argument("helmer: bad permutation");
  if (ebl < 0 || ebl >= n) ebl = 0;  // a block covering everything is plain DRL

  GB out;
  IntInput in;
  in.nvars = n;
  for (const Poly& f0 : F) {
    if (f0.ring().get() != R.get()) throw std::invalid_argument("helmer: generator in wrong ring");
    if (f0.is_zero()) continue;
    if (f0.is_constant()) {  // unit ideal
      out.polys = {Poly::constant(R, 1)};
      out.lead = {std::vector<int>(n, 0)};
      gb_stats().seconds +=
          std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
      return out;
    }
    Poly f = f0;  // fmpq_mpoly = content * zpoly with zpoly primitive integral
    fmpq_mpoly_struct* fr = f.raw();
    const long len = f.nterms();
    in.lens.push_back(static_cast<int32_t>(len));
    for (long i = 0; i < len; ++i) {
      std::vector<long> e = f.term_exp(i);
      for (int k = 0; k < n; ++k) in.exps.push_back(static_cast<int32_t>(e[perm[k]]));
      in.coeffs.emplace_back();
      fmpz_set(&in.coeffs.back().v, fmpq_mpoly_zpoly_term_coeff_ref(fr, i, R->ctx()));
    }
  }
  if (in.lens.empty()) {  // zero ideal
    gb_stats().seconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    return out;
  }

  std::map<std::vector<int32_t>, Accum> groups;
  ulong p = (1UL << 31);
  const Accum* winner = nullptr;
  for (int it = 0; it < gb_options().max_primes && !winner; ++it) {
    do { --p; } while (!n_is_prime(p));
    bool bad = false;  // p must not kill any input coefficient
    for (auto& c : in.coeffs)
      if (fmpz_fdiv_ui(&c.v, p) == 0) { bad = true; break; }
    if (bad) continue;
    ModResult r;
    if (!run_f4(in, static_cast<uint32_t>(p), ebl, r)) { gb_stats().unlucky++; continue; }
    Accum& A = groups[r.signature];
    if (A.has_candidate && A.confirms(r, static_cast<uint32_t>(p))) {
      winner = &A;
      break;
    }
    A.add(r, static_cast<uint32_t>(p));
    A.has_candidate = false;
    if (A.nprimes >= A.next_check) {
      A.next_check = A.nprimes + std::max(1, A.nprimes / 4);
      A.reconstruct();
    }
  }
  if (!winner) throw std::runtime_error("helmer: modular Gröbner basis did not stabilise");
  for (auto& [sig, A] : groups)
    if (&A != winner) gb_stats().unlucky += A.nprimes;

  // Build polynomials in ring order.
  const Accum& A = *winner;
  for (size_t g = 0; g < A.monos.size(); ++g) {
    std::vector<std::vector<long>> exps;
    std::vector<fmpq*> cs;
    for (size_t j = 0; j < A.monos[g].size(); ++j) {
      if (fmpq_is_zero(&A.cand[g][j].v)) continue;
      std::vector<long> e(n);
      for (int k = 0; k < n; ++k) e[perm[k]] = A.monos[g][j][k];
      exps.push_back(std::move(e));
      cs.push_back(const_cast<fmpq*>(&A.cand[g][j].v));
    }
    out.polys.push_back(Poly::from_terms(R, exps, cs));
    std::vector<int> lead(n);
    for (int k = 0; k < n; ++k) lead[perm[k]] = A.monos[g][0][k];  // lead term stored first
    out.lead.push_back(std::move(lead));
  }
  double el = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  gb_stats().seconds += el;
  if (std::getenv("HELMER_TRACE"))
    fprintf(stderr, "[gb] nv=%d ngens=%zu ebl=%d -> %zu elems, primes=%d, %.3fs\n", n, in.lens.size(), ebl,
            out.polys.size(), winner->nprimes + 1, el);
  return out;
}

GB groebner_drl(const RingPtr& R, const std::vector<Poly>& F) {
  std::vector<int> id(R->nvars());
  for (int i = 0; i < R->nvars(); ++i) id[i] = i;
  return groebner(R, F, id, 0);
}

}  // namespace helmer
