#include "helmer/certify.hpp"

#include <flint/fmpq.h>
#include <flint/fmpz.h>
#include <flint/nmod.h>
#include <flint/nmod_mat.h>
#include <flint/ulong_extras.h>

#include <algorithm>
#include <map>
#include <tuple>
#include <unordered_map>

namespace helmer {

namespace {

using Exp = std::vector<long>;

struct ExpHash {
  size_t operator()(const Exp& m) const {
    size_t h = 1469598103934665603ULL;
    for (long e : m) h = (h ^ static_cast<size_t>(e)) * 1099511628211ULL;
    return h;
  }
};

struct Zv {  // RAII fmpz
  fmpz v;
  Zv() { fmpz_init(&v); }
  Zv(const Zv& o) { fmpz_init(&v); fmpz_set(&v, &o.v); }
  Zv& operator=(const Zv& o) { fmpz_set(&v, &o.v); return *this; }
  ~Zv() { fmpz_clear(&v); }
};
struct Qv {  // RAII fmpq
  fmpq v;
  Qv() { fmpq_init(&v); }
  Qv(const Qv& o) { fmpq_init(&v); fmpq_set(&v, &o.v); }
  Qv& operator=(const Qv& o) { fmpq_set(&v, &o.v); return *this; }
  ~Qv() { fmpq_clear(&v); }
};

struct Terms {
  std::vector<Exp> exps;
  std::vector<Qv> coeffs;
};

Terms terms_of(const Poly& p) {
  Terms T;
  const long n = p.nterms();
  T.exps.reserve(n);
  T.coeffs.resize(n);
  for (long i = 0; i < n; ++i) {
    T.exps.push_back(p.term_exp(i));
    fmpq_mpoly_get_term_coeff_fmpq(&T.coeffs[i].v, p.raw(), i, p.ring()->ctx());
  }
  return T;
}

long block_degree(const Exp& e, const std::vector<char>& in_block, bool block) {
  long d = 0;
  for (size_t i = 0; i < e.size(); ++i)
    if (static_cast<bool>(in_block[i]) == block) d += e[i];
  return d;
}

long poly_block_degree(const Poly& p, const std::vector<char>& in_block, bool block) {
  long d = 0;
  for (long i = 0; i < p.nterms(); ++i) d = std::max(d, block_degree(p.term_exp(i), in_block, block));
  return d;
}

// All exponent vectors in the variables `vars` (of n) with total degree <= D.
void monomials(const std::vector<int>& vars, int n, int D, std::vector<Exp>& out) {
  out.assign(1, Exp(n, 0));
  for (int v : vars) {
    std::vector<Exp> next;
    for (const Exp& e : out) {
      long used = 0;
      for (int w : vars) used += e[w];
      for (long k = 0; used + k <= D; ++k) {
        Exp f = e;
        f[v] = k;
        next.push_back(std::move(f));
      }
    }
    out.swap(next);
  }
}

long support_size(int nblock, int nrest, int DA, int DB) {
  // C(DA + nblock, nblock) * C(DB + nrest, nrest), saturating
  auto binom = [](long a, long b) {
    long r = 1;
    for (long i = 1; i <= b; ++i) {
      r = r * (a - b + i) / i;
      if (r > (1L << 40)) return 1L << 40;
    }
    return r;
  };
  long s = binom(DA + nblock, nblock) * binom(DB + nrest, nrest);
  return std::min(s, 1L << 40);
}

ulong residue(const fmpq* c, ulong p, bool& bad) {
  ulong den = fmpz_fdiv_ui(fmpq_denref(c), p);
  if (den == 0) {
    bad = true;
    return 0;
  }
  ulong num = fmpz_fdiv_ui(fmpq_numref(c), p);
  nmod_t mod;
  nmod_init(&mod, p);
  return nmod_mul(num, n_invmod(den, p), mod);
}

}  // namespace

bool verify(const MembershipCertificate& C) {
  if (C.cofactors.size() != C.gens.size()) return false;
  const RingPtr& R = C.target.ring();
  Poly sum(R);
  for (size_t j = 0; j < C.gens.size(); ++j) {
    if (C.gens[j].ring() != R || C.cofactors[j].ring() != R) return false;
    sum = sum + C.cofactors[j] * C.gens[j];
  }
  return sum.equal(C.target);
}

std::optional<MembershipCertificate> certify_membership_bounded(
    const std::vector<Poly>& gens, const Poly& h, const std::vector<int>& block, int DA, int DB,
    int max_primes) {
  const RingPtr& R = h.ring();
  const int n = R->nvars();
  std::vector<char> in_block(n, 0);
  for (int v : block) in_block[v] = 1;
  std::vector<int> bvars, rvars;
  for (int v = 0; v < n; ++v) (in_block[v] ? bvars : rvars).push_back(v);

  // Support of every cofactor.
  std::vector<Exp> A, B, supp;
  monomials(bvars, n, DA, A);
  monomials(rvars, n, DB, B);
  for (const Exp& a : A)
    for (const Exp& b : B) {
      Exp e(n);
      for (int i = 0; i < n; ++i) e[i] = a[i] + b[i];
      supp.push_back(std::move(e));
    }
  const long m = static_cast<long>(gens.size());
  const long S = static_cast<long>(supp.size());
  const long ncols = m * S;

  // Rows: monomials of the products and of h. Entries reference (generator, term).
  std::vector<Terms> G;
  for (auto& g : gens) G.push_back(terms_of(g));
  Terms H = terms_of(h);
  std::unordered_map<Exp, long, ExpHash> row;
  struct Entry {
    long r, c;
    int j;
    long t;
  };
  std::vector<Entry> entries;
  for (long j = 0; j < m; ++j)
    for (long s = 0; s < S; ++s)
      for (long t = 0; t < static_cast<long>(G[j].exps.size()); ++t) {
        Exp e(n);
        for (int i = 0; i < n; ++i) e[i] = supp[s][i] + G[j].exps[t][i];
        auto it = row.emplace(std::move(e), static_cast<long>(row.size())).first;
        entries.push_back({it->second, j * S + s, static_cast<int>(j), t});
      }
  std::vector<long> hrow;
  for (auto& e : H.exps) hrow.push_back(row.emplace(e, static_cast<long>(row.size())).first->second);
  const long nrows = static_cast<long>(row.size());

  // Modular solves; keep primes with the best (rank, pivot pattern).
  long best_rank = -1;
  std::vector<long> best_piv;
  std::vector<Zv> crt;  // solution at pivot columns, CRT-combined
  Zv modulus;
  fmpz_one(&modulus.v);
  std::vector<Qv> last;  // previous successful reconstruction
  bool have_last = false;
  int inconsistent = 0, consistent = 0;
  flint_rand_t st;
  flint_randinit(st);
  ulong p = n_nextprime((UWORD(1) << 62) + n_randint(st, UWORD(1) << 40), 1);
  flint_randclear(st);

  for (int round = 0; round < max_primes; ++round, p = n_nextprime(p, 1)) {
    bool bad = false;
    nmod_mat_t M;
    nmod_mat_init(M, nrows, ncols + 1, p);
    for (const Entry& E : entries) {
      ulong c = residue(&G[E.j].coeffs[E.t].v, p, bad);
      nmod_mat_entry(M, E.r, E.c) = nmod_add(nmod_mat_entry(M, E.r, E.c), c, M->mod);
    }
    for (size_t t = 0; t < hrow.size(); ++t)
      nmod_mat_entry(M, hrow[t], ncols) = residue(&H.coeffs[t].v, p, bad);
    if (bad) {
      nmod_mat_clear(M);
      continue;
    }
    const long rank = nmod_mat_rref(M);
    std::vector<long> piv;
    bool incons = false;
    for (long i = 0; i < rank; ++i) {
      long c = 0;
      while (c <= ncols && nmod_mat_entry(M, i, c) == 0) ++c;
      if (c == ncols) incons = true;
      piv.push_back(c);
    }
    if (incons) {
      nmod_mat_clear(M);
      // Over QQ a consistent system stays consistent modulo every prime not dividing
      // a denominator of a solution; two misses with no hit: treat as no solution.
      if (++inconsistent >= 2 && consistent == 0) return std::nullopt;
      continue;
    }
    ++consistent;
    std::vector<ulong> sol(rank);
    for (long i = 0; i < rank; ++i) {
      ulong pv = nmod_mat_entry(M, i, piv[i]);
      sol[i] = nmod_mul(nmod_mat_entry(M, i, ncols), n_invmod(pv, p), M->mod);
    }
    nmod_mat_clear(M);

    if (rank > best_rank || (rank == best_rank && piv < best_piv)) {
      best_rank = rank;
      best_piv = piv;
      crt.assign(rank, Zv());
      for (long i = 0; i < rank; ++i) fmpz_set_ui(&crt[i].v, sol[i]);
      fmpz_set_ui(&modulus.v, p);
      have_last = false;
    } else if (rank == best_rank && piv == best_piv) {
      Zv pz;
      fmpz_set_ui(&pz.v, p);
      for (long i = 0; i < rank; ++i) {
        Zv r1;
        fmpz_set(&r1.v, &crt[i].v);
        fmpz_CRT_ui(&crt[i].v, &r1.v, &modulus.v, sol[i], p, 0);
      }
      fmpz_mul(&modulus.v, &modulus.v, &pz.v);
    } else {
      continue;  // unlucky prime
    }

    // Rational reconstruction; check exactly once it is stable.
    std::vector<Qv> q(best_rank);
    bool ok = true;
    for (long i = 0; i < best_rank && ok; ++i)
      ok = fmpq_reconstruct_fmpz(&q[i].v, &crt[i].v, &modulus.v);
    if (!ok) continue;
    bool stable = have_last;
    for (long i = 0; i < best_rank && stable; ++i) stable = fmpq_equal(&q[i].v, &last[i].v);
    last = q;
    have_last = true;
    if (!stable) continue;

    std::vector<std::vector<Exp>> ce(m);
    std::vector<std::vector<fmpq*>> cc(m);
    for (long i = 0; i < best_rank; ++i) {
      const long col = best_piv[i];
      if (fmpq_is_zero(&q[i].v)) continue;
      ce[col / S].push_back(supp[col % S]);
      cc[col / S].push_back(&q[i].v);
    }
    MembershipCertificate C{h, gens, {}, DA, DB};
    for (long j = 0; j < m; ++j) C.cofactors.push_back(Poly::from_terms(R, ce[j], cc[j]));
    if (verify(C)) return C;
  }
  return std::nullopt;
}

std::optional<MembershipCertificate> certify_membership(const std::vector<Poly>& gens,
                                                        const Poly& h,
                                                        const std::vector<int>& block,
                                                        const CertifyOptions& opt) {
  const RingPtr& R = h.ring();
  const int n = R->nvars();
  if (gens.empty()) return std::nullopt;
  if (h.is_zero()) {
    MembershipCertificate C{h, gens, {}, 0, 0};
    for (size_t j = 0; j < gens.size(); ++j) C.cofactors.emplace_back(R);
    return C;
  }
  std::vector<char> in_block(n, 0);
  for (int v : block) in_block[v] = 1;
  const int nb = static_cast<int>(block.size()), nr = n - nb;
  long dA = 0;
  for (auto& g : gens) dA = std::max(dA, poly_block_degree(g, in_block, true));
  const long dB = poly_block_degree(h, in_block, false);
  std::vector<std::tuple<long, int, int>> tries;  // (unknowns, DA, DB)
  for (long DA = std::max(0L, dA - 2); DA <= dA + opt.extra_block; ++DA)
    for (long DB = std::max(0L, dB - 1); DB <= dB + opt.extra_rest; ++DB) {
      const long u = support_size(nb, nr, static_cast<int>(DA), static_cast<int>(DB)) *
                     static_cast<long>(gens.size());
      if (u <= opt.max_unknowns) tries.emplace_back(u, static_cast<int>(DA), static_cast<int>(DB));
    }
  std::sort(tries.begin(), tries.end());
  for (auto& [u, DA, DB] : tries) {
    auto C = certify_membership_bounded(gens, h, block, DA, DB, opt.max_primes);
    if (C) return C;
  }
  return std::nullopt;
}

CertifiedElimination eliminate_certified(const Ideal& I, const std::vector<int>& elim,
                                         const CertifyOptions& opt) {
  CertifiedElimination out{eliminate(I, elim), {}};
  for (const Poly& P : out.ideal.gens())
    out.certificates.push_back(certify_membership(I.gens(), P, elim, opt));
  return out;
}

}  // namespace helmer
