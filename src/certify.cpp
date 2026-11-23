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
  // C(DA + nblock, nblock) * C(DB + nrest, nrest), saturating at 2^40
  constexpr long cap = 1L << 40;
  auto binom = [](long a, long b) {
    __int128 r = 1;
    for (long i = 1; i <= b; ++i) {
      r = r * (a - b + i) / i;
      if (r > cap) return cap;
    }
    return static_cast<long>(r);
  };
  const __int128 s = static_cast<__int128>(binom(DA + nblock, nblock)) * binom(DB + nrest, nrest);
  return s > cap ? cap : static_cast<long>(s);
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

bool verify(const Ideal& I, const MembershipCertificate& C) {
  const auto& G = I.gens();
  if (G.size() != C.gens.size()) return false;
  for (size_t j = 0; j < G.size(); ++j)
    if (G[j].ring() != C.gens[j].ring() || !G[j].equal(C.gens[j])) return false;
  return verify(C);
}

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

static std::optional<MembershipCertificate> solve_bounded(
    const std::vector<Poly>& gens, const Poly& h, const std::vector<int>& block, int DA, int DB,
    int max_primes, int first_attempt, int end_attempt) {
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

  // Integer form: scaling the columns of generator j by d_j (lcm of the denominators
  // of g_j) and h by d_h gives an integer system A' y = b' with x = D y / d_h.
  std::vector<Zv> dg(m);
  for (long j = 0; j < m; ++j) {
    fmpz_one(&dg[j].v);
    for (auto& c : G[j].coeffs) fmpz_lcm(&dg[j].v, &dg[j].v, fmpq_denref(&c.v));
  }
  Zv dh;
  fmpz_one(&dh.v);
  for (auto& c : H.coeffs) fmpz_lcm(&dh.v, &dh.v, fmpq_denref(&c.v));
  std::vector<Zv> ent(entries.size());  // integer entry of A'
  for (size_t k = 0; k < entries.size(); ++k) {
    const Entry& E = entries[k];
    fmpz_divexact(&ent[k].v, &dg[E.j].v, fmpq_denref(&G[E.j].coeffs[E.t].v));
    fmpz_mul(&ent[k].v, &ent[k].v, fmpq_numref(&G[E.j].coeffs[E.t].v));
  }
  std::vector<Zv> hb(H.coeffs.size());  // integer right-hand side b'
  for (size_t t = 0; t < H.coeffs.size(); ++t) {
    fmpz_divexact(&hb[t].v, &dh.v, fmpq_denref(&H.coeffs[t].v));
    fmpz_mul(&hb[t].v, &hb[t].v, fmpq_numref(&H.coeffs[t].v));
  }

  // Lifting budget: p-adic precision of 62 * max_primes bits (as many bits as
  // max_primes 62-bit primes would give).
  const long max_bits = 62L * std::max(1, max_primes);
  flint_rand_t st;
  flint_randinit(st);
  ulong p = n_nextprime((UWORD(1) << 25) + n_randint(st, UWORD(1) << 24), 1);
  flint_randclear(st);

  // 1. One elimination of the full system modulo a ~26-bit prime gives the rank, the
  //    pivot columns (the column rank profile, i.e. the RREF pivots) and r independent
  //    rows. A prime dividing a solution denominator can make a rationally
  //    consistent system inconsistent modulo p. Retry within the prime budget.
  // 2. The r x r subsystem on those rows and columns is nonsingular over QQ; its
  //    solution (free variables zero) is found by Dixon's p-adic lifting from one LU
  //    factorisation modulo p, stopping as soon as rational reconstruction gives
  //    cofactors that pass the exact check.
  // A pivot profile from an unlucky prime can only make that check fail; the whole
  // procedure is then repeated with another prime (three attempts in all).
  for (int i = 0; i < first_attempt; ++i) p = n_nextprime(p, 1);
  for (int attempt = first_attempt; attempt < end_attempt; ++attempt, p = n_nextprime(p, 1)) {
    nmod_t mod;
    nmod_init(&mod, p);
    auto red = [&](const fmpz* z) { return fmpz_fdiv_ui(z, p); };

    long r = 0;
    std::vector<long> piv, sel_rows;
    {
      nmod_mat_t M;
      nmod_mat_init(M, nrows, ncols + 1, p);
      for (size_t k = 0; k < entries.size(); ++k) {
        const Entry& E = entries[k];
        nmod_mat_entry(M, E.r, E.c) = nmod_add(nmod_mat_entry(M, E.r, E.c), red(&ent[k].v), mod);
      }
      for (size_t t = 0; t < hrow.size(); ++t) nmod_mat_entry(M, hrow[t], ncols) = red(&hb[t].v);
      std::vector<slong> P(nrows);
      r = nmod_mat_lu(P.data(), M, 0);
      bool incons = false;
      for (long i = 0; i < r; ++i) {
        long c = i;
        while (c <= ncols && nmod_mat_entry(M, i, c) == 0) ++c;
        if (c == ncols) incons = true;
        piv.push_back(c);
      }
      nmod_mat_clear(M);
      if (incons || r == 0) continue;
      sel_rows.assign(P.begin(), P.begin() + r);
    }

    // Square integer system: sparse entries (i, k, value) and right-hand side.
    std::vector<long> row_pos(nrows, -1), col_pos(ncols, -1);
    for (long i = 0; i < r; ++i) row_pos[sel_rows[i]] = i, col_pos[piv[i]] = i;
    struct SqEntry {
      long i, k;
      const fmpz* v;
    };
    std::vector<SqEntry> sq;
    for (size_t k = 0; k < entries.size(); ++k) {
      const long i = row_pos[entries[k].r], c = col_pos[entries[k].c];
      if (i >= 0 && c >= 0) sq.push_back({i, c, &ent[k].v});
    }
    std::vector<Zv> res(r);  // residual, starts at b'
    for (size_t t = 0; t < hrow.size(); ++t)
      if (row_pos[hrow[t]] >= 0) fmpz_set(&res[row_pos[hrow[t]]].v, &hb[t].v);

    nmod_mat_t LU, Bv, Yv, Zv_;
    nmod_mat_init(LU, r, r, p);
    for (auto& e : sq) nmod_mat_entry(LU, e.i, e.k) = nmod_add(nmod_mat_entry(LU, e.i, e.k), red(e.v), mod);
    std::vector<slong> Pq(r);
    if (nmod_mat_lu(Pq.data(), LU, 1) < r) {  // singular mod p: cannot happen for this p
      nmod_mat_clear(LU);
      continue;
    }
    nmod_mat_init(Bv, r, 1, p);
    nmod_mat_init(Yv, r, 1, p);
    nmod_mat_init(Zv_, r, 1, p);

    std::vector<Zv> y(r);  // p-adic accumulator, in [0, p^k)
    Zv pk;
    fmpz_one(&pk.v);
    std::vector<Qv> q(r);
    long next_try = 2;
    std::optional<MembershipCertificate> found;
    for (long k = 1; fmpz_bits(&pk.v) <= max_bits; ++k) {
      // z = A'^{-1} res mod p, from P A' = L U.
      for (long i = 0; i < r; ++i) nmod_mat_entry(Bv, i, 0) = red(&res[Pq[i]].v);
      nmod_mat_solve_tril(Yv, LU, Bv, 1);
      nmod_mat_solve_triu(Zv_, LU, Yv, 0);
      // y += z p^(k-1); res = (res - A' z) / p.
      for (long i = 0; i < r; ++i) fmpz_addmul_ui(&y[i].v, &pk.v, nmod_mat_entry(Zv_, i, 0));
      for (auto& e : sq) fmpz_submul_ui(&res[e.i].v, e.v, nmod_mat_entry(Zv_, e.k, 0));
      for (long i = 0; i < r; ++i) fmpz_divexact_ui(&res[i].v, &res[i].v, p);
      fmpz_mul_ui(&pk.v, &pk.v, p);

      bool exact = true;  // res == 0: y is the exact (integer) solution
      for (long i = 0; i < r && exact; ++i) exact = fmpz_is_zero(&res[i].v);
      if (!exact && k < next_try) continue;
      next_try = k + std::max(2L, k / 4);

      bool ok = true;
      for (long i = 0; i < r && ok; ++i) {
        if (exact) fmpq_set_fmpz(&q[i].v, &y[i].v);
        else ok = fmpq_reconstruct_fmpz(&q[i].v, &y[i].v, &pk.v);
      }
      if (!ok) continue;
      // x = D y / d_h
      std::vector<std::vector<Exp>> ce(m);
      std::vector<std::vector<fmpq*>> cc(m);
      for (long i = 0; i < r; ++i) {
        const long col = piv[i];
        if (fmpq_is_zero(&q[i].v)) continue;
        fmpq_mul_fmpz(&q[i].v, &q[i].v, &dg[col / S].v);
        fmpq_div_fmpz(&q[i].v, &q[i].v, &dh.v);
        ce[col / S].push_back(supp[col % S]);
        cc[col / S].push_back(&q[i].v);
      }
      MembershipCertificate C{h, gens, {}, DA, DB};
      for (long j = 0; j < m; ++j) C.cofactors.push_back(Poly::from_terms(R, ce[j], cc[j]));
      if (verify(C)) {
        found = std::move(C);
        break;
      }
      if (exact) break;
    }
    nmod_mat_clear(LU);
    nmod_mat_clear(Bv);
    nmod_mat_clear(Yv);
    nmod_mat_clear(Zv_);
    if (found) return found;
  }
  return std::nullopt;
}

std::optional<MembershipCertificate> certify_membership_bounded(
    const std::vector<Poly>& gens, const Poly& h, const std::vector<int>& block, int DA, int DB,
    int max_primes) {
  return solve_bounded(gens, h, block, DA, DB, max_primes, 0, 3);
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
  // Block bounds: from d_A - 2 up to max(d_A, deg_A h) + extra_block, where d_A is
  // the largest block degree of a generator. Other bounds: from d_B(h) - 1 up to
  // max(d_B(h), deg h - e) + extra_rest, where e is the smallest total degree of a
  // generator. For elimination (deg_A h = 0, deg h = d_B(h)) this is the window
  // d_A - 2..d_A + extra_block, d_B(h) - 1..d_B(h) + extra_rest. A bound in variables
  // that do not exist changes nothing, so it is pinned to 0.
  long dA = 0, emin = -1;
  for (auto& g : gens) {
    dA = std::max(dA, poly_block_degree(g, in_block, true));
    if (!g.is_zero()) emin = emin < 0 ? g.total_degree() : std::min(emin, g.total_degree());
  }
  if (emin < 0) return std::nullopt;  // all generators zero; h != 0 is not a member
  const long hA = poly_block_degree(h, in_block, true);
  const long dB = poly_block_degree(h, in_block, false);
  long loA = std::max(0L, dA - 2), hiA = std::max(dA, hA) + opt.extra_block;
  long loB = std::max(0L, dB - 1), hiB = std::max(dB, h.total_degree() - emin) + opt.extra_rest;
  if (nb == 0) loA = hiA = 0;
  if (nr == 0) loB = hiB = 0;
  std::vector<std::tuple<long, int, int>> tries;  // (unknowns, DA, DB)
  for (long DA = loA; DA <= hiA; ++DA)
    for (long DB = loB; DB <= hiB; ++DB) {
      const long u = support_size(nb, nr, static_cast<int>(DA), static_cast<int>(DB)) *
                     static_cast<long>(gens.size());
      if (u <= opt.max_unknowns) tries.emplace_back(u, static_cast<int>(DA), static_cast<int>(DB));
    }
  std::sort(tries.begin(), tries.end());
  // Most failed supports are genuinely too small. Try larger supports before
  // paying for a second elimination on each one. If none works, revisit all
  // eligible supports using the remaining two primes.
  for (int pass = 0; pass < 2; ++pass) {
    for (auto& [u, DA, DB] : tries) {
      auto C = solve_bounded(gens, h, block, DA, DB, opt.max_primes,
                             pass == 0 ? 0 : 1, pass == 0 ? 1 : 3);
      if (C) return C;
    }
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
