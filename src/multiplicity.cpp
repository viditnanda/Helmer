#include "helmer/multiplicity.hpp"

#include "helmer/modp.hpp"

#include <flint/nmod_poly.h>
#include <flint/nmod_poly_factor.h>

#include <algorithm>
#include <functional>
#include <map>
#include <tuple>

namespace helmer {

using modp::CtxPtr;
using modp::ModGB;
using modp::MPoly;

namespace {

// Repeat run(p) with fresh primes and random choices until `trials` results agree.
template <class T, class F>
T agree(const LocalOptions& o, const char* what, F&& run) {
  std::vector<std::pair<T, int>> seen;
  int runs = 0, failures = 0, pointless = 0;
  std::string last;
  while (runs < o.max_runs) {
    ++runs;
    try {
      T r = run(modp::random_prime());
      auto it = std::find_if(seen.begin(), seen.end(), [&](auto& s) { return s.first == r; });
      if (it == seen.end()) seen.emplace_back(r, 1), it = seen.end() - 1;
      else it->second++;
      if (it->second >= o.trials) return it->first;
    } catch (const modp::NoPoint& e) {
      --runs;  // a different prime is needed, not a different random choice
      if (++pointless > o.max_pointless_primes) {
        ++runs;
        ++failures;
        last = e.what();
        break;
      }
    } catch (const modp::Unlucky& e) {
      ++failures;
      last = e.what();
    } catch (const modp::Retry& e) {
      ++failures;
      last = e.what();
    }
  }
  throw std::runtime_error(std::string("helmer: ") + what + ": no " + std::to_string(o.trials) +
                           " of " + std::to_string(runs) + " randomised runs agreed (" +
                           std::to_string(failures) + " discarded" +
                           (last.empty() ? std::string() : "; last: " + last) + ")");
}

std::vector<MPoly> reduce_all(const Ideal& I, const CtxPtr& C) {
  std::vector<MPoly> out;
  for (auto& g : I.gens()) out.push_back(modp::reduce(g, C));
  return out;
}

std::vector<ulong> random_vector(ulong p, int n) {
  std::vector<ulong> v(n);
  for (auto& e : v) e = modp::random_elem(p);
  return v;
}

// The univariate polynomial in variable `v` of a block basis eliminating the others,
// with v moved last; or nothing.
bool eliminant(const CtxPtr& C, const std::vector<MPoly>& F, int v, nmod_poly_t P) {
  const int m = C->n;
  std::vector<MPoly> perm;  // swap v and m-1
  for (int i = 0; i < m; ++i) perm.push_back(MPoly::var(C, i == v ? m - 1 : i == m - 1 ? v : i));
  std::vector<MPoly> G;
  for (auto& f : F) G.push_back(f.compose(C, perm));
  ModGB B = modp::groebner(C, G, m - 1);
  if (B.unit) return false;
  for (auto& g : B.polys) {
    bool uni = true;
    for (int i = 0; i < m - 1 && uni; ++i) uni = g.degree_in(i) <= 0;
    if (!uni || g.total_degree() <= 0) continue;
    nmod_poly_zero(P);
    std::vector<ulong> e(m);
    for (long t = 0; t < g.nterms(); ++t) {
      nmod_mpoly_get_term_exp_ui(e.data(), g.raw(), t, C->ctx);
      nmod_poly_set_coeff_ui(P, e[m - 1],
                             nmod_mpoly_get_term_coeff_ui(g.raw(), t, C->ctx));
    }
    return true;
  }
  return false;
}

std::vector<ulong> roots(const nmod_poly_t P) {
  nmod_poly_factor_t fac;
  nmod_poly_factor_init(fac);
  nmod_poly_roots(fac, P, 0);
  std::vector<ulong> r;
  for (slong i = 0; i < fac->num; ++i) {  // factors x - r, monic
    ulong c0 = nmod_poly_get_coeff_ui(fac->p + i, 0);
    r.push_back(c0 == 0 ? 0 : P->mod.n - c0);
  }
  nmod_poly_factor_clear(fac);
  return r;
}

// A general F_p-point of V(F), F a Gröbner basis of an ideal of dimension k: slice by k
// random affine hyperplanes (parametrised: x = b + B u), take an F_p-root of the
// eliminant in the last slice coordinate, and solve for the others.
std::optional<std::vector<ulong>> general_point(const CtxPtr& C, const ModGB& GX, int k) {
  const int n = C->n;
  const ulong p = C->p;
  if (k < 0) throw std::invalid_argument("helmer: the subvariety is empty");
  if (k == n) return random_vector(p, n);
  const int m = n - k;
  auto Cu = modp::make_ctx(p, m);
  nmod_poly_t P;
  nmod_poly_init(P, p);
  for (int attempt = 0; attempt < (k == 0 ? 1 : 12); ++attempt) {
    auto img = modp::random_affine_map(Cu, random_vector(p, n), n);
    std::vector<MPoly> S;
    for (auto& g : GX.polys) S.push_back(g.compose(Cu, img));
    if (!eliminant(Cu, S, m - 1, P)) continue;
    auto rs = roots(P);
    if (rs.empty()) continue;
    const ulong r = rs[modp::random_elem(rs.size())];
    auto S2 = S;
    S2.push_back(MPoly::var(Cu, m - 1) - MPoly::constant(Cu, r));
    // Coordinates: from linear basis elements when the point is reduced, else from
    // the eliminant in each coordinate (a power of u_v - c).
    std::vector<long> u(m, -1);
    u[m - 1] = static_cast<long>(r);
    ModGB B = modp::groebner(Cu, S2);
    if (B.unit) continue;
    for (size_t g = 0; g < B.polys.size(); ++g) {
      const MPoly& f = B.polys[g];
      if (f.total_degree() != 1) continue;
      int v = -1;
      for (int i = 0; i < m; ++i)
        if (B.lead[g][i]) v = i;
      if (v < 0) continue;
      bool only = true;  // f = a u_v + c, after reduction
      for (int i = 0; i < m && only; ++i)
        if (i != v && f.degree_in(i) > 0) only = false;
      if (!only) continue;
      std::vector<ulong> zero(m, 0), unit(m, 0);
      unit[v] = 1;
      ulong c = f.eval(zero), a = (f.eval(unit) + p - c) % p;
      u[v] = static_cast<long>(n_mulmod2_preinv(p - c, n_invmod(a, p), p, n_preinvert_limb(p)));
    }
    bool ok = true;
    for (int v = 0; v < m - 1 && ok; ++v) {
      if (u[v] >= 0) continue;
      if (!eliminant(Cu, S2, v, P)) { ok = false; break; }
      auto rv = roots(P);
      if (rv.size() != 1) { ok = false; break; }
      u[v] = static_cast<long>(rv[0]);
    }
    if (!ok) continue;
    std::vector<ulong> uu(u.begin(), u.end()), x(n);
    for (int i = 0; i < n; ++i) x[i] = img[i].eval(uu);
    nmod_poly_clear(P);
    return x;
  }
  nmod_poly_clear(P);
  return std::nullopt;
}

// Multiplicity at the origin of V(G) : h^∞ in slice coordinates z (Cz), where the slice
// has dimension n - w + 1 (a general linear space through the point; all of k^n if w = 0)
// and w is the local dimension expected at the origin. 0 if the origin is not on a
// w-dimensional part.
long local_mult(const CtxPtr& Cz, const std::vector<MPoly>& G, const MPoly* h, int w) {
  const int s = Cz->n;
  auto at_origin = [](const std::vector<MPoly>& F) {
    for (auto& f : F)
      if (f.constant_term() != 0) return false;
    return true;
  };
  if (w == 0) {
    auto Cl = h ? modp::saturate(Cz, G, *h) : G;
    if (!at_origin(Cl)) return 0;
    // An isolated point: the length of O/(m^N) stabilises exactly at the local length.
    long prev = -1;
    for (int N = 1; N <= 64; ++N) {
      auto F = Cl;
      std::vector<int> e(s, 0);
      std::function<void(int, int)> mon = [&](int v, int left) {
        if (v == s - 1) {
          e[v] = left;
          MPoly t = MPoly::constant(Cz, 1);
          for (int i = 0; i < s; ++i)
            for (int j = 0; j < e[i]; ++j) t = t * MPoly::var(Cz, i);
          F.push_back(std::move(t));
          return;
        }
        for (int a = left; a >= 0; --a) {
          e[v] = a;
          mon(v + 1, left - a);
        }
      };
      mon(0, N);
      long L = modp::length(Cz, F);
      if (L == prev) return L;
      prev = L;
    }
    throw modp::Retry("origin is not an isolated point");
  }
  MPoly sat = MPoly::var(Cz, 0);
  if (h) sat = sat * *h;
  auto Cl = modp::saturate(Cz, G, sat);
  if (Cl.empty()) return s == 1 ? 1 : 0;
  if (!at_origin(Cl)) return 0;
  if (s == 1) return 1;
  auto Z = Cl;
  Z.push_back(MPoly::var(Cz, 1));
  const long tot = modp::length(Cz, Z);
  if (tot < 0) throw modp::Retry("hyperplane section is not zero-dimensional");
  const long off = modp::length_off(Cz, Z, MPoly::var(Cz, 0));
  if (off < 0) throw modp::Retry("hyperplane section is not zero-dimensional");
  return tot - off;
}

long multiplicity_run(const Ideal& X, const Ideal& Y, ulong p) {
  const int n = X.ring()->nvars();
  auto C = modp::make_ctx(p, n);
  ModGB GX = modp::groebner(C, reduce_all(X, C));
  const int k = modp::dim(GX, n);
  ModGB GY = modp::groebner(C, reduce_all(Y, C));
  if (GY.unit) return 0;
  const int m = modp::dim(GY, n);
  if (k > m) return 0;
  auto x = general_point(C, GX, k);
  if (!x) throw modp::NoPoint();
  for (auto& g : GY.polys)
    if (g.eval(*x) != 0) return 0;
  for (int w = m; w >= k; --w) {
    auto Cz = modp::make_ctx(p, w >= 1 ? n - w + 1 : n);
    auto img = modp::random_affine_map(Cz, *x, n);
    std::vector<MPoly> G;
    for (auto& g : GY.polys) G.push_back(g.compose(Cz, img));
    if (long e = local_mult(Cz, G, nullptr, w)) return e;
  }
  throw modp::Retry("no local dimension found");
}

// Entries ws (dimensions) of the polar sequence (all if empty); -1 for the others.
// The polar locus is saturated by h = det(JR * R2), a combination of c-minors vanishing
// on Sing Z; with it the minor equations behave like Lagrange equations grad = mu * a,
// which is much cheaper for F4 than saturating by low-degree equations of Sing Z.
std::vector<long> polar_run(const Ideal& X, const Ideal& Z, ulong p, const std::vector<int>& ws) {
  const int n = X.ring()->nvars();
  auto C = modp::make_ctx(p, n);
  ModGB GZ = modp::groebner(C, reduce_all(Z, C));
  if (GZ.unit) throw std::invalid_argument("helmer: polar_sequence: Z is empty");
  const int d = modp::dim(GZ, n), c = n - d;
  std::vector<long> seq(d + 1, -1);
  std::vector<char> want(d + 1, ws.empty());
  for (int w : ws)
    if (0 <= w && w <= d) want[w] = 1;
  for (int w = 0; w <= d; ++w)
    if (want[w]) seq[w] = 0;
  ModGB GX = modp::groebner(C, reduce_all(X, C));
  const int k = modp::dim(GX, n);
  if (k > d) return seq;
  auto x = general_point(C, GX, k);
  if (!x) throw modp::NoPoint();
  for (auto& g : GZ.polys)
    if (g.eval(*x) != 0) return seq;
  const auto& Zg = GZ.polys;
  long seq_d = -1;
  // c rows of random combinations of the Jacobian rows, and a random n x c matrix for
  // one c-minor-type element h = det(JR * R2) vanishing on Sing Z.
  std::vector<std::vector<MPoly>> JR(c, std::vector<MPoly>(n, MPoly(C)));
  for (int a = 0; a < c; ++a)
    for (auto& g : Zg) {
      const ulong r = modp::random_elem(p);
      for (int j = 0; j < n; ++j) JR[a][j] = JR[a][j] + g.derivative(j).scale(r);
    }
  std::vector<std::vector<ulong>> R2(n, std::vector<ulong>(c));
  for (auto& row : R2)
    for (auto& e : row) e = modp::random_elem(p);
  // Entries w < d with w <= k = dim X vanish: for w < k the polar variety is too small
  // to contain X; for w = k it is the closure of a w-dimensional subset of Z_reg, so it
  // contains X only as a component, impossible if X ⊆ Sing Z, while if X meets Z_reg a
  // generic projection is submersive at a general point of X. Likewise every entry
  // w < d vanishes when e_X(Z) = 1 (x a smooth point). So compute entry d first.
  bool need_lower = false;
  for (int w = k + 1; w < d; ++w) need_lower |= static_cast<bool>(want[w]);
  for (int l = d + 1; l >= 1; --l) {
    const int w = l - 1;
    if (w < d && (w <= k || !want[w])) continue;
    if (w == d && !want[d] && !need_lower) continue;
    if (w < d && seq_d == 1) break;
    auto Cz = modp::make_ctx(p, w >= 1 ? n - w + 1 : n);
    auto img = modp::random_affine_map(Cz, *x, n);
    std::vector<MPoly> G;
    for (auto& g : Zg) G.push_back(g.compose(Cz, img));
    std::vector<std::vector<MPoly>> JRz(c);
    for (int a = 0; a < c; ++a)
      for (int j = 0; j < n; ++j) JRz[a].push_back(JR[a][j].compose(Cz, img));
    std::optional<MPoly> h;
    if (c >= 1) {
      std::vector<std::vector<MPoly>> H(c, std::vector<MPoly>(c, MPoly(Cz)));
      for (int a = 0; a < c; ++a)
        for (int b = 0; b < c; ++b)
          for (int j = 0; j < n; ++j) H[a][b] = H[a][b] + JRz[a][j].scale(R2[j][b]);
      auto det = modp::maximal_minors(Cz, H);
      if (det.empty()) throw modp::Retry("degenerate singular-locus element");
      h = det[0];
    }
    if (c + l <= n) {  // rank [A_l ; JR] < c + l, constant rows first
      std::vector<std::vector<MPoly>> M;
      for (int r = 0; r < l; ++r) {
        M.emplace_back();
        for (int j = 0; j < n; ++j) M.back().push_back(MPoly::constant(Cz, modp::random_elem(p)));
      }
      for (auto& row : JRz) M.push_back(row);
      for (auto& f : modp::maximal_minors(Cz, M)) G.push_back(f);
    }
    const long e = local_mult(Cz, G, h ? &*h : nullptr, w);
    if (w == d) seq_d = e;
    if (want[w]) seq[w] = e;
  }
  return seq;
}

}  // namespace

long multiplicity(const Ideal& X, const Ideal& Y, const LocalOptions& o) {
  if (X.ring() != Y.ring()) throw std::invalid_argument("helmer: multiplicity: different rings");
  return agree<long>(o, "multiplicity", [&](ulong p) { return multiplicity_run(X, Y, p); });
}

bool contained(const Ideal& X, const Ideal& Y, const LocalOptions& o) {
  if (X.ring() != Y.ring()) throw std::invalid_argument("helmer: contained: different rings");
  return agree<bool>(o, "contained", [&](ulong p) {
    const int n = X.ring()->nvars();
    auto C = modp::make_ctx(p, n);
    auto Xp = reduce_all(X, C);
    MPoly h(C);
    for (auto& g : reduce_all(Y, C)) h = h + g.scale(modp::random_elem(p));
    return modp::in_radical(C, Xp, h);
  });
}

bool in_singular_locus(const Ideal& X, const Ideal& Y, const LocalOptions& o) {
  return multiplicity(X, Y, o) > 1;
}

std::vector<long> polar_sequence(const Ideal& X, const Ideal& Z, const LocalOptions& o) {
  if (X.ring() != Z.ring()) throw std::invalid_argument("helmer: polar_sequence: different rings");
  return agree<std::vector<long>>(o, "polar_sequence", [&](ulong p) { return polar_run(X, Z, p, {}); });
}

long polar_multiplicity(const Ideal& X, const Ideal& Z, int w, const LocalOptions& o) {
  if (X.ring() != Z.ring()) throw std::invalid_argument("helmer: polar_multiplicity: different rings");
  return agree<long>(o, "polar_multiplicity", [&](ulong p) {
    auto s = polar_run(X, Z, p, {w});
    return w >= 0 && w < static_cast<int>(s.size()) ? s[w] : 0L;
  });
}

Stratification min_coarsen(const Stratification& S, const LocalOptions& o,
                           const std::function<void(const std::string&)>& log) {
  struct Stratum {
    int dim;
    Ideal I;
  };
  std::vector<Stratum> all;
  for (auto& [d, L] : S.strata)
    for (auto& I : L) all.push_back({d, I});
  std::stable_sort(all.begin(), all.end(), [](auto& a, auto& b) { return a.dim > b.dim; });
  const size_t N = all.size();
  // above[i][j]: closure of i contains closure of j (ideal of i inside ideal of j: primes)
  std::map<std::pair<size_t, size_t>, bool> above_cache;
  auto above = [&](size_t i, size_t j) {
    auto key = std::make_pair(i, j);
    auto it = above_cache.find(key);
    if (it != above_cache.end()) return it->second;
    bool r = all[i].dim > all[j].dim && all[j].I.contains(all[i].I);
    above_cache[key] = r;
    return r;
  };
  std::map<std::tuple<size_t, size_t, int>, long> polar_cache;
  auto polar = [&](size_t x, size_t z, int w) {
    auto key = std::make_tuple(x, z, w);
    auto it = polar_cache.find(key);
    if (it != polar_cache.end()) return it->second;
    long r = agree<long>(o, "polar_multiplicity", [&](ulong p) {
      return polar_run(all[x].I, all[z].I, p, {w})[w];
    });
    polar_cache[key] = r;
    return r;
  };
  // Polar sequences of z along x and along y, entry by entry, cheapest-first order
  // (top dimension, i.e. the multiplicity, then upwards from points); stops at the
  // first difference.
  auto same_polar = [&](size_t x, size_t y, size_t z) {
    const int d = all[z].dim + (S.projective ? 1 : 0);
    std::vector<int> order{d};  // entries w <= dim X vanish for both (see polar_run)
    for (int w = std::max(all[x].dim + (S.projective ? 1 : 0) + 1, 1); w < d; ++w) order.push_back(w);
    for (int w : order)
      if (polar(x, z, w) != polar(y, z, w)) return false;
    return true;
  };
  std::vector<size_t> kept;
  for (size_t j = 0; j < N; ++j) {
    std::vector<size_t> cont;
    for (size_t i : kept)
      if (above(i, j)) cont.push_back(i);
    auto keep = [&](const std::string& why) {
      if (log) log("keep " + all[j].I.str() + ": " + why);
      kept.push_back(j);
    };
    if (cont.empty()) {
      kept.push_back(j);
      continue;
    }
    // the unique smallest stratum above j, below every other one above j
    size_t y = cont[0];
    for (size_t i : cont)
      if (all[i].dim < all[y].dim) y = i;
    bool unique = true;
    for (size_t i : cont)
      if (i != y && !above(i, y)) { unique = false; break; }
    if (!unique) { keep("in the closures of several strata"); continue; }
    if (multiplicity(all[j].I, all[y].I, o) != 1) { keep("singular point of " + all[y].I.str()); continue; }
    bool same = true;
    for (size_t z : cont) {
      if (z == y) continue;
      if (!same_polar(j, y, z)) {
        keep("polar multiplicities of " + all[z].I.str() + " jump");
        same = false;
        break;
      }
    }
    if (!same) continue;
    if (log) log("merge " + all[j].I.str() + " into " + all[y].I.str());
  }
  Stratification out{S.ring, S.projective, {}};
  for (auto& [d, L] : S.strata) out.strata[d];  // keep every dimension key
  for (size_t i : kept) out.strata[all[i].dim].push_back(all[i].I);
  return out;
}

}  // namespace helmer
