#include "helmer/whitney.hpp"

#include <algorithm>
#include <optional>

namespace helmer {

RingPtr conormal_ring(const RingPtr& R) {
  std::vector<std::string> names = R->names();
  const int n = R->nvars();
  for (int i = 0; i < n; ++i) {
    std::string v = "v_" + std::to_string(i);
    while (R->index_of(v) >= 0) v += "_";
    names.push_back(v);
  }
  return Ring::make(std::move(names));
}

static std::vector<int> identity(int n) {
  std::vector<int> m(n);
  for (int i = 0; i < n; ++i) m[i] = i;
  return m;
}

Ideal conormal(const Ideal& I, const RingPtr& S, const std::string& saturation) {
  const RingPtr& R = I.ring();
  const int n = R->nvars();
  const int c = I.codim();
  auto J1 = jacobian(I.gens());  // m x n over R
  auto idm = identity(n);
  // augmented matrix [v ; Jac] over S
  std::vector<std::vector<Poly>> M;
  std::vector<Poly> vrow;
  for (int i = 0; i < n; ++i) vrow.push_back(Poly::var(S, n + i));
  M.push_back(vrow);
  for (auto& row : J1) {
    std::vector<Poly> r;
    for (auto& p : row) r.push_back(p.embed(S, idm));
    M.push_back(std::move(r));
  }
  std::vector<Poly> gens;
  for (auto& g : I.gens()) gens.push_back(g.embed(S, idm));
  auto K = minors(c + 1, M);
  gens.insert(gens.end(), K.begin(), K.end());
  if (c == 0) return Ideal(S, gens);  // ambient space: the zero section, no saturation
  auto cminors = minors(c, J1);
  if (saturation == "minor") {
    // I prime: the conormal over the dense open set {g != 0} of X is smooth and
    // reduced, so K : g^oo is exactly I(Con X) for any c-minor g not in I.
    // Use the sparsest such minor.
    std::sort(cminors.begin(), cminors.end(),
              [](const Poly& a, const Poly& b) { return a.nterms() < b.nterms(); });
    for (auto& g : cminors)
      if (!I.contains(g)) return saturate(Ideal(S, gens), g.embed(S, idm));
    throw std::invalid_argument("helmer: conormal: every c-minor lies in I (not reduced?)");
  }
  std::vector<Poly> sat;
  for (auto& p : cminors) sat.push_back(p.embed(S, idm));
  return saturate(Ideal(S, gens), Ideal(S, sat), saturation);
}

std::vector<Ideal> polar_varieties(const Ideal& I, const Ideal& Jsing, long bound) {
  const RingPtr& R = I.ring();
  const int n = R->nvars();
  const int c = I.codim();
  const int d = n - c;
  auto Jac = jacobian(I.gens());
  std::vector<Ideal> out;
  for (int l = 1; l <= d + 1; ++l) {
    auto M = Jac;
    for (int r = 0; r < l; ++r) {
      std::vector<Poly> row;
      for (int i = 0; i < n; ++i) row.push_back(Poly::constant(R, random_coeff(bound)));
      M.push_back(std::move(row));
    }
    Ideal K = I + Ideal(R, minors(c + l, M));
    out.push_back(Jsing.is_unit() ? K : saturate(K, Jsing, "random"));
  }
  return out;
}

Ideal whitney_points(const Ideal& Z, const Ideal& ConX, const RingPtr& R) {
  const RingPtr& S = ConX.ring();
  const int n = R->nvars();
  auto idm = identity(n);
  std::vector<int> u = Z.independent_set();  // indices of x-variables, same in S
  Ideal CZ = ConX + embed(Z, S, idm);
  BlockData B = block_lc_lcm(CZ, u);
  if (B.h.is_constant()) return Ideal(R, {Poly::constant(R, 1)});
  std::vector<int> vs;
  for (int i = 0; i < n; ++i) vs.push_back(n + i);
  Ideal E = eliminate(Ideal(S, B.gb.polys) + B.h, vs);
  return restrict_to(E, R, idm);
}

namespace {

struct Store {
  int n_eff;
  std::map<int, std::vector<Ideal>> V;
  std::map<int, std::vector<std::string>> keys;

  int dim_of(const Ideal& K) const { return n_eff - K.codim(); }
  bool has(int d, const Ideal& K) const {
    auto it = keys.find(d);
    if (it == keys.end()) return false;
    auto k = K.key();
    return std::find(it->second.begin(), it->second.end(), k) != it->second.end();
  }
  // add K at its own dimension if new; returns true if added
  bool add(const Ideal& K) {
    int d = dim_of(K);
    if (d < 0) return false;
    if (has(d, K)) return false;
    V[d].push_back(K);
    keys[d].push_back(K.key());
    return true;
  }
};

}  // namespace

Stratification whitney_stratify(const Ideal& I0, const WhitneyOptions& opts) {
  auto log = [&](const std::string& s) {
    if (opts.log) opts.log(s);
  };
  const RingPtr& R = I0.ring();
  Ideal I = I0;
  const int nv = R->nvars();
  const int n = opts.projective ? nv - 1 : nv;
  Store st{n, {}, {}};
  const int dimI = n - I.codim();
  for (int i = 0; i <= std::max(dimI, 0); ++i) {
    st.V[i];
    st.keys[i];
  }
  Stratification out{R, opts.projective, {}};
  if (dimI < 0) return out;

  log("decompose X");
  auto comps = minimal_primes(I);
  bool equidim = true;
  for (auto& P : comps) {
    st.add(P);  // each component at its own dimension
    equidim &= st.dim_of(P) == dimI;
  }
  const int level = opts.level < 0 ? dimI : opts.level;
  auto finish = [&]() {
    for (auto& [d, L] : st.V) out.strata[d] = L;
    return out;
  };
  if (level == 0) return finish();
  // Intersections of distinct components lie in Sing X. (M2 gets this by filing every
  // component at the top dimension; we file components at their true dimension, so
  // intersect all pairs here, which also covers components of different dimensions.)
  if (!equidim)
    for (size_t i = 0; i < comps.size(); ++i)
      for (size_t j = i + 1; j < comps.size(); ++j)
        for (auto& th : minimal_primes(comps[i] + comps[j])) st.add(th);

  Ideal J = jacobian_ideal(I);
  if (equidim && dimI < 2) {
    for (auto& K : minimal_primes(J))
      if (st.dim_of(K) < dimI) st.add(K);
    return finish();
  }
  if (equidim && J.codim() == n)
    for (auto& K : minimal_primes(J)) st.add(K);
  if (equidim && J.codim() >= n) return finish();

  RingPtr S = conormal_ring(R);
  for (int a = 0; a < level; ++a) {
    const int k = dimI - a;
    std::vector<Ideal> cur = st.V[k];  // snapshot (additions below go to lower dims)
    log("level " + std::to_string(k) + ": " + std::to_string(cur.size()) + " pieces");
    // pairwise intersections
    for (size_t i = 0; i < cur.size(); ++i)
      for (size_t j = i + 1; j < cur.size(); ++j)
        for (auto& th : minimal_primes(cur[i] + cur[j])) st.add(th);
    for (auto& Ia : cur) {
      Ideal Ja = jacobian_ideal(Ia);
      const int mu = st.dim_of(Ja);  // dim of Sing(Ia); negative if smooth
      for (auto& K : minimal_primes(Ja))
        if (st.dim_of(K) < k) st.add(K);
      if (opts.method == "polar") {
        // Teissier: (B) fails along Y ⊆ Sing only where polar multiplicities jump;
        // such points lie in Sing ∩ (closures of polar varieties).
        if (mu > 0) {
          log("  polar varieties of a " + std::to_string(k) + "-dimensional piece");
          for (auto& P : polar_varieties(Ia, Ja, opts.polar_bound))
            for (auto& K : minimal_primes(Ja + P))
              if (st.dim_of(K) < k) st.add(K);
        }
        continue;
      }
      if (mu > dimI - level) {
        log("  conormal of a " + std::to_string(k) + "-dimensional piece");
        Ideal C = conormal(Ia, S, opts.saturation == "random" ? "minor" : opts.saturation);
        for (int ell = mu; ell >= 0; --ell) {
          // iterate by index: entries appended at dimension < ell do not affect this list
          for (size_t zi = 0; zi < st.V[ell].size(); ++zi) {
            Ideal Z = st.V[ell][zi];
            if (!Z.contains(Ia)) continue;  // need V(Z) ⊆ V(Ia)
            Ideal W = whitney_points(Z, C, R);
            if (W.is_unit()) continue;
            for (auto& K : minimal_primes(W)) st.add(K);
          }
        }
      }
    }
  }
  return finish();
}

}  // namespace helmer

// ======================================================================= maps
namespace helmer {

namespace {

struct Graph {
  RingPtr R, S, G;          // source, target, graph ring k[x, y]
  int n, m;
  std::vector<int> xmap, ymap;  // R -> G, S -> G
  Ideal Gam;                // (y_i - F_i)
  std::vector<int> xs, ys;  // variable indices in G
};

Graph make_graph(const std::vector<Poly>& F, const RingPtr& R, const RingPtr& S) {
  Graph g{R, S, nullptr, R->nvars(), S->nvars(), {}, {}, Ideal(R), {}, {}};
  if (static_cast<int>(F.size()) != g.m)
    throw std::invalid_argument("helmer: map has " + std::to_string(F.size()) +
                                " components but the target ring has " + std::to_string(g.m) + " variables");
  std::vector<std::string> names = R->names();
  for (int j = 0; j < g.m; ++j) names.push_back("helmer_y" + std::to_string(j));
  g.G = Ring::make(names);
  for (int i = 0; i < g.n; ++i) { g.xmap.push_back(i); g.xs.push_back(i); }
  for (int j = 0; j < g.m; ++j) { g.ymap.push_back(g.n + j); g.ys.push_back(g.n + j); }
  std::vector<Poly> gam;
  for (int j = 0; j < g.m; ++j) gam.push_back(Poly::var(g.G, g.n + j) - F[j].embed(g.G, g.xmap));
  g.Gam = Ideal(g.G, gam);
  return g;
}

// closure of f(K ∩ X) ⊆ k^m, plus IY
Ideal image_of(const Graph& g, const Ideal& K, const Ideal& IY) {
  Ideal E = eliminate(g.Gam + embed(K, g.G, g.xmap), g.xs);
  return restrict_to(E, g.S, g.ymap) + IY;
}
// f^{-1}(K) ∩ X
Ideal preimage_of(const Graph& g, const Ideal& K, const Ideal& IX) {
  Ideal E = eliminate(g.Gam + embed(K, g.G, g.ymap), g.ys);
  return restrict_to(E, g.R, g.xmap) + IX;
}

}  // namespace

Ideal non_proper_set(const std::vector<Poly>& F, const Ideal& IX, const Ideal& IY) {
  const RingPtr& R = IX.ring();
  const RingPtr& S = IY.ring();
  Graph g = make_graph(F, R, S);
  const int n = g.n, m = g.m;
  // Gröbner basis for [x] > [y]: refines x-degree, so x-homogenising it gives the
  // closure of the graph of f|X in P^n x k^m.
  std::vector<int> perm = g.xs;
  perm.insert(perm.end(), g.ys.begin(), g.ys.end());
  GB B = groebner(g.G, (g.Gam + embed(IX, g.G, g.xmap)).gens(), perm, n);
  // ring k[x, u, y]
  std::vector<std::string> names = R->names();
  names.push_back("helmer_u");
  for (int j = 0; j < m; ++j) names.push_back("helmer_y" + std::to_string(j));
  RingPtr H = Ring::make(names);
  const int u = n;
  std::vector<Poly> hom;
  for (auto& p : B.polys) {
    long d = 0;
    for (long t = 0; t < p.nterms(); ++t) {
      auto e = p.term_exp(t);
      long dx = 0;
      for (int i = 0; i < n; ++i) dx += e[i];
      d = std::max(d, dx);
    }
    std::vector<std::vector<long>> exps;
    std::vector<fmpq*> cs;
    std::vector<fmpq> store(p.nterms());
    for (long t = 0; t < p.nterms(); ++t) {
      auto e = p.term_exp(t);
      long dx = 0;
      for (int i = 0; i < n; ++i) dx += e[i];
      std::vector<long> f(n + 1 + m);
      for (int i = 0; i < n; ++i) f[i] = e[i];
      f[u] = d - dx;
      for (int j = 0; j < m; ++j) f[n + 1 + j] = e[n + j];
      fmpq_init(&store[t]);
      fmpq_mpoly_get_term_coeff_fmpq(&store[t], p.raw(), t, g.G->ctx());
      exps.push_back(f);
      cs.push_back(&store[t]);
    }
    hom.push_back(Poly::from_terms(H, exps, cs));
    for (auto* c : cs) fmpq_clear(c);
  }
  std::vector<int> elim;
  for (int i = 0; i <= n; ++i) elim.push_back(i);
  std::vector<int> ymapH;
  for (int j = 0; j < m; ++j) ymapH.push_back(n + 1 + j);
  std::optional<Ideal> acc;
  for (int i = 0; i < n; ++i) {  // chart x_i = 1 at infinity (u = 0)
    std::vector<Poly> vals;
    for (int k = 0; k < H->nvars(); ++k)
      vals.push_back(k == i ? Poly::constant(H, 1) : (k == u ? Poly(H) : Poly::var(H, k)));
    std::vector<Poly> gi;
    for (auto& p : hom) gi.push_back(p.substitute(H, vals));
    Ideal E = restrict_to(eliminate(Ideal(H, gi), elim), S, ymapH);
    acc = acc ? intersect(*acc, E) : E;
  }
  return *acc + IY;
}

MapStratification map_stratify(const std::vector<Poly>& F, const Ideal& IX, const Ideal& IY0,
                               const MapOptions& opts) {
  auto log = [&](const std::string& s) { if (opts.log) opts.log(s); };
  const RingPtr& R = IX.ring();
  const RingPtr& S = IY0.ring();
  Graph g = make_graph(F, R, S);
  const int n = g.n, m = g.m;
  RingPtr CR = conormal_ring(R), CS = conormal_ring(S);
  Ideal IY = IY0;

  Store V{n, {}, {}}, W{m, {}, {}};
  std::map<int, size_t> Vchk, Wchk;
  const int dimX = n - IX.codim();
  int dimY = m - IY.codim();
  Ideal fX = image_of(g, IX, IY);
  const int dimfX = m - fX.codim();
  log("decompose X and Y");
  for (auto& P : minimal_primes(IX)) V.add(P);
  for (auto& P : minimal_primes(IY)) W.add(P);
  if (!fX.equals(IY)) {
    IY = fX;
    dimY = dimfX;
    for (auto& P : minimal_primes(IY)) W.add(P);
  }
  int aX = 0, aY = 0;
  auto addV = [&](const Ideal& K) {
    if (V.add(K)) {
      int d = V.dim_of(K);
      if (d > dimX - aX) aX = dimX - d;  // unprocessed piece above: go back up
      return true;
    }
    return false;
  };
  auto addW = [&](const Ideal& K) {
    if (W.add(K)) {
      int d = W.dim_of(K);
      if (d > dimY - aY) aY = dimY - d;
      return true;
    }
    return false;
  };
  auto pull = [&](const Ideal& K) {
    for (auto& c : minimal_primes(preimage_of(g, K, IX))) addV(c);
  };
  auto push = [&](const Ideal& K) {
    for (auto& c : minimal_primes(image_of(g, K, IY))) addW(c);
  };

  if (!opts.proper) {
    log("Jelonek flag");
    std::vector<Ideal> JY{IY}, JX{IX};
    for (int i = 0; i < dimY; ++i) {
      JY.push_back(non_proper_set(F, JX.back(), IY));
      JX.push_back(preimage_of(g, JY.back(), IX));
    }
    for (auto& J : JY)
      for (auto& c : minimal_primes(J + IY)) W.add(c);
    for (auto& J : JX)
      for (auto& c : minimal_primes(J + IX)) V.add(c);
  }

  while (aX < dimX || aY < dimY) {
    if (dimY - aY > 0) {
      const int k = dimY - aY;
      auto& L = W.V[k];
      size_t t = Wchk[k];
      std::vector<Ideal> cur = L;
      for (size_t i = 0; i < cur.size(); ++i)
        for (size_t j = std::max(i + 1, t); j < cur.size(); ++j)
          for (auto& th : minimal_primes(cur[i] + cur[j])) {
            W.add(th);
            pull(th);
          }
      log("Y level " + std::to_string(k));
      for (size_t j2 = Wchk[k]; j2 < W.V[k].size(); ++j2) {
        Ideal Ya = W.V[k][j2];
        Wchk[k] = j2 + 1;
        Ideal JYa = jacobian_ideal(Ya);
        const int muY = W.dim_of(JYa);
        for (auto& K : minimal_primes(JYa)) {
          W.add(K);
          pull(K);
        }
        if (muY > 0 && opts.strata != "singular_only") {
          Ideal C = conormal(Ya, CS, "minor");
          for (int ell = muY; ell >= 0; --ell)
            for (size_t zi = 0; zi < W.V[ell].size(); ++zi) {
              Ideal Z = W.V[ell][zi];
              Ideal Wp = whitney_points(Z, C, S);
              if (Wp.is_unit()) continue;
              for (auto& K : minimal_primes(Wp)) {
                W.add(K);
                pull(K);
              }
            }
        }
      }
      aY++;
    }
    if (dimX - aX > 0) {
      const int k = dimX - aX;
      size_t t = Vchk[k];
      std::vector<Ideal> cur = V.V[k];
      for (size_t i = 0; i < cur.size(); ++i)
        for (size_t j = std::max(i + 1, t); j < cur.size(); ++j)
          for (auto& th : minimal_primes(cur[i] + cur[j])) {
            V.add(th);
            push(th);
          }
      log("X level " + std::to_string(k));
      for (size_t j2 = Vchk[k]; j2 < V.V[k].size(); ++j2) {
        Ideal Xa = V.V[k][j2];
        Vchk[k] = j2 + 1;
        Ideal JXa = jacobian_ideal(Xa);
        const int muX = V.dim_of(JXa);
        for (auto& K : minimal_primes(JXa)) V.add(K);
        // rank-drop loci of d(f|Xa): Xa + minors of [Jac Xa ; Jac F]
        std::vector<Poly> rowsrc = Xa.gens();
        rowsrc.insert(rowsrc.end(), F.begin(), F.end());
        auto AX = jacobian(rowsrc);
        const int r = std::min<int>(AX.size(), n);
        for (int i = 0; i <= r; ++i) {
          Ideal Ii = (r - i > 0) ? Xa + Ideal(R, minors(r - i, AX)) : Ideal(R, {Poly::constant(R, 1)});
          if (Ii.is_unit()) break;
          for (auto& K : minimal_primes(Ii)) {
            V.add(K);
            push(K);
          }
        }
        if (muX > 0 && opts.strata != "singular_only") {
          Ideal C = conormal(Xa, CR, "minor");
          for (int ell = muX; ell >= 0; --ell)
            for (size_t zi = 0; zi < V.V[ell].size(); ++zi) {
              Ideal Z = V.V[ell][zi];
              Ideal Wp = whitney_points(Z, C, R);
              if (Wp.is_unit()) continue;
              for (auto& K : minimal_primes(Wp))
                if (V.add(K)) push(K);
            }
        }
      }
      aX++;
    }
  }
  MapStratification out{{R, false, {}}, {S, false, {}}};
  for (int d = 0; d <= dimX; ++d) out.source.strata[d] = V.V[d];
  for (auto& [d, L] : W.V) out.target.strata[d] = L;
  for (int d = 0; d <= std::max(dimY, 0); ++d) out.target.strata[d];
  return out;
}

}  // namespace helmer
