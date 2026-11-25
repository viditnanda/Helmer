#include "helmer/characteristic.hpp"

#include "local_common.hpp"

namespace helmer {

using detail::agree;
using detail::random_vector;
using detail::reduce_all;
using modp::CtxPtr;
using modp::ModGB;
using modp::MPoly;

namespace {

// gamma_0, ..., gamma_d of V(F) for F generating a radical ideal of F_p[z] (empty if
// V(F) is empty). gamma_d: points of V(F) on d random affine hyperplanes. gamma_i, i < d:
// points of the i-dimensional polar locus, rank [A_{i+1} ; JR] < c + i + 1 with A random
// constant rows and JR = c random combinations of Jacobian rows, on i random affine
// hyperplanes, away from V(h), h = det(JR * R2) (vanishing on Sing and where the random
// combinations degenerate).
std::vector<long> global_polar_run(const CtxPtr& C, const std::vector<MPoly>& F) {
  const ulong p = C->p;
  const int n = C->n;
  ModGB G = modp::groebner(C, F);
  if (G.unit) return {};
  const int d = modp::dim(G, n), c = n - d;
  const auto& Zg = G.polys;
  std::vector<long> gamma(d + 1, 0);
  auto add_hyperplanes = [&](std::vector<MPoly>& K, int k) {
    for (int a = 0; a < k; ++a) {
      MPoly l = MPoly::constant(C, modp::random_elem(p));
      for (int j = 0; j < n; ++j) l = l + MPoly::var(C, j).scale(modp::random_elem(p));
      K.push_back(std::move(l));
    }
  };
  {
    auto K = Zg;
    add_hyperplanes(K, d);
    const long L = modp::length(C, K);
    if (L < 0) throw modp::Retry("linear section is not zero-dimensional");
    gamma[d] = L;
  }
  if (d == 0 || c == 0) return gamma;  // A^n: no critical points
  std::vector<std::vector<MPoly>> JR(c, std::vector<MPoly>(n, MPoly(C)));
  for (int a = 0; a < c; ++a)
    for (auto& g : Zg) {
      const ulong r = modp::random_elem(p);
      for (int j = 0; j < n; ++j) JR[a][j] = JR[a][j] + g.derivative(j).scale(r);
    }
  std::vector<std::vector<MPoly>> H(c, std::vector<MPoly>(c, MPoly(C)));
  for (int b = 0; b < c; ++b)
    for (int j = 0; j < n; ++j) {
      const ulong r = modp::random_elem(p);
      for (int a = 0; a < c; ++a) H[a][b] = H[a][b] + JR[a][j].scale(r);
    }
  auto hdet = modp::maximal_minors(C, H);
  if (hdet.empty()) throw modp::Retry("degenerate singular-locus element");
  const MPoly& h = hdet[0];
  for (int i = 0; i < d; ++i) {
    std::vector<std::vector<MPoly>> M;
    for (int r = 0; r <= i; ++r) {
      M.emplace_back();
      for (int j = 0; j < n; ++j) M.back().push_back(MPoly::constant(C, modp::random_elem(p)));
    }
    for (auto& row : JR) M.push_back(row);
    auto K = Zg;
    for (auto& f : modp::maximal_minors(C, M)) K.push_back(f);
    add_hyperplanes(K, i);
    const long L = modp::length_off(C, K, h);
    if (L < 0) throw modp::Retry("polar section is not zero-dimensional");
    gamma[i] = L;
  }
  return gamma;
}

long alternating(const std::vector<long>& gamma) {
  const int d = static_cast<int>(gamma.size()) - 1;
  long e = 0;
  for (int i = 0; i <= d; ++i) e += ((d - i) % 2 ? -1 : 1) * gamma[i];
  return e;
}

// Global Euler obstruction of V(I) ∩ L for a general affine-linear map
// u in F_p^m -> F_p^n (m = number of parameters); for homogeneous I and the chart of a
// general projective subspace this is the affine piece L_r \ L_{r+1}.
long sectional_global_eu(const Ideal& I, int m, ulong p) {
  const int n = I.ring()->nvars();
  auto C = modp::make_ctx(p, n);
  ModGB G = modp::groebner(C, reduce_all(I, C));
  if (G.unit) return 0;
  if (m == 0) {  // a general point
    auto x = random_vector(p, n);
    for (auto& g : G.polys)
      if (g.eval(x) != 0) return 0;
    return 1;
  }
  auto Cu = modp::make_ctx(p, m);
  auto img = modp::random_affine_map(Cu, random_vector(p, n), n);
  std::vector<MPoly> F;
  for (auto& g : G.polys) F.push_back(g.compose(Cu, img));
  auto gamma = global_polar_run(Cu, F);
  return gamma.empty() ? 0 : alternating(gamma);
}

// C(-r, m), the coefficient of h^m in (1 + h)^(-r).
long binom_neg(int r, int m) {
  if (m == 0) return 1;
  if (r == 0) return 0;
  long b = 1;  // C(r + m - 1, m)
  for (int k = 1; k <= m; ++k) b = b * (r + m - k) / k;
  return (m % 2 ? -1 : 1) * b;
}

// From s_r = integral of phi over a general codimension-r linear section (r = 0..N) to
// the degrees a_k of the dimension-k parts of c_*(phi) in A_*(P^N): by transversality,
// s_r = sum_{k >= r} a_k C(-r, k - r).
std::vector<long> classes_from_sections(const std::vector<long>& s) {
  const int N = static_cast<int>(s.size()) - 1;
  std::vector<long> a(N + 1, 0);
  for (int r = N; r >= 0; --r) {
    long v = s[r];
    for (int k = r + 1; k <= N; ++k) v -= a[k] * binom_neg(r, k - r);
    a[r] = v;
  }
  return a;
}

}  // namespace

long euler_obstruction(const Ideal& X, const Ideal& Z, const LocalOptions& o) {
  auto seq = polar_sequence(X, Z, o);
  const int d = static_cast<int>(seq.size()) - 1;
  if (d < 0 || seq[d] == 0) return 0;  // X not in V(Z)
  if (d == 0) return seq[0];
  long e = 0;  // Lê–Teissier
  for (int j = 1; j <= d; ++j) e += ((d - j) % 2 ? -1 : 1) * seq[j];
  return e;
}

std::vector<long> global_polar_degrees(const Ideal& Z, const LocalOptions& o) {
  return agree<std::vector<long>>(o, "global_polar_degrees", [&](ulong p) {
    const int n = Z.ring()->nvars();
    auto C = modp::make_ctx(p, n);
    return global_polar_run(C, reduce_all(Z, C));
  });
}

long global_euler_obstruction(const Ideal& Z, const LocalOptions& o) {
  auto g = global_polar_degrees(Z, o);
  return g.empty() ? 0 : alternating(g);
}

EulerData euler_obstruction_matrix(const Stratification& S, const LocalOptions& o,
                                   const std::function<void(const std::string&)>& log) {
  EulerData E{S.ring, S.projective, {}, {}, {}};
  for (auto& [d, L] : S.strata)  // increasing dimension
    for (auto& I : L) {
      E.dims.push_back(d);
      E.strata.push_back(I);
    }
  const size_t m = E.strata.size();
  E.eu.assign(m, std::vector<long>(m, 0));
  for (size_t t = 0; t < m; ++t) {
    E.eu[t][t] = 1;
    for (size_t s = 0; s < m; ++s) {
      if (E.dims[t] <= E.dims[s] || !E.strata[s].contains(E.strata[t])) continue;
      E.eu[s][t] = euler_obstruction(E.strata[s], E.strata[t], o);
      if (log) log("Eu of " + E.strata[t].str() + " along " + E.strata[s].str() + ": " +
                   std::to_string(E.eu[s][t]));
    }
  }
  return E;
}

std::vector<long> cc_multiplicities(const EulerData& E) {
  const size_t m = E.strata.size();
  int top = 0;
  for (int d : E.dims) top = std::max(top, d);
  auto sign = [&](size_t t) { return (top - E.dims[t]) % 2 ? -1L : 1L; };
  std::vector<long> mult(m, 0);
  for (size_t k = m; k-- > 0;) {  // decreasing dimension
    long rhs = 1;
    for (size_t t = 0; t < m; ++t)
      if (t != k && E.dims[t] > E.dims[k]) rhs -= E.eu[k][t] * sign(t) * mult[t];
    mult[k] = rhs * sign(k);
  }
  return mult;
}

ChernData chern_classes(const Stratification& S, const LocalOptions& o,
                        const std::function<void(const std::string&)>& log) {
  ChernData D;
  D.E = euler_obstruction_matrix(S, o, log);
  const auto& E = D.E;
  const size_t m = E.strata.size();
  const int n = S.ring->nvars();
  const int N = S.projective ? n - 1 : n;
  D.ambient = N;
  // chi_r(s): chi of open stratum s on the affine piece r (projective: L_r \ L_{r+1};
  // affine: a general codimension-r affine section), from the global Euler obstructions
  // of the closures by the unitriangular Euler-obstruction matrix.
  std::vector<std::vector<long>> piece(m, std::vector<long>(N + 1, 0));
  for (int r = 0; r <= N; ++r) {
    std::vector<long> g(m, 0);
    for (size_t t = 0; t < m; ++t) {
      if (E.dims[t] < r) continue;
      g[t] = agree<long>(o, "global Euler obstruction", [&](ulong p) {
        return sectional_global_eu(E.strata[t], N - r, p);
      });
    }
    for (size_t t = 0; t < m; ++t) {  // increasing dimension
      if (E.dims[t] < r) continue;
      long v = g[t];
      for (size_t s = 0; s < m; ++s)
        if (s != t && E.dims[s] < E.dims[t] && E.dims[s] >= r) v -= E.eu[s][t] * piece[s][r];
      piece[t][r] = v;
    }
    if (log) log("sections of codimension " + std::to_string(r) + " done");
  }
  D.sections.assign(m, std::vector<long>(N + 1, 0));
  for (size_t s = 0; s < m; ++s)
    for (int r = 0; r <= N; ++r) {
      if (S.projective)
        for (int j = r; j <= N; ++j) D.sections[s][r] += piece[s][j];
      else
        D.sections[s][r] = piece[s][r];
    }
  for (size_t s = 0; s < m; ++s) {
    D.chi.push_back(D.sections[s][0]);
    D.csm.push_back(classes_from_sections(D.sections[s]));
  }
  for (size_t t = 0; t < m; ++t) {
    std::vector<long> sec(N + 1, 0);
    for (size_t s = 0; s < m; ++s)
      if (E.eu[s][t])
        for (int r = 0; r <= N; ++r) sec[r] += E.eu[s][t] * D.sections[s][r];
    D.mather.push_back(classes_from_sections(sec));
  }
  return D;
}

}  // namespace helmer
