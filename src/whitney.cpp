#include "helmer/whitney.hpp"

#include <algorithm>
#include <optional>
#include <numeric>

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

static Stratification stratify_impl(const Ideal&, const WhitneyOptions&, bool defer_primes);

namespace {

// With a degree-compatible order, the linear elements of a reduced basis
// span all affine linear forms in the ideal.
std::optional<Poly> chart_form(const Ideal& X, const std::vector<Ideal>& singular) {
  if (singular.empty()) return {};
  auto R = X.ring();
  for (int i = R->nvars() - 1; i >= 0; --i) {
    auto x = Poly::var(R, i);
    if (std::all_of(singular.begin(), singular.end(),
                    [&](const Ideal& P) { return P.contains(x); }) && !X.contains(x))
      return x;
  }

  std::vector<Poly> forms;
  for (const auto& g : singular.front().gb().polys)
    if (g.total_degree() == 1) forms.push_back(g);
  for (size_t j = 1; j < singular.size() && !forms.empty(); ++j) {
    // Gaussian elimination on the images of the forms modulo the next prime.
    // Apply the same row operations to the forms to recover the kernel.
    std::vector<Poly> kernel;
    std::map<std::vector<long>, std::pair<Poly, Poly>> pivots;
    for (auto form : forms) {
      auto image = singular[j].normal_form(form);
      while (!image.is_zero()) {
        auto exp = image.term_exp(0);
        auto pos = pivots.find(exp);
        if (pos == pivots.end()) {
          pivots.emplace(exp, std::make_pair(image, form));
          break;
        }
        const auto& [other, lift] = pos->second;
        auto var = std::find(exp.begin(), exp.end(), 1L);
        // Images are affine linear. At a constant pivot the image is constant.
        auto coeff = [&](const Poly& p) {
          return var == exp.end() ? p : p.derivative(var - exp.begin());
        };
        auto a = coeff(other), b = coeff(image);
        image = a * image - b * other;
        form = a * form - b * lift;
      }
      if (image.is_zero()) kernel.push_back(form.primitive());
    }
    forms = std::move(kernel);
  }
  std::stable_sort(forms.begin(), forms.end(), [](const Poly& a, const Poly& b) {
    return a.nterms() < b.nterms();
  });
  for (const auto& form : forms)
    if (!X.contains(form)) return form.monic();
  return {};
}

// Express the equations in coordinates where the chosen monic form is x_p,
// homogenize, and use graded reverse lex with x_p last. In a homogeneous
// basis element its leading term has the smallest x_p exponent among its
// terms. Removing that power gives a basis of the saturation. Setting h=1
// and undoing the affine substitution recovers generators in the input ring.
Ideal chart_saturate(const Ideal& K, const Poly& form) {
  auto S = K.ring();
  const int n = S->nvars();
  auto lead = form.term_exp(0);
  const int pivot = std::find(lead.begin(), lead.end(), 1L) - lead.begin();
  auto names = S->names();
  std::string hname = "helmer_h";
  while (S->index_of(hname) >= 0) hname += "_";
  names.push_back(hname);
  auto H = Ring::make(names);
  auto ids = identity(n);
  std::vector<Poly> forward, back;
  for (int i = 0; i < n; ++i) {
    forward.push_back(Poly::var(H, i));
    back.push_back(Poly::var(S, i));
  }
  const bool tilted = form.nterms() > 1;
  if (tilted) {
    auto offset = form.embed(H, ids) - forward[pivot];  // form = x_p + offset
    forward[pivot] = forward[pivot] - offset;
    back[pivot] = form;
  }
  back.push_back(Poly::constant(S, 1));
  std::vector<Poly> homogeneous;
  for (const auto& f : K.gens()) {
    if (f.is_zero()) continue;
    auto p = tilted ? f.substitute(H, forward) : f.embed(H, ids);
    long degree = p.total_degree();
    for (long j = 0; j < p.nterms(); ++j) {
      auto e = p.term_exp(j);
      e[n] = degree - std::accumulate(e.begin(), e.end(), 0L);
      std::vector<ulong> exponents(e.begin(), e.end());
      fmpq_mpoly_set_term_exp_ui(p.raw(), j, exponents.data(), H->ctx());
    }
    fmpq_mpoly_sort_terms(p.raw(), H->ctx());
    homogeneous.push_back(std::move(p));
  }
  auto perm = identity(n + 1);
  perm.erase(perm.begin() + pivot);
  perm.push_back(pivot);
  auto G = groebner(H, homogeneous, perm, 0);
  std::vector<Poly> out;
  auto x = Poly::var(H, pivot);
  // After division by x_p the leading monomials need no longer be
  // minimal. A subset generating the same leading ideal is still a basis
  // of the saturation. Discard redundant elements before substitution.
  std::vector<size_t> order(G.polys.size());
  std::iota(order.begin(), order.end(), 0);
  std::vector<long> degrees;
  for (const auto& m : G.lead)
    degrees.push_back(std::accumulate(m.begin(), m.end(), 0L) - m[pivot]);
  std::stable_sort(order.begin(), order.end(),
                   [&](size_t a, size_t b) { return degrees[a] < degrees[b]; });
  std::vector<size_t> minimal;
  for (size_t i : order) {
    bool redundant = std::any_of(minimal.begin(), minimal.end(), [&](size_t j) {
      for (int k = 0; k <= n; ++k)
        if (k != pivot && G.lead[j][k] > G.lead[i][k]) return false;
      return true;
    });
    if (!redundant) minimal.push_back(i);
  }
  // Keep the original order for the downstream algebra.
  std::sort(minimal.begin(), minimal.end());
  for (size_t i : minimal) {
    auto& p = G.polys[i];
    if (int e = G.lead[i][pivot]) {
      Poly power(H);
      fmpq_mpoly_pow_ui(power.raw(), x.raw(), e, H->ctx());
      p = divexact(p, power);
    }
    out.push_back(p.substitute(S, back));
  }
  return Ideal(S, std::move(out));
}

}  // namespace

static Ideal conormal_impl(const Ideal& I, const RingPtr& S, const std::string& saturation,
                           const std::optional<Poly>& chart, bool grouped = false) {
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
  if (chart) return chart_saturate(Ideal(S, std::move(gens)), chart->embed(S, idm));
  auto cminors = minors(c, J1);
  if (saturation == "minor") {
    // I prime: the conormal over the dense open set {g != 0} of X is smooth and
    // reduced, so K : g^oo is exactly I(Con X) for any c-minor g not in I.
    // Use the sparsest such minor.
    std::sort(cminors.begin(), cminors.end(),
              [](const Poly& a, const Poly& b) { return a.nterms() < b.nterms(); });
    for (auto& g : cminors)
      if (grouped ? (I + g).dim() < I.dim() : !I.contains(g))
        return saturate(Ideal(S, gens), g.embed(S, idm));
    if (!grouped) throw std::invalid_argument("helmer: conormal: every c-minor lies in I (not reduced?)");
  }
  std::vector<Poly> sat;
  for (auto& p : cminors) sat.push_back(p.embed(S, idm));
  return saturate(Ideal(S, gens), Ideal(S, sat), grouped ? "exact" : saturation);
}

Ideal conormal(const Ideal& I, const RingPtr& S, const std::string& saturation) {
  return conormal_impl(I, S, saturation, {});
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

// Sparse column elimination over QQ. The first rows are polynomial images,
// with one scalar row per monomial. The trailing polynomials carry a lift
// through the same operations. A dependent column leaves its kernel vector
// in the trailing entries.
using ColumnPivots = std::map<std::pair<size_t, std::vector<long>>, std::vector<Poly>>;

bool independent_column(std::vector<Poly>& column, size_t images, ColumnPivots& pivots) {
  for (size_t i = 0; i < images; ++i) {
    while (!column[i].is_zero()) {
      const Poly& p = column[i];
      auto row = std::make_pair(i, p.term_exp(0));
      fmpq_t lc;
      fmpq_init(lc);
      fmpq_mpoly_get_term_coeff_fmpq(lc, p.raw(), 0, p.ring()->ctx());
      Poly coefficient(p.ring());
      fmpq_mpoly_set_fmpq(coefficient.raw(), lc, p.ring()->ctx());
      fmpq_clear(lc);
      auto it = pivots.find(row);
      if (it == pivots.end()) {
        for (size_t h = i; h < column.size(); ++h)
          column[h] = divexact(column[h], coefficient);
        pivots.emplace(std::move(row), std::move(column));
        return true;
      }
      for (size_t h = i; h < column.size(); ++h)
        column[h] = column[h] - coefficient * it->second[h];
    }
  }
  return false;
}

// Eliminate the derivative columns over QQ, carrying each column's direction
// with it. A dependency has the form e_j + sum_p a_p e_p, with p among the
// pivot columns. Thus x_p - sum_j a_p x_j are coordinates on the quotient.
std::optional<Stratification> affine_cylinder(const Ideal& I, const WhitneyOptions& opts, bool defer_primes) {
  const auto& R = I.ring();
  const int n = R->nvars();
  if (n < 2 || I.gens().empty()) return {};

  // The row index includes the generator: equal monomials in different
  // equations are independent constraints on a translation direction.
  ColumnPivots echelon;
  std::vector<int> pivots;
  std::vector<std::pair<int, Poly>> kernel;
  for (int j = 0; j < n; ++j) {
    std::vector<Poly> column;
    for (const auto& f : I.gens()) column.push_back(f.derivative(j));
    column.push_back(Poly::var(R, j));
    if (independent_column(column, I.gens().size(), echelon)) pivots.push_back(j);
    else kernel.emplace_back(j, std::move(column.back()));
  }
  if (pivots.empty() || kernel.empty()) return {};

  std::vector<std::string> names;
  for (int p : pivots) names.push_back(R->names()[p]);
  auto base = Ring::make(std::move(names));
  std::vector<Poly> section(n, Poly(base)), quotient;
  for (int i = 0; i < base->nvars(); ++i) {
    int p = pivots[i];
    section[p] = Poly::var(base, i);
    Poly form = Poly::var(R, p);
    for (const auto& [j, v] : kernel)
      form = form - v.derivative(p) * Poly::var(R, j);
    quotient.push_back(std::move(form));
  }

  const int shift = n - base->nvars();
  if (opts.log) opts.log("affine cylinder: removing " + std::to_string(shift) + (shift == 1 ? " free direction" : " free directions"));
  WhitneyOptions base_opts = opts;
  std::vector<Poly> equations;
  for (const auto& f : I.gens()) equations.push_back(f.substitute(base, section));
  Ideal base_ideal(base, std::move(equations));
  if (opts.level >= 0) base_opts.level = std::min(opts.level, base_ideal.dim());
  auto reduced = stratify_impl(base_ideal, base_opts, defer_primes);
  Stratification out{R, false, {}};
  if (reduced.strata.empty()) return out;  // empty variety
  for (int d = 0; d < shift; ++d) out.strata[d];
  for (const auto& [d, pieces] : reduced.strata) {
    auto& lifted = out.strata[d + shift];
    for (const auto& piece : pieces) {
      std::vector<Poly> gens;
      for (const auto& g : piece.gens()) gens.push_back(g.substitute(R, quotient));
      lifted.emplace_back(R, std::move(gens));
    }
  }
  return out;
}

// In dimension zero adjoining a squarefree annihilator for each coordinate
// gives a quotient of a finite etale QQ-algebra, hence a reduced algebra.
Ideal reduced_points(const Ideal& I) {
  Ideal result = I;
  auto vars = identity(I.ring()->nvars());
  for (int v : vars) {
    auto others = vars;
    others.erase(others.begin() + v);
    auto univariate = eliminate(I, others);
    const auto& basis = univariate.gb().polys;
    if (basis.size() != 1 || basis.front().is_constant())
      throw std::logic_error("expected a nonconstant coordinate annihilator");
    result = result + radical_part(basis.front());
  }
  result.canonicalise();
  return result;
}

// On a smooth curve, a nonzero tangent field has an open local orbit.
// If its flow preserves X and Z, Whitney (b) is invariant along this orbit.
// Generic Whitney regularity on each component then gives regularity there.
// The common zero set of all affine fields preserving the pair bounds the
// exceptional locus. A positive-dimensional bound gives no shortcut.
std::optional<Ideal> curve_orbit_bound(const Ideal& X, const Ideal& Z) {
  const auto& R = X.ring();
  const int n = R->nvars();
  const size_t constraints = X.gens().size() + Z.gens().size();
  ColumnPivots pivots;
  std::vector<Poly> zeros;
  for (int j = 0; j < n; ++j) {
    for (int k = -1; k < n; ++k) {
      Poly a = k < 0 ? Poly::constant(R, 1) : Poly::var(R, k);
      // Column for a*d/dx_j. Carry the field coefficients after its images.
      std::vector<Poly> column;
      for (const Ideal* L : {&X, &Z})
        for (const auto& f : L->gens())
          column.push_back(L->normal_form(a * f.derivative(j)));
      column.resize(constraints + n, Poly(R));
      column[constraints + j] = a;
      if (independent_column(column, constraints, pivots)) continue;
      for (size_t h = constraints; h < column.size(); ++h) {
        auto value = Z.normal_form(column[h]);
        if (value.is_constant() && !value.is_zero())
          return Ideal(R, {Poly::constant(R, 1)});
        if (!value.is_zero()) zeros.push_back(std::move(value));
      }
    }
  }
  if (zeros.empty()) return {};
  Ideal bound = Z + Ideal(R, std::move(zeros));
  if (bound.dim() > 0) return {};
  return bound;
}

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

static Stratification stratify_impl(const Ideal& I0, const WhitneyOptions& opts, bool defer_primes) {
  auto log = [&](const std::string& s) {
    if (opts.log) opts.log(s);
  };
  const RingPtr& R = I0.ring();
  if (!opts.projective) {
    if (auto cylinder = affine_cylinder(I0, opts, defer_primes)) return std::move(*cylinder);
  }
  const int nv = R->nvars();
  const int n = opts.projective ? nv - 1 : nv;
  Store st{n, {}, {}};
  const int dimI = n - I0.codim();
  for (int i = 0; i <= std::max(dimI, 0); ++i) {
    st.V[i];
    st.keys[i];
  }
  Stratification out{R, opts.projective, {}};
  if (dimI < 0) return out;

  // Height equal to the number of equations implies a complete intersection.
  // It is unmixed. Generic smoothness then establishes reducedness (R0 + S1).
  if (defer_primes) log("check top");
  Ideal candidate = I0;
  if (defer_primes && I0.gens().size() == 1)
    candidate = Ideal(R, {radical_part(I0.gens().front())});
  std::optional<Ideal> grouped_jacobian;
  if (defer_primes && candidate.codim() == static_cast<int>(candidate.gens().size())) {
    Ideal J = jacobian_ideal(candidate);
    if (J.dim() < candidate.dim()) grouped_jacobian = std::move(J);
  }
  log(grouped_jacobian ? "group top" : "decompose X");
  auto comps = grouped_jacobian ? std::vector<Ideal>{candidate} : minimal_primes(I0);
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

  auto lower_parts = [&](const Ideal& K) -> std::vector<Ideal> {
    if (defer_primes && K.is_unit()) return {};
    if (defer_primes && std::all_of(K.gb().polys.begin(), K.gb().polys.end(),
                                  [](const Poly& p) { return p.total_degree() <= 1; }))
      return {K};  // A proper affine-linear ideal is already prime.
    if (defer_primes && !opts.projective && K.dim() == 0) {
      log("group points");
      return {reduced_points(K)};
    }
    if (defer_primes) log("lower prime fallback");
    return minimal_primes(K);
  };

  if (equidim) {
    // The Jacobian criterion applies to the reduced equidimensional variety.
    if (comps.empty()) return finish();
    Ideal reduced = comps.front();
    for (size_t i = 1; i < comps.size(); ++i) reduced = intersect(reduced, comps[i]);
    Ideal J = grouped_jacobian ? *grouped_jacobian : jacobian_ideal(reduced);
    if (dimI < 2) {
      for (auto& K : lower_parts(J))
        if (st.dim_of(K) < dimI) st.add(K);
      return finish();
    }
    const int codimJ = J.codim();
    if (codimJ == n)
      for (auto& K : lower_parts(J)) st.add(K);
    if (codimJ >= n) return finish();
  }

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
      bool grouped = grouped_jacobian.has_value() && a == 0;
      Ideal Ja = grouped ? *grouped_jacobian : jacobian_ideal(Ia);
      const int mu = st.dim_of(Ja);  // dim of Sing(Ia); negative if smooth
      auto singular = lower_parts(Ja);
      for (auto& K : singular)
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
        std::optional<Ideal> C;
        for (int ell = mu; ell >= 0; --ell) {
          // Points have no smaller exceptional locus to add to a filtration.
          if (defer_primes && ell == 0) continue;
          // iterate by index: entries appended at dimension < ell do not affect this list
          for (size_t zi = 0; zi < st.V[ell].size(); ++zi) {
            Ideal Z = st.V[ell][zi];
            if (!Z.contains(Ia)) continue;  // need V(Z) ⊆ V(Ia)
            if (defer_primes && !opts.projective && ell == 1) {
              if (auto bound = curve_orbit_bound(Ia, Z)) {
                log("  curve orbit bound");
                for (auto& K : lower_parts(*bound)) st.add(K);
                continue;
              }
            }
            if (!C) {
              log("  conormal of a " + std::to_string(k) + "-dimensional piece");
              auto chart = opts.saturation == "random" ? chart_form(Ia, singular) : std::nullopt;
              // On a pure reduced union the open chart must meet every component.
              if (grouped && chart && (Ia + *chart).dim() >= Ia.dim()) chart.reset();
              if (chart) log("  conormal chart: " + chart->str());
              const auto saturation = opts.saturation == "random" ? "minor" : opts.saturation;
              C = conormal_impl(Ia, S, saturation, chart, grouped);
              if (defer_primes) log("conormal complete");
            }
            Ideal W = whitney_points(Z, *C, R);
            if (W.is_unit()) continue;
            for (auto& K : lower_parts(W)) st.add(K);
          }
        }
      }
    }
  }
  return finish();
}

Stratification whitney_stratify(const Ideal& I, const WhitneyOptions& opts) {
  return stratify_impl(I, opts, false);
}

std::map<int, Ideal> whitney_filtration(const Ideal& I, bool projective, bool defer_primes,
                                        std::function<void(const std::string&)> log) {
  WhitneyOptions opts;
  opts.projective = projective;
  opts.log = std::move(log);
  auto result = stratify_impl(I, opts, defer_primes);
  std::map<int, Ideal> filtration;
  Ideal acc(I.ring(), {Poly::constant(I.ring(), 1)});
  for (const auto& [d, parts] : result.strata) {
    for (const auto& part : parts) {
      if (acc.contains(part)) acc = part;
      else if (!part.contains(acc)) acc = intersect(acc, part);
    }
    filtration.emplace(d, acc);
  }
  return filtration;
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
