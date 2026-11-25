#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "helmer/certify.hpp"
#include "helmer/multiplicity.hpp"
#include "helmer/whitney.hpp"

namespace py = pybind11;
using namespace helmer;

namespace {

RingPtr ring_from(py::object vars) {
  std::vector<std::string> names;
  if (py::isinstance<py::str>(vars)) {
    std::string s = vars.cast<std::string>();
    std::string cur;
    for (char c : s) {
      if (c == ',' || c == ' ') {
        if (!cur.empty()) names.push_back(cur);
        cur.clear();
      } else {
        cur += c;
      }
    }
    if (!cur.empty()) names.push_back(cur);
  } else {
    names = vars.cast<std::vector<std::string>>();
  }
  return Ring::make(names);
}

Ideal ideal_from(const RingPtr& R, py::object gens) {
  if (py::isinstance<py::str>(gens)) return Ideal::parse(R, {gens.cast<std::string>()});
  std::vector<Poly> G;
  for (auto g : gens) {
    if (py::isinstance<Poly>(g)) G.push_back(g.cast<Poly>());
    else G.emplace_back(R, py::str(g).cast<std::string>());
  }
  return Ideal(R, G);
}

std::vector<int> var_indices(const RingPtr& R, const std::vector<std::string>& vars) {
  std::vector<int> idx;
  for (auto& v : vars) {
    int i = R->index_of(v);
    if (i < 0) throw std::invalid_argument("unknown variable " + v);
    idx.push_back(i);
  }
  return idx;
}

Poly poly_from(const RingPtr& R, py::handle f) {
  if (py::isinstance<Poly>(f)) return f.cast<Poly>();
  return Poly(R, py::str(f).cast<std::string>());
}

py::object certificate_obj(const std::optional<MembershipCertificate>& C) {
  if (!C) return py::none();
  py::dict d;
  d["target"] = C->target;
  d["generators"] = C->gens;
  d["cofactors"] = C->cofactors;
  d["bounds"] = py::make_tuple(C->bound_block, C->bound_rest);
  return d;
}

CertifyOptions certify_options(long max_unknowns, int max_primes) {
  CertifyOptions o;
  o.max_unknowns = max_unknowns;
  o.max_primes = max_primes;
  return o;
}

py::dict strata_dict(const Stratification& S) {
  py::dict d;
  for (auto& [dim, L] : S.strata) d[py::int_(dim)] = L;
  return d;
}

Stratification strata_from(py::dict d, bool projective) {
  Stratification S{nullptr, projective, {}};
  for (auto item : d) {
    const int dim = item.first.cast<int>();
    auto& L = S.strata[dim];
    for (auto I : item.second) {
      L.push_back(I.cast<Ideal>());
      if (!S.ring) S.ring = L.back().ring();
    }
  }
  if (!S.ring) throw std::invalid_argument("empty stratification");
  return S;
}

LocalOptions local_options(int trials) {
  LocalOptions o;
  o.trials = trials;
  o.max_runs = std::max(8, 4 * trials);
  return o;
}

std::function<void(const std::string&)> logger(py::object log) {
  if (log.is_none()) return {};
  auto f = log.cast<std::function<void(const std::string&)>>();
  return [f](const std::string& s) { py::gil_scoped_acquire g; f(s); };
}

}  // namespace

PYBIND11_MODULE(_helmer, m) {
  m.doc() = "Helmer: Whitney stratifications of algebraic varieties (C++ core)";

  py::class_<Ring, RingPtr>(m, "Ring")
      .def(py::init([](py::object vars) { return ring_from(vars); }), py::arg("variables"),
           "Polynomial ring over QQ, e.g. Ring('x,y,z') or Ring(['x','y','z']).")
      .def_property_readonly("variables", &Ring::names)
      .def_property_readonly("nvars", &Ring::nvars)
      .def("__repr__", [](const Ring& R) {
        std::string s = "QQ[";
        for (int i = 0; i < R.nvars(); ++i) s += (i ? ", " : "") + R.names()[i];
        return s + "]";
      });

  py::class_<Poly>(m, "Poly")
      .def(py::init<RingPtr, const std::string&>())
      .def("__str__", &Poly::str)
      .def("__repr__", &Poly::str)
      .def("__add__", &Poly::operator+)
      .def("__sub__", [](const Poly& a, const Poly& b) { return a - b; })
      .def("__mul__", &Poly::operator*)
      .def("__neg__", [](const Poly& a) { return -a; })
      .def("__eq__", &Poly::equal)
      .def_property_readonly("ring", &Poly::ring)
      .def("total_degree", &Poly::total_degree)
      .def("factor", [](const Poly& p) { return factor(p); });

  py::class_<Ideal>(m, "Ideal")
      .def(py::init([](RingPtr R, py::object gens) { return ideal_from(R, gens); }),
           py::arg("ring"), py::arg("generators"))
      .def_property_readonly("ring", &Ideal::ring)
      .def_property_readonly("gens", [](const Ideal& I) {
        std::vector<std::string> s;
        for (auto& g : I.gens()) s.push_back(g.str());
        return s;
      })
      .def("groebner_basis", [](const Ideal& I) {
        std::vector<std::string> s;
        for (auto& g : I.gb().polys) s.push_back(g.str());
        return s;
      }, "Reduced Gröbner basis (degree reverse lexicographic), as strings.")
      .def("dim", &Ideal::dim)
      .def("codim", &Ideal::codim)
      .def("is_unit", &Ideal::is_unit)
      .def("independent_set", [](const Ideal& I) {
        std::vector<std::string> s;
        for (int v : I.independent_set()) s.push_back(I.ring()->names()[v]);
        return s;
      })
      .def("contains", [](const Ideal& I, const Ideal& J) { return I.contains(J); })
      .def("contains_poly", [](const Ideal& I, const std::string& f) {
        return I.contains(Poly(I.ring(), f));
      })
      .def("__eq__", &Ideal::equals)
      .def("__add__", [](const Ideal& I, const Ideal& J) { return I + J; })
      .def("__str__", &Ideal::str)
      .def("__repr__", &Ideal::str)
      .def("__hash__", [](const Ideal& I) { return py::hash(py::str(I.key())); });

  m.def("minimal_primes", &minimal_primes, py::arg("ideal"),
        "Minimal primes over QQ (Macaulay2's `decompose`).");
  m.def("eliminate", [](const Ideal& I, const std::vector<std::string>& vars) {
    std::vector<int> idx;
    for (auto& v : vars) {
      int i = I.ring()->index_of(v);
      if (i < 0) throw std::invalid_argument("unknown variable " + v);
      idx.push_back(i);
    }
    return eliminate(I, idx);
  }, py::arg("ideal"), py::arg("variables"));
  m.def("saturate", [](const Ideal& I, const Ideal& J, const std::string& strategy) {
    return saturate(I, J, strategy);
  }, py::arg("ideal"), py::arg("by"), py::arg("strategy") = "random");
  m.def("intersect", &intersect);
  m.def("conormal", [](const Ideal& I, const std::string& sat) {
    return conormal(I, conormal_ring(I.ring()), sat);
  }, py::arg("ideal"), py::arg("saturation") = "random",
     "Conormal variety in k^n x P^(n-1) (variables v_0..v_{n-1}).");
  m.def("whitney_stratify", [](const Ideal& I, bool projective, int level,
                               const std::string& sat, py::object log,
                               const std::string& method, long polar_bound, bool minimal) {
    WhitneyOptions o;
    o.method = method;
    o.polar_bound = polar_bound;
    o.projective = projective;
    o.level = level;
    o.saturation = sat;
    if (!log.is_none()) {
      auto f = log.cast<std::function<void(const std::string&)>>();
      o.log = [f](const std::string& s) { py::gil_scoped_acquire g; f(s); };
    }
    auto S = whitney_stratify(I, o);
    if (minimal) S = min_coarsen(S, LocalOptions{}, o.log);
    return strata_dict(S);
  }, py::arg("ideal"), py::arg("projective") = false, py::arg("level") = -1,
     py::arg("saturation") = "random", py::arg("log") = py::none(),
     py::arg("method") = "conormal", py::arg("polar_bound") = 1000, py::arg("minimal") = false,
     "Whitney stratification. method='conormal' (Helmer–Mohr; default) or 'polar' "
     "(polar varieties; usually faster, randomised, generally non-minimal). "
     "minimal=True coarsens the result to the minimal Whitney stratification (min_coarsen).");
  m.def("min_coarsen", [](py::dict strata, bool projective, int trials, py::object log) {
    return strata_dict(min_coarsen(strata_from(strata, projective), local_options(trials), logger(log)));
  }, py::arg("strata"), py::arg("projective") = false, py::arg("trials") = 2,
     py::arg("log") = py::none(),
     "The minimal (coarsest) Whitney stratification refined by `strata` (a Whitney "
     "stratification, e.g. from whitney_stratify), by Teissier's criterion: a stratum is "
     "merged into the unique smallest stratum above it when that stratum's closure is "
     "smooth along it and every larger stratum has the same polar multiplicities along "
     "both (M2: minCoarsenWS). Monte Carlo; see multiplicity.");
  m.def("multiplicity", [](const Ideal& X, const Ideal& Y, int trials) {
    return multiplicity(X, Y, local_options(trials));
  }, py::arg("X"), py::arg("Y"), py::arg("trials") = 2,
     "Samuel multiplicity e_X(Y) of V(Y) along the irreducible subvariety V(X): the "
     "multiplicity of Y at a general point of X, counting the components of Y of largest "
     "dimension through X (0 if X is not in V(Y)). Homogeneous input is treated as affine "
     "cones (same multiplicities). Monte Carlo: computed at a general F_p-point for random "
     "31-bit primes p, repeated until `trials` independent runs agree.");
  m.def("contained", [](const Ideal& X, const Ideal& Y, int trials) {
    return contained(X, Y, local_options(trials));
  }, py::arg("X"), py::arg("Y"), py::arg("trials") = 2,
     "V(X) ⊆ V(Y), i.e. Y lies in the radical of X (Monte Carlo, modulo random primes). "
     "For containment of ideals use X.contains(Y).");
  m.def("in_singular_locus", [](const Ideal& X, const Ideal& Y, int trials) {
    return in_singular_locus(X, Y, local_options(trials));
  }, py::arg("X"), py::arg("Y"), py::arg("trials") = 2,
     "V(X) ⊆ Sing V(Y) for X irreducible and Y prime (e_X(Y) > 1); no singular-locus "
     "equations are computed.");
  m.def("polar_multiplicity", [](const Ideal& X, const Ideal& Z, int dim, int trials) {
    return polar_multiplicity(X, Z, dim, local_options(trials));
  }, py::arg("X"), py::arg("Z"), py::arg("dim"), py::arg("trials") = 2,
     "Entry `dim` of polar_sequence(X, Z), computed alone.");
  m.def("polar_sequence", [](const Ideal& X, const Ideal& Z, int trials) {
    return polar_sequence(X, Z, local_options(trials));
  }, py::arg("X"), py::arg("Z"), py::arg("trials") = 2,
     "Multiplicities at a general point of V(X) of the generic polar varieties of the "
     "prime Z, by dimension 0..dim Z (the last entry is e_X(Z)); 0 where X is not in the "
     "polar variety. Teissier: Whitney's conditions hold for (Z_reg, Y) at general points "
     "of Y iff this sequence is constant along Y.");
  m.def("polar_varieties", [](const Ideal& I, long bound) {
    return polar_varieties(I, jacobian_ideal(I), bound);
  }, py::arg("ideal"), py::arg("bound") = 1000,
     "Polar varieties of a prime ideal (M2: polarVars).");
  m.def("map_stratify", [](py::object F, const Ideal& X, const Ideal& Y, bool proper,
                           const std::string& strata, py::object log) {
    std::vector<Poly> Fp;
    for (auto f : F) {
      if (py::isinstance<Poly>(f)) Fp.push_back(f.cast<Poly>());
      else Fp.emplace_back(X.ring(), py::str(f).cast<std::string>());
    }
    MapOptions o;
    o.proper = proper;
    o.strata = strata;
    if (!log.is_none()) {
      auto fn = log.cast<std::function<void(const std::string&)>>();
      o.log = [fn](const std::string& s) { py::gil_scoped_acquire g; fn(s); };
    }
    auto r = map_stratify(Fp, X, Y, o);
    return py::make_tuple(strata_dict(r.source), strata_dict(r.target));
  }, py::arg("F"), py::arg("X"), py::arg("Y"), py::arg("proper") = true,
     py::arg("strata") = "all", py::arg("log") = py::none(),
     "Stratify f: X -> Y (F = components of f, as polynomials in X's ring). Returns "
     "(strata of X, strata of Y).");
  m.def("non_proper_set", [](py::object F, const Ideal& X, const Ideal& Y) {
    std::vector<Poly> Fp;
    for (auto f : F) Fp.emplace_back(X.ring(), py::str(f).cast<std::string>());
    return non_proper_set(Fp, X, Y);
  });
  m.def("certify_membership", [](const Ideal& I, py::object target, const std::vector<std::string>& block,
                                 long max_unknowns, int max_primes) {
    return certificate_obj(certify_membership(I.gens(), poly_from(I.ring(), target),
                                              var_indices(I.ring(), block),
                                              certify_options(max_unknowns, max_primes)));
  }, py::arg("ideal"), py::arg("target"), py::arg("block") = std::vector<std::string>(),
     py::arg("max_unknowns") = 4000, py::arg("max_primes") = 64,
     "A verified certificate target = sum_j c_j g_j over the generators g_j of `ideal`, as a "
     "dict (target, generators, cofactors, bounds), or None if none was found within the "
     "limits. Cofactor terms have degree <= D_A in the `block` variables and <= D_B in the "
     "others, for a small range of (D_A, D_B). The identity is checked exactly over QQ "
     "before returning; None proves nothing.");
  m.def("verify_membership", [](py::object gens, py::object target, py::object cofactors,
                                std::optional<RingPtr> ring) {
    const bool is_ideal = py::isinstance<Ideal>(gens);
    if (!is_ideal && !ring) throw std::invalid_argument("ring is required unless generators is an Ideal");
    const RingPtr R = is_ideal ? gens.cast<const Ideal&>().ring() : *ring;
    MembershipCertificate C{poly_from(R, target), {}, {}};
    for (auto c : cofactors) C.cofactors.push_back(poly_from(R, c));
    if (is_ideal) {
      const Ideal& I = gens.cast<const Ideal&>();
      C.gens = I.gens();
      return verify(I, C);
    }
    for (auto g : gens) C.gens.push_back(poly_from(R, g));
    return verify(C);
  }, py::arg("generators"), py::arg("target"), py::arg("cofactors"), py::arg("ring") = py::none(),
     "Exact check of target == sum_j cofactors[j] * g_j over QQ. `generators` is an Ideal "
     "(then g_j are its generators, the check proves target in it, and `ring` may be "
     "omitted) or a list of polynomials in `ring`.");
  m.def("eliminate_certified", [](const Ideal& I, const std::vector<std::string>& vars,
                                  long max_unknowns, int max_primes) {
    auto E = eliminate_certified(I, var_indices(I.ring(), vars),
                                 certify_options(max_unknowns, max_primes));
    py::list certs;
    for (auto& C : E.certificates) certs.append(certificate_obj(C));
    return py::make_tuple(E.ideal, certs);
  }, py::arg("ideal"), py::arg("variables"), py::arg("max_unknowns") = 4000,
     py::arg("max_primes") = 64,
     "eliminate(ideal, variables) together with, for each generator, a verified membership "
     "certificate in `ideal` (None where none was found). A certified generator vanishes "
     "on the zero set of `ideal`. Only a nonzero certified generator bounds anything: if "
     "the elimination ideal comes back zero, both lists are empty.");
  m.def("set_seed", &set_seed);
  m.def("stats", []() {
    py::dict d;
    d["gb_calls"] = gb_stats().calls;
    d["f4_runs"] = gb_stats().f4_runs;
    d["unlucky_primes"] = gb_stats().unlucky;
    d["gb_seconds"] = gb_stats().seconds;
    d["prime_calls"] = prime_stats().calls;
    d["factor_splits"] = prime_stats().splits;
    d["z_retries"] = prime_stats().z_retries;
    return d;
  });
}
