// Shared by the Monte Carlo local and global invariants (not installed).
#pragma once

#include "helmer/modp.hpp"
#include "helmer/multiplicity.hpp"

#include <algorithm>
#include <string>
#include <vector>

namespace helmer::detail {

using modp::CtxPtr;
using modp::ModGB;
using modp::MPoly;

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

inline std::vector<MPoly> reduce_all(const Ideal& I, const CtxPtr& C) {
  std::vector<MPoly> out;
  for (auto& g : I.gens()) out.push_back(modp::reduce(g, C));
  return out;
}

inline std::vector<ulong> random_vector(ulong p, int n) {
  std::vector<ulong> v(n);
  for (auto& e : v) e = modp::random_elem(p);
  return v;
}

}  // namespace helmer::detail
