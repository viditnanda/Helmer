#include "helmer/gb.hpp"
#include <flint/ulong_extras.h>
#include <iostream>
#include <set>

using namespace helmer;

// Reproduce the first draws on this platform, including Helmer's primality
// filter, so the tests do not depend on a library's uniform distribution.
static std::vector<ulong> primes(unsigned seed) {
  set_seed(seed);
  std::uniform_int_distribution<ulong> odd(1UL << 29, (1UL << 30) - 1);
  std::vector<ulong> result;
  std::set<ulong> used;
  while (result.size() < 3) {
    ulong p = 2 * odd(gb_prime_rng()) + 1;
    if (n_is_prime(p) && used.insert(p).second) result.push_back(p);
  }
  return result;
}

static int check_traces() {
  auto R = Ring::make({"x", "y", "z", "w", "u"});
  auto ps = primes(42);
  for (int scenario : {0, 1, 2, 3}) {
    int unlucky = scenario == 3 ? 2 : scenario;
    // At p, the S-pair loses the y^3 pivot: a+b^2 = p.
    auto b = Poly(R, "23456789012345678912343");
    auto a = Poly::constant(R, unlucky ? ps[unlucky] : 1) - b*b;
    // With u^2 present the row survives at the bad prime, with leading
    // monomial y*u^2. Without it the row reduces to zero.
    auto first = Poly(R, scenario == 3 ? "x^2+u^2" : "x^2");
    std::vector<Poly> input = {first + a * Poly(R, "y^2"),
        Poly(R, "x*y") + b * Poly(R, "y^2"), Poly(R, "z-1234567891231451235123"),
        Poly(R, "w-498537293872395723837491")};
    gb_options().trace_min_ms = -1;
    set_seed(42);
    GB expected = groebner_drl(R, input);
    gb_options().trace_min_ms = 0;
    set_seed(42);
    auto before = gb_stats();
    GB actual = groebner_drl(R, input);
    if (actual.lead != expected.lead || actual.polys.size() != expected.polys.size()) return 1;
    for (size_t i = 0; i < actual.polys.size(); ++i)
      if (!actual.polys[i].equal(expected.polys[i])) return 1;
    if (gb_stats().trace_learns != before.trace_learns + 1 ||
        gb_stats().trace_fallbacks != before.trace_fallbacks + (unlucky != 0) ||
        (unlucky != 1 && gb_stats().trace_replays == before.trace_replays)) return 1;
  }
  auto before = gb_stats();
  groebner_drl(R, {Poly(R, "x^2+y"), Poly(R, "x*y+1")});
  if (gb_stats().trace_learns != before.trace_learns) return 1;
  std::cout << "F4 trace selection and recovery passed\n";
  return 0;
}

int main() {
  int result = check_traces();
  flint_cleanup_master();
  return result;
}
