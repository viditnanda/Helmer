#include "helmer/primes.hpp"
#include <iostream>
#include <chrono>
using namespace helmer;
int main(int argc, char** argv){
  // usage: t_primes "x,y,z" "gen1" "gen2" ...
  std::vector<std::string> names; { std::string s=argv[1], cur; for(char c: s){ if(c==','){names.push_back(cur);cur.clear();} else cur+=c;} names.push_back(cur);} 
  auto R = Ring::make(names);
  std::vector<std::string> g(argv+2, argv+argc);
  Ideal I = Ideal::parse(R, g);
  auto t0=std::chrono::steady_clock::now();
  auto P = minimal_primes(I);
  double s=std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count();
  std::cout << "dim " << I.dim() << ", " << P.size() << " primes (" << s << "s; f4 runs " << gb_stats().f4_runs << ", z-retries " << prime_stats().z_retries << ")\n";
  for (auto& p: P) std::cout << "  [dim " << p.dim() << "] " << p.str() << "\n";
}
