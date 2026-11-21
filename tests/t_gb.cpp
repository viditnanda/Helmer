#include "helmer/gb.hpp"
#include <iostream>
using namespace helmer;
int main(){
  auto R = Ring::make({"x","y","z"});
  std::vector<Poly> F = {Poly(R,"x^2*y - 3/2*z^3 + 7"), Poly(R,"x*y^2 - z + 11/3"), Poly(R,"z^2*x - y")};
  GB G = groebner_drl(R, F);
  for (size_t i=0;i<G.polys.size();++i){ std::cout << G.polys[i].str() << "   lead:"; for(int e:G.lead[i]) std::cout<<" "<<e; std::cout << "  flint-first-term:"; for(long e: G.polys[i].term_exp(0)) std::cout<<" "<<e; std::cout<<"\n";}
  // elimination of x
  GB E = groebner(R, F, {0,1,2}, 1);
  std::cout << "elim x:\n"; for (auto&p:E.polys) if(!p.involves(0)) std::cout<<"  "<<p.str()<<"\n";
  std::cout << "runs " << gb_stats().f4_runs << " unlucky " << gb_stats().unlucky << "\n";
}
