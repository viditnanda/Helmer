// Whitney stratification from C++: build with the rest of Helmer (cmake -DHELMER_EXAMPLES=ON), then
//   ./stratify "x,y,z" "y^2*z - x^2"            (variables, then one or more equations)
//   ./stratify "x,y,z,w" "x*y*z + x*y*w + x*z*w + y*z*w" --projective
#include "helmer/whitney.hpp"

#include <cstring>
#include <iostream>
#include <sstream>

int main(int argc, char** argv) {
  if (argc < 3) {
    std::cerr << "usage: " << argv[0] << " \"x,y,z\" \"equation\" [\"equation\" ...] [--projective]\n";
    return 1;
  }
  std::vector<std::string> names;
  std::stringstream vars(argv[1]);
  for (std::string v; std::getline(vars, v, ',');) names.push_back(v);
  auto R = helmer::Ring::make(names);

  helmer::WhitneyOptions opts;
  std::vector<helmer::Poly> equations;
  for (int i = 2; i < argc; ++i) {
    if (!std::strcmp(argv[i], "--projective")) opts.projective = true;
    else equations.emplace_back(R, argv[i]);           // throws on a parse error
  }
  helmer::Ideal X(R, equations);

  helmer::Stratification S = helmer::whitney_stratify(X, opts);
  // S.strata[d]: the prime ideals whose zero sets are the closures of the d-dimensional strata.
  for (auto it = S.strata.rbegin(); it != S.strata.rend(); ++it)
    for (const helmer::Ideal& P : it->second) std::cout << "dim " << it->first << ": " << P.str() << "\n";
  return 0;
}
