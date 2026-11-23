"""Run one F4 trace benchmark in a fresh process and report wall time and counters."""
import argparse
import hashlib
import json
from time import perf_counter

import helmer


def katsura(n):
    xs = [f"x{i}" for i in range(n + 1)]
    equations = [xs[0] + "+2*(" + "+".join(xs[1:]) + ")-1"]
    for j in range(n):
        terms = [f"{xs[abs(i)]}*{xs[abs(j-i)]}" for i in range(-n, n+1) if abs(j-i) <= n]
        equations.append("+".join(terms) + "-" + xs[j])
    return helmer.Ideal(helmer.Ring(xs), equations)


def cyclic(n):
    xs = [f"x{i}" for i in range(n)]
    equations = []
    for degree in range(1, n):
        terms = ["*".join(xs[(i+j) % n] for j in range(degree)) for i in range(n)]
        equations.append("+".join(terms))
    return helmer.Ideal(helmer.Ring(xs), equations + ["*".join(xs) + "-1"])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", choices=["bs", "katsura6", "katsura7", "cyclic6", "cyclic7"])
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()
    helmer.set_seed(args.seed)
    if args.case == "bs":
        ideal = helmer.Ideal(helmer.Ring("x,y,z,t"), ["x^5+y^7*z+z^15+t*x*y^6"])
    else:
        constructor = katsura if args.case.startswith("katsura") else cyclic
        ideal = constructor(int(args.case[-1]))
    start = perf_counter()
    if args.case == "bs":
        result = {d: v.groebner_basis() for d, v in helmer.whitney_filtration(ideal).items()}
    else:
        result = ideal.groebner_basis()
    seconds = perf_counter() - start
    digest = hashlib.sha256(json.dumps(result, sort_keys=True).encode()).hexdigest()
    print(json.dumps(dict(case=args.case, seconds=seconds, sha256=digest, **helmer.stats())))
