"""Stratifying a polynomial map so that Thom's isotopy lemma applies.

map_stratify(F, X, Y) stratifies X and Y so that F maps each stratum of X submersively onto a stratum of Y.
Over each stratum of Y, the fibres of F are then topologically the same. Run:  python3 examples/02_maps.py
"""
import helmer

# The roots of a x^2 + b x + c, as the coefficients (a, b, c) vary: project the incidence variety
# X = {(a, b, c, x) : a x^2 + b x + c = 0} to Y = C^3 (the coefficients).
R, S = helmer.Ring("a, b, c, x"), helmer.Ring("a, b, c")
X = helmer.Ideal(R, ["a*x^2 + b*x + c"])
Y = helmer.Ideal(S, [])                                  # all of C^3
source, target = helmer.map_stratify(["a", "b", "c"], X, Y, proper=False)
print("strata of the coefficient space:")
for d in sorted(target, reverse=True):
    for P in target[d]:
        print(f"  dim {d}: V({', '.join(P.gens) or '0'})")
# The discriminant b^2 - 4ac (a double root) and a = 0 (a root escapes to infinity) are where the number
# of roots changes. proper=False adds the set over which the map is not proper (here a = 0); for a proper
# map it would be empty.
print("not proper over:", helmer.non_proper_set(["a", "b", "c"], X, Y))

# The source strata: each lies over one target stratum.
print("source strata by dimension:", {d: len(L) for d, L in sorted(source.items())})
