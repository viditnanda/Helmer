"""The commutative algebra underneath: Groebner bases, dimension, elimination, saturation, minimal primes.

Everything Macaulay2's top-level packages provided to WhitneyStratifications is implemented in Helmer and
available on its own. Run:  python3 examples/03_algebra.py
"""
import helmer

R = helmer.Ring("x, y, z")
I = helmer.Ideal(R, ["x^2 - y*z", "y^2 - x*z"])
print("I =", I, "| dim", I.dim(), "| codim", I.codim())
print("Groebner basis (degree reverse lexicographic):", I.groebner_basis())

# Minimal primes over Q (Macaulay2's `decompose`): here three components.
for P in helmer.minimal_primes(I):
    print("  minimal prime:", P, "of dimension", P.dim())

# Elimination: the image of the twisted cubic t -> (t, t^2, t^3) is cut out by eliminating t.
T = helmer.Ring("t, x, y, z")
graph = helmer.Ideal(T, ["x - t", "y - t^2", "z - t^3"])
print("twisted cubic:", helmer.eliminate(graph, ["t"]))

# Saturation removes components: (I : J^infinity) is the closure of V(I) minus V(J).
print("I : (x)^infinity =", helmer.saturate(I, helmer.Ideal(R, ["x"])))

# Intersections, sums, membership.
A, B = helmer.Ideal(R, ["x"]), helmer.Ideal(R, ["y"])
print("(x) ∩ (y) =", helmer.intersect(A, B), "| (x) + (y) =", A + B)
print("x*y in (x)?", A.contains_poly("x*y"), "| (x) contains (x^2, x*y)?", A.contains(helmer.Ideal(R, ["x^2", "x*y"])))

# Conormal varieties, in C^n x P^(n-1) with dual coordinates v_0, ..., v_{n-1}:
print("conormal of the cone x^2 + y^2 - z^2:", helmer.conormal(helmer.Ideal(R, ["x^2 + y^2 - z^2"])))
