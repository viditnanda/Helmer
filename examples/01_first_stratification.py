"""A first Whitney stratification, and how to read it.

Run:  python3 examples/01_first_stratification.py
"""
import helmer

# A polynomial ring over Q, an ideal, and its Whitney stratification.
R = helmer.Ring("x, y, z")
X = helmer.Ideal(R, ["y^2*z - x^2"])           # the Whitney umbrella
WS = helmer.whitney_stratify(X)

# WS maps each dimension d to a list of prime ideals: the closures of the d-dimensional strata.
for d in sorted(WS, reverse=True):
    for P in WS[d]:
        print(f"dim {d}: V({', '.join(P.gens)})")
# The open strata are differences: here the surface minus the z-axis (dim 2), the z-axis minus the origin
# (dim 1), and the origin (dim 0). The singular locus is the whole z-axis, but the origin is a stratum of its
# own: Whitney's condition (b) fails there, because the "handle" of the umbrella (z < 0) meets the axis.

# Over the reals the same polynomials give a valid stratification: some strata may have no real points, or
# real points of lower dimension. Here, real points with z < 0 on the surface force x = y = 0: the handle.

# A curve with a cusp and a node, a surface with a line of nodes, a cone:
for gens in (["y^2 - x^3"], ["y^2 - x^2*(x + 1)"], ["x*y"], ["x^2 + y^2 - z^2"]):
    S = helmer.whitney_stratify(helmer.Ideal(R if "z" in gens[0] else helmer.Ring("x, y"), gens))
    print(gens[0], "->", {d: [", ".join(P.gens) for P in L] for d, L in sorted(S.items()) if L})

# Varieties given by several equations: the curve where a sphere meets the two planes x = 0 and y = 0.
S = helmer.whitney_stratify(helmer.Ideal(R, ["x^2 + y^2 + z^2 - 1", "x*y"]))
print("sphere and xy = 0:", {d: [", ".join(P.gens) for P in L] for d, L in sorted(S.items()) if L})
# Two great circles meeting in the two points (0, 0, 1) and (0, 0, -1). Strata are prime ideals over Q, so a
# single one may have several points, conjugate or real (as V(x, y, z^2 - 2) would).

# Projective varieties: homogeneous equations and projective=True. The Cayley cubic has four nodes.
P3 = helmer.Ring("x, y, z, w")
S = helmer.whitney_stratify(helmer.Ideal(P3, ["x*y*z + x*y*w + x*z*w + y*z*w"]), projective=True)
print("Cayley cubic in P^3:", len(S[0]), "singular points,", len(S[2]), "surface")

# Progress messages, for long computations:
helmer.whitney_stratify(X, log=lambda message: print("  log:", message))
