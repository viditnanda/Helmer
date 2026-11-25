"""Singular configurations of a mechanism: the planar four-bar linkage.

Two cranks of lengths a and c, pivoted at (0, 0) and (d, 0), are joined by a coupler of length b. The
configuration space is the variety of the joint positions (x1, y1) and (x2, y2). Its singular strata are the
configurations where the mechanism can change mode or lock. Run:  python3 examples/05_linkages.py
"""
import helmer

R = helmer.Ring("x1, y1, x2, y2")


def four_bar(a, b, c, d):
    return helmer.Ideal(R, [f"x1^2 + y1^2 - {a * a}", f"(x2 - {d})^2 + y2^2 - {c * c}",
                            f"(x1 - x2)^2 + (y1 - y2)^2 - {b * b}"])


for lengths, name in [((1, 3, 2, 3), "generic (a, b, c, d) = (1, 3, 2, 3)"),
                      ((1, 3, 2, 2), "change-point (a + b = c + d): (1, 3, 2, 2)"),
                      ((1, 2, 1, 2), "parallelogram (a = c, b = d): (1, 2, 1, 2)")]:
    S = helmer.whitney_stratify(four_bar(*lengths))
    print(name)
    print("  curves (configuration space components):", len(S.get(1, [])))
    for P in S.get(0, []):
        print("  singular configuration: V(" + ", ".join(P.gens) + ")")
# The generic linkage is a smooth curve. At a change point the coupler and cranks can line up flat, and the
# configuration space acquires a singular point (1, 0, 4, 0). The parallelogram linkage has two modes
# (parallel and anti-parallel), two curves that cross at the two flat configurations.

# Equations in the differences of joint coordinates are invariant under translating every joint: Helmer
# detects such free directions and stratifies a slice (log message "affine cylinder").
T = helmer.Ring("x1, y1, x2, y2")
I = helmer.Ideal(T, ["(x1 - x2)^2 + (y1 - y2)^2 - 1", "(x1 - x2)*(y1 - y2)"])
helmer.whitney_stratify(I, log=lambda m: print("  log:", m) if "cylinder" in m else None)
