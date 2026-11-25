"""Two algorithms: the conormal route (default) and the polar-variety route.

Both are ports of Martin Helmer's Macaulay2 package. The conormal route is deterministic up to the modular
lift of Groebner bases; the polar route draws random projections and returns a valid, generally non-minimal
stratification. Run:  python3 examples/06_conormal_or_polar.py
"""
import time

import helmer

# Example X_1 of Helmer-Mohr (a complete intersection in C^4): Macaulay2 takes 210 s on the conormal route.
R = helmer.Ring("x_1, x_2, x_3, x_4")
X1 = helmer.Ideal(R, ["x_1^2*x_3 - x_2^2", "x_2^4 - x_1*x_2^2 - x_3*x_4^2", "x_1^2*x_2^2 - x_1^3 - x_4^2"])

for method in ("conormal", "polar"):
    helmer.set_seed(1)                       # reproducible random choices
    t = time.perf_counter()
    S = helmer.whitney_stratify(X1, method=method)
    print(f"{method:9s} {time.perf_counter() - t:6.2f} s  strata per dimension:", {d: len(L) for d, L in sorted(S.items())})

print("work done by the last call:", helmer.stats())
# The polar route's extra strata are artefacts of its random projections; both are Whitney stratifications.
