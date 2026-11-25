"""Certified elimination: polynomial identities that prove where the critical values of a function lie.

eliminate_certified returns eliminate(I, variables) together with, for each generator P, cofactors c_j such
that P = sum_j c_j g_j over the generators g_j of I. The identity is checked exactly over Q, so it proves
P lies in I, whatever happened in the modular arithmetic. Run:  python3 examples/04_certified.py
"""
import helmer

# The critical values of f = (x^2 + y^2 - 1)(x^2 + y^2 - 4) on the plane: eliminate x, y from
# (df/dx, df/dy, f - epsilon). Its critical points are the origin and a whole circle, not isolated.
R = helmer.Ring("x, y, epsilon")
f = "(x^2+y^2-1)*(x^2+y^2-4)"
I = helmer.Ideal(R, ["2*x*(2*x^2+2*y^2-5)", "2*y*(2*x^2+2*y^2-5)", f + " - epsilon"])
J, certificates = helmer.eliminate_certified(I, ["x", "y"])
print("elimination ideal:", J)
print("every critical value of f is a root of", J.gens[0], "(4 at the origin, -9/4 on a circle)")

cert = certificates[0]
print("certificate: target", cert["target"])
print("  cofactors:", cert["cofactors"])
print("  re-checked exactly:", helmer.verify_membership(I, cert["target"], cert["cofactors"]))

# Membership certificates for any polynomial: x^2 + y^2 - 5/2 vanishes on the critical circle.
C = helmer.Ring("x, y")
grad = helmer.Ideal(C, ["2*x*(2*x^2+2*y^2-5)", "2*y*(2*x^2+2*y^2-5)"])
c = helmer.certify_membership(grad, "x*(x^2+y^2-5/2)")
print("\nx (x^2 + y^2 - 5/2) =", " + ".join(f"({a})*({g})" for a, g in zip(c["cofactors"], c["generators"])))
