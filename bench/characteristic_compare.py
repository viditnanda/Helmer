"""Helmer side of bench/m2_characteristic_checks.m2 (same inputs, same order)."""
import helmer
helmer.set_seed(1)
PROJ = [["D4_P", "x,y,z,w", ["x^2*w-y^3+y*z^2"]], ["cayley_P", "x,y,z,w", ["x*y*z+x*y*w+x*z*w+y*z*w"]], ["quartic_P", "a_0,a_1,a_2,a_3,a_4", ["a_0^2*a_4-a_1*a_2^2+a_3^3-a_3*a_0^2-a_4^2*a_3"]], ["umbrella_P", "x,y,z,w", ["y^2*z-x^2*w"]], ["nodal_cubic", "x,y,z", ["y^2*z-x^3-x^2*z"]], ["cuspidal_cubic", "x,y,z", ["y^2*z-x^3"]], ["cayley_cubic", "x,y,z,w", ["x*y*z+x*y*w+x*z*w+y*z*w"]]]
AFF = [["umbrella", "x,y,z", ["x^2-y^2*z"]], ["cubic_cone", "x,y,z", ["x^3+y^3+z^3"]], ["bs_small", "x,y,z,t", ["x^3+y^4*z+z^9+t*x*y^3"]], ["HM_X2", "x_1,x_2,x_3,x_4", ["x_1^6+x_2^6+x_1^4*x_3*x_4+x_3^3"]], ["HM_ex1.1", "x_1,x_2,x_3,x_4", ["x_1^7-2*x_1^5*x_4+x_1^3*x_4^2+x_1^2*x_2^2-x_1*x_2^2*x_3+x_2^3"]]]
for n, v, g in PROJ:
    I = helmer.Ideal(helmer.Ring(v), g)
    c = helmer.csm_class(I, projective=True)
    print(f"{n}  CSM {c}   Euler {c[0]}")
for n, v, g in AFF:
    S = helmer.whitney_stratify(helmer.Ideal(helmer.Ring(v), g))
    D = helmer.cc_multiplicities(S)
    print(f"{n}  multCC", [(m, d, str(P)) for d in sorted(S) for P, m in zip(S[d], D.get(d, []))])
