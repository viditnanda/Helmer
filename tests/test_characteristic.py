"""Euler obstructions, characteristic cycles, Euler characteristics and Chern classes.

Closed forms (cones over plane curves, smooth hypersurfaces, singular plane cubics, the
Cayley cubic), the global Euler obstruction of curves, and an independent consistency
check: chi(X-bar) = chi(X ∩ chart) + chi(X-bar ∩ hyperplane), computed by the projective
and the affine code paths.
"""
import pytest
import helmer
from helmer import Ideal, Ring
from test_helmer import FIX


def ideal(R, *g):
    return Ideal(R, list(g))


@pytest.fixture(autouse=True)
def seed():
    helmer.set_seed(1)


@pytest.mark.parametrize("d", [2, 3, 4])
def test_cone_over_smooth_plane_curve(d):
    R = Ring("x,y,z")
    Z = ideal(R, f"x^{d}+y^{d}+z^{d}")
    assert helmer.euler_obstruction(ideal(R, "x", "y", "z"), Z) == 2 * d - d * d
    S = helmer.whitney_stratify(Z)
    D = helmer.chern_classes(S)
    assert D["cc"] == {0: [(d - 1) ** 2], 2: [1]}  # the Milnor number of the vertex
    assert D["euler"] == {0: [1], 2: [0]}
    # sections: chi(Eu) = 2d - d^2; a general plane meets the cone in the plane curve
    # minus its d points at infinity (chi 3d - d^2 - d); a general line in d points.
    assert D["chern_mather"][2] == [[2 * d - d * d, 3 * d - d * d, d, 0]]


@pytest.mark.parametrize("gens,gamma", [
    (["x"], [0, 1]),
    (["x*y-1"], [2, 2]),            # C*: chi 0
    (["y-x^2"], [1, 2]),            # C: chi 1
    (["y^2-x^3-x^2"], [2, 3]),      # nodal cubic: global Eu 1 = chi(C - node) + 2
    (["x^2+y^2-1"], [2, 2]),
])
def test_global_polar_degrees(gens, gamma):
    R = Ring("x,y")
    assert helmer.global_polar_degrees(ideal(R, *gens)) == gamma
    assert helmer.global_euler_obstruction(ideal(R, *gens)) == gamma[1] - gamma[0]


@pytest.mark.parametrize("variables,gens,chi", [
    ("x,y", ["y^2-x^3-x^2"], 0),
    ("x,y", ["x*y-1"], 0),
    ("x,y", ["y^2-x^3"], 1),
    ("x,y", ["x*y"], 1),
    ("x,y,z", ["x^2+y^2-z^2"], 1),
    ("x,y,z", ["x^2-y^2*z"], 1),             # Whitney umbrella (contractible)
    ("x,y,z", ["x^2+y^2+z^2-1"], 2),         # T*S^2
    ("x,y,z", ["(x^2+y^2-z^2)*((x-1)^2+y^2-z^2)"], 2),
])
def test_affine_euler_characteristic(variables, gens, chi):
    assert helmer.euler_characteristic(ideal(Ring(variables), *gens)) == chi


def smooth_hypersurface_csm(n, d):
    """[a_0..a_n] for a smooth degree-d hypersurface in P^n: (1+h)^(n+1) d h / (1 + d h)."""
    from math import comb
    c = [0] * (n + 1)  # coefficients of h^j
    for j in range(1, n + 1):
        c[j] = sum(comb(n + 1, j - 1 - k) * d * (-d) ** k for k in range(j))
    return [c[n - k] for k in range(n + 1)]


@pytest.mark.parametrize("n,d", [(2, 1), (2, 2), (2, 3), (2, 4), (3, 2), (3, 3), (3, 4)])
def test_smooth_hypersurface_csm(n, d):
    R = Ring(",".join(f"x{i}" for i in range(n + 1)))
    F = ideal(R, "+".join(f"x{i}^{d}" for i in range(n + 1)))
    assert helmer.csm_class(F, projective=True) == smooth_hypersurface_csm(n, d)


@pytest.mark.parametrize("variables,gens,csm", [
    ("x,y,z", [], [3, 3, 1]),                                   # P^2
    ("x,y,z", ["y^2*z-x^3-x^2*z"], [1, 3, 0]),                  # nodal cubic
    ("x,y,z", ["y^2*z-x^3"], [2, 3, 0]),                        # cuspidal cubic
    ("x,y,z,w", ["x*y*z+x*y*w+x*z*w+y*z*w"], [5, 3, 3, 0]),     # Cayley cubic: 9 - 4 nodes
])
def test_singular_projective_csm(variables, gens, csm):
    assert helmer.csm_class(ideal(Ring(variables), *gens), projective=True) == csm


@pytest.mark.parametrize("name", sorted(n for n in FIX if FIX[n]["projective"]))
def test_projective_equals_chart_plus_hyperplane(name):
    ex = FIX[name]
    v = [s.strip() for s in ex["vars"].split(",")]
    sub = lambda g, val: g.replace(v[-1], f"({val})")
    X = ideal(Ring(ex["vars"]), *ex["gens"])
    chart = ideal(Ring(",".join(v[:-1])), *[sub(g, 1) for g in ex["gens"]])
    infinity = ideal(Ring(",".join(v[:-1])), *[sub(g, 0) for g in ex["gens"]])
    assert helmer.euler_characteristic(X, projective=True) == (
        helmer.euler_characteristic(chart) + helmer.euler_characteristic(infinity, projective=True))


def test_briancon_speder_small_origin_invisible_to_cc():
    # The origin is a Whitney stratum, but the family is mu-constant there (topologically
    # trivial along the t-axis), so the characteristic cycle of the constant function has
    # multiplicity 0 over it; at the 12 points where mu jumps it is 1.
    R = Ring("x,y,z,t")
    S = helmer.whitney_stratify(ideal(R, "x^3+y^4*z+z^9+t*x*y^3"), minimal=True)
    D = helmer.chern_classes(S)
    origin = [i for i, P in enumerate(S[0]) if P == ideal(R, "x", "y", "z", "t")][0]
    assert D["cc"][0][origin] == 0 and D["cc"][0][1 - origin] == 1
    assert D["euler"][0][origin] == 1 and D["euler"][0][1 - origin] == 12
    assert sum(sum(v) for v in D["euler"].values()) == 1  # contracts onto the t-axis


def test_per_stratum_classes_add_up():
    R = Ring("x,y,z,w")
    X = ideal(R, "x*y*z+x*y*w+x*z*w+y*z*w")
    S = helmer.whitney_stratify(X, projective=True)
    D = helmer.chern_classes(S, projective=True)
    total = [sum(c[k] for L in D["csm"].values() for c in L) for k in range(4)]
    assert total == helmer.csm_class(X, projective=True)
    assert D["euler"][0] == [1, 1, 1, 1]  # the four nodes
    assert D["chern_mather"][2][0][2] == 3  # c_M has the degree as its top part
