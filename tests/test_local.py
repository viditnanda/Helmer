"""Multiplicities, containment, polar sequences and minimal Whitney stratifications.

Closed-form values, Macaulay2's minCoarsenWS where recorded, and two consistency
checks: the randomised polar route (which adds random non-minimal strata) must coarsen
to the same stratification as the conormal route, and a projective variety must coarsen
consistently with its affine cone.
"""
import pytest
import helmer
from helmer import Ideal, Ring
from test_helmer import FIX, as_strata, same


def ideal(R, *g):
    return Ideal(R, list(g))


@pytest.fixture(autouse=True)
def seed():
    helmer.set_seed(1)


@pytest.mark.parametrize("variables,X,Y,e", [
    ("x,y", ["x", "y"], ["y^2-x^3"], 2),                       # cusp
    ("x,y", ["x", "y"], ["x*y"], 2),                           # node
    ("x,y", ["x-1", "y"], ["y^2-x^3+1"], 1),                   # smooth point
    ("x,y", ["x", "y"], ["x^2", "x*y"], 1),                    # embedded point: line x = 0 counts
    ("x,y", ["x"], ["x^2"], 2),                                # double line
    ("x,y,z", ["x", "y", "z"], ["x^2+y^2-z^2"], 2),            # quadric cone
    ("x,y,z", ["x", "y", "z"], ["x^3+y^3+z^3"], 3),            # cubic cone
    ("x,y,z", ["x", "y"], ["x^2-y^3"], 2),                     # along a curve
    ("x,y,z", ["x", "y"], ["x^2-y^2*z"], 2),                   # umbrella along its handle
    ("x,y,z", ["x-1", "y", "z"], ["x^2+y^2-z^2"], 0),          # not contained
    # two planes meeting in a point (not Cohen-Macaulay: a general linear section through
    # the point has length 3, the multiplicity is 2)
    ("a,b,c,d", ["a", "b", "c", "d"], ["a*c", "a*d", "b*c", "b*d"], 2),
    ("a,b,c,d", ["a", "b", "c", "d"], ["a*c-b^2", "b*d-c^2", "a*d-b*c"], 3),  # twisted cubic cone
    # nodes at the conjugate points (±sqrt 2, 0): no F_p-point unless 2 is a square mod p
    ("x,y", ["x^2-2", "y"], ["y^2-(x^2-2)^2"], 2),
])
def test_multiplicity(variables, X, Y, e):
    R = Ring(variables)
    assert helmer.multiplicity(ideal(R, *X), ideal(R, *Y)) == e


def test_multiplicity_projective_cone():
    R = Ring("x,y,z,w")
    # a point of the projective curve y^2 z = x^3 (the cusp) and a smooth point
    C = ideal(R, "y^2*z-x^3", "w")
    assert helmer.multiplicity(ideal(R, "x", "y", "w"), C) == 2
    assert helmer.multiplicity(ideal(R, "x-z", "y-z", "w"), C) == 1


def test_contained_and_singular():
    R = Ring("x,y")
    assert helmer.contained(ideal(R, "x", "y"), ideal(R, "x*y"))
    assert helmer.contained(ideal(R, "x^3"), ideal(R, "x*y"))  # radical containment
    assert not helmer.contained(ideal(R, "x-1"), ideal(R, "x*y"))
    S = Ring("x,y,z")
    cone = ideal(S, "x^2+y^2-z^2")
    assert helmer.in_singular_locus(ideal(S, "x", "y", "z"), cone)
    assert not helmer.in_singular_locus(ideal(S, "x-3", "y-4", "z-5"), cone)


@pytest.mark.parametrize("d", [2, 3, 4])
def test_polar_sequence_cone_over_smooth_curve(d):
    # Cone over a smooth plane curve of degree d: the polar curve is d(d-1) lines (the
    # class), the surface has multiplicity d; Euler obstruction 2d - d^2 at the vertex.
    R = Ring("x,y,z")
    Z = ideal(R, f"x^{d}+y^{d}+z^{d}")
    s = helmer.polar_sequence(ideal(R, "x", "y", "z"), Z)
    assert s == [0, d * (d - 1), d]
    assert s[2] - s[1] + s[0] == 2 * d - d * d


def test_polar_sequence_umbrella():
    R = Ring("x,y,z")
    U = ideal(R, "x^2-y^2*z")
    handle = helmer.polar_sequence(ideal(R, "x", "y"), U)
    origin = helmer.polar_sequence(ideal(R, "x", "y", "z"), U)
    assert handle == [0, 0, 2] and origin == [0, 1, 2]  # the polar curve passes through 0
    assert helmer.polar_multiplicity(ideal(R, "x", "y", "z"), U, 1) == 1


@pytest.mark.parametrize("name", sorted(n for n in FIX if "m2_minimal_strata" in FIX[n]))
def test_min_coarsen_matches_m2(name):
    ex = FIX[name]
    R = Ring(ex["vars"])
    M = helmer.whitney_stratify(Ideal(R, ex["gens"]), projective=ex["projective"], minimal=True)
    assert same(M, as_strata(R, ex["m2_minimal_strata"]))


@pytest.mark.parametrize("name", sorted(n for n in FIX if not FIX[n]["projective"]))
def test_polar_route_coarsens_to_the_same(name):
    ex = FIX[name]
    R = Ring(ex["vars"])
    I = Ideal(R, ex["gens"])
    A = helmer.min_coarsen(helmer.whitney_stratify(I))
    helmer.set_seed(7)
    B = helmer.min_coarsen(helmer.whitney_stratify(I, method="polar", polar_bound=100))
    assert same(A, B)


@pytest.mark.parametrize("name", sorted(n for n in FIX if FIX[n]["projective"]))
def test_projective_consistent_with_cone(name):
    ex = FIX[name]
    R = Ring(ex["vars"])
    I = Ideal(R, ex["gens"])
    P = helmer.min_coarsen(helmer.whitney_stratify(I, projective=True), projective=True)
    A = helmer.min_coarsen(helmer.whitney_stratify(I))
    cone = {d + 1: L for d, L in P.items()}
    cone.setdefault(0, []).append(Ideal(R, [v.strip() for v in ex["vars"].split(",")]))
    assert same(A, cone)


def test_point_on_no_curve_stratum_is_kept():
    # Cone ∪ plane z = 5: the vertex lies on no curve stratum while curve strata exist.
    # (Macaulay2's minCoarsenWS 2.23 drops it: with no containing stratum at the next
    # level, the point is never added.)
    R = Ring("x,y,z")
    M = helmer.whitney_stratify(ideal(R, "(x^2+y^2-z^2)*(z-5)"), minimal=True)
    assert any(P == ideal(R, "x", "y", "z") for P in M[0])


def test_polar_test_uses_the_containing_stratum():
    # Umbrella ∪ plane x = 5: the origin must stay (polar curve of the umbrella). In
    # minCoarsenWS 2.23 the polar test indexes the top level by position (Zlist_b, not
    # Zlist_(contY_b)), so if the plane is listed first it compares the plane's (zero)
    # polar sequences and merges the origin into the handle.
    R = Ring("x,y,z")
    M = helmer.whitney_stratify(ideal(R, "(x^2-y^2*z)*(x-5)"), minimal=True)
    assert any(P == ideal(R, "x", "y", "z") for P in M[0])
    assert any(P == ideal(R, "x", "y") for P in M[1])
