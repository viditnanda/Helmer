import pytest
import helmer
from test_helmer import FIX


def test_pinned_polygon_filtration():
    R = helmer.Ring("x0,y0,x1,y1,x2,y2")
    equations = [f"x{i}^2+y{i}^2-1" for i in range(3)]
    equations.append("(x0+x1+x2-2)^2+(y0+y1+y2)^2-1")
    I = helmer.Ideal(R, equations)
    log = []
    before = helmer.stats()["prime_calls"]
    F = helmer.whitney_filtration(I, log=log.append)
    assert helmer.stats()["prime_calls"] == before
    assert "group top" in log and "group points" in log
    # Rank drops exactly when all three moving bars and the closing bar
    # are parallel. Their sum fixes the horizontal direction.
    points = helmer.Ideal(R, ["1"])
    for signs in [(1, 1, 1), (-1, 1, 1), (1, -1, 1), (1, 1, -1)]:
        point = helmer.Ideal(R, [f"x{i}-({s})" for i, s in enumerate(signs)]
                              + [f"y{i}" for i in range(3)])
        points = helmer.intersect(points, point)
    assert F == {0: points, 1: points, 2: I}


@pytest.mark.parametrize("gens", [
    ["(x*y)^2"],             # a nonreduced hypersurface
    ["x*z", "y*z"],          # components of different dimensions
    ["x^2", "y^2"],          # generically nonreduced complete intersection
    [], ["1"],               # ambient space and empty variety
])
def test_filtration_fallback_and_cylinders(gens):
    R = helmer.Ring("x,y,z")
    I = helmer.Ideal(R, gens)
    eager = helmer.whitney_filtration(I, defer=False)
    assert helmer.whitney_filtration(I) == eager


def test_two_cones_filtration():
    R = helmer.Ring("x,y,z")
    I = helmer.Ideal(R, ["(x^2+y^2-z^2)*((x-1)^2+y^2-z^2)"])
    eager = helmer.whitney_filtration(I, defer=False)
    grouped = helmer.whitney_filtration(I)
    assert grouped[2] == eager[2] and grouped[1] == eager[1]
    assert eager[0].contains(grouped[0])
    assert grouped[0].dim() == 0


def test_conormal_of_union_with_no_dense_jacobian_minor():
    R = helmer.Ring("x,y,z")
    I = helmer.Ideal(R, ["x*y*z"])
    F = helmer.whitney_filtration(I)
    assert F == {0: helmer.Ideal(R, ["x", "y", "z"]),
                 1: helmer.Ideal(R, ["x*y", "x*z", "y*z"]), 2: I}


def test_rotating_polygon_has_no_exceptional_points():
    R = helmer.Ring("x0,y0,x1,y1,x2,y2")
    I = helmer.Ideal(R, ["x0^2+y0^2-1", "x1^2+y1^2-4", "x2^2+y2^2-9",
                         "(x0+x1+x2)^2+(y0+y1+y2)^2-4"])
    # The singular configurations are all rotations of (1,-2,3) along one axis.
    circle = helmer.Ideal(R, ["3*x0-x2", "3*y0-y2", "3*x1+2*x2", "3*y1+2*y2",
                              "x2^2+y2^2-9"])
    log = []
    F = helmer.whitney_filtration(I, log=log.append)
    assert F == {0: helmer.Ideal(R, ["1"]), 1: circle, 2: I}
    assert "  curve orbit bound" in log and "conormal complete" not in log


def test_quartic_orbit_bound_after_translation_and_shear():
    R = helmer.Ring("x,y,z,w,t")
    a = "(x+t-2)"
    I = helmer.Ideal(R, [f"{a}^4+y^3*z+z^4+(t-2)*{a}*y^2+w^4"])
    axis = helmer.Ideal(R, [a, "y", "z", "w"])
    vertex = helmer.Ideal(R, ["x", "y", "z", "w", "t-2"])
    log = []
    F = helmer.whitney_filtration(I, log=log.append)
    assert F == {0: vertex, 1: axis, 2: axis, 3: axis, 4: I}
    assert "  curve orbit bound" in log and "conormal complete" not in log


def test_symmetry_fixing_the_curve_requires_conormal():
    R = helmer.Ring("x,y,t")
    I = helmer.Ideal(R, ["x^4+y^4+(t^2+t+1)*x^2*y^2"])
    log = []
    F = helmer.whitney_filtration(I, log=log.append)
    # Scaling x and y fixes the t-axis. Two branches coalesce for c(t)=±2.
    assert F[0] == helmer.Ideal(R, ["x", "y", "(t^2+t-1)*(t^2+t+3)"])
    assert "conormal complete" in log and "  curve orbit bound" not in log


@pytest.mark.parametrize("name", sorted(FIX))
def test_filtration_contains_the_minimal_stratification(name):
    # A filtration may be finer than the minimal Whitney stratification, never coarser.
    ex = FIX[name]
    R = helmer.Ring(ex["vars"])
    I = helmer.Ideal(R, ex["gens"])
    F = helmer.whitney_filtration(I, projective=ex["projective"])
    M = helmer.whitney_stratify(I, projective=ex["projective"], minimal=True)
    for d, strata in M.items():
        for S in strata:
            assert d in F and helmer.contained(S, F[d])


def test_shifted_quartic_filtration_is_minimal():
    # Eager stratification takes minutes; the polar route then min_coarsen gives the same.
    R = helmer.Ring("x,y,z,w,t")
    I = helmer.Ideal(R, ["(x+y)^4+y^3*z+z^4+(t-2)*(x+y)*y^2+w^4"])
    F = helmer.whitney_filtration(I)
    assert F[0] == helmer.Ideal(R, ["x", "y", "z", "w", "t-2"])
    assert F[1] == F[3] == helmer.Ideal(R, ["x", "y", "z", "w"])
