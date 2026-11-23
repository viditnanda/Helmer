import pytest
import helmer


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
