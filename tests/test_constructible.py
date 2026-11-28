"""Constructible functions: arithmetic, Euler integrals, pushforward, pullback, Euler
obstructions and characteristic cycles. Every pushforward here has a closed form."""
import pytest
import helmer
from helmer import Ideal, Ring, indicator, ConstructibleFunction


@pytest.fixture(autouse=True)
def seed():
    helmer.set_seed(1)


def cf(R, terms):
    return ConstructibleFunction(R, {Ideal(R, g): c for g, c in terms})


def test_indicator_inclusion_exclusion():
    R = Ring("x,y")
    assert indicator(Ideal(R, ["x*y"])) == cf(R, [(["x"], 1), (["y"], 1), (["x", "y"], -1)])
    assert indicator(Ideal(R, ["x*y"])).integral() == 1
    assert (indicator(Ideal(R, ["x"])) * indicator(Ideal(R, ["y"]))) == cf(R, [(["x", "y"], 1)])


def test_push_finite_maps():
    S = Ring("s")
    assert indicator(Ideal(Ring("x"), [])).push(["x^2"], S) == cf(S, [([], 2), (["s"], -1)])
    R = Ring("x,y")
    circle = indicator(Ideal(R, ["x^2+y^2-1"]))
    assert circle.push(["x"], S) == cf(S, [([], 2), (["s-1"], -1), (["s+1"], -1)])


def test_push_non_proper_finite_jelonek_set():
    # x: (xy = 1) -> C is proper away from 0, where the fibre escapes to infinity
    R, S = Ring("x,y"), Ring("s")
    assert indicator(Ideal(R, ["x*y-1"])).push(["x"], S) == cf(S, [([], 1), (["s"], -1)])


def test_push_broughton_atypical_value_at_infinity():
    # x(xy - 1): C^2 -> C has no critical points, yet the fibre over 0 has chi 1 and all
    # others chi 0; non-proper over every value, so this needs the charts at infinity.
    R, S = Ring("x,y"), Ring("s")
    assert indicator(Ideal(R, [])).push(["x*(x*y-1)"], S) == cf(S, [(["s"], 1)])


def test_push_family_of_plane_cubics():
    # general fibre: affine elliptic curve (chi -1); over 0: the cusp (chi 1)
    R, S = Ring("x,y,t"), Ring("s")
    assert indicator(Ideal(R, ["y^2-x^3-t"])).push(["t"], S) == cf(S, [([], -1), (["s"], 2)])


def test_push_roots_of_a_quadratic():
    # number of roots of a x^2 + b x + c (and chi of the solution set, a line, at 0)
    R, T = Ring("a,b,c,x"), Ring("a,b,c")
    X = indicator(Ideal(R, ["a*x^2+b*x+c"]))
    r = X.push(["a", "b", "c"], T)
    assert r == cf(T, [([], 2), (["a"], -1), (["b^2-4*a*c"], -1), (["a", "b", "c"], 1)])
    assert r.integral() == X.integral() == 1
    assert r.value(Ideal(T, ["a", "b", "c-1"])) == 0  # 0 = c != 0: no roots


def test_push_two_link_arm():
    # inverse kinematics of a planar arm with unit links: 2 configurations generically,
    # 1 on the reach circle, a circle (chi 0) over the origin, none over the rest of the
    # isotropic cone x^2 + y^2 = 0
    C, W = Ring("a,b,c,d"), Ring("x,y")
    conf = indicator(Ideal(C, ["a^2+b^2-1", "c^2+d^2-1"]))
    ik = conf.push(["a+c", "b+d"], W)
    assert ik == cf(W, [([], 2), (["x^2+y^2-4"], -1), (["x^2+y^2"], -2)])
    assert ik.integral() == conf.integral() == 0


def test_pullback():
    R, U = Ring("x,y"), Ring("u")
    assert indicator(Ideal(U, ["u"])).pull(["x*y"], Ideal(R, [])) == indicator(Ideal(R, ["x*y"]))


def test_euler_obstruction_and_characteristic_cycle():
    R = Ring("x,y,z")
    cone = Ideal(R, ["x^2+y^2-z^2"])
    assert indicator(cone).eu_decomposition() == {cone: 1, Ideal(R, ["x", "y", "z"]): 1}
    assert indicator(cone).characteristic_cycle() == {cone: 1, Ideal(R, ["x", "y", "z"]): 1}
    eu = helmer.euler_obstruction_function(Ideal(R, ["x^3+y^3+z^3"]))
    assert eu.value(Ideal(R, ["x", "y", "z"])) == -3
    assert eu.value(Ideal(R, ["x^3+y^3+z^3"])) == 1


def test_from_strata_round_trip():
    R = Ring("x,y,z")
    S = helmer.whitney_stratify(Ideal(R, ["x^2-y^2*z"]))
    values = {d: [d + 1 for _ in L] for d, L in S.items()}
    assert helmer.from_strata(S, values).on(S) == values
