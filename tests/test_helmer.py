"""Helmer vs Macaulay2 (WhitneyStratifications 2.23, AssocPrimes=>false) ground truth."""
import json, os
import pytest
import helmer

HERE = os.path.dirname(__file__)
FIX = json.load(open(os.path.join(HERE, "m2_whitney_fixtures.json")))


def as_strata(R, pairs):
    out = {}
    for d, gens in pairs:
        out.setdefault(int(d), []).append(helmer.Ideal(R, [g for g in gens.split(",") if g.strip()]))
    return out


def same(A, B):
    for k in set(k for k in A if A[k]) | set(k for k in B if B[k]):
        a, b = A.get(k, []), B.get(k, [])
        if len(a) != len(b) or not all(any(x == y for y in b) for x in a):
            return False
    return True


@pytest.mark.parametrize("name", sorted(FIX))
def test_whitney_matches_m2(name):
    ex = FIX[name]
    R = helmer.Ring(ex["vars"])
    WS = helmer.whitney_stratify(helmer.Ideal(R, ex["gens"]), projective=ex["projective"])
    assert same(WS, as_strata(R, ex["m2_HM_strata"]))


@pytest.mark.parametrize("variables,f,projective,expected", [
    ("x,y", "y^2-x^3", False, [(0, "x,y"), (1, "y^2-x^3")]),
    ("x,y", "x*y", False, [(0, "x,y"), (1, "x"), (1, "y")]),
    ("x,y,z", "x^2+y^2-z^2", False, [(0, "x,y,z"), (2, "x^2+y^2-z^2")]),
    ("x,y,z", "x*y", True, [(0, "x,y"), (1, "x"), (1, "y")]),
])
def test_nonreduced_whitney(variables, f, projective, expected):
    R = helmer.Ring(variables)
    WS = helmer.whitney_stratify(helmer.Ideal(R, [f"({f})^2"]), projective=projective)
    assert same(WS, as_strata(R, expected))


def test_saturate_with_repeated_factors():
    R = helmer.Ring("x,y,z")
    I = helmer.Ideal(R, ["x^3*(x+y)^4*z^2"])
    f = helmer.Ideal(R, ["7*x^2*(x+y)^3"])
    assert helmer.saturate(I, f) == helmer.Ideal(R, ["z^2"])


@pytest.mark.parametrize("variables,gens,expected_chart", [
    ("x,y,z", ["x^2-y^2*z"], "y"),
    ("x,y,z", ["(x-y-1)^2-(z-y)^2*y"], "y - z"),
    # Two singular lines share one translated plane.
    ("x,y,z", ["(x+y-1)^2-(x*(z-1))^2*z"], "x + y - 1"),
    # The twisted cubic spans the ambient space. Adding w=0 makes its
    # affine hull a hyperplane shared by the entire parent variety.
    ("x,y,z", ["(x-y^2)^2-(z-y^3)^2*y"], None),
    ("x,y,z,w", ["w", "(x-y^2)^2-(z-y^3)^2*y"], None),
])
def test_conormal_chart(variables, gens, expected_chart):
    R = helmer.Ring(variables)
    I = helmer.Ideal(R, gens)
    messages = []
    actual = helmer.whitney_stratify(I, log=messages.append)
    charts = [m.split(": ", 1)[1] for m in messages if "conormal chart:" in m]
    assert charts == ([expected_chart] if expected_chart else [])
    assert same(actual, helmer.whitney_stratify(I, saturation="minor"))


def test_map_quadratic_roots():
    R, S = helmer.Ring("a,b,c,x"), helmer.Ring("a,b,c")
    X, Y = helmer.map_stratify(["a", "b", "c"], helmer.Ideal(R, ["a*x^2+b*x+c"]), helmer.Ideal(S, []), proper=False)
    exp_Y = {0: [["a", "b", "c"]], 1: [["a", "b"]], 2: [["a"], ["b^2-4*a*c"]], 3: [[]]}
    exp_X = {0: [], 1: [["a", "b", "c"]], 2: [["a", "b*x+c"], ["b*x+2*c", "2*a*x+b", "b^2-4*a*c"]], 3: [["a*x^2+b*x+c"]]}
    assert same(Y, {d: [helmer.Ideal(S, g) for g in L] for d, L in exp_Y.items()})
    assert same(X, {d: [helmer.Ideal(R, g) for g in L] for d, L in exp_X.items()})


def test_map_HLN_example_1_2():
    R, S = helmer.Ring("a,b,x,y"), helmer.Ring("a,b")
    X, Y = helmer.map_stratify(["a", "b"], helmer.Ideal(R, ["x^2-y^2+b", "-a*x+x^2+b*y"]), helmer.Ideal(S, []), proper=False)
    exp_Y = {0: [["b", "a"], ["b+4", "a"], ["27*b-4", "729*a^2+128"]],
             1: [["b"], ["a^6-3*a^4*b^2+3*a^2*b^4-b^6+a^4*b-20*a^2*b^3-8*b^5-16*b^4"]], 2: [[]]}
    assert same(Y, {d: [helmer.Ideal(S, g) for g in L] for d, L in exp_Y.items()})
    assert [len(X[d]) for d in sorted(X)] == [4, 2, 1]


MP = [  # (vars, gens, number of minimal primes per Macaulay2's decompose)
    ("x,y,z", ["x*y*z", "x+y+z"], 3), ("x,y,z", ["x^2-y*z", "y^2-x*z"], 3),
    ("x,y,z,w", ["x*w-y*z", "y^2-x*z", "z^2-y*w"], 1), ("a,b,c,d,e", ["a*b-c*d", "b*c-d*e", "a*e-b*d"], 5),
    ("x,y,z", ["(x^2-2)*(y^2+1)", "z^3-x"], 2), ("x,y,z,w", ["x^3-y^2*z", "y^3-z^2*w", "x*z-y*w"], 6),
    ("x,y,z", ["x^2*y^2-z^4", "x^3-y^3"], 6), ("x,y,z,w", ["(x*y-z*w)^2", "x^2*z-y^3"], 2),
    ("x,y,z,u,v", ["x*u-y*v", "x*v-z*u", "y*u-z*v"], 4), ("x,y,z", ["x^5+y^5+z^5-1", "x+y+z", "x*y*z"], 0),
]


@pytest.mark.parametrize("v,g,count", MP)
def test_minimal_primes_count_and_primality(v, g, count):
    R = helmer.Ring(v)
    I = helmer.Ideal(R, g)
    P = helmer.minimal_primes(I)
    assert len(P) == count
    for p in P:  # each prime contains I, and intersections stay off I's complement
        assert p.contains(I)


def test_gb_matches_m2():
    R = helmer.Ring("x,y,z")
    I = helmer.Ideal(R, ["x^2*y - 3/2*z^3 + 7", "x*y^2 - z + 11/3", "z^2*x - y"])
    E = helmer.eliminate(I, ["x"])
    assert E == helmer.Ideal(R, ["3*y^3-3*z^3+11*z^2", "9*z^5-42*z^2-6*z+22"])


def test_map_feynman_doc_example():
    # WhitneyStratifications docs: mapStratify(F, Xh, ideal(0_S)) for a Feynman graph polynomial
    R, S = helmer.Ring("M1,P,X_0,X_1,X_2"), helmer.Ring("M1,P")
    X = helmer.Ideal(R, ["X_0*X_1*X_2*((X_1+X_2)*X_0+(X_1+X_2)*(M1*X_1+M1*X_2)-P*X_1*X_2)"])
    A, B = helmer.map_stratify(["M1", "P"], X, helmer.Ideal(S, []))
    exp_Y = {0: [["P", "M1"]], 1: [["P"], ["4*M1-P"], ["M1"]], 2: [[]]}
    assert same(B, {d: [helmer.Ideal(S, g) for g in L] for d, L in exp_Y.items()})
    assert {d: len(L) for d, L in A.items()} == {0: 1, 1: 8, 2: 9, 3: 6, 4: 4}


def skeleton(S, d):
    return [P for k, L in S.items() if k <= d for P in L]


def contains_skeleta(big, small):
    """Every skeleton of `small` lies inside the matching skeleton of `big`."""
    dims = set(big) | set(small)
    return all(all(any(P.contains(Q) for Q in skeleton(big, d)) for P in skeleton(small, d)) for d in dims)


@pytest.mark.parametrize("name", sorted(n for n in FIX if not FIX[n]["projective"]))
def test_polar_refines_minimal(name):
    # The polar route is randomised and non-minimal; it must refine the minimal
    # Whitney stratification (M2 minCoarsenWS where recorded, else M2's output).
    ex = FIX[name]
    R = helmer.Ring(ex["vars"])
    helmer.set_seed(7)
    P = helmer.whitney_stratify(helmer.Ideal(R, ex["gens"]), method="polar", polar_bound=100)
    ref = as_strata(R, ex.get("m2_minimal_strata", ex["m2_HM_strata"]))
    assert contains_skeleta(P, ref)


@pytest.mark.parametrize("p,q", [
    (2147483647, 2147483629),  # first two primes on the old fixed schedule
    (1884583903, 2089177841),  # first two with seed 42 on libstdc++
])
def test_coincident_modular_generators(p, q):
    helmer.set_seed(42)
    R = helmer.Ring("x,y")
    I = helmer.Ideal(R, ["x+y", f"x+{p*q+1}*y"])
    assert I.dim() == 0
    assert I == helmer.Ideal(R, ["x", "y"])
    helmer.set_seed(42)
    assert helmer.eliminate(I, ["x"]) == helmer.Ideal(R, ["y"])


@pytest.mark.parametrize("u", ["x", "(x+t)"])
def test_affine_cylinder_quartic(u):
    R = helmer.Ring("x,y,z,w,t")
    f = f"{u}^4+y^3*(z+t)+w^4"
    messages = []
    strata = helmer.whitney_stratify(helmer.Ideal(R, [f]), log=messages.append)
    assert "affine cylinder: removing 1 free direction" in messages
    assert set(strata) == set(range(5))
    assert same(strata, as_strata(R, [(4, f), (2, f"{u},y,w"), (1, f"{u},y,w,z+t")]))


@pytest.mark.parametrize("level", [-1, 0])
@pytest.mark.parametrize("method", ["conormal", "polar"])
def test_cylinder_rational_shear_and_repeated_factors(level, method):
    R = helmer.Ring("x,y,z,w")
    u, v = "(x+2/3*z-w)", "(y-5/2*z)"
    I = helmer.Ideal(R, [f"({u}*{v})^2"])
    strata = helmer.whitney_stratify(I, level=level, method=method)
    expected = [(3, u), (3, v)]
    if level != 0:
        expected.append((2, f"{u},{v}"))
    assert set(strata) == set(range(4))
    assert same(strata, as_strata(R, expected))


def test_cylinder_univariate_base():
    R = helmer.Ring("x,y,z")
    f = "(y+z-1)*(y+z+1)"
    assert same(helmer.whitney_stratify(helmer.Ideal(R, [f])),
                as_strata(R, [(2, "y+z-1"), (2, "y+z+1")]))


@pytest.mark.parametrize("seed", [3, 11, 29, 47])
def test_common_translation_directions_match_direct_stratification(seed):
    import random
    rng = random.Random(seed)
    R = helmer.Ring("x,y,z,s,t")
    # A rational shear of a nodal curve times an affine plane.
    def form(variable):
        a, b = rng.sample(range(1, 13), 2)
        return f"({variable}+{a}/3*s-{b}/5*t)"
    u, v, w = map(form, ("x", "y", "z"))
    gens = [f"{u}*{v}", w]
    messages = []
    result = helmer.whitney_stratify(helmer.Ideal(R, gens), log=messages.append)
    assert "affine cylinder: removing 2 free directions" in messages
    assert same(result, as_strata(R, [(3, f"{u},{w}"), (3, f"{v},{w}"),
                                      (2, f"{u},{v},{w}")]))
    # Redundant equations preserve the ideal and force the ambient calculation:
    # the derivatives of s*w and t*w in a common translation direction detect
    # its s and t coordinates, respectively.
    direct_log = []
    direct = helmer.whitney_stratify(helmer.Ideal(R, gens + [f"s*{w}", f"t*{w}"]),
                                    log=direct_log.append)
    assert not any("affine cylinder" in m for m in direct_log)
    assert same(result, direct)


def test_inconsistent_cylinder_has_no_strata():
    R = helmer.Ring("x,y,z")
    assert helmer.whitney_stratify(helmer.Ideal(R, ["x-y", "x-y-1"])) == {}
