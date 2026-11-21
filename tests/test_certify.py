"""Certified membership: certificates verify exactly, corrupted ones are rejected,
and non-members never receive a certificate."""
import pytest
import helmer


def critical_ideal(f, grads, names):
    R = helmer.Ring(",".join(names + ["epsilon"]))
    return R, helmer.Ideal(R, grads + [f"{f}-epsilon"])


CASES = {
    "circle": ("x^2+y^2-1", ["2*x", "2*y"], ["x", "y"]),
    "two circles": ("(x^2+y^2-1)*((x-3)^2+y^2-1)",
                    ["2*x*((x-3)^2+y^2-1)+2*(x-3)*(x^2+y^2-1)", "2*y*((x-3)^2+y^2-1)+2*y*(x^2+y^2-1)"],
                    ["x", "y"]),
    # critical points of f form whole circles here
    "concentric": ("(x^2+y^2-1)*(x^2+y^2-4)*(x^2+y^2-9)",
                   ["2*x*((x^2+y^2-4)*(x^2+y^2-9)+(x^2+y^2-1)*(x^2+y^2-9)+(x^2+y^2-1)*(x^2+y^2-4))",
                    "2*y*((x^2+y^2-4)*(x^2+y^2-9)+(x^2+y^2-1)*(x^2+y^2-9)+(x^2+y^2-1)*(x^2+y^2-4))"],
                   ["x", "y"]),
    "torus": ("(x^2+y^2+z^2+3)^2-16*(x^2+y^2)",
              ["4*x*(x^2+y^2+z^2+3)-32*x", "4*y*(x^2+y^2+z^2+3)-32*y", "4*z*(x^2+y^2+z^2+3)"],
              ["x", "y", "z"]),
}


@pytest.mark.parametrize("name", sorted(CASES))
def test_certified_elimination(name):
    f, grads, names = CASES[name]
    R, I = critical_ideal(f, grads, names)
    J, certs = helmer.eliminate_certified(I, names)
    assert J == helmer.eliminate(I, names)
    assert len(certs) == len(J.gens) and all(c is not None for c in certs)
    for g, c in zip(J.gens, certs):
        assert c["target"] == helmer.Poly(R, g)
        assert helmer.verify_membership(c["generators"], c["target"], c["cofactors"], R)


def test_corrupted_certificates_are_rejected():
    f, grads, names = CASES["two circles"]
    R, I = critical_ideal(f, grads, names)
    _, (c,) = helmer.eliminate_certified(I, names)
    gens, target, cof = c["generators"], c["target"], list(c["cofactors"])
    bump = helmer.Poly(R, "1/7*x*y")
    for j in range(len(cof)):
        bad = list(cof)
        bad[j] = bad[j] + bump
        assert not helmer.verify_membership(gens, target, bad, R)
    assert not helmer.verify_membership(gens, target + helmer.Poly(R, "1"), cof, R)
    assert not helmer.verify_membership(gens, target, cof[:-1], R)


def test_members_and_non_members():
    R = helmer.Ring("x,y")
    I = helmer.Ideal(R, ["x^2-y", "y^3-1"])
    c = helmer.certify_membership(I, "x^6-1", ["x"])
    assert c is not None and helmer.verify_membership(c["generators"], c["target"], c["cofactors"], R)
    assert helmer.certify_membership(I, "x-1", ["x"]) is None
    assert helmer.certify_membership(helmer.Ideal(R, ["x^2"]), "x", ["x"]) is None
    z = helmer.certify_membership(I, "0", ["x"])
    assert z is not None and all(str(q) == "0" for q in z["cofactors"])
