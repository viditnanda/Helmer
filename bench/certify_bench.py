"""Time eliminate_certified on critical-value ideals. usage: python bench/certify_bench.py [case ...] (needs sympy)"""
import sys, time, json

import helmer, sympy as sp

def crit(f, names):
    syms = sp.symbols(names)
    fe = sp.expand(f(*syms))
    R = helmer.Ring(",".join(names + ["epsilon"]))
    gens = [str(sp.expand(sp.diff(fe, s))).replace("**", "^") for s in syms] + [str(fe).replace("**", "^") + "-epsilon"]
    return R, helmer.Ideal(R, gens), names

def prod(it):
    r = 1
    for t in it: r = r * t
    return r

CASES = {
  "torus": (lambda x,y,z: (x**2+y**2+z**2+3)**2-16*(x**2+y**2), ["x","y","z"]),
  "two spheres": (lambda x,y,z: (x**2+y**2+z**2-1)*((x-1)**2+y**2+z**2-1), ["x","y","z"]),
  "clifford F": (lambda x,y,z,w: (x**2+y**2-sp.Rational(1,2))**2+(z**2+w**2-sp.Rational(1,2))**2, ["x","y","z","w"]),
  "concentric 5": (lambda x,y: prod(x**2+y**2-(i+1)**2 for i in range(5)), ["x","y"]),
  "separated 5": (lambda x,y: prod((x-3*i)**2+y**2-1 for i in range(5)), ["x","y"]),
  "separated 6": (lambda x,y: prod((x-3*i)**2+y**2-1 for i in range(6)), ["x","y"]),
  "concentric 8": (lambda x,y: prod(x**2+y**2-(i+1)**2 for i in range(8)), ["x","y"]),
}
names = sys.argv[1:] or list(CASES)
out = {}
for k in names:
    R, I, vs = crit(*CASES[k])
    t0 = time.time(); J = helmer.eliminate(I, vs); t1 = time.time()
    J2, certs = helmer.eliminate_certified(I, vs); t2 = time.time()
    ok = all(c is not None and helmer.verify_membership(I, c["target"], c["cofactors"]) for c in certs) and len(certs) > 0
    bounds = [c["bounds"] for c in certs if c]
    out[k] = dict(cert_s=round((t2-t1)-(t1-t0), 3), ok=ok, bounds=bounds)
    print(k, out[k], flush=True)
