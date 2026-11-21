import sys, json, time, subprocess, tempfile, os, re
sys.path.insert(0, '.')
import helmer
from examples import EXAMPLES
from compare import strata_from_m2, skeleton, var_subset
fx = json.load(open('../tests/m2_whitney_fixtures.json'))
PKG = os.environ.get("WHITSTRAT_M2", "WhitneyStratifications.m2")
def m2_pol(v, g):
    src = f'needsPackage("WhitneyStratifications", FileName=>"{PKG}");\nR=QQ[{v}]; I=ideal({",".join(g)});\nt0=cpuTime(); WS=whitneyStratifyPol(I); print("TIME "|toString(cpuTime()-t0));\n'
    fd, fn = tempfile.mkstemp(suffix='.m2'); os.write(fd, src.encode()); os.close(fd)
    try:
        out = subprocess.run(['M2', '--script', fn], capture_output=True, text=True, timeout=int(os.environ.get("M2_TIMEOUT", "900"))).stdout
        m = re.search(r'TIME (\S+)', out); return float(m.group(1)) if m else None
    except subprocess.TimeoutExpired:
        return None
names = sys.argv[1:] or sorted(n for n in fx if not fx[n]["projective"])
print("example | M2 conormal | M2 polar | helmer conormal | helmer polar | polar valid | extra strata")
for name in names:
    ex = fx[name]; R = helmer.Ring(ex["vars"]); I = helmer.Ideal(R, ex["gens"])
    t = time.time(); C = helmer.whitney_stratify(I); tc = time.time() - t
    t = time.time(); P = helmer.whitney_stratify(I, method="polar"); tp = time.time() - t
    M = strata_from_m2(R, ex["m2_HM_strata"])
    dims = sorted(set(C) | set(P))
    valid = all(var_subset(skeleton(C, d), skeleton(P, d)) for d in dims)
    extra = sum(len(L) for L in P.values()) - sum(len(L) for L in C.values())
    tm = m2_pol(ex["vars"], ex["gens"])
    print(f"{name} | {ex['m2_seconds']:.2f} | {tm if tm is None else round(tm,2)} | {tc:.2f} | {tp:.2f} | {valid} | {extra:+d}", flush=True)
