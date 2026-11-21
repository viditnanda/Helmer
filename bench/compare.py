import sys, json, time
sys.path.insert(0, '.')
import helmer
from examples import EXAMPLES

def strata_from_m2(R, lst):
    out = {}
    for d, gens in lst:
        out.setdefault(d, []).append(helmer.Ideal(R, [g.strip() for g in gens.split(',')] if gens.strip() else []))
    return out

def skeleton(S, d):
    return [P for k, L in S.items() if k <= d for P in L]

def var_subset(A, B):
    """V(A) ⊆ V(B) for lists of primes."""
    return all(any(P.contains(Q) for Q in B) for P in A)

def same(S1, S2):
    ks = set(k for k in S1 if S1[k]) | set(k for k in S2 if S2[k])
    for k in ks:
        a, b = S1.get(k, []), S2.get(k, [])
        if len(a) != len(b) or not all(any(x == y for y in b) for x in a): return False
    return True

def run(names, m2file):
    m2 = json.load(open(m2file)) if m2file else {}
    rows = []
    for name in names:
        v, g, p = EXAMPLES[name]
        R = helmer.Ring(v); I = helmer.Ideal(R, g)
        t = time.time(); WS = helmer.whitney_stratify(I, projective=p); t = time.time() - t
        counts = {d: len(L) for d, L in sorted(WS.items())}
        row = [name, f"{t:.3f}", str(counts)]
        for alg in ("HM", "HN"):
            r = m2.get(name, {}).get(alg)
            if not r or r['time'] is None:
                row += ["-", "-"]; continue
            S = strata_from_m2(R, r['strata'])
            dims = sorted(set(S) | set(WS))
            eq = same(WS, S)
            # skeleton comparison: is M2's skeleton contained in ours / ours in M2's?
            sub = all(var_subset(skeleton(S, d), skeleton(WS, d)) for d in dims)
            sup = all(var_subset(skeleton(WS, d), skeleton(S, d)) for d in dims)
            rel = "equal" if eq else ("M2⊂ours" if sub and not sup else "ours⊂M2" if sup and not sub else "same skel" if sub and sup else "differ")
            row += [f"{r['time']:.3f}", rel]
        rows.append(row)
        print(" | ".join(row), flush=True)
    return rows

if __name__ == "__main__":
    run(sys.argv[2:], sys.argv[1])
