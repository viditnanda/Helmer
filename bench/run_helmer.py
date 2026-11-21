import sys, time, json
sys.path.insert(0,'.')
import helmer
from examples import EXAMPLES
name = sys.argv[1]
v,g,p = EXAMPLES[name]
R = helmer.Ring(v); I = helmer.Ideal(R, g)
t0=time.time()
WS = helmer.whitney_stratify(I, projective=p, log=lambda s: print(f"  [{time.time()-t0:7.2f}s] {s}", flush=True))
print(name, f"{time.time()-t0:.2f}s", {d: len(L) for d,L in sorted(WS.items())})
for d in sorted(WS):
    for P in WS[d]: print("  ", d, str(P)[:200])
print(helmer.stats())
