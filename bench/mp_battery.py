import sys, subprocess, time, re
import helmer
CASES = [
 ("x,y,z", ["x*y*z", "x+y+z"]),
 ("x,y,z", ["x^2-y*z", "y^2-x*z"]),
 ("x,y,z,w", ["x*w-y*z", "y^2-x*z", "z^2-y*w"]),               # twisted cubic cone
 ("x,y,z,w", ["x^2*y-z^3", "x*w^2-y^2"]),
 ("x,y,z", ["x^4+y^4-z^4", "x*y*z-1"]),
 ("a,b,c,d,e", ["a*b-c*d", "b*c-d*e", "a*e-b*d"]),
 ("x,y,z", ["(x^2-2)*(y^2+1)", "z^3-x"]),
 ("x,y,z,w", ["x^3-y^2*z", "y^3-z^2*w", "x*z-y*w"]),
 ("x,y,z", ["x^2*y^2-z^4", "x^3-y^3"]),
 ("x,y,z,w", ["(x*y-z*w)^2", "x^2*z-y^3"]),
 ("x,y,z,u,v", ["x*u-y*v", "x*v-z*u", "y*u-z*v"]),
 ("x,y,z", ["x^5+y^5+z^5-1", "x+y+z", "x*y*z"]),
]
m2 = ['R0=QQ[x];']
for i,(v,g) in enumerate(CASES):
    m2.append(f'R=QQ[{v}]; L=decompose ideal({",".join(g)}); print("CASE {i} "|toString(#L)); for P in L do print("P {i} "|toString(flatten entries gens gb P));')
open('/tmp/mp.m2','w').write("\n".join(m2))
t=time.time(); out = subprocess.run(['M2','--script','/tmp/mp.m2'],capture_output=True,text=True,timeout=1200).stdout; tm2=time.time()-t
bad=0; th=0
for i,(v,g) in enumerate(CASES):
    R=helmer.Ring(v)
    m2P=[helmer.Ideal(R,[x for x in m.group(1).split(',')]) for m in re.finditer(rf'^P {i} \{{(.*)\}}$', out, re.M)]
    t=time.time(); P=helmer.minimal_primes(helmer.Ideal(R,g)); th+=time.time()-t
    ok = len(P)==len(m2P) and all(any(p==q for q in m2P) for p in P)
    bad += not ok
    print(f"case {i:2d}: helmer {len(P)} primes, M2 {len(m2P)}  {'OK' if ok else 'MISMATCH'}")
print(f"{len(CASES)-bad}/{len(CASES)} agree; helmer total {th:.2f}s, M2 process {tm2:.1f}s")
