import sys, subprocess, json, time, re
sys.path.insert(0, '.')
from examples import EXAMPLES
PKG = __import__("os").environ.get("WHITSTRAT_M2", "WhitneyStratifications.m2")  # path to the M2 package file
def m2_script(vars_, gens, proj, assoc):
    return f'''needsPackage("WhitneyStratifications", FileName=>"{PKG}");
R=QQ[{vars_}]; I=ideal({",".join(gens)});
t0=cpuTime(); WS=whitneyStratify(I, Projective=>{str(proj).lower()}, AssocPrimes=>{str(assoc).lower()}); t1=cpuTime();
print("TIME "|toString(t1-t0));
for k in sort keys WS do for P in toList(WS#k) do print("STRATUM "|toString(k)|" "|toString(flatten entries gens P));
'''
out = {}
for name in sys.argv[2:] or EXAMPLES:
    v, g, p = EXAMPLES[name]
    res = {}
    for assoc in [a == "HN" for a in __import__("os").environ.get("ALGS","HM,HN").split(",")]:
        import tempfile, os
        fd, fn = tempfile.mkstemp(suffix='.m2'); os.write(fd, m2_script(v, g, p, assoc).encode()); os.close(fd)
        try:
            r = subprocess.run(['M2','--script',fn], capture_output=True, text=True, timeout=int(sys.argv[1]))
            txt = r.stdout + r.stderr
            tm = re.search(r'TIME (\S+)', txt)
            strata = [(int(m.group(1)), m.group(2)) for m in re.finditer(r'STRATUM (\d+) \{(.*)\}', txt)]
            res['HN' if assoc else 'HM'] = {'time': float(tm.group(1)) if tm else None, 'strata': strata, 'err': None if tm else txt[-400:]}
        except subprocess.TimeoutExpired:
            res['HN' if assoc else 'HM'] = {'time': None, 'strata': [], 'err': 'timeout'}
        print(name, 'HN' if assoc else 'HM', res['HN' if assoc else 'HM']['time'], res['HN' if assoc else 'HM']['err'] and res['HN' if assoc else 'HM']['err'][:200], flush=True)
    out[name] = res
    json.dump(out, open('m2_results.json' if len(sys.argv) <= 2 else f'm2_{"_".join(sys.argv[2:])}.json', 'w'), indent=1)
