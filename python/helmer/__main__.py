"""Command line: Whitney-stratify the variety of some polynomial equations.

    helmer "y^2*z - x^2"                              # the Whitney umbrella
    helmer "x^2 + y^2 + z^2 - 1" "x*y" --vars x,y,z   # several equations, variables in a chosen order
    helmer --projective "x*y*z + x*y*w + x*z*w + y*z*w" --vars x,y,z,w
    helmer --json "y^2 - x^3"

(or `python -m helmer ...`). Coefficients are rational: write 3/2, not 1.5.
"""
import argparse
import json
import re
import sys
import time

from . import Ideal, Ring, __version__, set_seed, whitney_stratify

_NAME = re.compile(r"[A-Za-z_][A-Za-z0-9_]*")


def infer_variables(equations):
    """The variable names in the equations, sorted naturally (x, y, z; x_2 before x_10)."""
    names = {v for e in equations for v in _NAME.findall(e)}
    key = lambda v: [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", v)]
    return sorted(names, key=key)


def main(argv=None):
    p = argparse.ArgumentParser(
        prog="helmer",
        description="Whitney stratification of the variety defined by polynomial equations over Q. "
                    "Prints, for each dimension, the prime ideals whose zero sets are the closures of the strata.",
        epilog="The open strata of dimension d are the connected components of V(P) minus the strata of lower "
               "dimension, for P listed under d. Over R the same polynomials give a valid stratification.")
    p.add_argument("equations", nargs="+", help="polynomials, e.g. 'y^2*z - x^2' (several: their common zeros)")
    p.add_argument("--vars", help="comma-separated variables, in order (default: the names in the equations, sorted)")
    p.add_argument("--projective", action="store_true", help="homogeneous equations in projective space")
    p.add_argument("--method", choices=["conormal", "polar"], default="conormal",
                   help="conormal (default; deterministic up to the modular lift) or polar (randomised, non-minimal)")
    p.add_argument("--seed", type=int, help="seed the random choices, for reproducible runs")
    p.add_argument("--json", action="store_true", help="print {dimension: [[generators], ...]} as JSON")
    p.add_argument("--verbose", "-v", action="store_true", help="progress messages on stderr")
    p.add_argument("--version", action="version", version="helmer " + __version__)
    a = p.parse_args(argv)

    names = [v.strip() for v in a.vars.split(",")] if a.vars else infer_variables(a.equations)
    if not names:
        p.error("no variables: give --vars, or equations that contain some")
    if a.seed is not None:
        set_seed(a.seed)
    R = Ring(",".join(names))
    try:
        I = Ideal(R, a.equations)
    except Exception as err:  # a parse error names the offending polynomial
        p.error(str(err))
    t0 = time.perf_counter()
    log = (lambda s: print(f"[{time.perf_counter() - t0:8.2f} s] {s}", file=sys.stderr)) if a.verbose else None
    S = whitney_stratify(I, projective=a.projective, method=a.method, log=log)
    dt = time.perf_counter() - t0

    if a.json:
        print(json.dumps({str(d): [P.gens for P in S[d]] for d in sorted(S, reverse=True) if S[d]}, indent=1))
        return 0
    space = f"P^{len(names) - 1}" if a.projective else f"C^{len(names)}"
    print(f"Whitney stratification in {space} with coordinates {', '.join(names)} ({dt:.2f} s)")
    if not any(S.values()):
        print("  the variety is empty")
    for d in sorted(S, reverse=True):
        for P in S[d]:
            print(f"  dim {d}:  " + ", ".join(P.gens))
    return 0


if __name__ == "__main__":
    sys.exit(main())
