"""Constructible functions on complex affine varieties, and Euler calculus.

A constructible function is stored as a finite Z-linear combination of indicator
functions 1_V of irreducible closed subvarieties V (prime ideals). These indicators are a
free basis of the group of constructible functions, so the representation is canonical
and functions built from different stratifications add without a common refinement.

Pushforward along a polynomial map f uses map_stratify(..., proper=False), whose Jelonek
flag makes f restricted to each preimage proper over the target strata, so by Thom's
isotopy lemma chi(V ∩ f^{-1}(y)) is constant on each target stratum T. The value on T is
read off the preimage of deg T general points of T (a variety over QQ):
chi(V ∩ f^{-1}(T ∩ L)) / deg T. Euler characteristics come from chern_classes, so
everything here is Monte Carlo in the same sense; pushforward values are moreover computed
for two independent random slices and must agree.
"""
from itertools import combinations
import random as _random

from . import _helmer as _h

__all__ = ["ConstructibleFunction", "indicator", "from_strata", "euler_obstruction_function"]


class ConstructibleFunction:
    """sum_V c_V 1_V over prime ideals V of one polynomial ring (affine space)."""

    def __init__(self, ring, terms=None):
        self.ring = ring
        self.terms = {}
        for V, c in (terms or {}).items():
            self._add(V, c)

    # -- arithmetic --------------------------------------------------------------------
    def _add(self, V, c):
        if c:
            self.terms[V] = self.terms.get(V, 0) + c
            if self.terms[V] == 0:
                del self.terms[V]

    def copy(self):
        return ConstructibleFunction(self.ring, dict(self.terms))

    def __add__(self, other):
        r = self.copy()
        for V, c in other.terms.items():
            r._add(V, c)
        return r

    def __neg__(self):
        return ConstructibleFunction(self.ring, {V: -c for V, c in self.terms.items()})

    def __sub__(self, other):
        return self + (-other)

    def __mul__(self, other):
        if isinstance(other, int):
            return ConstructibleFunction(self.ring, {V: other * c for V, c in self.terms.items()})
        r = ConstructibleFunction(self.ring)
        for V, c in self.terms.items():
            for W, d in other.terms.items():
                r = r + indicator(V + W) * (c * d)
        return r

    __rmul__ = __mul__

    def __eq__(self, other):
        return isinstance(other, ConstructibleFunction) and self.terms == other.terms

    def __bool__(self):
        return bool(self.terms)

    def __repr__(self):
        if not self.terms:
            return "0"
        parts = []
        for V, c in sorted(self.terms.items(), key=lambda t: (-t[0].dim(), str(t[0]))):
            parts.append(f"{c:+d}·1[{V}]")
        return " ".join(parts)

    # -- evaluation ----------------------------------------------------------------------
    def value(self, X):
        """The value at a general point of the irreducible V(X)."""
        return sum(c for V, c in self.terms.items() if X.contains(V))

    def on(self, strata):
        """{dim: [value on each stratum]} for a stratification {dim: [closures]}."""
        return {d: [self.value(S) for S in L] for d, L in strata.items()}

    def integral(self, trials=2):
        """The Euler integral: sum_V c_V chi(V)."""
        return sum(c * _chi(V, trials) for V, c in self.terms.items())

    euler_characteristic = integral

    # -- Euler obstructions and characteristic cycles -------------------------------------
    def eu_decomposition(self, trials=2):
        """{W: e_W} with self = sum_W e_W Eu_W (W prime)."""
        out = ConstructibleFunction(self.ring)
        for V, c in self.terms.items():
            S = _h.whitney_stratify(V)
            m = _h.cc_multiplicities(S, trials=trials)
            dV = V.dim()
            for d, L in S.items():
                for W, mult in zip(L, m.get(d, [])):
                    out._add(W, c * (-1) ** (dV - d) * mult)
        return dict(out.terms)

    def characteristic_cycle(self, trials=2):
        """{W: multiplicity of the conormal variety of W in CC(self)}: (-1)^dim W e_W."""
        return {W: (-1) ** W.dim() * e for W, e in self.eu_decomposition(trials).items()}

    # -- functoriality --------------------------------------------------------------------
    def pull(self, F, source):
        """f^* of a function on the target, along f = F: V(source) -> target (F: polynomials
        in source's ring, one per target variable)."""
        Fp = _polys(F, source.ring)
        r = ConstructibleFunction(source.ring)
        for V, c in self.terms.items():
            r = r + indicator(source + _pullback_ideal(V, Fp)) * c
        return r

    def push(self, F, target_ring, target=None, trials=2):
        """f_* along f = F: (affine source) -> V(target) ⊆ target_ring's affine space:
        (f_* phi)(y) = integral of phi over f^{-1}(y)."""
        Y = target if target is not None else _h.Ideal(target_ring, [])
        Fp = _polys(F, self.ring)
        r = ConstructibleFunction(target_ring)
        for V, c in self.terms.items():
            r = r + _push_indicator(V, Fp, Y, trials) * c
        return r


# ------------------------------------------------------------------------------ helpers

def _polys(F, ring):
    return [f if isinstance(f, _h.Poly) else _h.Poly(ring, str(f)) for f in F]


def _pullback_ideal(W, Fp):
    return _h.Ideal(Fp[0].ring, [_h.Poly(W.ring, str(g)).compose(Fp) for g in W.groebner_basis()])


def _chi(V, trials=2):
    if V.is_unit():
        return 0
    D = _h.chern_classes(_h.whitney_stratify(V), trials=trials)
    return sum(sum(v) for v in D["euler"].values())


def _random_affine_forms(ring, k, rng):
    names = _variables(ring)
    out = []
    for _ in range(k):
        terms = [f"{rng.randint(-997, 997)}*{v}" for v in names] + [str(rng.randint(-997, 997))]
        out.append("+".join(terms).replace("+-", "-"))
    return out


def _variables(ring):
    return list(ring.variables)


# Random slices for pushforward; reseeded by helmer.set_seed.
_rng = _random.Random()


def _push_indicator(V, Fp, Y, trials):
    if len(Fp) == 1 and all(g.strip() == "0" for g in Y.gens):
        return _push_to_line(V, Fp[0], Y.ring, trials)
    _check_jelonek_flag(V, Fp, Y)
    _, tgt = _h.map_stratify(Fp, V, Y, proper=False)
    rng = _rng
    values = {}
    for d, L in tgt.items():
        values[d] = []
        for T in L:
            seen = []
            while True:
                pts = _h.Ideal(T.ring, [str(g) for g in T.groebner_basis()] +
                               _random_affine_forms(T.ring, d, rng))
                deg = _chi(pts, trials)
                if deg <= 0:
                    continue  # slice missed T (cannot happen for general slices)
                chi = _chi(V + _pullback_ideal(pts, Fp), trials)
                if chi % deg:
                    raise RuntimeError(f"push: fibre Euler characteristic {chi}/{deg} over {T}")
                seen.append(chi // deg)
                if seen.count(seen[-1]) >= trials:
                    values[d].append(seen[-1])
                    break
                if len(seen) > 4 * trials:
                    raise RuntimeError(f"push: slices over {T} disagree: {seen}")
    return from_strata(tgt, values)


def _check_jelonek_flag(V, Fp, Y):
    """map_stratify(proper=False) makes f proper over the target strata only if every map
    in the Jelonek flag is proper over a dense subset of its image."""
    X = V
    for _ in range(len(Fp) + 1):
        if X.is_unit():
            return
        J = _h.non_proper_set(Fp, X, Y)
        for Q in _h.minimal_primes(_h.image_closure(Fp, X, Y)):
            if Q.dim() > 0 and _h.contained(Q, J):  # over a point there is nothing to vary
                raise NotImplementedError(
                    "push: the map is non-proper over a dense part of the image of "
                    f"{X} (so the Jelonek flag does not control the fibres); pushforward "
                    "is implemented for such maps only when the target is a line")
        X = V + _pullback_ideal(J, Fp)


def _push_to_line(V, f, T, trials):
    """f_* 1_V for f: V -> line. Bifurcation values: critical values of f on the strata of
    V, and of the projection on the strata of (closure of the graph in P^n x C) ∪ (H x C),
    H the hyperplane at infinity, chart by chart (a superset: any Whitney refinement only
    adds critical values, and the minimal stratification of that union is global and
    contains the boundary as a union of strata). Outside them the projection is a proper
    stratified submersion, so the fibres of V are topologically constant (Thom)."""
    R = V.ring
    xs = list(R.variables)
    (svar,) = T.variables
    zero = _h.Ideal(T, [])
    cands = []

    def add(cv):
        for Q in _h.minimal_primes(cv):
            if Q.is_unit():
                continue
            if Q.dim() != 0:
                raise RuntimeError(f"push: critical values {cv} are not finite")
            if Q not in cands:
                cands.append(Q)

    for d, L in _h.whitney_stratify(V).items():
        for P in L:
            add(_h.critical_values([f], P, zero))
    # The Jelonek set of a map to a line is finite or everything. If finite, f is proper
    # over its complement and its points are candidates; the charts at infinity are only
    # needed when f is non-proper over every value (Broughton's polynomial, families of
    # affine hypersurfaces), and can be expensive.
    J = _h.non_proper_set([f], V, zero)
    if not all(g.strip() == "0" for g in J.gens):
        add(J)
        return _line_values(V, f, T, cands, trials)
    gens = [_h.Poly(R, str(g)) for g in V.groebner_basis()]
    for i, xi in enumerate(xs):
        names = ["helmer_u0"] + [x for x in xs if x != xi] + [svar, "helmer_v"]
        Rc = _h.Ring(names)
        v = _h.Poly(Rc, "helmer_v")
        images = [v if x == xi else _h.Poly(Rc, x) * v for x in xs]
        G = [g.compose(images) for g in gens]
        G.append(_h.Poly(Rc, svar) - f.compose(images))
        G.append(_h.Poly(Rc, "helmer_u0") * v - _h.Poly(Rc, "1"))
        closure = _h.eliminate(_h.Ideal(Rc, G), ["helmer_v"])
        Ri = _h.Ring(names[:-1])
        down = [_h.Poly(Ri, x) for x in names[:-1]] + [_h.Poly(Ri, "0")]
        Gi = _h.Ideal(Ri, [_h.Poly(Rc, str(g)).compose(down) for g in closure.groebner_basis()])
        U = _h.intersect(Gi, _h.Ideal(Ri, ["helmer_u0"]))
        s_i = _h.Poly(Ri, svar)
        for d, L in _h.whitney_stratify(U).items():
            for P in L:
                add(_h.critical_values([s_i], P, zero))
    return _line_values(V, f, T, cands, trials)


def _line_values(V, f, T, cands, trials):
    """f_* 1_V on the line, given every value over which the fibre may change: the value
    at two independent general points (which must agree) and over each candidate
    (Galois-conjugate values share one computation)."""
    (svar,) = T.variables
    zero = _h.Ideal(T, [])

    def at(Q):
        pts = Q if isinstance(Q, _h.Ideal) else _h.Ideal(T, [Q])
        deg = _chi(pts, trials)
        chi = _chi(V + _pullback_ideal(pts, [f]), trials)
        if chi % deg:
            raise RuntimeError(f"push: fibre Euler characteristic {chi}/{deg} over {pts}")
        return chi // deg
    generic = []
    while True:
        y0 = _rng.randint(-10**6, 10**6)
        if any(_h.Ideal(T, [f"{svar}-({y0})"]).contains(Q) for Q in cands):
            continue
        generic.append(at(f"{svar}-({y0})"))
        if generic.count(generic[-1]) >= trials:
            break
        if len(generic) > 4 * trials:
            raise RuntimeError(f"push: general fibres disagree: {generic}")
    g = generic[-1]
    out = ConstructibleFunction(T, {zero: g})
    for Q in cands:
        out._add(Q, at(Q) - g)
    return out


# ------------------------------------------------------------------------------ constructors

def indicator(I):
    """1_{V(I)} for an ideal I (any, not necessarily radical or irreducible)."""
    R = I.ring
    if I.is_unit():
        return ConstructibleFunction(R)
    P = _h.minimal_primes(I)
    r = ConstructibleFunction(R)
    if len(P) == 1:
        r._add(P[0], 1)
        return r
    # inclusion-exclusion over the components; intersections recurse
    for k in range(1, len(P) + 1):
        for sub in combinations(P, k):
            J = sub[0]
            for Q in sub[1:]:
                J = J + Q
            r = r + indicator(J) * ((-1) ** (k + 1))
    return r


def from_strata(strata, values):
    """The function equal to values[d][i] on the open stratum strata[d][i] (a Whitney
    stratification, or any stratification by prime closures with the frontier property)."""
    flat = [(d, S, values[d][i]) for d, L in strata.items() for i, S in enumerate(L)]
    if not flat:
        raise ValueError("from_strata: empty stratification")
    R = flat[0][1].ring
    flat.sort(key=lambda t: -t[0])
    coeff = []
    for d, S, v in flat:  # top down: v = sum of coefficients of the closures containing S
        above = sum(c for (e, W, c) in coeff if e > d and S.contains(W))
        coeff.append((d, S, v - above))
    return ConstructibleFunction(R, {S: c for _, S, c in coeff if c})


def euler_obstruction_function(V, trials=2):
    """Eu_V as a constructible function, for V prime."""
    S = _h.whitney_stratify(V)
    labels, M = _h.euler_obstruction_matrix(S, trials=trials)
    t = next(i for i, (_, W) in enumerate(labels) if W == V)
    values = {d: [] for d in S}
    for i, (d, W) in enumerate(labels):
        values.setdefault(d, [])
    pos = {(d, i): k for k, (d, W) in enumerate(labels) for i, X in enumerate(S[d]) if X == W}
    for d, L in S.items():
        values[d] = [M[pos[(d, i)]][t] for i in range(len(L))]
    return from_strata(S, values)
