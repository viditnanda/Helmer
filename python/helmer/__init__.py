"""Helmer: Whitney stratifications of real and complex algebraic varieties.

Named in memory of Martin Helmer, whose Macaulay2 package
WhitneyStratifications this library reimplements in C++.
"""
from ._helmer import (Ring, Poly, Ideal, minimal_primes, eliminate, saturate,
                      intersect, conormal, whitney_stratify, whitney_filtration, map_stratify,
                      non_proper_set, polar_varieties, set_seed, stats,
                      certify_membership, verify_membership, eliminate_certified,
                      multiplicity, contained, in_singular_locus, polar_sequence,
                      polar_multiplicity, min_coarsen,
                      euler_obstruction, global_polar_degrees, global_euler_obstruction,
                      euler_obstruction_matrix, cc_multiplicities, chern_classes)

__all__ = ["Ring", "Poly", "Ideal", "minimal_primes", "eliminate", "saturate",
           "intersect", "conormal", "whitney_stratify", "whitney_filtration", "map_stratify",
           "non_proper_set", "polar_varieties", "set_seed", "stats",
           "certify_membership", "verify_membership", "eliminate_certified",
           "multiplicity", "contained", "in_singular_locus", "polar_sequence",
           "polar_multiplicity", "min_coarsen",
           "euler_obstruction", "global_polar_degrees", "global_euler_obstruction",
           "euler_obstruction_matrix", "cc_multiplicities", "chern_classes",
           "euler_characteristic", "csm_class"]
__version__ = "1.0.0"


def _strata(ideal, projective):
    return whitney_stratify(ideal, projective=projective)


def euler_characteristic(ideal, projective=False, trials=2):
    """chi of V(ideal) (affine or projective), from a Whitney stratification."""
    D = chern_classes(_strata(ideal, projective), projective=projective, trials=trials)
    return sum(sum(v) for v in D["euler"].values())


def csm_class(ideal, projective=False, trials=2):
    """The Chern-Schwartz-MacPherson class of V(ideal) in A_*(P^N) as [a_0, ..., a_N]
    (a_k = degree of the dimension-k part); an affine variety is taken in A^N ⊂ P^N."""
    D = chern_classes(_strata(ideal, projective), projective=projective, trials=trials)
    total = [0] * (D["ambient"] + 1)
    for classes in D["csm"].values():
        for c in classes:
            total = [a + b for a, b in zip(total, c)]
    return total


from . import constructible
from .constructible import (ConstructibleFunction, indicator, from_strata,
                            euler_obstruction_function)

_set_seed_compiled = set_seed


def set_seed(seed):
    """Seed Helmer's random choices (Gröbner primes, slices, and the random linear
    sections used by constructible-function pushforward)."""
    _set_seed_compiled(seed)
    constructible._rng.seed(seed)


__all__ += ["ConstructibleFunction", "indicator", "from_strata", "euler_obstruction_function"]
