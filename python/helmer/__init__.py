"""Helmer: Whitney stratifications of real and complex algebraic varieties.

Named in memory of Martin Helmer, whose Macaulay2 package
WhitneyStratifications this library reimplements in C++.
"""
from ._helmer import (Ring, Poly, Ideal, minimal_primes, eliminate, saturate,
                      intersect, conormal, whitney_stratify, map_stratify,
                      non_proper_set, polar_varieties, set_seed, stats,
                      certify_membership, verify_membership, eliminate_certified,
                      multiplicity, contained, in_singular_locus, polar_sequence,
                      polar_multiplicity, min_coarsen)

__all__ = ["Ring", "Poly", "Ideal", "minimal_primes", "eliminate", "saturate",
           "intersect", "conormal", "whitney_stratify", "map_stratify",
           "non_proper_set", "polar_varieties", "set_seed", "stats",
           "certify_membership", "verify_membership", "eliminate_certified",
           "multiplicity", "contained", "in_singular_locus", "polar_sequence",
           "polar_multiplicity", "min_coarsen"]
__version__ = "1.0.0"
