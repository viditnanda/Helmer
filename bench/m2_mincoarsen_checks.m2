-- Two cases where WhitneyStratifications 2.23's minCoarsenWS returns a stratification that
-- is not Whitney (see README, "Differences from the Macaulay2 package"). Helmer's
-- min_coarsen keeps the point in both (tests/test_local.py).
--   M2 --script bench/m2_mincoarsen_checks.m2
needsPackage "WhitneyStratifications";
R = QQ[x,y,z];

-- (c) cone ∪ plane z = 5: the cone's vertex lies on no curve stratum. Expect the vertex in
-- dimension 0; minCoarsenWS drops it.
V = whitneyStratify(ideal((x^2+y^2-z^2)*(z-5)), AssocPrimes=>false);
print("cone + plane, whitneyStratify: " | toString(peek V));
print("cone + plane, minCoarsenWS:    " | toString(peek minCoarsenWS V));

-- (b) umbrella ∪ plane x = 5: the origin (the umbrella's polar curve passes through it)
-- must stay. The polar test reads the top level by position, so the answer depends on the
-- order of the two surfaces.
strata = L -> (H := new MutableHashTable; scan(L, p -> H#(first p) = last p); H);
pts = {ideal(x, y, z)};
curves = {ideal(x-5, y^2*z-25), ideal(x, y)};
A = strata {(2, {ideal(x^2-y^2*z), ideal(x-5)}), (1, curves), (0, pts)};
B = strata {(2, {ideal(x-5), ideal(x^2-y^2*z)}), (1, curves), (0, pts)};
print("umbrella first: " | toString(peek minCoarsenWS A));
print("plane first:    " | toString(peek minCoarsenWS B));
