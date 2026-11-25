-- Macaulay2 side of bench/characteristic_compare.py.  M2 --script bench/m2_characteristic_checks.m2
-- CSM I = sum_k a_k H^(N-k); Helmer prints [a_0, ..., a_N].
needsPackage "CharacteristicClasses"; needsPackage "WhitneyStratifications";
R = QQ[x,y,z,w]; I = ideal(x^2*w-y^3+y*z^2);
print("D4_P  CSM " | toString CSM I | "   Euler " | toString Euler I);
R = QQ[x,y,z,w]; I = ideal(x*y*z+x*y*w+x*z*w+y*z*w);
print("cayley_P  CSM " | toString CSM I | "   Euler " | toString Euler I);
R = QQ[a_0,a_1,a_2,a_3,a_4]; I = ideal(a_0^2*a_4-a_1*a_2^2+a_3^3-a_3*a_0^2-a_4^2*a_3);
print("quartic_P  CSM " | toString CSM I | "   Euler " | toString Euler I);
R = QQ[x,y,z,w]; I = ideal(y^2*z-x^2*w);
print("umbrella_P  CSM " | toString CSM I | "   Euler " | toString Euler I);
R = QQ[x,y,z]; I = ideal(y^2*z-x^3-x^2*z);
print("nodal_cubic  CSM " | toString CSM I | "   Euler " | toString Euler I);
R = QQ[x,y,z]; I = ideal(y^2*z-x^3);
print("cuspidal_cubic  CSM " | toString CSM I | "   Euler " | toString Euler I);
R = QQ[x,y,z,w]; I = ideal(x*y*z+x*y*w+x*z*w+y*z*w);
print("cayley_cubic  CSM " | toString CSM I | "   Euler " | toString Euler I);
R = QQ[x,y,z]; V = whitneyStratify(ideal(x^2-y^2*z), AssocPrimes=>false);
C = multCC V; print("umbrella  multCC " | toString last C);
R = QQ[x,y,z]; V = whitneyStratify(ideal(x^3+y^3+z^3), AssocPrimes=>false);
C = multCC V; print("cubic_cone  multCC " | toString last C);
R = QQ[x,y,z,t]; V = whitneyStratify(ideal(x^3+y^4*z+z^9+t*x*y^3), AssocPrimes=>false);
C = multCC V; print("bs_small  multCC " | toString last C);
R = QQ[x_1,x_2,x_3,x_4]; V = whitneyStratify(ideal(x_1^6+x_2^6+x_1^4*x_3*x_4+x_3^3), AssocPrimes=>false);
C = multCC V; print("HM_X2  multCC " | toString last C);
R = QQ[x_1,x_2,x_3,x_4]; V = whitneyStratify(ideal(x_1^7-2*x_1^5*x_4+x_1^3*x_4^2+x_1^2*x_2^2-x_1*x_2^2*x_3+x_2^3), AssocPrimes=>false);
C = multCC V; print("HM_ex1.1  multCC " | toString last C);
