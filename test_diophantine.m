% Test Diophantine_equation.m

clear all

C = [1];
D = [1 -1];
[F, M]=Diophantine_equation(C, D, 1);
assert(F == 1)
assert(M == 1)

C = [0];
D = [1 -0.8];
[F, M]=Diophantine_equation(C, D, 1);
assert(F == 0)
assert(M == 0)

C = [4];
D = [1 -0.8];
[F, M]=Diophantine_equation(C, D, 1);
assert(F == 4)
assert(M == 3.2)

C = [1 -0.2];
D = [1 -0.8];
[F, M]=Diophantine_equation(C, D, 2);
assert(all(abs(F - [1.0 0.0; 1.0 0.6]) < 1e-10, [1 2]))
assert(all(abs(M - [0.60; 0.48]) < 1e-10))
