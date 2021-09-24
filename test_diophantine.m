% Test diophantine_recursive.m

clear all


%% Test sub_poly.m

assert(isequal(sub_poly(8, 5), 3))
a = [1 2 3];
b = [1 2];
assert(isequal(sub_poly(a, a), [0 0 0]))
assert(isequal(sub_poly(b, b), [0 0]))
assert(isequal(sub_poly(a, b), [0 0 3]))
assert(isequal(sub_poly(b, a), [0 0 -3]))

a = [0.5 -0.6 0.7 -0.8];
assert(isequal(sub_poly(a, a), [0 0 0 0]))
assert(isequal(sub_poly(a, b), [-0.5 -2.6  0.7 -0.8]))
assert(isequal(sub_poly(b, a), [ 0.5 2.6  -0.7 0.8]))

% Example from docstring
p = sub_poly([1 0.2], [1 0.1 0.8]);
assert(isequal(p, [0 0.1 -0.8]))

% Edge cases
p = sub_poly([1 0.2], []);
assert(isequal(p, [1 0.2]))
p = sub_poly([], [1 0.2]);
assert(isequal(p, [-1 -0.2]))


%% Test diophantine.m

assert(diophantine(3, 4) == 0.75)
A = [0.18 1.14 1.43 0.48 0.49];
[quotient, remainder] = diophantine(A, A);
assert(quotient == 1)
assert(isequal(remainder, zeros(1, length(A) - 1)))

C = 0; D = [1 -0.8];
[quotient, remainder] = diophantine(C, D);
assert(quotient == 0.0)
assert(remainder == 0.0)

C = 4; D = [1 -0.8];
[quotient, remainder] = diophantine(C, D);
assert(quotient == 4)
assert(remainder == 3.2)


%% Test diophantine_recursive.m

C = [1];
D = [1 -1];
[F, M] = diophantine_recursive(C, D, 1);
assert(F == 1)
assert(M == 1)

C = [0];
D = [1 -0.8];
[F, M] = diophantine_recursive(C, D, 1);
assert(F == 0)
assert(M == 0)

C = [4];
D = [1 -0.8];
[F, M] = diophantine_recursive(C, D, 1);
assert(F == 4)
assert(M == 3.2)

C = [1 -0.2];
D = [1 -0.8];
[F, M] = diophantine_recursive(C, D, 2);
assert(all(abs(F - [1.0 0.0; 1.0 0.6]) < 1e-10, [1 2]))
assert(all(abs(M - [0.60; 0.48]) < 1e-10))

C = [1];
D = [1 -1.8 0.8];
[F, M] = diophantine_recursive(C, D, 3);
F_test = [ ...
    1.0000         0         0;
    1.0000    1.8000         0;
    1.0000    1.8000    2.4400 ...
];
M_test = [ ...
    1.8000   -0.8000;
    2.4400   -1.4400;
    2.9520   -1.9520 ...
];
assert(all(abs(F - F_test) < 1e-10, [1 2]))
assert(all(abs(M - M_test) < 1e-10, [1 2]))

C = [1 0.6 -0.2];
D = [1 0.6 -0.2];
[F, M] = diophantine_recursive(C, D, 3);
F_test = [ ...
    1    0     0;
    1    0     0;
    1    0     0 ...
];
M_test = zeros(3, 2);
assert(all(abs(F - F_test) < 1e-10, [1 2]))
assert(all(abs(M - M_test) < 1e-10, [1 2]))
