% Test predictive control functions

clear all; clc;
rng(0)

%% 1-state SISO State space model
A = -0.9;
B = 0.1;
C = 1;
D = 0;

Hp = 1; Hc = 1;

[E1,F1] = calcef(A,B,C,D,Hp);
[E2,F2] = pred_mats(A,B,C,D,Hp,Hc);

assert(all(E1 == E2, [2 1]))
assert(all(F1 == F2, [2 1]))  % Different for some reason

Hp = 3; Hc = 3;

[E1,F1] = calcef(A,B,C,D,Hp);
[E2,F2] = pred_mats(A,B,C,D,Hp,Hc);

assert(all(E1 == E2, [2 1]))
assert(all(F1 == F2, [2 1]))

Hp = 5;

[E1,F1] = calcef(A,B,C,D,Hp);
% Note: calcef does not account for short control horizons
[E2,F2] = pred_mats(A,B,C,D,Hp,Hp);

assert(all(abs(E1 - E2) < 1e-14, [2 1]))
assert(all(abs(F1 - F2) < 1e-14, [2 1]))


%% Multi-variable MIMO State space model

% Example from GEL-7029 exercise 13.1
% See also predmat_methods.m

A = [0.6807         0;
          0    1.0000];
B = [0.5000;
          0];
C = [0.7152    1.0000];
D = 0;
Hp = 5;
Hc = 1;

% André's method for F (which is E here):
f=C*A;
F=[];
for i=1:Hp
  F=[F; f];
  f=f*A;
end

% André's method for G (which is F here):
g=C*B;
G=[];
for i=1:Hp
    G=[G;g];
    g=g+C*A^i*B;
end
E1 = F;
F1 = G;

[E2,F2] = pred_mats(A,B,C,D,Hp,Hc);

assert(all(abs(F1 - F2) < 1e-14, [2 1]))
assert(all(abs(F1 - F2) < 1e-14, [2 1]))


% Another example
A = [0.7 1  0;
     0   1  1;
     0   0  1];
B = [1 0;
     0 0;
     0 1];
C = [0.3 0 0];
D = zeros(1, 2);

Hp = 1;

[E1,F1] = calcef(A,B,C,D,Hp);
[E2,F2] = pred_mats(A,B,C,D,Hp,Hp);

assert(all(E1 == E2, [2 1]))
assert(all(abs(F1 - F2) < 1e-14, [2 1]))

Hp = 3; Hc = 3;

[E1,F1] = calcef(A,B,C,D,Hp);
[E2,F2] = pred_mats(A,B,C,D,Hp,Hc);

assert(all(abs(E1 - E2) < 1e-14, [2 1]))
assert(all(abs(F1 - F2) < 1e-14, [2 1]))


%% Random multivariable state-space model

n = 5; ny = 2; nu = 3;
sys = drss(n,ny,nu);
A = sys.A;
B = sys.B;
C = sys.C;
D = sys.D;

Hp = 3; Hc = 2;

[E1, F1] = calcef(A,B,C,D,Hp);
[E2, F2] = pred_mats(A,B,C,D,Hp,Hc);

E3 = [    C*A;
        C*A^2;  
        C*A^3];
Z = zeros(2,3);
F3 = [    C*B      D    Z  Z;
        C*A*B    C*B    D  Z;
      C*A^2*B  C*A*B  C*B  D];
% Concise version (not used here):
% F = [    C*B            D;
%        C*A*B        C*B+D;
%      C*A^2*B  C*B+C*A*B+D]

assert(all(abs(E1 - E2) < 1e-12, [2 1]))
assert(all(abs(E1 - E3) < 1e-12, [2 1]))
%assert(all(abs(F1 - F2) < 1e-12, [2 1]))

%% Test prediction equations in simulation

% 1 state system
A = -0.9;
B = 0.1;
C = 1;
D = 0;
n = size(A,1);
nu = size(B,2);
ny = size(C,1);
sys = ss(A,B,C,D,1)

nT = 10;
t = (0:nT)';
U = idinput([nT+1 nu]);
x0 = 0.1;
[Y, t, X] = lsim(sys,U,t,x0);

% 1-step ahead prediction
Hp = 1; Hc = 1;
[E, F] = pred_mats(A,B,C,D,Hp,Hc);

x = X(1, :)';
u = U(1:Hc, :);
Y_p = E*x + F*u;
Y_pred = [nan; Y_p; nan(nT-Hp,1)];
table(t, U, X, Y, Y_pred)
assert(all(abs(Y_p - Y(2:1+Hp)) < 1e-12))

% 2-step ahead prediction
Hp = 2; Hc = 2;
[E, F] = pred_mats(A,B,C,D,Hp,Hc);

x = X(1, :)';
u = U(1:Hc, :);
Y_p = E*x + F*u;
Y_pred = [nan; Y_p; nan(nT-Hp,1)];
table(t, U, X, Y, Y_pred)
assert(all(abs(Y_p - Y(2:1+Hp)) < 1e-12))

% 5-step ahead prediction with Hc=3
Hp = 5; Hc = 3;
[E, F] = pred_mats(A,B,C,D,Hp,Hc);

U(Hc:end, :) = U(Hc, :);
x0 = 0.1;
[Y, t, X] = lsim(sys,U,t,x0);

x = X(1, :)';
u = U(1:Hc, :);
Y_p = E*x + F*u;
Y_pred = [nan; Y_p; nan(nT-Hp,1)];
table(t, U, X, Y, Y_pred)
assert(all(abs(Y_p - Y(2:1+Hp)) < 1e-12))

% 10-step ahead prediction with Hc=1
Hp = 10; Hc = 1;
[E, F] = pred_mats(A,B,C,D,Hp,Hc);

U(Hc:end, :) = U(Hc, :);
x0 = 0.1;
[Y, t, X] = lsim(sys,U,t,x0);

x = X(1, :)';
u = U(1:Hc, :);
Y_p = E*x + F*u;
Y_pred = [nan; Y_p; nan(nT-Hp,1)];
table(t, U, X, Y, Y_pred)
assert(all(abs(Y_p - Y(2:1+Hp)) < 1e-12))

%% Random multivariable state-space model

rng(0)
n = 5; ny = 2; nu = 3;
sys = drss(n,ny,nu);
sys.Ts = 1;
A = sys.A;
B = sys.B;
C = sys.C;
D = sys.D;

nT = 10;
t = (0:nT)';

% 1-step ahead prediction
Hp = 1; Hc = 1;
U = idinput([nT+1 nu]);
U(Hc:end, :) = repmat(U(Hc, :), nT-Hc+2, 1);
x0 = 0.1*ones(n,1);
[Y, t, X] = lsim(sys,U,t,x0);
table(t, U, X, Y)

[E, F] = pred_mats(A,B,C,D,Hp,Hc);

x = X(1, :)';
u = U(1:Hc, :)';
Y_p = E*x + F*u;
Y_pred = [nan(1,ny); Y_p'; nan(nT-Hp,ny)];
table(t, U, X, Y, Y_pred)
assert(all(abs(Y_p' - Y(2:1+Hp,:)) < 1e-12))

% 2-step ahead prediction
Hp = 2; Hc = 2;
U = idinput([nT+1 nu]);
U(Hc:end, :) = repmat(U(Hc, :), nT-Hc+2, 1);
x0 = 0.1*ones(n,1);
[Y, t, X] = lsim(sys,U,t,x0);
table(t, U, X, Y)

[E, F] = pred_mats(A,B,C,D,Hp,Hc);

x = X(1, :)';
u = reshape(U(1:Hc, :)', [Hc*nu 1]);
Y_p = E*x + F*u;
Y_p = reshape(Y_p, [ny, Hp])';
Y_pred = [nan(1,ny); Y_p; nan(nT-Hp,ny)];
table(t, U, X, Y, Y_pred)
assert(all(abs(Y_p - Y(2:1+Hp,:)) < 1e-12, [2 1]))

% 5-step ahead prediction with Hc=3
Hp = 5; Hc = 3;
U = randn([nT+1 nu]);
U(Hc:end, :) = repmat(U(Hc, :), nT-Hc+2, 1);
x0 = 0.1*ones(n,1);
[Y, t, X] = lsim(sys,U,t,x0);
table(t, U, X, Y)

[E, F] = pred_mats(A,B,C,D,Hp,Hc);

x = X(1, :)';
u = reshape(U(1:Hc, :)', [Hc*nu 1]);
Y_p = E*x + F*u;
Y_p = reshape(Y_p, [ny, Hp])';
Y_pred = [nan(1,ny); Y_p; nan(nT-Hp,ny)];
table(t, U, X, Y, Y_pred)
assert(all(abs(Y_p - Y(2:1+Hp,:)) < 1e-12, [2 1]))

% 10-step ahead prediction with Hc=1
Hp = 10; Hc = 1;
U = idinput([nT+1 nu]);
U(Hc:end, :) = repmat(U(Hc, :), nT-Hc+2, 1);
x0 = 0.1*ones(n,1);
[Y, t, X] = lsim(sys,U,t,x0);
table(t, U, X, Y)

[E, F] = pred_mats(A,B,C,D,Hp,Hc);

x = X(1, :)';
u = reshape(U(1:Hc, :)', [Hc*nu 1]);
Y_p = E*x + F*u;
Y_p = reshape(Y_p, [ny, Hp])';
Y_pred = [nan(1,ny); Y_p; nan(nT-Hp,ny)];
table(t, U, X, Y, Y_pred)
assert(all(abs(Y_p - Y(2:1+Hp,:)) < 1e-12, [2 1]))

disp("Tests complete")