function [M, G2, K1] = gpc_arrays(A, B, C, D, Hp, Hc, lambda)
% [M, G2, K1] = gpc_arrays(A, B, C, D, Hp, Hc, lambda)
% computes arrays needed to implement GPC controller.
% See script test_gpc_eqn.m for usage.

    % Solve diophantine equations
    [F, M] = Diophantine_equation(C,D,Hp)

    % Prepare G matrix
    G = tril(B*toeplitz(F(Hp,:)));

    % Split into free and forced components:
    G1 = G(:,2:end)
    G2 = G(:,1)

    % Results from Homework 5:
    G_test = [
        0     0     0;
        0.4   0     0;
        0.72  0.4   0;
        0.976 0.72  0.4;
    ];
    assert(all(abs(G1 - G_test) < 1e-10, [1,2]))

    M_test = [
        1.8     -0.8;
        2.44    -1.44;
        2.952   -1.952;
        3.3616  -2.3616;
    ];
    assert(all(abs(M - M_test) < 1e-10, [1,2]))

    % Truncate G to match the control horizon
    Gc = G1(:,1:Hc);

    % GPC control law
    K1 = (Gc'*Gc + lambda)^-1*Gc';

end