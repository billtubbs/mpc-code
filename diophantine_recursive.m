function [F, M] = diophantine_recursive(C, D, hp)
% Recursive solution of the Diophantine equation.
% 
% Returns the polynomials `F` and `M` that satisfy
% 
%     C / D = F + z^(-j) * M / D
% 
% for j = 1 to `hp`, where `C`, `D`, `F`, and `M` are
% polynomials represented as lists or arrays containing
% the coefficients in decreasing powers of z.
% 
% Arguments:
%     C : array
%         Coefficients of `C`.
%     D : array
%         Coefficients of `D`.
%     hp : int
%         Number of iterations of recursion.
% 
% Returns:
%     F : array (2-D)
%         Coefficients of `F` in rows for j = 1 to `hp`.
%     M : array (2-D)
%         Coefficients of `M` in rows for j = 1 to `hp`.
% 
% Example:
% >> # Solve for 2-steps ahead:
% >> # (1 - 0.2z^-1) / (1 - 0.8z^-1) = 1 + 0.6z^-1 + 0.48z^-2 / (1 - 0.8z^-1)
% >> [F, M] = diophantine_recursive([1 -0.2], [1 -0.8], 2)
% 
% F =
% 
%     1.0000         0
%     1.0000    0.6000
% 
% 
% M =
% 
%     0.6000
%     0.4800
%
    F = zeros(hp, hp);
    M = zeros(hp, max(length(C), length(D)) - 1);
    quotients = nan(1, hp);
    numerator = C;
    for j = 1:hp
        [q, r] = diophantine(numerator, D);
        quotients(j) = q;
        F(j, 1:j) = quotients(1:j);
        M(j, :) = r;
        numerator = r;
    end