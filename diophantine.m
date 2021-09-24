function [quotient, remainder] = diophantine(C, D)
% [quotient, remainder] = diophantine(C, D)
% Solves the following diophantine equation:
% 
%     C / D = quotient + remainder / D
% 
% where `C`, `D`, `quotient` and `remainder` are
% polynomials represented as lists or arrays containing
% the coefficients in decreasing powers of z.
%
    quotient = C(1) / D(1);
    remainder = sub_poly(C, quotient .* D);
    assert(remainder(1) == 0)
    remainder = remainder(2:end);
    
end