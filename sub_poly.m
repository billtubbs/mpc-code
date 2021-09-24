function p = sub_poly(a, b)
% p = sub_poly(a, b)
% Returns the coefficients of the polynomial
% P(z^-1) resulting from the subtraction of two 
% polynomials A(z^-1) and B(z^-1).
% 
%     P(z^-1) = A(z^-1) - B(z^-1)
% 
% where the polynomial coefficients `p`, `a`, 
% and `b` are represented as arrays containing 
% the coefficients in decreasinghelp_ powers of z.
% 
% `a` and `b` do not need to be the same length.
% 
% Parameters:
%     a : Coefficients of A(z^-1)
%     b : Coefficients of B(z^-1)
% 
% Returns:
%     p : array
%         Coefficients of P(z^-1)
% 
% Example:
% >> sub_poly([1 0.2], [1 0.1 0.8])
% [ 0.0  0.1 -0.8]

    la = length(a);
    lb = length(b);

    if la == lb
        p = a - b;
    else
        p = zeros(1, max(la, lb));
        p(1:la) = a;
        p(1:lb) = p(1:lb) - b;
    end

end