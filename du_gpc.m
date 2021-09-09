function duk = du_gpc(dukm1, y_past, R, M, G, K1)
% duk = du_gpc(dukm1, y_past, R, M, G, K1) computes the change
% in the control action for a GPC controller.
% See script test_gpc_eqn.m for usage.
% 
% Arguments
%   dukm1 : change in control action at time k-1
%   y_past : vector of past outputs y(k-i) for i = 0, 1, ...
%   R : vector of reference values r(k+i) for i = 1:Hp
%   M : polynomial coefficients
%   G : G matrix from GPC control law
%   K1 : K1 vector from GPC control law
%
    f = G * dukm1 + M * y_past;
    duk = K1 * (R - f);
end