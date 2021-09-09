function uk = solve_mpc_SISO(xk,ukm1,R,E,F,lambda)
% U = solve_mpc_SISO(xk,ukm1,R,E,F,lambda) solves the model-based
% predictive control problem in the case of a SISO linear system.
%
% Control criterion which is minimized:
%
%  J = sum((r(k+i) - y(k+i))^2) for i = 1:Hp
%      + lambda * sum((u(k+i-1) - u(k+i-2))^2) for i = 1:Hc
%
%  uk = argmin(J(uk))
%
% See p. 28 of GlobPC course notes (GEL-7029): 
%  - MPC_A_general_framework_Desbiens_etal_2016.pdf
%    (Desbiens, Hodouin, Milot, August 2016)
%
    Hc = size(F,2);
    H = R - E*xk;
    P = diag(-ones(1, Hc-1), -1) + eye(Hc);
    Q = zeros(Hc, 1); Q(1,1) = 1;
    uk = (F'*F + lambda*(P'*P)) \ (F'*H + lambda*(P'*Q*ukm1));
end