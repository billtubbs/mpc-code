function uk = solve_mpc_MIMO(xk,ukm1,R,E,F,psi,lambda)
% U = solve_mpc_MIMO(xk,ukm1,R,E,F,lambda) solves the model-based
% predictive control problem in the case of a MIMO linear system.
%
% Control criterion which is minimized:
%
%  J = (sum(psi(m) * sum((r(k+i) - y(k+i))^2) for i = 1:Hp)
%       for m=1:ny)
%      + (sum(lambda(m) * sum((u(k+i-1) - u(k+i-2))^2) for i = 1:Hc)
%         for m=1:nu)
%
%  uk = argmin(J(uk))
%
% See p. 31 of GlobPC course notes (GEL-7029): 
%  - MPC_A_general_framework_Desbiens_etal_2016.pdf
%    (Desbiens, Hodouin, Milot, August 2016)
%
    % Control horizon length
    nu = length(lambda);
    Hc = size(F,2) / nu;
    % Prediction horizon length
    ny = length(psi);
    Hp = size(E,1) / ny;
    H = R - E*xk;
    C = {-ones(3), zeros(3), ones(3)};  % block matrices
    P = diag(-ones(1, Hc-1), -1) + eye(Hc);
    P = cell2mat(C(P + 2));
    Q = zeros(Hc, 1); Q(1,1) = 1;
    Q = cell2mat(C(Q + 2)');
    Psi = diag(repmat(psi,1,Hp));
    Lambda = diag(repmat(lambda,1,Hc));
    uk = (F'*Psi*F + P'*Lambda*P) \ (F'*Psi*H + P'*Lambda*Q*ukm1);
end