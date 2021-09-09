function [E,F] = pred_mats(A,B,C,D,Hp,Hc)
% [E,F] = pred_mats(A,B,C,D,Hp,Hc) computes prediction
% matrices E, F in the prediction equation for Y(k+Hp):
%
%   Y(k+Hp) = E x(k) + F U(k)
%
%  where
%
%   Y(k+Hp) = [y(k+1) y(k+2) ... y(k+Hp)]'
%      y(k) = [y_1(k) y_2(k) ... y_ny(k)]'
%      x(k) = [x_1(k) x_2(k) ... x_n(k)]'
%      U(k) = [u(k) u(k+1) ... u(k+Hc)]'
%      u(k) = [u_1(k) u_2(k) ... u_nu(k)]'
% 
% Example 1: for Hp=2, Hc=1:
% E = [C*A;  C*A^2];
% F = [C*B+D;
%      C*A*B+C*B+D];
%
% Example 2: for Hp=3, Hc=2:
% E = [C*A;   C*A^2;  C*A^3];
% F = [C*B      D;
%      C*A*B    C*B+D;
%      C*A^2*B  C*B+C*A*B+D];

    % System dimensions
    n = size(A,1);  % number of states
    nu = size(B,2);  % number of inputs
    ny = size(C,1);  % number of outputs
    if n == 1
        % System with only 1 state
        E = C*A.^(1:Hp)';
        F = toeplitz(C*(A.^(0:Hp-1))*B, [C*B D zeros(1, Hp-2)]);
        % Add cols Hc+1:Hp+1 to col Hc
        F(:, Hc) = sum(F(:, Hc:end), 2);
        % Drop cols Hc+1:Hp+1
        F = F(:, 1:Hc);
    else
        A_blocks = cellfun(@(i) A^i, num2cell(0:Hp), 'UniformOutput', false);
        E_blocks = cellfun(@(X) C*X, A_blocks(2:end), 'UniformOutput', false);
        E = cell2mat(E_blocks(1:Hp)');
        F_blocks = [cellfun(@(X) C*X*B, A_blocks(1:Hp), 'UniformOutput', false) {D} {zeros(ny,nu)}];
        F = cell2mat(F_blocks(toeplitz(1:Hp, [1 Hp+1 (Hp+2)*ones(1,Hp-1)])));
        % Add blocks in cols Hc+1:Hp+1 to blocks in col Hc
        for i = 1:nu
            F(:, (Hc-1)*nu + i) = sum(F(:, (Hc-1:Hp)*nu+i), 2);
        end
        % Drop cols Hc+1:Hp+1
        F = F(:, 1:Hc*nu);
    end
end
    