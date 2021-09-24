% Test solve_mpc_SISO

clear all
rng(0)

% Multi-variable SISO State space model

% Example from GEL-7029 exercise 13.1
% See also predmat_methods.m

A = [0.6807         0;
          0    1.0000];
B = [0.5000;
          0];
C = [0.7152    1.0000];
D = 0;
n = size(A,1);
nu = size(B,2);
ny = size(C,1);
sys = ss(A,B,C,D,1)

% MPC settings
Hp = 10;
Hc = 5;
lambda = 0;

% Reference (future setpoints)
R = zeros(Hp,1);
R(2:3) = [1.5; 1];
R(4:end) = 0.25;

ukm1 = 0;  % previous input
xk = 0.1*randi(5,[n 1]);  % initial condition

% Prediction matrices
[E,F] = pred_mats(A,B,C,D,Hp,Hc);

% MPC control actions
uk = solve_mpc_SISO(xk,ukm1,R,E,F,lambda)

% Output prediction
y_pred = E*xk + F*uk;

% Predicted output error
errors = R - y_pred;
J = sum(errors.^2) + lambda * sum(diff([ukm1; uk]).^2)

% Check it is a local minimum
d = 0.01;
uk1 = uk - d;
y_pred1 = E*xk + F*uk1;
J1 = sum((R - y_pred1).^2) + lambda * sum(diff([ukm1; uk1]).^2);
uk2 = uk + d;
y_pred2 = E*xk + F*uk2;
J2 = sum((R - y_pred1).^2) + lambda * sum(diff([ukm1; uk2]).^2);
assert(all([J1 J2] > J))

% Simulate system under control actions
k_sim = (0:Hp)';
u_sim = [uk; repmat(uk(end), Hp-Hc+1, 1)];
[y_sim, t, x_sim] = lsim(sys, u_sim, k_sim, xk);

% Plot reference, predicted and simulated output

R_sim = [0; R];
y_pred_sim = [0; y_pred];

figure(1); clf

subplot(2,1,1)
stairs(k_sim, R_sim, 'k--'); hold on
stairs(k_sim, y_sim, 'o');
stairs(k_sim, y_pred_sim, 'Linewidth', 2);
xlabel('k')
ylabel('y_pred(k)', 'Interpreter', 'none')
grid on
legend('r(k)', 'y(k)', 'y_pred(k)', 'Interpreter', 'none')

subplot(2,1,2)
stairs(k_sim, u_sim, 'Linewidth', 2)
xlabel('k')
ylabel('u(k)')
grid on