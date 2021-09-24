% Test solve_mpc_MIMO.m

clear all

% Random multivariable state-space model

rng(0)
n = 5; ny = 2; nu = 3;
sys = drss(n,ny,nu);
sys.Ts = 1;
A = sys.A;
B = sys.B;
C = sys.C;
D = sys.D;

% MPC settings
Hp = 15;
Hc = 3;
psi = [5 1];
lambda = [0 0 0];

% Reference (future setpoints)
R = repmat((1:ny)', [Hp, 1]);
%R(2:end, :) = ones(Hp-1,ny);

ukm1 = zeros(nu, 1);  % previous input
xk = 0.1*randi(5,[n 1]);  % initial condition

% Prediction matrices
[E,F] = pred_mats(A,B,C,D,Hp,Hc);

% MPC control actions
uk = solve_mpc_MIMO(xk,ukm1,R,E,F,psi,lambda)

% Output prediction
y_pred = E*xk + F*uk;

% Predicted output error
errors = R - y_pred;
J = sum(psi * reshape(errors.^2, [ny, Hp]))

% Check it is a local minimum
d = 1e-5;
for i=1:nu
    uk1 = uk; uk1(i) = uk1(i) - d;
    y_pred1 = E*xk + F*uk1;
    J1 = sum(psi * reshape((R - y_pred1).^2, [ny, Hp]));
    uk2 = uk; uk2(i) = uk2(i) - d;
    y_pred2 = E*xk + F*uk2;
    J2 = sum(psi * reshape((R - y_pred2).^2, [ny, Hp]));
    assert(all([J2 J2] > J))
end

% Simulate system under control actions
k_sim = (0:Hp)';
U = reshape(uk, [nu, Hc])';
u_sim = [U; repmat(U(end,:), Hp-Hc+1, 1)];
[y_sim, t, x_sim] = lsim(sys, u_sim, k_sim, xk);


% Plot reference, predicted and simulated output
R_sim = [zeros(1,ny); reshape(R, [ny, Hp])'];
y_pred_sim = [zeros(1,ny); reshape(y_pred, [ny, Hp])'];

figure(1); clf

subplot(2,1,1)
stairs(k_sim, R_sim, 'k--'); hold on
stairs(k_sim, y_sim, 'o');
stairs(k_sim, y_pred_sim, 'Linewidth', 2);
xlabel('k')
ylabel('y_pred(k)', 'Interpreter', 'none')
grid on
labels = cell(1, ny*3);
for i=1:ny
    labels(i) = {sprintf('$r_%d(k)$', i)};
    labels(ny+i) = {sprintf('$y_%d(k)$', i)};
    labels(ny*2+i) = {strcat('$\hat{y}', sprintf('_%d(k)$', i))}; 
end
legend(labels, 'Interpreter', 'Latex')

subplot(2,1,2)
stairs(k_sim, u_sim, 'Linewidth', 2)
xlabel('k')
ylabel('u(k)')
grid on
labels = cell(1, nu);
for i=1:nu
    labels(i) = {sprintf('$u_%d(k)$', i)};
end
legend(labels, 'Interpreter', 'Latex')