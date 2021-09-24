% Test GPC control law calculation
% This is based on GEL-7029 homework 6, Q2

rng(1)  % seed for RNG

% CARIMA model
a = 0.8;
b = 0.4;
d = 2;
c = 0;
A = [1 -a];
B = b;
C = [1];
D = conv(A,[1 -1]);
d = 2;

% GPC parameters
Hp = 4; 
Hc = 1;
lambda = 0.1;

% Solve diophantine equations
[F, M] = diophantine_recursive(C,D,Hp)

% Prepare G matrix
G = tril(B*toeplitz(F(Hp,:)));

% Split into free and forced components:
G1 = G(:,2:end)
G2 = G(:,1)

% Results from Homework 5:
G1_test = [
    0     0     0;
    0.4   0     0;
    0.72  0.4   0;
    0.976 0.72  0.4;
];
assert(all(abs(G1 - G1_test) < 1e-10, [1,2]))

G2_test = [
    0.4;
    0.72;
    0.976;
    1.1808;
];
assert(all(abs(G2 - G2_test) < 1e-10))

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
K1 = (Gc'*Gc + lambda)^-1*Gc'
K1_test = [0 0.4 0.72 0.976] ./ (lambda + 1.63098);
assert(all(abs(K1 - K1_test) < 1e-5))

% Calculate using external function
[M, G2, K1] = gpc_arrays(A, B, C, D, Hp, Hc, lambda);
assert(all(abs(M - M_test) < 1e-10, [1,2]))
assert(all(abs(G2 - G2_test) < 1e-10, [1,2]))
assert(all(abs(K1 - K1_test) < 1e-5))

% Simulation setup
nT = 20;
Vq = 0.001;
eint = false;

% Run simulation
sim_data = simulate(nT,a,b,c,d,Hp,lambda,M,G2,K1,Vq,eint)

% Make plots
annotation = "\lambda="+lambda;
fig = plot_results(sim_data,annotation);

% Check simulation output with MATLAB sim function
sys = idpoly(A,B,C,1,1,1,1,'IODelay',d);
opt = simOptions('AddNoise', true, 'NoiseData', sim_data.e + sim_data.q);
y_sim = sim(sys, sim_data.u, opt);
assert(all(abs(y_sim - sim_data.y) < 1e-10))

% Check outputs match original from homework
% See file gel-7029/homework/hw6/hw6_p2.mlx
y_test = [
   -0.0205    0.0209    0.5074    0.6814    0.8384 ...
    0.9429    0.9865    1.0189    1.0282    1.0251 ...
    1.4980    1.8936    1.4707    1.2082    1.1532 ...
    1.0456    0.9672    0.9179    1.0600    1.1602 ...
    1.0478 ...
]';
assert(all(abs(y_test - sim_data.y) < 1e-4))


%% Simulation and plotting functions

function results = simulate(nT,a,b,c,d,Hp,lambda,M,G,K1,Vq,eint)
    % nT : number of timesteps
    % a,b,c : process model coefficients
    % d  : process time delay
    % Hp : Prediction horizon
    % lambda : GPC parameter
    % Vq : noise variance
    % eint : integrate noise if true
    
    k_ind = [-d:nT+Hp]';  % timesteps index (starts at k=-d)
    
    % Temporary arrays to store data
    r = zeros(size(k_ind));
    u = zeros(size(k_ind));
    du = zeros(size(k_ind));
    e = zeros(size(k_ind));
    q = zeros(size(k_ind));
    y = zeros(size(k_ind));
    
    % Reference signal
    r(k_ind >= 0) = ones(nT+Hp+1,1);  % step
    
    % Noise signal
    e(k_ind >= 0) = sqrt(Vq)*randn(nT+Hp+1,1);
    
    if eint == true
        % Integrate noise
        e = cumsum(e);
    end
    
    % Add a step disturbance
    q(k_ind >= 10) = 0.5;
    
    i_ind = [];  % array indices
    for k=0:nT
    
        % Array index for timestep k
        i = find(k_ind == k);

        % Process model
        y(i) = a*y(i-1) + b*u(i-d) + e(i) + c*e(i-1);       
        y(i) = y(i) + q(i);  % add output disturbance

        % GPC law - derived by hand in homework
%         du(i) = 1/(lambda + 1.631) * (0.4*r(i+2) + 0.72 * r(i+3) ...
%                    + 0.976 * r(i+4) - 2.1432 * du(i-1) - 6.3824 * y(i) ...
%                    + 4.2864 * y(i-1));
        % GPC law
        y_past = [y(i); y(i-1)];
        r_fut = r(i+1:i+Hp);
        du(i) = du_gpc(du(i-1), y_past, r_fut, M, G, K1);
        u(i) = u(i-1) + du(i);

        i_ind = [i_ind i];  % store array indices
    
    end
    
    % Combine results
    results = table( ...
        k_ind(i_ind), r(i_ind), u(i_ind), du(i_ind), ...
        e(i_ind), q(i_ind), y(i_ind), ...
        'VariableNames', {'k','r','u','du','e','q','y'} ...
    );

end

function fig = plot_results(sim_data,annotation)
    fig = figure;
    subplot(2,1,1)
    plot(sim_data.k, sim_data.y, 'o-', 'LineWidth',2); hold on
    plot(sim_data.k, sim_data.r, 'k--', 'LineWidth',1);
    ylim([-0.5 2.5])
    grid on
    title('Output vs Reference')
    legend('y(k)','r(k)','Location','best')
    text(1,1.9,annotation)
    subplot(2,1,2)
    stairs(sim_data.k, sim_data.u, 'LineWidth',2);
    grid on
    title('Input')
    ylim([-5 5])
    legend('u(k)','Location','best')
    xlabel('t')
end
