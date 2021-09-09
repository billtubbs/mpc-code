


%% 1-state SISO State space model
A = -0.9;
B = 0.1;
C = 1;
D = 0;

Hp = 1; Hc = 1;

[E1,F1] = calcEF1(A,B,C,D,Hp);
[E2,F2] = calcEF2(A,B,C,D,Hp);

assert(all(E1 == E2, [2 1]))
assert(all(F1 == F2, [2 1]))  % Different for some reason


%% Speed tests of both methods

n = 14; ny = 3; nu = 6;
sys = drss(n,ny,nu);
A = sys.A;
B = sys.B;
C = sys.C;
D = sys.D;

Hp = 30; Hc = 2;

% André's function from GlobPC
fprintf("calcEF1: ")
tic
[E1,F1] = calcEF1(A,B,C,D,Hp);
toc

% My code
fprintf("calcEF2: ")
tic
[E2,F2] = calcEF2(A,B,C,D,Hp);
toc

% My code
fprintf("calcEF3: ")
tic
[E3,F3] = calcEF3(A,B,C,D,Hp);
toc

assert(all(abs(E1 - E2) < 1e-14, [2 1]))
assert(all(abs(F1 - F2) < 1e-12, [2 1]))
assert(all(abs(E1 - E3) < 1e-14, [2 1]))
assert(all(abs(F1 - F3) < 1e-12, [2 1]))

% Run repeated tests
f1 = @() calcEF1(A,B,C,D,Hp);
f2 = @() calcEF2(A,B,C,D,Hp);
f3 = @() calcEF3(A,B,C,D,Hp);
results = zeros(1, 3);
n = 5;
for i=1:n
    t1 = timeit(f1); 
    t2 = timeit(f2); 
    t3 = timeit(f2);
    fprintf("%f %f %f\n", t1, t2, t3)
    results = results + [t1 t2 t3];
end
results = results ./ n;
fprintf("Avgs: %f %f %f\n", results)


%% Functions to compare

function [E,F] = calcEF1(A,B,C,D,j)
% Calculate the prediction matrices E et F

    % System dimensions
    cA = size(A,2);
    cB = size(B,2);
    lC = size(C,1);
    cD = size(D,2);

	E = zeros(lC*j,cA);
	F = zeros(lC*j,cB*j+cD);
	E(1:lC,:) = C*A;
	F(1:lC,1:cB+cD) = [C*B D];
    for i=1:(j-1)
		E(lC*i+1:lC*(1+i),:) = E(lC*(i-1)+1:lC*i,:)*A;
		F(lC*i+1:lC*(i+1),1:j*cB+cD) = [E(lC*(i-1)+1:lC*i,:)*B ...
            F(lC*(i-1)+1:lC*i,1:(j-1)*cB+cD)];
    end

end


function [E,F] = calcEF2(A,B,C,D,j)
% Calculate the prediction matrices E et F

    % System dimensions
    n = size(A,1);  % number of states
    nu = size(B,2);  % number of inputs
    ny = size(C,1);  % number of outputs

    A_blocks = cellfun(@(i) A^i, num2cell(1:j), 'UniformOutput', false);
    E_blocks = cellfun(@(X) C*X, A_blocks, 'UniformOutput', false);
    E = cell2mat(E_blocks(1:j)');
    F_blocks = [cellfun(@(X) C*X*B, [{eye(n)} A_blocks(1:j-1)], ...
        'UniformOutput', false) {D} {zeros(ny,nu)}];
    F = cell2mat(F_blocks(toeplitz(1:j, [1 j+1 (j+2)*ones(1,j-1)])));

end


function [E,F] = calcEF3(A,B,C,D,j)
% Calculate the prediction matrices E et F

    % System dimensions
    n = size(A,1);  % number of states
    nu = size(B,2);  % number of inputs
    ny = size(C,1);  % number of outputs

    A_blocks = cellfun(@(i) A^i, num2cell(0:j), 'UniformOutput', false);
    E_blocks = cellfun(@(X) C*X, A_blocks(2:end), 'UniformOutput', false);
    E = cell2mat(E_blocks(1:j)');
    F_blocks = [cellfun(@(X) C*X*B, A_blocks(1:j), 'UniformOutput', ...
        false) {D} {zeros(ny,nu)}];
    F = cell2mat(F_blocks(toeplitz(1:j, [1 j+1 (j+2)*ones(1,j-1)])));

end