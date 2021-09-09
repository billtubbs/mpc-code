%% Speed tests

n = 14; ny = 3; nu = 6;
sys = drss(n,ny,nu);
A = sys.A;
B = sys.B;
C = sys.C;
D = sys.D;

Hp = 30; Hc = 2;

%% Calculation of E matrix

% Method 1 - fastest accurate method
tic
E1 = zeros(Hp*ny, n);
for i=1:Hp
    E1(1+ny*(i-1):ny*i, :) = C*A^i;
end
toc

% Method 2 - fastest but numerically less accurate
tic
E2 = zeros(Hp*ny, n);
X = C*A;
for i=1:Hp
    E2(1+ny*(i-1):ny*i, :) = X;
    X = X*A;
end
toc

% Method 3 - no for loop but slower
tic
blocks = [eye(n) cellfun(@(i) C*A^i, num2cell(1:Hp), 'UniformOutput', false) zeros(n)];
E3 = cell2mat(blocks(2:Hp+1)');
toc

assert(all(abs(E1 - E2) < 1e-14, [2 1]))
assert(all(abs(E1 - E3) < 1e-14, [2 1]))


%% Calculation of both matrices using cellfun and cell2mat

% André's function from GlobPC
fprintf("   calcef: ")
tic
[E1,F1] = calcef(A,B,C,D,Hp);
toc

% Note this method assumes Hc=Hp
% size(F1) : [ ny*Hp nu*Hp]
% Convert to same output as pred_mats
for i = 1:nu
	F1(:, (Hc-1)*nu + i) = sum(F1(:, (Hc-1:Hp)*nu+i), 2);
end
% Drop blocks in cols Hc+1:Hp+1
F1 = F1(:, 1:Hc*nu);

% My code
fprintf("pred_mats: ")
tic
[E2,F2] = pred_mats(A,B,C,D,Hp,Hc);
toc
% size(F2) : [ ny*Hp nu*Hc]

assert(all(abs(E1 - E2) < 1e-14, [2 1]))
assert(all(abs(F1 - F2) < 1e-12, [2 1]))

