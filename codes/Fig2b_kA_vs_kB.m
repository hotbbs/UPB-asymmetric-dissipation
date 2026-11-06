%% Optimal ratio kappaB = 2*kappaA: (kappaA, kappaB) plane
clear; tic

% ---------------- System parameters ----------------
Delta = 2.0;
g_ans = 0.15;          % fixed nonlinear coupling
Fa    = 0.05;          % fixed drive of mode a
Fb    = -2*Fa*Fa/g_ans;         % fixed drive of mode b
N1    = 4;  N2 = 4;

% ---------------- Operators ----------------
a  = tensor(destroy(N1), identity(N2));
b  = tensor(identity(N1), destroy(N2));
ad = a';  bd = b';

aDa   = ad*a;
bDb   = bd*b;
ad2a2 = ad*ad*a*a;

% ---------------- Hamiltonian (independent of kappaA, kappaB) -----
H = Delta*aDa + 2*Delta*bDb + g_ans*(ad*ad*b + bd*a*a) ...
    + Fa*(ad + a) + Fb*(bd + b);
LH = -1i*(spre(H) - spost(H));

% Dissipation superoperators
Sa = 0.5*(2*spre(a)*spost(ad) - spre(aDa) - spost(aDa));
Sb = 0.5*(2*spre(b)*spost(bd) - spre(bDb) - spost(bDb));

% ---------------- 2D scan: x = kappaA, y = kappaB ----------------
M     = 6*(100+1);
kAset = linspace(0.01, 10.0, M);   % x axis
kBset = linspace(0.01, 10.0, M);   % y axis

g02 = zeros(M, M);               % rows -> kB (y), cols -> kA (x)

dq = parallel.pool.DataQueue;
afterEach(dq, @(k) fprintf('Progress: %5.2f%%\n', 100*k/M));

for kj = 1:M                  % outer loop: y axis (kB)
    kB  = kBset(kj);
    row = zeros(1, M);
    parfor ki = 1:M                 % inner loop: x axis (kA)
        kA = kAset(ki);
        L  = LH + kA*Sa + kB*Sb;
        rhoss = steady(L);

        na      = real(trace(aDa*rhoss));
        row(ki) = real(trace(ad2a2*rhoss)) / max(na^2, 1e-30);
    end
    g02(kj, :) = row;
    send(dq, kj);
end

% ---------------- Plot ----------------
figure
imagesc(kAset, kBset, log10(abs(g02)))   % x = kA, y = kB
axis xy
xlabel('$\kappa_A$','Interpreter','latex','FontSize',20)
ylabel('$\kappa_B$','Interpreter','latex','FontSize',20)
title('$\lg[g^{(2)}(0)]$','Interpreter','latex','FontSize',20)
colormap('jet'); colorbar; hold on

% Optimal line: kappaB = 2*kappaA
kA_line = linspace(min(kAset), max(kAset), 400);
kB_line = 2*kA_line;
mask = (kB_line >= min(kBset)) & (kB_line <= max(kBset)) & ...
       (kA_line >= min(kAset)) & (kA_line <= max(kAset));
plot(kA_line(mask), kB_line(mask), 'w--', 'LineWidth', 2.5)

save('Fig_kB_vs_kA.mat')
print(gcf, 'Fig_kB_vs_kA.png', '-dpng', '-r300')
toc