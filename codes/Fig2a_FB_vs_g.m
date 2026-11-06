%% g vs Fb scan (Fa fixed)
clear; tic

% ---------------- System parameters ----------------
Delta = 0.0;
kkA   = 1.0;
kkB   = 1.0;

Fa    = 0.05;         % fixed drive amplitude of mode a (key parameter)
N1    = 4;  N2 = 4;  % photon-number truncation for modes a and b

% ---------------- Operators ----------------
a  = tensor(destroy(N1), identity(N2));
b  = tensor(identity(N1), destroy(N2));
ad = a';   bd = b';

aDa   = ad*a;
bDb   = bd*b;
ad2a2 = ad*ad*a*a;   % a^dagger^2 a^2

% ---------------- Superoperators (built once) ----------------
% Terms independent of Fb and g are precomputed to speed up the scan
La = 0.5*kkA*(2*spre(a)*spost(ad) - spre(aDa) - spost(aDa));
Lb = 0.5*kkB*(2*spre(b)*spost(bd) - spre(bDb) - spost(bDb));

Hbase = Delta*aDa + 2*Delta*bDb + Fa*(ad + a);   % contains fixed Fa drive
Hnl   = ad*ad*b + bd*a*a;                        % nonlinear coupling term
Hbd   = bd + b;

Lbase = -1i*(spre(Hbase) - spost(Hbase)) + La + Lb;
LFb   = -1i*(spre(Hbd)   - spost(Hbd));
Lg    = -1i*(spre(Hnl)   - spost(Hnl));

% ---------------- Scan grid ----------------
M     = 6*(100+1);
gset  = linspace(-0.5, 0.5, M);   % y axis: g
Fbset = linspace(-0.5, 0.5, M);   % x axis: Fb

g02 = zeros(M, M);

% Progress reporting that is safe under parfor
dq = parallel.pool.DataQueue;
afterEach(dq, @(k) fprintf('Progress: %5.2f%%  (%d/%d)\n', 100*k/M, k, M));

% ---------------- Main loop: outer g, inner Fb ----------------
for kg = 1:M
    g   = gset(kg);
    row = zeros(1, M);
    parfor kf = 1:M
        Fb = Fbset(kf);
        L  = Lbase + Fb*LFb + g*Lg;
        rhoss = steady(L);

        nb      = real(trace(aDa*rhoss));
        row(kf) = real(trace(ad2a2*rhoss)) / max(nb^2, 1e-30);
    end
    g02(kg, :) = row;
    send(dq, kg);
end

% ---------------- Plot ----------------
figure
imagesc(Fbset, gset, log10(abs(g02)))
axis xy
xlabel('$F_B$','Interpreter','latex','FontSize',20)
ylabel('$g$','Interpreter','latex','FontSize',20)
title('$\lg[g^{(2)}(0)]$','Interpreter','latex','FontSize',20)
colormap('jet'); colorbar; hold on

% ---- Analytic optimal condition: g = -Fa^2 * kB / (Fb * kA) ----
Fb_ana = linspace(min(Fbset), max(Fbset), 2000);   % 点多一点，曲线更平滑
g_ana  = -(Fa^2 * kkB) ./ (Fb_ana * kkA);

% Set out-of-window and near-singular points to NaN (do NOT delete them!)
g_ana(g_ana < min(gset) | g_ana > max(gset)) = NaN;
g_ana(abs(Fb_ana) < 1e-6)                    = NaN;   % guard the Fb = 0 singularity

plot(Fb_ana, g_ana, 'w--', 'LineWidth', 2.5)

save('Fig_g_vs_Fb_Fa0.1.mat')
print(gcf, 'Fig_g_vs_Fb_Fa0.1.png', '-dpng', '-r300')
toc