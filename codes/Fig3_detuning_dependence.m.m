%% Compare two solution families vs detuning (updated parameters)
clear; tic

% ---------------- Common parameters ----------------
kA    = 1.0;
g_ans = 0.15;  
Fa    = 0.05;  
N1 = 4; 
N2 = 4;

% ---------------- Operators ----------------
a  = tensor(destroy(N1), identity(N2));
b  = tensor(identity(N1), destroy(N2));
ad = a';  bd = b';
aDa = ad*a;  bDb = bd*b;  ad2a2 = ad*ad*a*a;

Sa = 0.5*(2*spre(a)*spost(ad) - spre(aDa) - spost(aDa));
Sb = 0.5*(2*spre(b)*spost(bd) - spre(bDb) - spost(bDb));

% ---------------- Two solution families ----------------
fam = struct();

% Family 1: Delta = 0, kB = kA, Fb = -Fa^2*kB/(g*kA)
fam(1).kB    = kA;
fam(1).Fb    = -Fa^2*fam(1).kB/(g_ans*kA);
fam(1).label = '$\kappa_B=\kappa_A$';
fam(1).color = [0.10 0.30 0.85];
fam(1).style = '--';

% Family 2: kB = 2*kA, Fb = -2*Fa^2/g, Delta arbitrary
fam(2).kB    = 2*kA;
fam(2).Fb    = -2*Fa^2/g_ans;
fam(2).label = '$\kappa_B=2\kappa_A$';
fam(2).color = [0.85 0.10 0.10];
fam(2).style = '-';

fprintf('Family 1: kB=%.2f, Fb=%.5f\n', fam(1).kB, fam(1).Fb);
fprintf('Family 2: kB=%.2f, Fb=%.5f\n', fam(2).kB, fam(2).Fb);

% ---------------- Detuning scan ----------------
M    = 201;
Dset = linspace(-3.0, 3.0, M);

g2_data = zeros(2, M);
na_data = zeros(2, M);

for ic = 1:2
    kB = fam(ic).kB;
    Fb = fam(ic).Fb;
    for kd = 1:M
        Delta = Dset(kd);
        H  = Delta*aDa + 2*Delta*bDb + g_ans*(ad*ad*b + bd*a*a) ...
             + Fa*(ad + a) + Fb*(bd + b);
        LH = -1i*(spre(H) - spost(H));
        L  = LH + kA*Sa + kB*Sb;

        rhoss = steady(L);   % tighter tolerance for weak drive

        na = real(trace(aDa*rhoss));
        g2_data(ic, kd) = real(trace(ad2a2*rhoss)) / max(na^2, 1e-300);
        na_data(ic, kd) = na;
    end
    fprintf('Family %d done.\n', ic);
end

% ---------------- Plot ----------------
figure('Position',[100 100 1100 460])

% (a) log10 g2(0) vs Delta
subplot(1,2,1); hold on
for ic = 1:2
    plot(Dset, log10(max(g2_data(ic,:),1e-300)), fam(ic).style, ...
        'Color', fam(ic).color, 'LineWidth', 2.0, 'DisplayName', fam(ic).label);
end

xlabel('$\Delta/\kappa_A$','Interpreter','latex','FontSize',16)
ylabel('$\lg\,g^{(2)}(0)$','Interpreter','latex','FontSize',16)
title('(a) Photon antibunching','Interpreter','latex','FontSize',16)
legend('Location','best','Interpreter','latex','FontSize',11)
box on

% ---- (b) log10 <n_a> vs Delta: solid line first, then dashed ----
subplot(1,2,2); hold on

h_b = gobjects(2,1);                 % store handles indexed by family
order = [2 1];                       % draw Family 2 (solid) first, then Family 1 (dashed)
for k = 1:numel(order)
    ic = order(k);
    h_b(ic) = plot(Dset, log10(max(na_data(ic,:),1e-300)), fam(ic).style, ...
        'Color', fam(ic).color, 'LineWidth', 2.0);
end

xlabel('$\Delta/\kappa_A$','Interpreter','latex','FontSize',16)
ylabel('$\lg\,\langle n_a\rangle$','Interpreter','latex','FontSize',16)
title('(b) Mean photon number','Interpreter','latex','FontSize',16)
legend(h_b, {fam(1).label, fam(2).label}, ...
       'Location','best','Interpreter','latex','FontSize',11)
box on

sgtitle('Comparison of two solution families','Interpreter','latex','FontSize',16)

save('Fig_two_families_g015_Fa005.mat')
print(gcf, 'Fig_two_families_g015_Fa005.png', '-dpng', '-r300')
toc