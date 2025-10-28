function [] = ResonantCheck()
%% Physical constants (for conductance units)
e = 1.602176634e-19;      % C
h = 6.62607015e-34;       % J s
G0 = 2*e^2/h;             % conductance quantum (S) including spin
%% temperature
T_K = 4.0;                % Kelvin for finite-T conductance
kB_ev = 8.617333262145e-5; % eV/K
kBT = kB_ev * T_K;
%% setup
[H] = setupEMvsBW();

%% calc
[Elist, T_NEGF, T_BW, GammaL_eff, GammaR_eff, Delta_eff, SigmaL, SigmaR, GammaL_mat, GammaR_mat] = calcEMvsBW(H, kBT, G0);

plotEMvsBW(Elist, T_NEGF, T_BW, GammaL_eff, GammaR_eff, Delta_eff)
sweepEMvsGW(H, SigmaL, SigmaR, GammaL_mat, GammaR_mat, G0)
disp('Finished')
end

function [H] = setupEMvsBW()
%% -------------------- Model parameters --------------------
t_lead = 1.0;             % lead hopping (energy units)
epsilon_lead = 0.0;       % lead onsite energy
t_c = 0.2;                % coupling dot <-> first lead site
epsilon_dot = 0.0;        % dot onsite energy (resonant level)
N_lead_each = 30;         % number of lead sites included on each side (extended molecule)

%% -------------------- Build Extended Molecule Hamiltonian (1D chain) --------------------
N = 2*N_lead_each + 1;    % total EM sites
H = zeros(N,N);

idx_dot = N_lead_each + 1;

% Left lead sites (1..N_lead_each)
for n = 1:N_lead_each
    H(n,n) = epsilon_lead;
    if n < N_lead_each
        H(n, n+1) = -t_lead;
        H(n+1, n) = -t_lead;
    end
end

% Right lead sites (idx_dot+1 .. N)
for n = idx_dot+1:N
    H(n,n) = epsilon_lead;
end
for n = idx_dot+1:(N-1)
    H(n, n+1) = -t_lead;
    H(n+1, n) = -t_lead;
end

% Connect dot to first lead sites on each side
H(idx_dot, idx_dot-1) = -t_c;
H(idx_dot-1, idx_dot) = -t_c;
H(idx_dot, idx_dot+1) = -t_c;
H(idx_dot+1, idx_dot) = -t_c;

% Dot onsite
H(idx_dot, idx_dot) = epsilon_dot;
end

function [Elist, T_NEGF, T_BW, GammaL_eff, GammaR_eff, Delta_eff, SigmaL, SigmaR, GammaL_mat, GammaR_mat] = calcEMvsBW(H, kBT, G0)
%% Energy grid
t_lead = 1.0;             % lead hopping (energy units)
Emin = -3*t_lead;
Emax = 3*t_lead;
nE = 2001;
Elist = linspace(Emin, Emax, nE);

%% -------------------- Absorbing self-energies on outer EM sites --------------------
N_lead_each = 30;         % number of lead sites included on each side (extended molecule)
N = 2*N_lead_each + 1;    % total EM sites

SigmaL = zeros(N,N);
SigmaR = zeros(N,N);

N_absorb_layers = 4;      % number of outer-most sites per side with absorbing Sigma
left_abs_idx = 1 : min(N_absorb_layers, N_lead_each);
right_abs_idx = N - (0:(min(N_absorb_layers, N_lead_each)-1));

eta = 0.06;               % absorbing strength (positive) -> Sigma = -i*eta on outer sites
SigmaL(sub2ind([N,N], left_abs_idx, left_abs_idx)) = -1i * eta;
SigmaR(sub2ind([N,N], right_abs_idx, right_abs_idx)) = -1i * eta;

%% calc
[T_NEGF, GammaL_mat, GammaR_mat] = doEM(Elist, H, SigmaL, SigmaR);
[T_BW, GammaL_eff, GammaR_eff, Delta_eff] = doBW(Elist, H, SigmaL, SigmaR);

%% -------------------- Conductances (finite-T Landauer) --------------------
Conductance(Elist, T_NEGF, T_BW, kBT, G0)
end

function [T_NEGF, GammaL_mat, GammaR_mat] = doEM(Elist, H, SigmaL, SigmaR)
N_lead_each = 30;         % number of lead sites included on each side (extended molecule)
N = 2*N_lead_each + 1;    % total EM sites
I_N = eye(N);

% Gamma matrices for the outer absorbers (used in NEGF transmission)
GammaL_mat = 1i * (SigmaL - SigmaL');
GammaR_mat = 1i * (SigmaR - SigmaR');

T_NEGF = zeros(size(Elist));    % NEGF full transmission
%% -------------------- Energy loop: compute GR (full) & projections --------------------
for ii = 1:length(Elist)
    E = Elist(ii);
    % Full Green's function with both absorbers
    num_eta = 1e-12;          % tiny numerical broadening for matrix inversion
    GR_full = ((E + 1i*num_eta) * I_N - H - SigmaL - SigmaR) \ I_N;
    GA_full = GR_full';
    
    % NEGF transmission (matrix formula)
    T_NEGF(ii) = real(trace(GammaL_mat * GR_full * GammaR_mat * GA_full));
    if T_NEGF(ii) < 0 && T_NEGF(ii) > -1e-14, T_NEGF(ii) = 0; end
end
end

function [T_BW, GammaL_eff, GammaR_eff, Delta_eff] = doBW(Elist, H, SigmaL, SigmaR)
N_lead_each = 30;         % number of lead sites included on each side (extended molecule)
N = 2*N_lead_each + 1;    % total EM sites
I_N = eye(N);

%% -------------------- Preallocations --------------------
T_BW   = zeros(size(Elist));    % Breit-Wigner using EM-extracted Gammas
GammaL_eff = zeros(size(Elist));
GammaR_eff = zeros(size(Elist));
Delta_eff = zeros(size(Elist)); % real level shift
A_dot = zeros(size(Elist));     % spectral function at dot from full GR (A = -2 Im G_dd)

%% -------------------- Energy loop: compute GR (full) & projections --------------------    
for ii = 1:length(Elist)
    E = Elist(ii);
    % Full Green's function with both absorbers
    num_eta = 1e-12;          % tiny numerical broadening for matrix inversion
    GR_full = ((E + 1i*num_eta) * I_N - H - SigmaL - SigmaR) \ I_N;

    % Dot Green's function element (dd)
    idx_dot = N_lead_each + 1;
    Gdd_full = GR_full(idx_dot, idx_dot);
    A_dot(ii) = -2 * imag(Gdd_full);
    
    epsilon_dot = 0.0;        % dot onsite energy (resonant level)

    left = true;
    if left == true
    % Now compute dot Green's function for left-only absorber (SigmaR = 0)
    GR_Lonly = ((E + 1i*num_eta) * I_N - H - SigmaL) \ I_N;
    Gdd_Lonly = GR_Lonly(idx_dot, idx_dot);
    % effective self-energy due to left lead as seen by dot:
    Sigma_eff_L = E - epsilon_dot - 1.0 / Gdd_Lonly;
    % Extract Gammas and shift (energy dependent)
    GammaL_eff(ii) = -2 * imag(Sigma_eff_L);
    end
    
    right = true;
    if right == true
    % Now compute dot Green's function for right-only absorber (SigmaL = 0)
    GR_Ronly = ((E + 1i*num_eta) * I_N - H - SigmaR) \ I_N;
    Gdd_Ronly = GR_Ronly(idx_dot, idx_dot);
    % effective self-energy due to left lead as seen by dot:
    Sigma_eff_R = E - epsilon_dot - 1.0 / Gdd_Ronly;
    % Extract Gammas and shift (energy dependent)
    GammaR_eff(ii) = -2 * imag(Sigma_eff_R);
    end
    
    % Effective total self-energy from both sides (optionally we can also do Sigma_eff_total = E - eps - 1/Gdd_full)
    Sigma_eff_total = Sigma_eff_L + Sigma_eff_R;  % should be close to E - eps - 1/Gdd_full
    
    % Extract Gammas (energy dependent)
    Gamma_tot = GammaL_eff(ii) + GammaR_eff(ii);

    % Extract shift (energy dependent)
    Delta_eff(ii) = real(Sigma_eff_total);  % level shift = Re Sigma_eff_total
    
    % Breit-Wigner transmission using extracted Gammas and shift
    numerator = GammaL_eff(ii) * GammaR_eff(ii);
    denom = (E - epsilon_dot - Delta_eff(ii))^2 + (Gamma_tot/2)^2;
    if denom <= 0
        T_BW(ii) = 0;
    else
        T_BW(ii) = numerator / denom;
    end
end
end

function [] = Conductance(Elist, T_NEGF, T_BW, kBT, G0)
%% -------------------- Conductances (finite-T Landauer) --------------------
% Fermi derivative
if kBT > 0
    EF = 0.0;
    dfdE = -1 ./ (4 * kBT * (cosh((Elist - EF) ./ (2*kBT))).^2);
else
    % approximate delta -> pick value of T at EF
    dfdE = zeros(size(Elist));
    % put a narrow spike near EF
    EF = 0.0;
    [~, idxEF] = min(abs(Elist - EF));
    dfdE(idxEF) = -1;
end

G_NEGF = G0 * trapz(Elist, -dfdE .* T_NEGF);
G_BW   = G0 * trapz(Elist, -dfdE .* T_BW);

T_K = 4.0;                % Kelvin for finite-T conductance
fprintf('Finite-T linear conductance at T=%.2f K (EF=%.3g):\n', T_K, EF);
fprintf('  G (NEGF, EM) = %.6g S  (%.6g G0)\n', G_NEGF, G_NEGF / G0);
fprintf('  G (Breit-Wigner using EM Gammas) = %.6g S  (%.6g G0)\n\n', G_BW, G_BW / G0);
end

function [] = plotEMvsBW(Elist, T_NEGF, T_BW, GammaL_eff, GammaR_eff, Delta_eff)
%% -------------------- Plots --------------------
figure('Name','Transmission comparison','NumberTitle','off','Position',[100 100 900 700]);
subplot(3,1,1);
plot(Elist, T_NEGF, 'LineWidth', 1.5); hold on;
plot(Elist, T_BW, '--', 'LineWidth', 1.4);
xlabel('Energy (units of t)');
ylabel('T(E)');
legend('T_{NEGF}','T_{BW (EM)}','Location','Best');
title('Transmission: full NEGF vs Breit–Wigner using EM-extracted \Gamma_{L,R}');
grid on;

subplot(3,1,2);
plot(Elist, GammaL_eff, 'LineWidth', 1.2); hold on;
plot(Elist, GammaR_eff, 'LineWidth', 1.2);
plot(Elist, GammaL_eff + GammaR_eff, ':k','LineWidth',1);
xlabel('Energy');
ylabel('\Gamma(E)');
legend('\Gamma_L','\Gamma_R','\Gamma_{tot}','Location','Best');
title('Effective broadenings (projected to dot) computed from EM');
grid on;

subplot(3,1,3);
plot(Elist, Delta_eff, 'LineWidth', 1.2);
xlabel('Energy');
ylabel('\Delta(E) = Re[\Sigma^{eff}_L + \Sigma^{eff}_R]');
title('Effective level shift');
grid on;

N_lead_each = 30;         % number of lead sites included on each side (extended molecule)
eta = 0.06;               % absorbing strength (positive) -> Sigma = -i*eta on outer sites
N_absorb_layers = 4;      % number of outer-most sites per side with absorbing Sigma
sgtitle(sprintf('EM size: %d sites each lead, eta = %.3g, N_absorb = %d', N_lead_each, eta, N_absorb_layers));
end

function [] = sweepEMvsGW(H, SigmaL, SigmaR, GammaL_mat, GammaR_mat, G0)
%% -------------------- Gate sweep (zero-T approx using T(EF)) --------------------
doGateSweep = true;
if doGateSweep
    neps = 201;
    epsilon_dot = 0.0;        % dot onsite energy (resonant level)
    eps_range = linspace(epsilon_dot - 2.5, epsilon_dot + 2.5, neps);
    G_vs_eps_NEGF = zeros(size(eps_range));
    G_vs_eps_BW = zeros(size(eps_range));
    % For speed, evaluate approximate zero-T conductance by using T(EF) with shifted eps
    for k = 1:neps
        eps_tmp = eps_range(k);
        % Update dot onsite in H and recompute dot GR at EF for NEGF (full)
        N_lead_each = 30;         % number of lead sites included on each side (extended molecule)
        idx_dot = N_lead_each + 1;
        H_tmp = H;
        H_tmp(idx_dot, idx_dot) = eps_tmp;

        N = 2*N_lead_each + 1;    % total EM sites
        I_N = eye(N);
        EF = 0.0;
        num_eta = 1e-12;          % tiny numerical broadening for matrix inversion
        GR_full_tmp = ((EF + 1i*num_eta) * I_N - H_tmp - SigmaL - SigmaR) \ I_N;
        GA_full_tmp = GR_full_tmp';
        T_NEGF_EF = real(trace(GammaL_mat * GR_full_tmp * GammaR_mat * GA_full_tmp));
        if T_NEGF_EF < 0 && T_NEGF_EF > -1e-14, T_NEGF_EF = 0; end
        G_vs_eps_NEGF(k) = G0 * T_NEGF_EF;
        
        % For BW: compute GR_dd with left-only and right-only to get Gammas at EF
        GR_L_tmp = ((EF + 1i*num_eta) * I_N - H_tmp - SigmaL) \ I_N;
        Gdd_L = GR_L_tmp(idx_dot, idx_dot);
        Sigma_eff_L_tmp = EF - eps_tmp - 1.0 / Gdd_L;
        GR_R_tmp = ((EF + 1i*num_eta) * I_N - H_tmp - SigmaR) \ I_N;
        Gdd_R = GR_R_tmp(idx_dot, idx_dot);
        Sigma_eff_R_tmp = EF - eps_tmp - 1.0 / Gdd_R;
        GammaL_tmp = -2 * imag(Sigma_eff_L_tmp);
        GammaR_tmp = -2 * imag(Sigma_eff_R_tmp);
        Delta_tmp = real(Sigma_eff_L_tmp + Sigma_eff_R_tmp);
        denom = (EF - eps_tmp - Delta_tmp)^2 + ( (GammaL_tmp + GammaR_tmp)/2 )^2;
        if denom <= 0
            Tbw_tmp = 0;
        else
            Tbw_tmp = (GammaL_tmp * GammaR_tmp) / denom;
        end
        G_vs_eps_BW(k) = G0 * Tbw_tmp;
    end
    
    figure('Name','Gate sweep (zero-T approx)','NumberTitle','off');
    plot(eps_range, G_vs_eps_NEGF / G0, 'LineWidth', 1.4); hold on;
    plot(eps_range, G_vs_eps_BW / G0, '--', 'LineWidth', 1.4);
    xlabel('\epsilon_d');
    ylabel('G / G_0');
    legend('G_{NEGF}','G_{BW}','Location','Best');
    title('Gate sweep: NEGF (EM) vs Breit–Wigner using EM-extracted \Gamma');
    grid on;
end

%% -------------------- Final suggestions & checks --------------------
disp('--- Suggested convergence checks ---');
disp('- Increase N_lead_each until T_NEGF and T_BW (and G) stop changing significantly.');
disp('- Sweep eta across a few decades (e.g. 1e-3 .. 1e-0) and find an absorbing plateau.');
disp('- Ensure N_absorb_layers is large enough to prevent reflections from the outer boundary.');
disp('- If you need energy-independent Gamma, evaluate GammaL_eff/ GammaR_eff at E=EF and use constants in BW formula.');
end

function [] = sweepEM()
    Nd = 121;
    epsilons = linspace(-1.5, 1.5, Nd);
    G_vs_eps = zeros(size(epsilons));
    % reuse T_of_E calculation but easier: compute zero-T conductance approx by T at E=EF
    % We'll compute T(EF) for each gate (shift dot onsite)
    for k = 1:length(epsilons)
        H_tmp = H;
        H_tmp(idx_dot, idx_dot) = epsilons(k);
        % compute GR at EF
        GR_tmp = ((EF + 1i*num_eta) * I_N - H_tmp - SigmaL - SigmaR) \ I_N;
        GA_tmp = GR_tmp';
        GammaL = 1i * (SigmaL - SigmaL');
        GammaR = 1i * (SigmaR - SigmaR');
        T0 = real(trace(GammaL * GR_tmp * GammaR * GA_tmp));
        if T0 < 0 && T0 > -1e-12, T0 = 0; end
        G_vs_eps(k) = G0 * T0;  % zero-T linear conductance in S (approx)
    end
    
    figure;
    plot(epsilons, G_vs_eps / G0, 'LineWidth', 1.5);
    xlabel('\epsilon_{d}');
    ylabel('G / G_0');
    title('Zero-T (approx) Conductance vs Dot Level (gate sweep)');
    grid on;
end

function [] = sweepBW()
    eps_range = linspace(eps_d - 2*Gamma - 1.0, eps_d + 2*Gamma + 1.0, 401);
    G_vs_eps = zeros(size(eps_range));
    for ii = 1:length(eps_range)
        eps_tmp = eps_range(ii);
        T_tmp = (GammaL * GammaR) ./ ( (Gamma/2).^2 + (EF - eps_tmp).^2 );
        if useFiniteT
            % finite T convolution but since T(E) is analytic we evaluate T(E) for existing E grid shifted by eps_tmp
            T_shifted = (GammaL .* GammaR) ./ ( (Gamma/2).^2 + (E - eps_tmp).^2 );
            integrand = -dfdE .* T_shifted;
            G_vs_eps(ii) = G0 * trapz(E, integrand);
        else
            G_vs_eps(ii) = G0 * T_tmp;
        end
    end

    figure('Name','Gate sweep: G vs \epsilon_d','NumberTitle','off');
    plot(eps_range, G_vs_eps / G0, 'LineWidth', 1.5);
    xlabel('\epsilon_d (gate)');
    ylabel('G / G_0');
    title('Conductance vs dot level position');
    grid on;
end