function [] = ResonantCheck()
N_lead_each = 30;         % number of lead sites included on each side (extended molecule)
N = 2*N_lead_each + 1;    % total EM sites

N_absorb_layers = 4;      % number of outer-most sites per side with absorbing Sigma
left_abs_idx = 1 : min(N_absorb_layers, N_lead_each);
right_abs_idx = N - (0:(min(N_absorb_layers, N_lead_each)-1));

% Energy grid
t_lead = 1.0;             % lead hopping (energy units)
Emin = -3*t_lead;
Emax = 3*t_lead;
nE = 2001;
Elist = linspace(Emin, Emax, nE);

% -------------------- Absorbing self-energies on outer EM sites --------------------
eta = 0.06;               % absorbing strength (positive) -> Sigma = -i*eta on outer sites

SigmaL = zeros(N,N);
SigmaL(sub2ind([N,N], left_abs_idx, left_abs_idx)) = -1i * eta;

SigmaR = zeros(N,N);
SigmaR(sub2ind([N,N], right_abs_idx, right_abs_idx)) = -1i * eta;

%% setup
[H] = setupH();

%% calc
[T_NEGF] = doEM(Elist, H, SigmaL, SigmaR);
[T_BW, GammaL_eff, GammaR_eff, Delta_eff] = doBW(Elist, H, SigmaL, SigmaR);

%%
sweepEMvsGW(H, SigmaL, SigmaR)

plotEMvsBW(Elist, T_NEGF, T_BW, GammaL_eff, GammaR_eff, Delta_eff)

if false
% -------------------- Conductances (finite-T Landauer) --------------------
Conductance(Elist, T_NEGF, T_BW)
end

%% -------------------- Final suggestions & checks --------------------
disp('--- Suggested convergence checks ---');
disp('- Increase N_lead_each until T_NEGF and T_BW (and G) stop changing significantly.');
disp('- Sweep eta across a few decades (e.g. 1e-3 .. 1e-0) and find an absorbing plateau.');
disp('- Ensure N_absorb_layers is large enough to prevent reflections from the outer boundary.');
disp('- If you need energy-independent Gamma, evaluate GammaL_eff/ GammaR_eff at E=EF and use constants in BW formula.');
end

function [H] = setupH()
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

function [] = sweepEMvsGW(H, SigmaL, SigmaR)
N_lead_each = 30;         % number of lead sites included on each side (extended molecule)
idx_dot = N_lead_each + 1;

% -------------------- Gate sweep (zero-T approx using T(EF)) --------------------
    neps = 201;
    epsilon_dot = 0.0;        % dot onsite energy (resonant level)
    eps_range = linspace(epsilon_dot - 2.5, epsilon_dot + 2.5, neps);
    % For speed, evaluate approximate zero-T conductance by using T(EF) with shifted eps

    [G_NEGF] = sweepEM(eps_range, H, SigmaL, SigmaR, idx_dot);

    [G_BW] = sweepBW(eps_range, H, SigmaL, SigmaR, idx_dot);

    figure('Name','Gate sweep (zero-T approx)','NumberTitle','off');
    plot(eps_range, G_NEGF);
    hold on;
    plot(eps_range, G_BW, '--');
    xlabel('\epsilon_d');
    ylabel('G / G_0');
    legend('G_{NEGF}','G_{BW}','Location','Best');
    title('Gate sweep: NEGF (EM) vs Breit–Wigner using EM-extracted \Gamma');
    grid on;
end

%% -------------------- Other functions: EM -------------------- 
function [T_NEGF] = doEM(Elist, H, SigmaL, SigmaR)
num_eta = 1e-12;          % tiny numerical broadening for matrix inversion

% Gamma matrices for the outer absorbers (used in NEGF transmission)
GammaL_mat = 1i * (SigmaL - SigmaL');
GammaR_mat = 1i * (SigmaR - SigmaR');

% -------------------- Energy loop: compute GR (full) & projections --------------------
T_NEGF = zeros(size(Elist));    % NEGF full transmission
I_N = eye(size(H));
for ii = 1:length(Elist)
    E = Elist(ii);

    % Full Green's function with both absorbers
    GR_full = ((E + 1i*num_eta) * I_N - H - SigmaL - SigmaR) \ I_N;
    GA_full = GR_full';
    
    % NEGF transmission (matrix formula)
    T_NEGF(ii) = real(trace(GammaL_mat * GR_full * GammaR_mat * GA_full));
end
end

function [G_NEGF] = sweepEM(eps_range, H, SigmaL, SigmaR, idx_dot)
num_eta = 1e-12;          % tiny numerical broadening for matrix inversion
EF = 0.0;

% Gamma matrices for the outer absorbers (used in NEGF transmission)
GammaL_mat = 1i * (SigmaL - SigmaL');
GammaR_mat = 1i * (SigmaR - SigmaR');

% reuse T_of_E calculation but easier: compute zero-T conductance approx by T at E=EF
% We'll compute T(EF) for each gate (shift dot onsite)
G_NEGF = zeros(size(eps_range));
I_N = eye(size(H));
for k = 1:length(eps_range)
    eps_tmp = eps_range(k);
    % Update dot onsite in H and recompute dot GR at EF for NEGF (full)
    H_tmp = H;
    H_tmp(idx_dot, idx_dot) = eps_tmp;

    % compute GR at EF
    GR_full_tmp = ((EF + 1i*num_eta) * I_N - H_tmp - SigmaL - SigmaR) \ I_N;
    GA_full_tmp = GR_full_tmp';

    G_NEGF(k) = real(trace(GammaL_mat * GR_full_tmp * GammaR_mat * GA_full_tmp));
    % zero-T linear conductance in S (approx)
end
end

%% -------------------- Other functions: BW -------------------- 
function [T_BW, GammaL_eff, GammaR_eff, Delta_eff] = doBW(Elist, H, SigmaL, SigmaR)
% -------------------- Energy loop: compute GR (full) & projections --------------------
N_lead_each = 30;         % number of lead sites included on each side (extended molecule)
idx_dot = N_lead_each + 1;

[SigmaL_eff, GammaL_eff] = doBWside(Elist, H, SigmaL, idx_dot);

[SigmaR_eff, GammaR_eff] = doBWside(Elist, H, SigmaR, idx_dot);

[T_BW, Delta_eff] = doBWboth(Elist, SigmaL_eff, SigmaR_eff, GammaL_eff, GammaR_eff);

if false
[A_dot, Gdd_full] = doBWfull(Elist, H, SigmaL, SigmaR);
end
end

function [G_BW] = sweepBW(eps_range, H, SigmaL, SigmaR, idx_dot)
    % reuse T_of_E calculation but easier: compute zero-T conductance approx by T at E=EF
    % We'll compute T(EF) for each gate (shift dot onsite)

    % For BW: compute GR_dd with left-only to get GammaL at EF
    [SigmaL_tmp, GammaL] = sweepBWside(eps_range, H, SigmaL, idx_dot);

    % For BW: compute GR_dd with right-only to get GammaR at EF
    [SigmaR_tmp, GammaR] = sweepBWside(eps_range, H, SigmaR, idx_dot);
    
    [G_BW] = sweepBWboth(eps_range, SigmaL_tmp, SigmaR_tmp, GammaL, GammaR);
end

function [Sigma_eff, Gamma_eff] = doBWside(Elist, H, Sigma, idx_dot)
num_eta = 1e-12;          % tiny numerical broadening for matrix inversion
epsilon_dot = 0.0;        % dot onsite energy (resonant level)

I_N = eye(size(H));
Sigma_eff = zeros(size(Elist));
Gamma_eff = zeros(size(Elist));
for ii = 1:length(Elist)
    E = Elist(ii);

    % Now compute dot Green's function for left-only absorber (SigmaR = 0)
    GR = ((E + 1i*num_eta) * I_N - H - Sigma) \ I_N;
    Gdd = GR(idx_dot, idx_dot);

    % effective self-energy due to left lead as seen by dot
    Sigma_eff(ii) = E - epsilon_dot - 1.0 / Gdd;

    % Extract Gammas and shift (energy dependent)
    Gamma_eff(ii) = -2 * imag(Sigma_eff(ii));
end
end

function [Sigma_tmp, Gamma_tmp] = sweepBWside(eps_range, H, Sigma, idx_dot)
num_eta = 1e-12;          % tiny numerical broadening for matrix inversion
EF = 0.0;

I_N = eye(size(H));
Sigma_tmp = zeros(size(eps_range));
Gamma_tmp = zeros(size(eps_range));
for k = 1:length(eps_range)
    eps_tmp = eps_range(k);
        
    % Update dot onsite in H and recompute dot GR at EF for NEGF (full)
    H(idx_dot, idx_dot) = eps_tmp;
        
    % For BW: compute GR_dd with left-only to get GammaL at EF
    GR = ((EF + 1i*num_eta) * I_N - H - Sigma) \ I_N;
    Gdd = GR(idx_dot, idx_dot);
        
    % effective self-energy due to left lead as seen by dot
    Sigma_tmp(k) = EF - eps_tmp - 1.0 / Gdd;

    % Extract Gammas and shift (energy dependent)
    Gamma_tmp(k) = -2 * imag(Sigma_tmp(k));
end
end

function [T_BW, Delta] = doBWboth(Elist, SigmaL, SigmaR, GammaL, GammaR)
epsilon_dot = 0.0;        % dot onsite energy (resonant level)
    
T_BW   = zeros(size(Elist));    % Breit-Wigner using EM-extracted Gammas
Delta = zeros(size(Elist)); % real level shift
for ii = 1:length(Elist)
    E = Elist(ii);
    % Effective total self-energy from both sides (optionally we can also do Sigma_eff_total = E - eps - 1/Gdd_full)
    Sigma = SigmaL(ii) + SigmaR(ii);  % should be close to E - eps - 1/Gdd_full
    % Extract Gammas (energy dependent)
    Gamma = GammaL(ii) + GammaR(ii);
    % Extract shift (energy dependent)
    Delta(ii) = real(Sigma);  % level shift = Re Sigma_eff_total
    
    % Breit-Wigner transmission using extracted Gammas and shift
    T_BW(ii) = (GammaL(ii) * GammaR(ii)) / ((E - epsilon_dot - Delta(ii))^2 + (Gamma/2)^2);
end
end

function [G_BW] = sweepBWboth(eps_range, SigmaL, SigmaR, GammaL, GammaR)
EF = 0.0;
    
G_BW = zeros(size(eps_range));
Delta = zeros(size(eps_range));
for k = 1:length(eps_range)
    eps_tmp = eps_range(k);
    % Effective total self-energy from both sides
    Sigma = SigmaL(k) + SigmaR(k);
    % Extract Gammas (energy dependent)
    Gamma = GammaL(k) + GammaR(k);
    % Extract shift (energy dependent)
    Delta(k) = real(Sigma);

    % Breit-Wigner conductance using extracted Gammas and shift
    G_BW(k) = (GammaL(k) * GammaR(k)) / ((EF - eps_tmp - Delta(k))^2 + (Gamma/2)^2);
end
end

%% -------------------- Plotting functions -------------------- 
function [] = plotEMvsBW(Elist, T_NEGF, T_BW, GammaL_eff, GammaR_eff, Delta_eff)
%% -------------------- Plots --------------------
figure('Name','Transmission comparison','NumberTitle','off','Position',[100 100 900 700]);

% 1st subplot
subplot(3,1,1);
plot(Elist, T_NEGF);
hold on;
plot(Elist, T_BW, '--');
xlabel('Energy (units of t)');
ylabel('T(E)');
legend('T_{NEGF}','T_{BW (EM)}','Location','Best');
title('Transmission: full NEGF vs Breit–Wigner using EM-extracted \Gamma_{L,R}');
grid on;

% 2nd subplot
subplot(3,1,2);
plot(Elist, GammaL_eff);
hold on;
plot(Elist, GammaR_eff);
plot(Elist, GammaL_eff + GammaR_eff, ':k');
xlabel('Energy');
ylabel('\Gamma(E)');
legend('\Gamma_L','\Gamma_R','\Gamma_{tot}','Location','Best');
title('Effective broadenings (projected to dot) computed from EM');
grid on;

% 3rd subplot
subplot(3,1,3);
plot(Elist, Delta_eff);
xlabel('Energy');
ylabel('\Delta(E) = Re[\Sigma^{eff}_L + \Sigma^{eff}_R]');
title('Effective level shift');
grid on;

% Title
N_lead_each = 30;         % number of lead sites included on each side (extended molecule)
eta = 0.06;               % absorbing strength (positive) -> Sigma = -i*eta on outer sites
N_absorb_layers = 4;      % number of outer-most sites per side with absorbing Sigma
sgtitle(sprintf('EM size: %d sites each lead, eta = %.3g, N_absorb = %d', N_lead_each, eta, N_absorb_layers));
end

%% -------------------- Other functions -------------------- 
function [] = Conductance(Elist, T_NEGF, T_BW)
% Physical constants (for conductance units)
e = 1.602176634e-19;      % C
h = 6.62607015e-34;       % J s
G0 = 2*e^2/h;             % conductance quantum (S) including spin

% temperature
T_K = 4.0;                % Kelvin for finite-T conductance
kB_ev = 8.617333262145e-5; % eV/K
kBT = kB_ev * T_K;

% -------------------- Conductances (finite-T Landauer) --------------------
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

function [A_dot, Gdd_full] = doBWfull(Elist, H, SigmaL, SigmaR)
N_lead_each = 30;         % number of lead sites included on each side (extended molecule)
idx_dot = N_lead_each + 1;

I_N = eye(size(H));
A_dot = zeros(size(Elist));     % spectral function at dot from full GR (A = -2 Im G_dd)
for ii = 1:length(Elist)
    E = Elist(ii);
    num_eta = 1e-12;          % tiny numerical broadening for matrix inversion
    % Full Green's function with both absorbers
    GR_full = ((E + 1i*num_eta) * I_N - H - SigmaL - SigmaR) \ I_N;

    % Dot Green's function element (dd)
    Gdd_full = GR_full(idx_dot, idx_dot);
    A_dot(ii) = -2 * imag(Gdd_full);
end
end