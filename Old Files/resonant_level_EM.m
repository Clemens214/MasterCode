% resonant_level_EM.m
% Calculate resonant-level conductance using the extended-molecule formalism
% (local absorbing self-energies on outer surface atoms).
%
% Outputs:
%  - Transmission T(E) vs energy
%  - Linear conductance G (in units of G0 = 2 e^2 / h) at specified temperature
%
% Author: ChatGPT (template)
% Date: 2025-10-28

clear; close all; clc;

%% Physical constants
e = 1.602176634e-19;      % elementary charge (C)
h = 6.62607015e-34;       % Planck constant (J s)
G0 = 2*e^2/h;             % conductance quantum (S) including spin degeneracy

%% Model parameters (change as needed)
t_lead = 1.0;             % lead hopping (energy units)
epsilon_lead = 0.0;       % lead onsite energy (sets band center)
t_c = 0.2;                % coupling between dot and first lead site
epsilon_dot = 0.0;        % dot onsite energy (resonant level) - can be swept
N_lead_each = 40;         % number of lead sites included on each side (extended molecule)
N_absorb_layers = 4;      % number of outer-most sites on each side to apply absorbing Sigma
eta = 0.05;               % absorbing strength (positive real). We'll apply Sigma = -i*eta on those sites
                          % recommended: test several eta values to find plateau (e.g. 1e-3 .. 0.5)
                          
%% Energy grid for T(E)
Emin = -3*t_lead;         % energy grid limits (choose wide enough to include lead bands)
Emax = 3*t_lead;
nE = 2001;
Elist = linspace(Emin, Emax, nE);

%% Temperature and Fermi energy (for conductance)
kB = 8.617333262145e-5;   % eV/K (Boltzmann constant)
T_K = 10;                 % temperature in Kelvin (set to small value for near-zero T)
kBT = kB * T_K;           % in eV (consistent with energies assumed)
EF = 0.0;                 % Fermi energy (same energy units as Hamiltonian)

%% Build extended-molecule Hamiltonian (1D chain: left lead, dot, right lead)
N = N_lead_each*2 + 1;    % total sites (left lead sites + dot + right lead sites)
H = zeros(N,N);

% Indexing:
% 1 .. N_lead_each         = left lead sites (1 is outer-most left surface)
% N_lead_each + 1         = dot site
% N_lead_each + 2 .. N    = right lead sites (N is outer-most right surface)

idx_dot = N_lead_each + 1;

% Left lead block (we connect them in-site order from surface inward)
for n = 1:N_lead_each
    H(n,n) = epsilon_lead;
    if n < N_lead_each
        H(n, n+1) = -t_lead;
        H(n+1, n) = -t_lead;
    end
end

% Right lead block
offsetR = N_lead_each + 1;
for n = (idx_dot+1):N
    H(n,n) = epsilon_lead;
end
for n = (idx_dot+1):(N-1)
    H(n, n+1) = -t_lead;
    H(n+1, n) = -t_lead;
end

% Connect dot to first sites of left and right lead (these are sites N_lead_each and idx_dot+1)
if N_lead_each >= 1
    H(idx_dot, idx_dot-1) = -t_c;   % dot to left-first site
    H(idx_dot-1, idx_dot) = -t_c;
    H(idx_dot, idx_dot+1) = -t_c;   % dot to right-first site
    H(idx_dot+1, idx_dot) = -t_c;
else
    error('N_lead_each must be at least 1.');
end

% Also put dot onsite energy
H(idx_dot, idx_dot) = epsilon_dot;

% Optionally show Hamiltonian sparsity
% spy(H)

%% Construct absorbing self-energies on outer sites
% We'll apply Sigma_L = -i*eta on the first N_absorb_layers sites (left outermost),
% and Sigma_R = -i*eta on the last N_absorb_layers sites (right outermost).
SigmaL = zeros(N,N);
SigmaR = zeros(N,N);

% Left absorbing region indices (1 is outermost)
left_abs_idx = 1 : min(N_absorb_layers, N_lead_each);
% Right absorbing region indices (outermost right)
right_abs_idx = N - (0:(min(N_absorb_layers, N_lead_each)-1));

SigmaL(sub2ind([N N], left_abs_idx, left_abs_idx)) = -1i * eta;
SigmaR(sub2ind([N N], right_abs_idx, right_abs_idx)) = -1i * eta;

% total Sigma (left + right)
% Note: these are energy-independent approximations (local imaginary potentials).
%% Preallocate
T_of_E = zeros(size(Elist));

I_N = eye(N);

% small numerical broadening to ensure retarded inversion well-defined
num_eta = 1e-12;

% Compute T(E) loop
for ii = 1:length(Elist)
    E = Elist(ii);
    % Retarded Green's function: G^R = inv((E + i0+)I - H - SigmaL - SigmaR)
    GR = ((E + 1i*num_eta) * I_N - H - SigmaL - SigmaR) \ I_N;
    GA = GR';  % hermitian conjugate (since GR is matrix, GA = GR^\dagger)
    
    % Compute Gamma matrices: Gamma = i (Sigma - Sigma^\dagger)
    % With Sigma = -i*eta_diag, Gamma becomes positive 2*eta on those diag sites
    GammaL = 1i * (SigmaL - SigmaL');
    GammaR = 1i * (SigmaR - SigmaR');
    
    % Transmission: T = trace(GammaL * GR * GammaR * GA)
    % ensure real small negative numerical values are clipped
    Tval = real(trace(GammaL * GR * GammaR * GA));
    if Tval < 0 && Tval > -1e-12
        Tval = 0;
    end
    T_of_E(ii) = Tval;
end

%% Compute linear conductance at finite temperature
% G = G0 * integral [-df/dE * T(E) dE]
% df/dE = -1/(4 kBT cosh^2((E-EF)/(2 kBT))) for Fermi function, but we'll compute numerically
f = @(E) 1./(1 + exp((E - EF)/kBT));
dfdE_numeric = - ( f(Elist + 1e-6) - f(Elist - 1e-6) ) / (2e-6);  % numerical derivative (small step)
% better: analytic derivative:
dfdE_analytic = @(E) -1./(4 * kBT * (cosh((E - EF)/(2*kBT))).^2);
dfdE = dfdE_analytic(Elist);

% integrate using trapezoidal rule
G_linear = G0 * trapz(Elist, -dfdE .* T_of_E);  % -dfdE is positive
G_in_G0 = G_linear / G0;   % conductance in units of G0

%% Plots
figure;
subplot(2,1,1);
plot(Elist, T_of_E, 'LineWidth', 1.5);
xlabel('Energy (units of t)');
ylabel('T(E)');
title('Transmission vs Energy');
grid on;

subplot(2,1,2);
plot(Elist, -dfdE, 'LineWidth', 1);
xlabel('Energy (units of t)');
ylabel('-df/dE');
title(sprintf('Thermal window (-df/dE) at T = %.2f K (kBT = %.3e eV)', T_K, kBT));
grid on;

sgtitle(sprintf('Resonant level: epsilon_d = %.3g, eta = %.3g, Nlead_each = %d', epsilon_dot, eta, N_lead_each));

fprintf('Linear conductance (finite T = %.2f K): G = %.6g S  (%.6g G0)\n', T_K, G_linear, G_in_G0);
fprintf('Peak transmission (max T) = %.6g\n', max(T_of_E));

%% Optional: sweep dot energy (gate) and compute G0 at T=0 (approx)
% This shows resonance line-shape vs epsilon_dot
doGateSweep = true;
if doGateSweep
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

%% Recommendations & checks (printed)
disp('--- Recommendations & checks ---');
disp('- Vary N_lead_each (increase) to confirm convergence of T(E) and G.');
disp('- Vary eta across a few decades (e.g. 1e-3 .. 1.0) and look for a plateau in results.');
disp('- If resonance is shifted relative to expected position, consider increasing the extended-molecule size or include a real part of Sigma.');
disp('- For more realistic leads use surface Green''s functions instead of local -i*eta if energy dependence is important.');

