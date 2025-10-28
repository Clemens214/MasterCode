% resonant_level_BW.m
% Compute resonant-level conductance using Breit-Wigner transmission (from the Quantum Transport notes).
% Uses both finite-T Landauer integral and the T=0 analytic expression.
%
% References:
% Donarini & Grifoni, "Quantum Transport in Interacting Nanojunctions"
% Eqs. (1.88)-(1.90), (3.71). See uploaded PDF. 
%
% Author: ChatGPT
% Date: 2025-10-28

clear; close all; clc;

%% Physical constants
e = 1.602176634e-19;   % elementary charge (C)
h = 6.62607015e-34;    % Planck constant (J s)
G0 = 2*e^2/h;          % conductance quantum (S), includes spin degeneracy

%% Model parameters (user-changeable)
GammaL = 0.05;         % Gamma_L (energy units, e.g. eV)
GammaR = 0.05;         % Gamma_R
eps_d  = 0.0;          % resonant level energy (same units as Gamma)
EF     = 0.0;          % Fermi energy
kB_ev  = 8.617333262145e-5; % Boltzmann constant in eV/K

%% Temperature and numeric control
T_K   = 4.0;           % temperature (Kelvin)
kBT   = kB_ev * T_K;   % in eV
useFiniteT = true;     % true -> evaluate finite-T integral; false -> use zero-T formula
Ewindow = 10 * max(1, (GammaL+GammaR));  % energy window around eps_d for integration (in same units)
nE = 20001;            % number of energy points for numerical integration (odd for symmetry)

%% Energy grid
Emin = eps_d - Ewindow;
Emax = eps_d + Ewindow;
E = linspace(Emin, Emax, nE);

%% Transmission (Breit-Wigner) (Eq. (1.88) in the notes)
Gamma = GammaL + GammaR;
T_E = (GammaL .* GammaR) ./ ( (Gamma/2).^2 + (E - eps_d).^2 );

%% Finite-T conductance (Landauer) using analytic derivative of Fermi function
if useFiniteT
    % Fermi function: f(E) = 1/(1+exp((E-EF)/kBT))
    % derivative: df/dE = -1/(4 kBT cosh^2((E-EF)/(2 kBT)))
    if kBT == 0
        error('kBT is zero but useFiniteT==true. Set T_K>0 or use useFiniteT=false');
    end
    dfdE = -1 ./ (4 * kBT * (cosh((E - EF) ./ (2*kBT))).^2);
    integrand = -dfdE .* T_E; % -df/dE is positive, so integrand >= 0
    G = G0 * trapz(E, integrand);   % numerical integration (Trapezoid)
    % also compute T(EF) for comparison
    T_at_EF = (GammaL * GammaR) / ( (Gamma/2)^2 + (EF - eps_d)^2 );
    fprintf('Finite-T conductance at T=%.3g K: G = %.6g S (%.6g G0)\n', T_K, G, G / G0);
    fprintf('Transmission at EF: T(EF) = %.6g\n', T_at_EF);
else
    % Zero-temperature analytic expression (Eq. (1.90))
    G = G0 * (GammaL * GammaR) / ( (Gamma/2)^2 + (EF - eps_d)^2 );
    fprintf('T=0 conductance (analytic): G = %.6g S (%.6g G0)\n', G, G / G0);
end

%% Plot results
figure('Name','Resonant level conductance (Breit-Wigner)','NumberTitle','off','Position',[200 200 700 500]);
subplot(2,1,1);
plot(E, T_E, 'LineWidth', 1.4);
xlabel('Energy (same units as \Gamma, \epsilon_d)');
ylabel('T(E)');
title('Breit–Wigner transmission T(E) = \Gamma_L\Gamma_R / [(\Gamma/2)^2 + (E-\epsilon_d)^2]');
grid on;
xlim([Emin, Emax]);

subplot(2,1,2);
if useFiniteT
    plot(E, -dfdE, 'LineWidth', 1.2);
    hold on;
    yyaxis right
    plot(E, T_E, '--', 'LineWidth', 1.1);
    ylabel('T(E) (dashed)');
    yyaxis left
    ylabel('-\partial f / \partial E');
    legend({'-df/dE','T(E)'}, 'Location','best');
    title(sprintf('Thermal window (-df/dE) at T=%.2f K and T(E) (dashed)', T_K));
else
    plot(E, T_E, 'LineWidth', 1.4);
    xlabel('Energy');
    ylabel('T(E)');
    title('Transmission (T=0 result uses T(EF) only)');
end
grid on;

%% Gate sweep: conductance vs epsilon_d (zero-T approx using T(EF) or finite-T convolution)
doGateSweep = true;
if doGateSweep
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

%% Final remarks printed (with citation)
fprintf('\nNotes:\n - This script uses the Breit–Wigner transmission and Landauer integral\n');
fprintf('  (see Eqs. (1.88)-(1.90), (3.71) in Donarini & Grifoni, Quantum Transport). \n');
fprintf('  Adjust GammaL/GammaR, temperature, and energy window for convergence.\n');
fprintf('Reference: Quantum Transport in Interacting Nanojunctions (uploaded PDF). \n');
