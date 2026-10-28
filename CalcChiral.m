%% CalcChiral.m
% Torque in a junction whose molecule carries an axially chiral end group.
%
% Geometry
%   The molecule is rotated rigidly by theta against a fixed stator.
%   The contact angles in the molecular frame are
%       Theta_L = theta - phiL ,   Theta_R = theta - phiR ,
%   so the relative twist Delta theta = Theta_L - Theta_R = phiR - phiL stays
%   fixed while theta runs over a full revolution.
%   The end group is fixed in the molecular frame (makeSample).
%
% Torques (all energy-resolved traces include the 1/(2 pi) of the thesis)
%   net     : tau     = -dH/dtheta = tau_L + tau_R  (rotation of the molecule)
%   contact : tau_L   = -dH/dTheta_L                 (torsion at the left contact)
%   s / a   : Gamma_L + Gamma_R  (equilibrium-like)  /  Gamma_L - Gamma_R (bias)
%
% Cycle quantities (motor criterion)
%   C_s(omega) = oint dtheta <tau>_s(omega, theta)    (= 0 exactly, conservative)
%   C_a(omega) = oint dtheta <tau>_a(omega, theta)
%   W(V)       = 1/2 int_{mu0 - V/2}^{mu0 + V/2} d omega  C_a(omega)
%              = work done on the rotor per revolution at T = 0.
%
% Requires the updated makeSample.m, makeSystemSI.m and TorqueCalc.m.

clear; clc;

%% Variables

% variables for the sample
sizeSample = 5;
orderSample = 2;
energySample = 0;
hopping = 1;
hoppingsSample = hopping*eye(orderSample);
sampleVals = struct('size', sizeSample, 'order', orderSample, 'energy', energySample, 'hopping', hoppingsSample);

% axially chiral end group (see makeSample.m)
endGroup = struct( ...
    'chiralSites', 3, ...          % number of chiral bonds (hopping) or sites (onsite)
    'anisotropy',  0.3, ...        % kappa; 0 restores the isotropic chain
    'twist',       pi/6, ...       % director rotation per bond; -pi/6 = enantiomer, 0 = achiral
    'phase',       0, ...          % director angle of the first bond
    'endType',     "hopping", ...  % "hopping" (particle-hole symmetric) or "onsite"
    'shift',       0, ...          % on-site shift of the end-group sites
    'side',        "left");

% variables for the leads
sizeLead = 1;
energyLead = energySample;
hoppingLead = hopping;
leadVals = struct('size', sizeLead, 'energy', energyLead, 'hopping', hoppingLead);

% stator geometry: relative twists Delta theta = phiR - phiL
relTwists = [0, pi/4, pi/3];
phiR = 0;
phiL = phiR - relTwists;

% rigid rotation of the molecule over one revolution (periodic grid)
nTheta = 72;
thetas = 2*pi*(0:nTheta-1)/nTheta;

% energies (symmetric grid, so that omega -> -omega can be checked)
EnergyMax = 2.5;
EnergyStep = 0.005;
Energies = -EnergyMax : EnergyStep : EnergyMax;

% bias windows for the work per cycle
voltages = 0 : 0.02 : 2;
gates = [0, 0.5];            % mu0: centre of the bias window

%% Sample
egArgs = namedargs2cell(endGroup);
sample = makeSample(energySample, hoppingsSample, sizeSample, orderSample, egArgs{:});

%% Calculation
nTw = numel(relTwists);
nE  = numel(Energies);
netS = zeros(nTheta, nE, nTw);   netA = zeros(nTheta, nE, nTw);
conS = zeros(nTheta, nE, nTw);   conA = zeros(nTheta, nE, nTw);

for a = 1:nTw
    for k = 1:nTheta
        ThL = thetas(k) - phiL(a);
        ThR = thetas(k) - phiR;
        hoppingsInter = hopping*[cos(ThL), sin(ThL); cos(ThR), sin(ThR)];
        dL = hopping*[-sin(ThL), cos(ThL)];      % dt_L/dtheta
        dR = hopping*[-sin(ThR), cos(ThR)];      % dt_R/dtheta
        derivNet = [dL; dR];                     % rigid rotation of the molecule
        derivL   = [dL; 0, 0];                   % left contact only

        netS(k, :, a) = TorqueCalc(sample, Energies, sampleVals, leadVals, hoppingsInter, derivNet, ...
                                   linearResponse=true, conservative=true)/(2*pi);
        netA(k, :, a) = TorqueCalc(sample, Energies, sampleVals, leadVals, hoppingsInter, derivNet, ...
                                   linearResponse=true, nonconservative=true)/(2*pi);
        conS(k, :, a) = TorqueCalc(sample, Energies, sampleVals, leadVals, hoppingsInter, derivL, ...
                                   linearResponse=true, conservative=true)/(2*pi);
        conA(k, :, a) = TorqueCalc(sample, Energies, sampleVals, leadVals, hoppingsInter, derivL, ...
                                   linearResponse=true, nonconservative=true)/(2*pi);
    end
    disp(['Relative twist: ', num2str(relTwists(a)/pi), ' pi done'])
end

%% Cycle integrals and work per revolution
dTheta = 2*pi/nTheta;
cycS = squeeze(sum(netS, 1))*dTheta;   % nE x nTw, conservative: should vanish
cycA = squeeze(sum(netA, 1))*dTheta;   % nE x nTw
if nTw == 1
    cycS = cycS(:);  cycA = cycA(:);
end

W = zeros(numel(voltages), numel(gates), nTw);
for a = 1:nTw
    for g = 1:numel(gates)
        for v = 1:numel(voltages)
            W(v, g, a) = workPerCycle(Energies, cycA(:, a), gates(g), voltages(v));
        end
    end
end

%% Checks
fprintf('\nEnd group: %s, kappa = %.2f, twist = %.3f pi, %d sites\n', ...
        endGroup.endType, endGroup.anisotropy, endGroup.twist/pi, endGroup.chiralSites);
for a = 1:nTw
    flipS = netS(:, end:-1:1, a);  flipA = netA(:, end:-1:1, a);
    fprintf('Delta theta = %.3f pi\n', relTwists(a)/pi);
    fprintf('  max |<tau>_s|                  : %.3e\n', max(abs(netS(:, :, a)), [], 'all'));
    fprintf('  max |<tau>_a|                  : %.3e\n', max(abs(netA(:, :, a)), [], 'all'));
    fprintf('  max |oint <tau>_s dtheta|      : %.3e   (conservative -> 0)\n', max(abs(cycS(:, a))));
    fprintf('  max |oint <tau>_a dtheta|      : %.3e\n', max(abs(cycA(:, a))));
    fprintf('  max |<tau>_s(w) + <tau>_s(-w)| : %.3e   (0 if particle-hole symmetric)\n', ...
            max(abs(netS(:, :, a) + flipS), [], 'all'));
    fprintf('  max |<tau>_a(w) + <tau>_a(-w)| : %.3e\n', max(abs(netA(:, :, a) + flipA), [], 'all'));
    for g = 1:numel(gates)
        fprintf('  W(V = %.2f, mu0 = %.2f)        : %+.3e\n', voltages(end), gates(g), W(end, g, a));
    end
end

%% Plots
a = 1;   % twist shown in the maps
figure('Name', 'Net torque maps', 'Color', 'w');
subplot(1, 2, 1);
imagesc(thetas/pi, Energies, netS(:, :, a).'); axis xy; colorbar;
xlabel('$\theta/\pi$', 'Interpreter', 'latex'); ylabel('$\omega/|t|$', 'Interpreter', 'latex');
title(sprintf('$\\langle\\hat\\tau\\rangle_s$, $\\Delta\\theta = %.2f\\pi$', relTwists(a)/pi), 'Interpreter', 'latex');
subplot(1, 2, 2);
imagesc(thetas/pi, Energies, netA(:, :, a).'); axis xy; colorbar;
xlabel('$\theta/\pi$', 'Interpreter', 'latex'); ylabel('$\omega/|t|$', 'Interpreter', 'latex');
title(sprintf('$\\langle\\hat\\tau\\rangle_a$, $\\Delta\\theta = %.2f\\pi$', relTwists(a)/pi), 'Interpreter', 'latex');

legTw = arrayfun(@(x) sprintf('\\Delta\\theta = %.3f\\pi', x/pi), relTwists, 'UniformOutput', false);
figure('Name', 'Cycle integrals', 'Color', 'w');
subplot(2, 1, 1); plot(Energies, cycS, 'LineWidth', 1.2); grid on; legend(legTw);
ylabel('$\oint \langle\hat\tau\rangle_s\, d\theta$', 'Interpreter', 'latex');
subplot(2, 1, 2); plot(Energies, cycA, 'LineWidth', 1.2); grid on;
ylabel('$\oint \langle\hat\tau\rangle_a\, d\theta$', 'Interpreter', 'latex');
xlabel('$\omega/|t|$', 'Interpreter', 'latex');

figure('Name', 'Work per revolution', 'Color', 'w'); hold on;
for g = 1:numel(gates)
    plot(voltages, squeeze(W(:, g, :)), 'LineWidth', 1.2);
end
hold off; grid on;
xlabel('$V/|t|$', 'Interpreter', 'latex'); ylabel('$W$ per revolution', 'Interpreter', 'latex');
legW = {};
for g = 1:numel(gates)
    for k = 1:nTw
        legW{end+1} = sprintf('\\mu_0 = %.2f, \\Delta\\theta = %.2f\\pi', gates(g), relTwists(k)/pi); %#ok<SAGROW>
    end
end
legend(legW);

%% helping functions
function [W] = workPerCycle(Energies, cycA, mu0, V)
% W = 1/2 int_{mu0-V/2}^{mu0+V/2} d omega  oint <tau>_a d theta   (T = 0)
    lo = mu0 - V/2;
    hi = mu0 + V/2;
    if hi <= lo
        W = 0;
        return
    end
    nPts = max(2, ceil((hi - lo)/(Energies(2) - Energies(1))) + 1);
    w = linspace(lo, hi, nPts);
    W = 0.5*trapz(w, interp1(Energies, cycA, w, 'linear', 0));
end
