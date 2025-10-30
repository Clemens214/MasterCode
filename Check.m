%% Variables

% variables for the sample
sizeSample = 1;
orderSample = 1;
eigenenergy = 0;
hopping = 1;
hoppingsSample = hopping*eye(orderSample);

% variables for the leads
sizeLead = 104;
[leadVals, derivVals] = calcVals(maxVal = 1, decay = 0.2, offset = 32);
hoppingLead = hopping;

% variables for the coupling
valInter = 0.2;
hoppingsInter = [valInter*hopping; valInter*hopping];
hoppingsDeriv = [0; 0];

%variables for the calculation of the current
voltageMax = 5;
voltageStep = 0.5;
voltages = makeList(voltageMax, voltageStep);
chemPots = setupPots(voltages);
Energies = getEnergies(chemPots);
dotEnergy = 0;
EnergiesDot = Energies + dotEnergy;

%% Calculation
% compute the Hamiltonian of the Sample
sample = makeSample(eigenenergy, hoppingsSample, sizeSample,  orderSample);
    
% preparing the Extended Molecule Hamiltonian
[totalSystem, GammaL, GammaR, SigmaL, SigmaR] = makeSystemEM(sample, sizeSample, orderSample, sizeLead, hoppingLead, hoppingsInter, leadVals);

checkMatrix(totalSystem);

%compute the Eigenvectors and the Eigenvalues of the system
[Eigenvals, leftEVs, rightEVs, Product] = eigenvectors(totalSystem);%, checkMore=true);

%% OLD
    N_lead_each = 30;         % number of lead sites included on each side (extended molecule)
    idx_dot = N_lead_each + 1;
    N = 2*N_lead_each + 1;    % total EM sites
    
    N_absorb_layers = 4;      % number of outer-most sites per side with absorbing Sigma
    left_abs_idx = 1 : min(N_absorb_layers, N_lead_each);
    right_abs_idx = N - (0:(min(N_absorb_layers, N_lead_each)-1));
    
    % Energy grid
    t_lead = 1.0;             % lead hopping (energy units)
    Emin = -3*t_lead;
    Emax = 3*t_lead;
    nE = 2001;
    Energies = linspace(Emin, Emax, nE);
    
    % -------------------- Gate sweep (zero-T approx using T(EF)) --------------------
    neps = 201;
    epsilon_dot = 0.0;        % dot onsite energy (resonant level)
    EnergiesDot = linspace(epsilon_dot - 2.5, epsilon_dot + 2.5, neps);
    
    % -------------------- Absorbing self-energies on outer EM sites --------------------
    eta = 0.06;               % absorbing strength (positive) -> Sigma = -i*eta on outer sites
    
    SigmaL = zeros(N,N);
    SigmaL(sub2ind([N,N], left_abs_idx, left_abs_idx)) = -1i * eta;
    
    SigmaR = zeros(N,N);
    SigmaR(sub2ind([N,N], right_abs_idx, right_abs_idx)) = -1i * eta;

    % Gamma matrices for the outer absorbers (used in NEGF transmission)
    GammaL = 1i * (SigmaL - SigmaL');
    GammaR = 1i * (SigmaR - SigmaR');
    
    % setup
    H = setupH();
    totalSystem = H + SigmaL + SigmaR;

[G_NEGF, G_BW, T_NEGF, T_BW, GammaL_eff, GammaR_eff, Delta] = ResonantCheck(Energies, EnergiesDot, totalSystem, GammaL, GammaR, SigmaL, SigmaR);

%% plot
plotLoretzian(EnergiesDot, G_NEGF, G_BW)

plotTransport(Energies, T_NEGF, T_BW, GammaL_eff, GammaR_eff, Delta)

%% chemPots
function [chemPots] = setupPots(voltages)
    chemPots = struct('left', [], 'right', []);
    for j = 1:length(voltages)
        chemPotL = voltages(j)/2;
        chemPotR = -1*voltages(j)/2;
        chemPots(j) = struct('left', chemPotL, 'right', chemPotR);
    end
end

function [totalSysDeriv] = makeDeriv(sizeSample, orderSample, sizeLead, hoppingsDeriv, derivVals)
    sampleDeriv = zeros(sizeSample*orderSample, sizeSample*orderSample);
    hoppingDeriv = 0;
    [totalSysDeriv, ~, ~] = makeSystemEM(sampleDeriv, sizeSample, orderSample, sizeLead, hoppingDeriv, hoppingsDeriv, derivVals, check=false);
end

function [H] = setupH()
    % -------------------- Model parameters --------------------
    t_lead = 1.0;             % lead hopping (energy units)
    epsilon_lead = 0.0;       % lead onsite energy
    t_c = 0.2;                % coupling dot <-> first lead site
    epsilon_dot = 0.0;        % dot onsite energy (resonant level)
    N_lead_each = 30;         % number of lead sites included on each side (extended molecule)
    
    % -------------------- Build Extended Molecule Hamiltonian (1D chain) --------------------
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

%% plotting functions
function [] = plotLoretzian(eps_range, G_NEGF, G_BW)
    figure('Name','Gate sweep (zero-T approx)','NumberTitle','off');
    % plot
    plot(eps_range, G_NEGF);
    hold on;
    plot(eps_range, G_BW, '--');
    % labels
    xlabel('\epsilon_d');
    ylabel('G / G_0');
    legend('G_{NEGF}','G_{BW}','Location','Best');
    title('Gate sweep: NEGF (EM) vs Breit–Wigner using EM-extracted \Gamma');
    grid on;
end

function [] = plotTransport(Elist, T_NEGF, T_BW, GammaL_eff, GammaR_eff, Delta_eff)
    % -------------------- Plots --------------------
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

%% helping functions
function [values] = makeList(maxVal, stepVal, options)
    arguments
        maxVal 
        stepVal 
        options.full = false
    end
    if options.full == false
        minVal = 0;
    else
        minVal = -1*maxVal;
    end
    numVal = (maxVal-minVal)/stepVal+1;
    values = linspace(minVal, maxVal, numVal);
end

function [leadVals, derivVals] = calcVals(opt)%(maxVal, decay, offset)
    arguments
        opt.maxVal = 1
        opt.decay = 0.3
        % offset should be at most half the length of the leads
        opt.offset = 32
        % normal: 32
    end
    maxVal = opt.maxVal;
    decay = opt.decay;
    offset = opt.offset;
    leadVals = {maxVal, decay, offset};
    derivVals = {0, decay, offset};
end

function [Filtered] = getEnergies(chemPots)
    Energies = zeros(1, length(chemPots)*2);
    for i = 1:length(chemPots)
        Energies(2*i-1) = chemPots(i).left;
        Energies(2*i) = chemPots(i).right;
    end
    Sorted = sort(Energies);
    Filtered = unique(Sorted);
end

function [] = saveVar(var, order)
    filename = append('Indices', int2str(order), '.mat');
    save(filename, "var")
end

function [] = checkMatrix(totalSystem)
    % Example matrix (replace with your A)
    A = totalSystem;
    % 1) Compute right eigenvectors and eigenvalues
    [V, ~] = eig(A);      % A * V = V * D
    % Conditioning of eigenvector matrix
    condV = cond(V);
    % Rank of eigenvector matrix
    rankV = rank(V);
    % Check if matrix is defective
    if rankV < size(A,1)
        disp('Matrix appears defective (not diagonalizable).');
    else
        disp('Matrix is diagonalizable but ill-conditioned.');
    end
    fprintf('Condition number of V: %g\n', condV);
end