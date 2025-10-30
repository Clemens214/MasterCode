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
EnergyMax = 3;
EnergyStep = 0.005;
Energies = makeList(EnergyMax, -1*EnergyMax, EnergyStep);
dotEnergy = 0;
EnergiesDot = Energies + dotEnergy;

%% Calculation
% compute the Hamiltonian of the Sample
sample = makeSample(eigenenergy, hoppingsSample, sizeSample,  orderSample);

% preparing the Extended Molecule Hamiltonian
[totalSystem, GammaL, GammaR, SigmaL, SigmaR] = makeSystemEM(sample, sizeSample, orderSample, sizeLead, hoppingLead, hoppingsInter, leadVals);

% compute the resonance curves
%[G_NEGF, G_BW, T_NEGF, T_BW, GammaL_eff, GammaR_eff, Delta]
[G_NEGF, G_BW] = ResonantCheck(Energies, EnergiesDot, totalSystem, GammaL, GammaR, SigmaL, SigmaR);



%% plot
plotLoretzian(EnergiesDot, G_NEGF, G_BW)
%plotTransport(Energies, T_NEGF, T_BW, GammaL_eff, GammaR_eff, Delta)

%% plotting functions
function [] = plotLoretzian(Energies, G_NEGF, G_BW, yFit)
    figure('Name','Gate sweep (zero-T approx)','NumberTitle','off');
    % plot
    plot(Energies, G_NEGF);
    hold on;
    plot(Energies, G_BW);% '--');
    plot(Energies, yFit, '--');
    % labels
    xlabel('\epsilon_d');
    ylabel('G / G_0');
    legend('G_{NEGF}','G_{BW}','Fit','Location','Best');
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
    plot(Elist, T_BW);% '--');
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
    sgtitle(sprintf('EM size: %d sites each lead', sizeLead));
end

%% helping functions
function [values] = makeList(maxVal, minVal, stepVal)
    arguments
        maxVal
        minVal
        stepVal
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