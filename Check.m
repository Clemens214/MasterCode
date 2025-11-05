%% Variables

% variables for the sample
orderSample = 1;
eigenenergy = 0;
dotEnergy = 0;
hopping = 1;
hoppingsSample = hopping*eye(orderSample);

% variables for the leads
sizeLead = 104;
[leadVals, derivVals] = calcVals(maxVal = 1, decay = 0.2, offset = 32);
hoppingLead = hopping;

% variables for the calculation of the current
EnergyMax = 3;
EnergyStep = 0.005;
Energies = makeList(EnergyMax, -1*EnergyMax, EnergyStep);

%% Check the Conductance
sizeSample = 48;
% variables for the coupling
valInter = 1;
hoppingsInter = [valInter*hopping; valInter*hopping];

% compute the Hamiltonian of the Sample
sample = makeSample(eigenenergy, hoppingsSample, sizeSample,  orderSample);

% preparing the Extended Molecule Hamiltonian
[totalSystem, GammaL, GammaR] = makeSystemEM(sample, sizeSample, orderSample, sizeLead, hoppingLead, hoppingsInter, leadVals);

% calculate the Transmission
Conductance = TransCheck(Energies, totalSystem, GammaL, GammaR);

%% Check the Spectrum
sizeSample = 1;
% variables for the coupling
valInter = 0.2;
hoppingsInter = [valInter*hopping; valInter*hopping];

% compute the Hamiltonian of the Sample
sample = makeSample(eigenenergy, hoppingsSample, sizeSample,  orderSample);

% preparing the Extended Molecule Hamiltonian
[totalSystem, GammaL, GammaR, SigmaL, SigmaR] = makeSystemEM(sample, sizeSample, orderSample, sizeLead, hoppingLead, hoppingsInter, leadVals);

% compute the resonance curves
[G_NEGF, G_BW] = ResonantCheck(Energies, totalSystem, GammaL, GammaR, SigmaL, SigmaR, eigEnergy=dotEnergy);
%[G_NEGF, G_BW, T_NEGF, T_BW, GammaL_eff, GammaR_eff, Delta]


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

function [yData] = Lorentzian(xData, max, gamma, x0)
    yData = zeros(size(xData));
    for i = 1:length(xData)
        % max       -> maximum value of the function
        % 2*gamma   -> half-width at half-maximum (HWHM)
        % x0        -> shifts the peak of the distribution
        yData(i) = max/pi * gamma/((xData-x0)^2 + gamma^2);
    end
end