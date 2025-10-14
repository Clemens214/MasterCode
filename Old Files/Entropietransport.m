%% Variables

% variables for the sample
sizeSample = 48;
orderSample = 1;
eigenenergy = 0;
hopping = 1;
hoppingsSample = hopping*eye(orderSample);

% variables for the leads
sizeLead = 104;
[leadVals, derivVals] = calcVals(maxVal = 1, decay = 0.2, offset = 32);
hoppingLead = hopping;
hoppingsInter = [hopping; hopping];

%variables for the calculation of the current
chemPotMax = 0;%1;
chemPotStep = 1;
chemPots = makeList(chemPotMax, chemPotStep, full=true);

%variables for the calculation of the transmission
omegaVal = 2.5;
omegaMax = omegaVal*hopping;
omegaStep = 0.05;
omegas = makeList(omegaMax, omegaStep, full=true);

%% Calculation

for i = 1:2
if i == 1
    orderSample = 1;
    hoppingsInter = [hopping; hopping];
elseif i == 2
    orderSample = 2;
    hoppingsInter = [hopping, 0; hopping, 0];
end
hoppingsSample = hopping*eye(orderSample);

% compute the Hamiltonian of the Sample
sample = makeSample(eigenenergy, hoppingsSample, sizeSample,  orderSample);

% preparing the Extended Molecule Hamiltonian
[totalSystem, gammaL, gammaR] = makeSystemEM(sample, sizeSample, orderSample, sizeLead, hoppingLead, hoppingsInter, leadVals);

%compute the Eigenvectors and the Eigenvalues of the system
disp('Starting calculation of the Eigenvectors.')
[Eigenvals, leftEVs, rightEVs] = eigenvectors(totalSystem);
disp('Finished calculation of the Eigenvectors.')

Transmission(1:length(chemPots)) = {zeros(1,length(omegas))};
for j = 1:length(chemPots)
    for k = 1:length(omegas)
        TransmissionResult = transmission(totalSystem, gammaL, gammaR, Eigenvals, leftEVs, rightEVs, omegas(k));
        Transmission{j}(k) = TransmissionResult;
        
        disp(['chemPot: ', num2str(chemPots(j)), ', Energy: ', num2str(omegas(k))])
    end
end

if i == 1
    totalSystem1D = totalSystem;
    Transmission1D = Transmission;
elseif i == 2
    totalSystem2D = totalSystem;
    Transmission2D = Transmission;
end
end

Difference = {Transmission1D{1}-Transmission2D{1}};

%% plot
hold on
plotGraph (1, 'Transmission: 1D', omegas, Transmission1D, chemPots)
plotGraph (1, 'Transmission: 1D', [0], {[0]}, chemPots)

plotGraph (2, 'Transmission: 2D', omegas, Transmission2D, chemPots)
plotGraph (2, 'Transmission: 2D', [0], {[0]}, chemPots)

plotGraph (3, 'Transmission: 1D', omegas, Transmission1D, chemPots)
plotGraph (3, 'Transmission: 2D', omegas, Transmission2D, chemPots)
plotGraph (3, 'Both: 1->1D, 2->2D', [0], {[0]}, chemPots)

plotGraph (4, 'Difference: (1D-2D)', omegas, Difference, chemPots)
plotGraph (4, 'Difference: (1D-2D)', [0], {[0]}, chemPots)

%% plotting functions
function [] = plotGraph (value, Title, omegas, Vals, chemPots)
    figure(value);
    title(Title);
    for i = 1:length(chemPots)
        plot(omegas, Vals{i})
        hold on
    end
    labels = strcat('chemPot = ',cellstr(num2str(chemPots.')));
    legend(labels)
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

function [] = saveVar(var, order)
    filename = append('Indices', int2str(order), '.mat');
    save(filename, "var")
end