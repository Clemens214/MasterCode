%% Variables

% variables for the sample
lengthSample = 48;
orderSample = 1;
hopping = 1;

% variables for the leads
lengthLead = 104;
maxVal = 1;
decay = 0.2; %0.3
offset = 32; %32
%offset should be at most half the length of the leads; normal: 32
hoppingLead = hopping;
hoppingInter = hopping;

lengthTotal = lengthSample+ 2*lengthLead; %256

%variables for the calculation of the current
TempMax = 2; %2;
TempStep = 0.05; %0.05;
TempNum = TempMax/TempStep+1;
Temps = linspace(0, TempMax, TempNum);

chemPotMax = 1; %1;
chemPotStep = 1;
chemPotNum = 2*chemPotMax/chemPotStep+1;
chemPots = linspace(-chemPotMax, chemPotMax, chemPotNum);

%variables for the calculation of the transmission
omegaVal = 2;
omegaMax = omegaVal*hopping;
omegaStep = 0.005;
omegaNum = 2*omegaMax/omegaStep+1;
omegas = linspace(-omegaMax, omegaMax, omegaNum);

%% Calculation
averageTimes = 1;
AllParticle(1:averageTimes) = {cell([1, length(chemPots)])};
for i = 1:averageTimes
    % compute the Hamiltonian of the Sample
    sample = makeHamiltonian(randomNum(disorderStrength, lengthSample), hopping, lengthSample, 'top left');
    
    % preparing the Extended Molecule Hamiltonian
    [totalSystemEM, gammaL_EM, gammaR_EM] = prepareEM(sample, hopping, lengthSample, lengthTotal, maxVal, decay, offset);
    
    % compute the Transmissions
    %plotTransmission (omegas, sample, totalSystemEM, gammaL_EM, gammaR_EM, hoppingInter, hoppingLead, lengthSample, lengthLead);

    %compute the Eigenvectors and the Eigenvalues of the system
    disp('Starting calculation of the Eigenvectors.')
    [Eigenvals, leftEVs, rightEVs, ControlEV, MatchLeft, MatchRight, DiffLeft, DiffRight] = eigenvectors(totalSystemEM);
    disp('Finished calculation of the Eigenvectors.')
    
    PotsParticle(1:length(chemPots)) = {zeros(1,length(Temps))};
    for j = 1:length(chemPots)
        currentsParticle = zeros(1,length(Temps));
        for k = 1:length(Temps)
            [particleResult] = currentQuick(totalSystemEM, gammaL_EM, gammaR_EM, Eigenvals, leftEVs, rightEVs, Temps(k), chemPots(j), lengthLead);
            currentsParticle(k) = particleResult;
            disp(['Average: ', num2str(i), ', chemPot: ', num2str(chemPots(j)), ', Temp: ', num2str(Temps(k))])
        end
        PotsParticle{j} = currentsParticle;
    end
    AllParticle{i} = PotsParticle;
end

[AvgParticle, StdParticle] = Average(AllParticle, chemPots, Temps, averageTimes);

plotGraph (1, 'Particle', Temps, AvgParticle, StdParticle, chemPots)

%%
function [] = plotGraph (value, Title, Temps, AvgVals, StdVals, chemPots)
    figure(value);
    title(Title);
    for i = 1:length(chemPots)
        errorbar(Temps, AvgVals{i}, StdVals{i})
        hold on
    end
    labels = strcat('chemPot = ',cellstr(num2str(chemPots.')));
    legend(labels)
end

%%
function [AvgVals, StdVals] = Average (Data, chemPots, Temps, averageTimes)
    AvgVals(1:length(chemPots)) = {zeros(1,length(Temps))};
    StdVals(1:length(chemPots)) = {zeros(1,length(Temps))};
    for i = 1:length(chemPots)
        avgVals = zeros(1,length(Temps));
        stdVals = zeros(1,length(Temps));
        for j = 1:length(Temps)
            allVals = zeros(1,length(averageTimes));
            for k = 1:averageTimes
                allVals(k) = Data{k}{i}(j);
            end
            avgVals(j) = mean(allVals);
            stdVals(j) = std(allVals);
        end
        AvgVals{i} = avgVals;
        StdVals{i} = stdVals;
    end
end

%%
function [] = plotTransmission (omegas, sample, totalSystemEM, gammaL_EM, gammaR_EM, hoppingInter, hoppingLead, lengthSample, lengthLead)
    [TransmissionsEM, TransmissionsSI] = Transmission(omegas, sample, totalSystemEM, gammaL_EM, gammaR_EM, hoppingInter, hoppingLead, lengthSample, lengthLead);
    plot(omegas, [TransmissionsEM, TransmissionsSI])
    yscale log
end