function [Results, varargout] = TransCalc(sample, voltages, sampleVals, leadVals, hoppingsInter, mode, options)
% calculate the transmission through a molecule for zero temperature
arguments
    sample
    voltages
    sampleVals
    leadVals
    hoppingsInter
    mode.EM = false
    mode.SI = true
    options.linearResponse = true
end
    if mode.EM == true
        disp('Using the extended molecule formalism.')
        sizeSample = sampleVals.size;
        orderSample = sampleVals.order;
        sizeLead = leadVals.size;
        hoppingLead = leadVals.hopping;
        [totalSystem, gammaL, gammaR] = makeSystemEM(sample, sizeSample, orderSample, sizeLead, hoppingLead, hoppingsInter);
    elseif mode.SI == true
        disp('Using semi-infinite leads.')
        totalSystem = sample;
        gammaL = zeros(size(sample));
        gammaR = zeros(size(sample));
        mode.energy = sampleVals.energy;
        mode.hopping = leadVals.hopping;
    end
    %disp('Starting calculation of the angular momentum.')
    if options.linearResponse == true
        Energies = voltages;
        Results = Transmission(Energies, totalSystem, gammaL, gammaR, hoppingsInter, mode);
    elseif options.linearResponse == false
        chemPots = setupPots(voltages);
        [Results, values] = integrate(chemPots, totalSystem, gammaL, gammaR, hoppingsInter, mode);
        varargout{1} = values;
    end
    %disp('Starting calculation of the angular momentum.')
end

function [chemPots] = setupPots(voltages)
    chemPots = struct('left', [], 'right', []);
    for j = 1:length(voltages)
        chemPotL = voltages(j)/2;
        chemPotR = -1*voltages(j)/2;
        chemPots(j) = struct('left', chemPotL, 'right', chemPotR);
    end
end

%% integrate the transmission
function [Results, varargout] = integrate(chemPots, totalSystem, gammaL, gammaR, hoppingsInter, mode, options)
arguments
    chemPots
    totalSystem
    gammaL
    gammaR
    hoppingsInter
    mode
    options.stepMult = 10
    options.stepMin = 0.05
end
    % get the step size
    Energies = getEnergies(chemPots);
    Diffs = zeros(1, length(Energies)-1);
    for i = 2:length(Energies)
        Diffs(i) = Energies(i) - Energies(i-1); 
    end
    LCD = lcd(Diffs);
    stepSize = min((1/LCD) / options.stepMult, options.stepMin);

    % calculate the transmissions
    evalPoints = makeList(max(Energies), min(Energies), stepSize);
    values = Transmission(evalPoints, totalSystem, gammaL, gammaR, hoppingsInter, mode);

    % calculate the integrals
    Results = zeros(1, length(chemPots));
    for i = 1:length(chemPots)
        fermiFunc = getFermiFunc(evalPoints, chemPots(i).left) - getFermiFunc(evalPoints, chemPots(i).right);
        % calculate the Result
        yData = fermiFunc .* values;
        if length(evalPoints) > 1
            Results(i) = trapz(evalPoints, yData);
        elseif isscalar(evalPoints)
            Results(i) = 0;
        end
        %disp(['Voltage: ', num2str(chemPots(i).left-chemPots(i).right), ', j=', num2str(i)])
    end
    varargout{1} = values;
end

function [fermiFunc] = getFermiFunc(evalPoints, chemPot, Temp)
arguments
    evalPoints
    chemPot
    Temp = 0
end
    fermiFunc = zeros(size(evalPoints));
    for i = 1:length(evalPoints)
        E = evalPoints(i);
        if Temp ~= 0
            fermiFunc(i) = 1/(exp((E-chemPot)/Temp)+1);
        elseif Temp == 0
            if E <= chemPot
                fermiFunc(i) = 1;
            else
                fermiFunc(i) = 0;
            end
        end
    end
end

%% total transmission in the linear transport approximation
function [Results] = Transmission(Energies, sample, gammaL_EM, gammaR_EM, hoppingsInter, mode)
    %calculates the transport through a molecule in the linear transport approximation
    arguments
        Energies
        sample
        gammaL_EM
        gammaR_EM
        hoppingsInter
        mode
    end
    % calculate the transport matrix and the trace
    Traces = zeros(1, length(Energies));
    parfor i = 1:length(Energies)
        if mode.SI == true
            eigenenergy = mode.energy;
            hoppingLead = mode.hopping;
            [totalSystem, gammaL, gammaR] = makeSystemSI(Energies(i), sample, eigenenergy, hoppingLead, hoppingsInter);
        elseif mode.EM == true
            totalSystem = sample;
            gammaL = gammaL_EM;
            gammaR = gammaR_EM;
        end
        Matrix = TransmissionAlt(Energies(i), totalSystem, gammaL, gammaR);
        Traces(i) = trace(real(Matrix));
    end
    % return the results
    Results = Traces;
end

function [Result] = TransmissionAlt(Energy, totalSystem, gammaL, gammaR)
    % G*B*Gt*C
    % F = decomposition(GreensInv,'lu');
    % Y = F \ B;
    % W = C * Y;
    % Z = F' \ W;
    % t = trace(Z);
    
    % GreensFunc * gammaL * GreensFunc' * gammaR
    GreensInv = Energy*eye(length(totalSystem)) - totalSystem;
    F = decomposition(GreensInv,'lu');    % create reusable LU object (works for sparse/dense)
    
    Y = F \ gammaL;
    W = gammaR * Y;
    Z = F' \ W;                   % uses transpose of factorization
    % t = trace(Z);
    Result = Z;
end

%% helping functions
function [Filtered] = getEnergies(chemPots)
    Energies = zeros(1, length(chemPots)*2);
    for i = 1:length(chemPots)
        Energies(2*i-1) = chemPots(i).left;
        Energies(2*i) = chemPots(i).right;
    end
    Sorted = sort(Energies);
    Filtered = unique(Sorted);
end

function [values] = makeList(maxVal, minVal, stepVal)
    arguments
        maxVal
        minVal
        stepVal
    end
    numVal = (maxVal-minVal)/stepVal+1;
    values = linspace(minVal, maxVal, numVal);
end