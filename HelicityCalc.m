function [Results, varargout] = HelicityCalc(sample, voltages, sampleVals, leadVals, hoppingsInter, choice, mode, options)
% calculate the helicity of a molecule for zero temperature
arguments
    sample
    voltages
    sampleVals
    leadVals
    hoppingsInter
    choice.conservative = true
    choice.nonconservative = false
    choice.left = false
    choice.right = false
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
        operator = HelicityOperator(length(totalSystem), sizeLead);
    elseif mode.SI == true
        disp('Using semi-infinite leads.')
        totalSystem = sample;
        operator = zeros(size(sample));
        gammaL = zeros(size(sample));
        gammaR = zeros(size(sample));
        mode.energy = sampleVals.energy;
        mode.hopping = leadVals.hopping;
    end
    %disp('Starting calculation of the helicity.')
    if options.linearResponse == true
        Energies = voltages;
        if choice.conservative == true || choice.nonconservative == true || choice.left == true || choice.right == true
            Results = Helicity(Energies, totalSystem, operator, gammaL, gammaR, hoppingsInter, choice, mode);
        else
            choiceL = choice;
            choiceL.left = true;
            ResultsL =  Helicity(Energies, totalSystem, operator, gammaL, gammaR, hoppingsInter, choiceL, mode);
            choiceR = choice;
            choiceR.right = true;
            ResultsR =  Helicity(Energies, totalSystem, operator, gammaL, gammaR, hoppingsInter, choiceR, mode);
            Results = ResultsL + ResultsR;
        end
    elseif options.linearResponse == false
        chemPots = setupPots(voltages);
        if choice.conservative == true || choice.nonconservative == true || choice.left == true || choice.right == true
            [Results, values, Energies] = integrate(chemPots, totalSystem, operator, gammaL, gammaR, hoppingsInter, choice, mode);
        else
            choiceL = choice;
            choiceL.left = true;
            [ResultsL, valuesL, EnergiesL] = integrate(chemPots, totalSystem, operator, gammaL, gammaR, hoppingsInter, choiceL, mode);
            choiceR = choice;
            choiceR.right = true;
            [ResultsR, valuesR, EnergiesR] = integrate(chemPots, totalSystem, operator, gammaL, gammaR, hoppingsInter, choiceR, mode);
            Results = ResultsL + ResultsR;
            [values, Energies] = combine(valuesL, valuesR, EnergiesL, EnergiesR);
        end
        varargout{1} = values;
        varargout{2} = Energies;
    end
    %disp('Finished calculation of the helicity.')
end

function [chemPots] = setupPots(voltages)
    chemPots = struct('left', [], 'right', []);
    for j = 1:length(voltages)
        chemPotL = voltages(j)/2;
        chemPotR = -1*voltages(j)/2;
        chemPots(j) = struct('left', chemPotL, 'right', chemPotR);
    end
end

function [operator] = HelicityOperator(sizeTotal, sizeLead, order)
arguments
    sizeTotal
    sizeLead
    order = 2
end
    Pauli = [0, -1j; 1j, 0];
    % generate the helicity operator
    sizeCenter = sizeTotal - 2*sizeLead;
    sizeSample = sizeCenter/order;
    small = zeros(sizeSample, sizeSample);
    for row = 1:sizeSample
        for column = 1:sizeSample
            if column == row+1
                small(row, column) = 1;
            elseif row == column+1
                small(row, column) = -1;
            end
        end
    end
    sample = zeros(sizeCenter, sizeCenter);
    for row = 1:sizeSample
        for column = 1:sizeSample
            rangeRow = order*(row-1)+1 : order*row;
            rangeColumn = order*(column-1)+1 : order*column;
            if small(row, column) ~= 0
                sample(rangeRow, rangeColumn) = small(row, column) * Pauli;
            end
        end
    end
    % embed the operator in total Hamiltonian
    operator = zeros(sizeTotal, sizeTotal);
    mid = sizeLead+1 : sizeTotal-sizeLead;
    operator(mid, mid) = sample;
end

function [values, Energies] = combine(valuesL, valuesR, EnergiesL, EnergiesR)
    idx = 1;
    jdx = 1;
    size = max(length(EnergiesL), length(EnergiesR));
    values = zeros(1, size);
    Energies = zeros(1, size);
    for i = 1:size
        if jdx > length(EnergiesR) || EnergiesL(idx) < EnergiesR(jdx)
            idx = idx+1;
            values(i) = valuesL(idx);
            Energies(i) = EnergiesL(idx);
        elseif idx > length(EnergiesR) || EnergiesL(idx) > EnergiesR(jdx)
            jdx = jdx+1;
            values(i) = valuesR(idx);
            Energies(i) = EnergiesR(idx);
        elseif EnergiesL(idx) == EnergiesR(jdx)
            values(i) = valuesL(idx) + valuesR(jdx);
            Energies(i) = EnergiesL(idx);
            idx = idx+1;
            jdx = jdx+1;
        end
    end
end

%% integrate the helicity
function [Results, varargout] = integrate(chemPots, totalSystem, operator, gammaL, gammaR, hoppingsInter, choice, mode, options)
arguments
    chemPots
    totalSystem
    operator
    gammaL
    gammaR
    hoppingsInter
    choice
    mode
    options.stepMult = 10
    options.stepMin = 0.05
    options.minVal = -3
    options.check = false
end
    % get the bounds
    [maxPoint, minPoint] = choiceBounds(chemPots, choice);

    % get the step size
    Energies = getEnergies(chemPots);
    Diffs = zeros(1, length(Energies)-1);
    for i = 2:length(Energies)
        Diffs(i) = Energies(i) - Energies(i-1); 
    end
    LCD = lcd(Diffs);
    stepSize = min((1/LCD) / options.stepMult, options.stepMin);

    % calculate the transmissions
    evalPoints = makeList(maxPoint, minPoint, stepSize);
    values = Helicity(evalPoints, totalSystem, operator, gammaL, gammaR, hoppingsInter, choice, mode);
    
    % calculate the integrals
    Results = zeros(1, length(chemPots));
    for i = 1:length(chemPots)
        fermiFunc = choiceFermiFunc(evalPoints, chemPots(i).left, chemPots(i).right, choice);
        if options.check == true
            testFermi(evalPoints, fermiFunc, choice);
        end
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
    varargout{2} = evalPoints;
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
            if E < chemPot
                fermiFunc(i) = 1;
            elseif E > chemPot
                fermiFunc(i) = 0;
            elseif E == chemPot
                fermiFunc(i) = 0.5;
            end
        end
    end
end

%% total helicity in the linear transport approximation
function [Results] = Helicity(Energies, sample, operator_EM, gammaL_EM, gammaR_EM, hoppingsInter, choice, mode)
arguments
    Energies
    sample
    operator_EM
    gammaL_EM
    gammaR_EM
    hoppingsInter
    choice
    mode
end
    % calculate the transport matrix and the trace
    Traces = zeros(1, length(Energies));
    parfor i = 1:length(Energies)
        if mode.SI == true
            eigenenergy = mode.energy;
            hoppingLead = mode.hopping;
            [totalSystem, gammaL, gammaR] = makeSystemSI(Energies(i), sample, eigenenergy, hoppingLead, hoppingsInter);
            sizeExtra = (size(totalSystem) - size(sample))/2;
            operator = HelicityOperator(length(totalSystem), sizeExtra(1));
        elseif mode.EM == true
            totalSystem = sample;
            operator = operator_EM;
            gammaL = gammaL_EM;
            gammaR = gammaR_EM;
        end
        Matrix = choiceLin(Energies(i), operator, totalSystem, gammaL, gammaR, choice);
        Traces(i) = trace(real(Matrix));
    end
    % return the results
    Results = Traces;
end

function [Result] = HelicityAlt(Energy, operator, totalSystem, midFactor, options)
    % A*G*B*Gt
    % F = decomposition(GreensInv,'lu');
    % Y = F \ B;
    % T = A * Y;
    % Z = F' \ T;
    % t = trace(Z);
arguments
    Energy
    operator
    totalSystem
    midFactor
    options.eta = 1E-12
end
    eta = 1j*options.eta;
    % operator * GreensFunc * midFactor * GreensFunc'
    GreensInv = (Energy+eta)*eye(length(totalSystem)) - totalSystem;
    F = decomposition(GreensInv,'lu');    % create reusable factorization object
    
    Y = F \ midFactor;                     % solves Aw * Y = B
    T = operator * Y;
    Z = F' \ T;                    % solves Aw' * Z = T
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

%% choosing functions
function [maxPoint, minPoint] = choiceBounds(chemPots, choice, options)
arguments
    chemPots
    choice
    options.minVal = 0%-3
    options.print = false
end
    % get the bounds
    if choice.conservative == true
        strDisp = 'Conservative';
        maxPoint = max([[chemPots.left], [chemPots.right]]);
        minPoint = min([[chemPots.left], [chemPots.right], options.minVal]);
    elseif choice.nonconservative == true
        strDisp = 'Nonconservative';
        maxPoint = max([[chemPots.left], [chemPots.right]]);
        minPoint = min([[chemPots.left], [chemPots.right]]);
    elseif choice.left == true
        strDisp = 'Left';
        maxPoint = max([chemPots.left]);
        minPoint = min([[chemPots.left], options.minVal]);
    elseif choice.right == true
        strDisp = 'Right';
        maxPoint = max([chemPots.right]);
        minPoint = min([[chemPots.right], options.minVal]);
    end
    if options.print == true
        disp([strDisp, '; Maximum: ',num2str(maxPoint), ', Minimum: ',num2str(minPoint)])
    end
end

function [fermiFunc] = choiceFermiFunc(evalPoints, chemPotL, chemPotR, choice)
    if choice.conservative == true
        fermiFuncL = getFermiFunc(evalPoints, chemPotL);
        fermiFuncR = getFermiFunc(evalPoints, chemPotR);
        fermiFunc = 0.5 * (fermiFuncL + fermiFuncR);
    elseif choice.nonconservative == true
        fermiFuncL = getFermiFunc(evalPoints, chemPotL);
        fermiFuncR = getFermiFunc(evalPoints, chemPotR);
        fermiFunc = 0.5 * (fermiFuncL - fermiFuncR);
    elseif choice.left == true
        fermiFunc = getFermiFunc(evalPoints, chemPotL);
    elseif choice.right == true
        fermiFunc = getFermiFunc(evalPoints, chemPotR);
    end
end

function [TotalResult] = choiceLin(Energy, operator, totalSystem, gammaL, gammaR, choice)
    if choice.conservative == true || choice.nonconservative == true || choice.left == true || choice.right == true
        if choice.conservative == true
            midFactor = gammaL + gammaR;
        elseif choice.nonconservative == true
            midFactor = gammaL - gammaR;
        elseif choice.left == true
            midFactor = gammaL;
        elseif choice.right == true
            midFactor = gammaR;
        end
        TotalResult = HelicityAlt(Energy, operator, totalSystem, midFactor);
    else
        ResultL = HelicityAlt(Energy, operator, totalSystem, gammaL);
        ResultR = HelicityAlt(Energy, operator, totalSystem, gammaR);
        TotalResult = ResultL + ResultR;
    end
end

function [] = testFermi(evalPoints, fermiFunc, choice)
    evalFunc = evalPoints(fermiFunc~=0);
    maxPoint = max(evalFunc);
    minPoint = min(evalFunc);
    if choice.conservative == true
        strDisp = 'Conservative';
    elseif choice.nonconservative == true
        strDisp = 'Nonconservative';
    elseif choice.left == true
        strDisp = 'Left';
    elseif choice.right == true
        strDisp = 'Right';
    end
    disp(['Fermi: ', strDisp, '; Maximum: ',num2str(maxPoint), ', Minimum: ',num2str(minPoint)])
end
