function [Results] = AngularCalc(totalSystem, gammaL, gammaR, voltages, sizeLead, choice, options)
% calculate the angular momentum experienced by a molecule for zero temperature
arguments
    totalSystem
    gammaL
    gammaR
    voltages
    sizeLead
    choice.conservative = false
    choice.nonconservative = false
    choice.left = false
    choice.right = false
    options.linearResponse = true
end
    operator = AngularOperator(length(totalSystem), sizeLead);
    %disp('Starting calculation of the angular momentum.')
    if options.linearResponse == true
        Energies = voltages;
        if choice.conservative == true || choice.nonconservative == true || choice.left == true || choice.right == true
            Results = AngularMomentum(Energies, operator, totalSystem, gammaL, gammaR, choice);
        else
            choiceL = choice;
            choiceL.left = true;
            ResultsL =  AngularMomentum(Energies, operator, totalSystem, gammaL, gammaR, choiceL);
            choiceR = choice;
            choiceR.right = true;
            ResultsR =  AngularMomentum(Energies, operator, totalSystem, gammaL, gammaR, choiceR);
            Results = ResultsL + ResultsR;
        end
    elseif options.linearResponse == false
        chemPots = setupPots(voltages);
        if choice.conservative == true || choice.nonconservative == true || choice.left == true || choice.right == true
            Results = integrate(chemPots, operator, totalSystem, gammaL, gammaR, choice);
        else
            choiceL = choice;
            choiceL.left = true;
            ResultsL = integrate(chemPots, operator, totalSystem, gammaL, gammaR, choiceL);
            choiceR = choice;
            choiceR.right = true;
            ResultsR = integrate(chemPots, operator, totalSystem, gammaL, gammaR, choiceR);
            Results = ResultsL + ResultsR;
        end
    end
    %disp('Finished calculation of the angular momentum.')
end

function [chemPots] = setupPots(voltages)
    chemPots = struct('left', [], 'right', []);
    for j = 1:length(voltages)
        chemPotL = voltages(j)/2;
        chemPotR = -1*voltages(j)/2;
        chemPots(j) = struct('left', chemPotL, 'right', chemPotR);
    end
end

function [operator] = AngularOperator(sizeTotal, sizeLead, order)
arguments
    sizeTotal
    sizeLead
    order = 2
end
    Pauli = [0, -1j; 1j, 0];
    % generate the angular momentum operator
    sizeCenter = sizeTotal - 2*sizeLead;
    sizeSample = sizeCenter/order;
    small = eye(sizeSample, sizeSample);
    sample = zeros(sizeCenter, sizeCenter);
    for row = 1:sizeSample
        for column = 1:sizeSample
            if small(row, column) ~= 0
                range = order*(row-1)+1 : order*row;
                sample(range, range) = small(row, column) * Pauli;
            end
        end
    end
    % embed the operator in total Hamiltonian
    operator = zeros(sizeTotal, sizeTotal);
    mid = sizeLead+1 : sizeTotal-sizeLead;
    operator(mid, mid) = sample;
end

%% integrate the angular momentum
function [Results] = integrate(chemPots, operator, totalSystem, gammaL, gammaR, choice, options)
arguments
    chemPots
    operator
    totalSystem
    gammaL
    gammaR
    choice
    options.stepMult = 10
    options.minVal = -3
    options.print = false
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
    stepSize = (1/LCD) / options.stepMult;

    % calculate the transmissions
    evalPoints = makeList(maxPoint, minPoint, stepSize);
    values = AngularMomentum(evalPoints, operator, totalSystem, gammaL, gammaR, choice);
    
    % calculate the integrals
    Results = zeros(1, length(chemPots));
    for i = 1:length(chemPots)
        fermiFunc = choiceFermiFunc(evalPoints, chemPots(i).left, chemPots(i).right, choice);
        %testFermi(evalPoints, fermiFunc, choice);
        % calculate the Result
        yData = fermiFunc .* values;
        if length(evalPoints) > 1
            Results(i) = trapz(evalPoints, yData);
        elseif isscalar(evalPoints)
            Results(i) = 0;
        end
        %disp(['Voltage: ', num2str(chemPots(i).left-chemPots(i).right), ', j=', num2str(i)])
    end
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

%% total angular momentum in the linear transport approximation
function [Results] = AngularMomentum(Energies, operator, totalSystem, gammaL, gammaR, choice)
arguments
    Energies
    operator
    totalSystem
    gammaL
    gammaR
    choice
end
    % calculate the transport matrix and the trace
    Traces = zeros(1, length(Energies));
    parfor i = 1:length(Energies)
        Matrix = choiceLin(Energies(i), operator, totalSystem, gammaL, gammaR, choice);
        Traces(i) = trace(real(Matrix));
    end
    % return the results
    Results = Traces;
end

function [Result] = AngularAlt(Energy, operator, totalSystem, midFactor)
    % A*G*B*Gt
    % F = decomposition(GreensInv,'lu');
    % Y = F \ B;
    % T = A * Y;
    % Z = F' \ T;
    % t = trace(Z);

    % operator * GreensFunc * midFactor * GreensFunc'
    GreensInv = Energy*eye(length(totalSystem)) - totalSystem;
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
        TotalResult = AngularAlt(Energy, operator, totalSystem, midFactor);
    else
        ResultL = AngularAlt(Energy, operator, totalSystem, gammaL);
        ResultR = AngularAlt(Energy, operator, totalSystem, gammaR);
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