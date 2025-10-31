function [Results] = TorqueCalc(totalSystem, totalSysDeriv, gammaL, gammaR, chemPots, choice, options)
% calculate the torque through a molecule for zero temperature
arguments
    totalSystem
    totalSysDeriv
    gammaL
    gammaR
    chemPots
    choice.conservative = false
    choice.nonconservative = false
    choice.left = false
    choice.right = false
    options.linearResponse = false
end
    %disp('Starting calculation of the torque.')
    if options.linearResponse == true
        Energies = getEnergies(chemPots);
        if choice.conservative == true || choice.nonconservative == true || choice.left == true || choice.right == true
            Results = Torque(Energies, totalSystem, totalSysDeriv, gammaL, gammaR, choice);
        else
            choiceL = choice;
            choiceL.left = true;
            ResultsL = Torque(Energies, totalSystem, totalSysDeriv, gammaL, gammaR, choiceL);
            choiceR = choice;
            choiceR.right = true;
            ResultsR = Torque(Energies, totalSystem, totalSysDeriv, gammaL, gammaR, choiceR);
            Results = ResultsL + ResultsR;
        end
    elseif options.linearResponse == false
        if choice.conservative == true || choice.nonconservative == true || choice.left == true || choice.right == true
            Results = integrate(chemPots, totalSystem, totalSysDeriv, gammaL, gammaR, choice);
        else
            choiceL = choice;
            choiceL.left = true;
            ResultsL = integrate(chemPots, totalSystem, totalSysDeriv, gammaL, gammaR, choiceL);
            choiceR = choice;
            choiceR.right = true;
            ResultsR = integrate(chemPots, totalSystem, totalSysDeriv, gammaL, gammaR, choiceR);
            Results = ResultsL + ResultsR;
        end
    end
    %disp('Finished calculation of the torque.')
end

%% integrate the torque
function [Results] = integrate(chemPots, totalSystem, totalSysDeriv, gammaL, gammaR, choice, options)
arguments
    chemPots
    totalSystem
    totalSysDeriv
    gammaL
    gammaR
    choice
    options.stepMult = 10
    options.minVal = -3
end
    % get the bounds
    if choice.conservative == true
        maxPoint = max([[chemPots.left], [chemPots.right]]);
        minPoint = min([[chemPots.left], [chemPots.right], options.minVal]);
    elseif choice.nonconservative == true
        maxPoint = max([[chemPots.left], [chemPots.right]]);
        minPoint = min([[chemPots.left], [chemPots.right]]);
    elseif choice.left == true
        maxPoint = max([chemPots.left]);
        minPoint = min([[chemPots.left], options.minVal]);
    elseif choice.right == true
        maxPoint = max([chemPots.right]);
        minPoint = min([[chemPots.right], options.minVal]);
    end
    % get the step size
    Energies = getEnergies(chemPots);
    Diffs = zeros(1, length(Energies)-1);
    for i = 2:length(Energies)
        Diffs(i) = Energies(i) - Energies(i-1); 
    end
    lcd = LowestCommonDenominator(Diffs);
    stepSize = (1/lcd) / options.stepMult;

    % calculate the transmissions
    evalPoints = makeList(maxPoint, minPoint, stepSize);
    values = Torque(evalPoints, totalSystem, totalSysDeriv, gammaL, gammaR, choice);
    
    % calculate the integrals
    Results = zeros(1, length(chemPots));
    for i = 1:length(chemPots)
        fermiFunc = choiceFermiFunc(evalPoints, chemPots(i).left, chemPots(i).right, choice);
        % calculate the Result
        yData = fermiFunc .* values;
        if length(evalPoints) > 1
            Results(i) = trapz(evalPoints, yData);
        elseif isscalar(evalPoints)
            Results(i) = 0;
        end
        disp(['Voltage: ', num2str(chemPots(i).left-chemPots(i).right), ', j=', num2str(i)])
    end
end

function [fermiFunc] = choiceFermiFunc(evalPoints, chemPotL, chemPotR, choice)
    if choice.conservative == true
        fermiFuncL = getFermiFunc(evalPoints, chemPotL);
        fermiFuncR = getFermiFunc(evalPoints, chemPotR);
        fermiFunc = fermiFuncL + fermiFuncR;
    elseif choice.nonconservative == true
        fermiFuncL = getFermiFunc(evalPoints, chemPotL);
        fermiFuncR = getFermiFunc(evalPoints, chemPotR);
        fermiFunc = fermiFuncL - fermiFuncR;
    elseif choice.left == true
        fermiFunc = getFermiFunc(evalPoints, chemPotL);
    elseif choice.right == true
        fermiFunc = getFermiFunc(evalPoints, chemPotR);
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
            if E <= chemPot
                fermiFunc(i) = 1;
            else
                fermiFunc(i) = 0;
            end
        end
    end
end

%% total torque in the linear transport approximation
function [Results] = Torque(Energies, totalSystem, totalSysDeriv, gammaL, gammaR, choice)
    %calculates the transport through a molecule in the linear transport approximation
    arguments
        Energies
        totalSystem
        totalSysDeriv
        gammaL
        gammaR
        choice
    end
    % calculate the transport matrix and the trace
    Traces = zeros(1, length(Energies));
    parfor i = 1:length(Energies)
        Matrix = choiceLin(Energies(i), totalSystem, totalSysDeriv, gammaL, gammaR, choice);
        Traces(i) = trace(real(Matrix));
    end
    % return the results
    Results = Traces;
end

function [TotalResult] = choiceLin(Energy, totalSystem, totalSysDeriv, gammaL, gammaR, choice)
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
        TotalResult = TorqueAlt(Energy, totalSystem, totalSysDeriv, midFactor);
    else
        ResultL = TorqueAlt(Energy, totalSystem, totalSysDeriv, gammaL);
        ResultR = TorqueAlt(Energy, totalSystem, totalSysDeriv, gammaR);
        TotalResult = ResultL + ResultR;
    end
end

function [Result] = TorqueAlt(Energy, totalSystem, totalSysDeriv, midFactor)
    % A*G*B*Gt
    % F = decomposition(GreensInv,'lu');
    % Y = F \ B;
    % T = A * Y;
    % Z = F' \ T;
    % t = trace(Z);

    % totalSysDeriv * GreensFunc * midFactor * GreensFunc'
    GreensInv = Energy*eye(length(totalSystem)) - totalSystem;
    F = decomposition(GreensInv,'lu');    % create reusable factorization object
    
    Y = F \ midFactor;                     % solves Aw * Y = B
    T = totalSysDeriv * Y;
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