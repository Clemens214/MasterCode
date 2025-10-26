function [Results] = TorqueCalc(totalSystem, totalSysDeriv, gammaL, gammaR, Eigenvals, leftEVs, rightEVs, chemPots, options)
% calculate the torque through a molecule for zero temperature
arguments
    totalSystem
    totalSysDeriv
    gammaL
    gammaR
    Eigenvals
    leftEVs
    rightEVs
    chemPots
    options.linearResponse = false
    options.conservative = false
    options.nonconservative = false
    options.left = false
    options.right = false
end
    choice = struct('conservative', options.conservative, ...
                    'nonconservative', options.nonconservative, ...
                    'left', options.left, 'right', options.right);
    %disp('Starting calculation of the torque.')
    if options.linearResponse == true
        Energies = getEnergies(chemPots);
        Results = Torque(Energies, totalSystem, totalSysDeriv, gammaL, gammaR, choice);
    elseif options.linearResponse == false
        Results = zeros(1, length(chemPots));
        for i = 1:length(chemPots)
            chemPotL = chemPots(i).left;
            chemPotR = chemPots(i).right;
            TotalResult = TorqueChoice(Eigenvals, leftEVs, rightEVs, totalSysDeriv, gammaL, gammaR, chemPotL, chemPotR, choice);
            Results(i) = real(trace(TotalResult));

            voltage = chemPotL - chemPotR;
            disp(['Voltage: ', num2str(voltage), ', j=', num2str(i)])
        end
    end
    %disp('Finished calculation of the torque.')
end

%% total torque for finite voltages
function [TotalResult] = TorqueChoice(Eigenvals, leftEVs, rightEVs, totalSysDeriv, gammaL, gammaR, chemPotL, chemPotR, choice)
    if choice.conservative == true
        midFactor = gammaL + gammaR;
        TotalResult = TorqueMatrix(Eigenvals, leftEVs, rightEVs, totalSysDeriv, midFactor, chemPotL, chemPotR, choice);
    elseif choice.nonconservative == true
        midFactor = gammaL + gammaR;
        TotalResult = TorqueMatrix(Eigenvals, leftEVs, rightEVs, totalSysDeriv, midFactor, chemPotL, chemPotR, choice);
    elseif choice.left == true
        midFactor = gammaL;
        TotalResult = TorqueMatrix(Eigenvals, leftEVs, rightEVs, totalSysDeriv, midFactor, chemPotL, chemPotR, choice);
    elseif choice.right == true
        midFactor = gammaR;
        TotalResult = TorqueMatrix(Eigenvals, leftEVs, rightEVs, totalSysDeriv, midFactor, chemPotL, chemPotR, choice);
    else
        choiceL = choice;
        choiceL.left = true;
        ResultL = TorqueMatrix(Eigenvals, leftEVs, rightEVs, totalSysDeriv, gammaL, chemPotL, chemPotR, choiceL);
        choiceR = choice;
        choiceR.right = true;
        ResultR = TorqueMatrix(Eigenvals, leftEVs, rightEVs, totalSysDeriv, gammaR, chemPotL, chemPotR, choiceR);
        TotalResult = ResultL + ResultR;
    end
end

function [Result] = TorqueMatrix(Eigenvals, leftEVs, rightEVs, totalSysDeriv, midFactor, chemPotL, chemPotR, choice)
    index = struct('i', [], 'j', [], ...
                    'Eigenval', [], 'EigenvalD', [], ...
                    'leftEV', [], 'leftEVD', [], ...
                    'rightEV', [], 'rightEVD', []);
    for i = 1:length(Eigenvals)
        for j = 1:length(Eigenvals)
            idx = (i-1)*length(Eigenvals) + j;
            index(idx).i = i;
            index(idx).j = j;
            index(idx).Eigenval = Eigenvals(i,i);
            index(idx).leftEV = leftEVs(:,i)';
            index(idx).rightEV = rightEVs(:,i);
            index(idx).EigenvalD = Eigenvals(j,j)';
            index(idx).leftEVD = leftEVs(:,j);
            index(idx).rightEVD = rightEVs(:,j)';
        end
    end

    %disp('Starting calculation of the torque element.')
    Result = 0;
    parfor idx = 1:length(index)
        % get the normal left and right Eigenvectors
        EigVal = index(idx).Eigenval;
        leftEV = index(idx).leftEV;
        rightEV = index(idx).rightEV;
        
        % get the daggered Eigenvectors
        EigValDagger = index(idx).EigenvalD;
        leftEVdagger = index(idx).leftEVD;
        rightEVdagger = index(idx).rightEVD;
            
        % compute the matrix element for chosen i and j
        ProductLeft = totalSysDeriv * rightEV;
        ProductMid = leftEV * midFactor * leftEVdagger;
        ProductRight = rightEVdagger;
        
        Product = ProductLeft * ProductMid * ProductRight;
        
        % compute the additional matrix element
        factor = choiceFactor(EigVal, EigValDagger, chemPotL, chemPotR, choice);
        
        Result = Result + Product*factor;
    end
    %disp('Finished calculation of the torque element.')
end

%% calculate the factor
function [Factor] = choiceFactor(EigVal, EigValDagger, chemPotL, chemPotR, choice)
    if choice.conservative == true
        Factor = factorElement(EigVal, EigValDagger, chemPotL) + factorElement(EigVal, EigValDagger, chemPotR);
    elseif choice.nonconservative == true
        Factor = factorElement(EigVal, EigValDagger, chemPotL) - factorElement(EigVal, EigValDagger, chemPotR);
    elseif choice.left == true
        Factor = factorElement(EigVal, EigValDagger, chemPotL);
    elseif choice.right == true
        Factor = factorElement(EigVal, EigValDagger, chemPotR);
    end
end

function [result] = factorElement(eig1, eig2, chemPot)
    %pot = chemPot(1);
    factor = 1/(eig1 -eig2);
    element1 = log(chemPot - eig1);
    element2 = log(chemPot - eig2);
    result = factor*(element1 - element2);
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
    for i = 1:length(Energies)
        Matrix = choiceLin(Energies(i), totalSystem, totalSysDeriv, gammaL, gammaR, choice);
        Traces(i) = trace(real(Matrix));
        
        disp(['Energy: ', num2str(Energies(i)), ', j=', num2str(i)])
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
        %TotalResult = TorqueZeroTemp(Energy, totalSystem, totalSysDeriv, midFactor);
        TotalResult = TorqueAlt(Energy, totalSystem, totalSysDeriv, midFactor);
    else
        %ResultL = TorqueZeroTemp(Energy, totalSystem, totalSysDeriv, gammaL);
        ResultL = TorqueAlt(Energy, totalSystem, totalSysDeriv, gammaL);
        %ResultR = TorqueZeroTemp(Energy, totalSystem, totalSysDeriv, gammaR);
        ResultR = TorqueAlt(Energy, totalSystem, totalSysDeriv, gammaR);
        TotalResult = ResultL + ResultR;
    end
end

function [Result] = TorqueZeroTemp(Energy, totalSystem, totalSysDeriv, midFactor)
    % calculate the Greens Function
    GreensFuncInv = Energy*eye(length(totalSystem)) - totalSystem;
    GreensFunc = inv(GreensFuncInv);
    
    % calculate the matrix product
    Result = totalSysDeriv * GreensFunc * midFactor * GreensFunc';
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