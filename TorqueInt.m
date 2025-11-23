function [Results] = TorqueInt(totalSystem, totalSysDeriv, gammaL, gammaR, voltages, choice, options)
% calculate the torque through a molecule for zero temperature
arguments
    totalSystem
    totalSysDeriv
    gammaL
    gammaR
    voltages
    choice.conservative = false
    choice.nonconservative = false
    choice.left = false
    choice.right = false
    options.Schur = false
    options.linearResponse = true
    options.print = false
end
    %disp('Starting calculation of the torque.')
    if options.linearResponse == false
        if options.Schur == false
            % compute the Eigenvectors and the Eigenvalues of the system
            [Eigenvals, leftEVs, rightEVs] = getEigenvectors(totalSystem);%, checkMore=true);
        elseif options.Schur == true
            % compute the Schur decomposition of the System's pseudo Hamiltonian
            [Diag, upperTriag, SchurVec] = getSchur(totalSystem, options);
        end
        chemPots = setupPots(voltages);
        Results = zeros(1, length(chemPots));
        for i = 1:length(chemPots)
            chemPotL = chemPots(i).left;
            chemPotR = chemPots(i).right;
            if options.Schur == true
                Matrix = TorqueChoice(Diag, upperTriag, SchurVec, totalSysDeriv, gammaL, gammaR, chemPotL, chemPotR, choice);
            elseif options.Schur == false
                Matrix = TorqueChoice(Eigenvals, leftEVs, rightEVs, totalSysDeriv, gammaL, gammaR, chemPotL, chemPotR, choice);
            end
            Results(i) = trace(real(Matrix));
            if options.print == true
                disp(['Voltage: ', num2str(chemPotL - chemPotR), ', j=', num2str(i)])
            end
        end
    elseif options.linearResponse == true
        Energies = voltages;
        Results = zeros(1, length(Energies));
        for i = 1:length(Energies)
            Matrix = TorqueLin(Energies(i), totalSystem, totalSysDeriv, gammaL, gammaR, choice);
            Results(i) = trace(real(Matrix));
            if options.print == true
                disp(['Energy: ', num2str(Energies(i)), ', j=', num2str(i)])
            end
        end
    end
    %disp('Finished calculation of the torque.')
end

%% total torque for finite voltages
function [TotalResult] = TorqueChoice(Val1, Val2, Val3, Val4, gammaL, gammaR, chemPotL, chemPotR, choice, options)
arguments
    Val1
    Val2
    Val3
    Val4
    gammaL
    gammaR
    chemPotL
    chemPotR
    choice
    options.Schur = true
end
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
        if options.Schur == true
            TotalResult = TorqueMatrixSchur(Val1, Val2, Val3, Val4, midFactor, chemPotL, chemPotR, choice);
        else
            TotalResult = TorqueMatrixEV(Val1, Val2, Val3, Val4, midFactor, chemPotL, chemPotR, choice);
        end
    else
        choiceL = choice;
        choiceL.left = true;
        choiceR = choice;
        choiceR.right = true;
        if options.Schur == true
            ResultL = TorqueMatrixSchur(Val1, Val2, Val3, Val4, gammaL, chemPotL, chemPotR, choiceL);
            ResultR = TorqueMatrixSchur(Val1, Val2, Val3, Val4, gammaR, chemPotL, chemPotR, choiceR);
        else
            ResultL = TorqueMatrixEV(Val1, Val2, Val3, Val4, gammaL, chemPotL, chemPotR, choiceL);
            ResultR = TorqueMatrixEV(Val1, Val2, Val3, Val4, gammaR, chemPotL, chemPotR, choiceR);
        end
        TotalResult = ResultL + ResultR;
    end
end

function [Result] = TorqueMatrixSchur(Eigenvals, leftEVs, rightEVs, totalSysDeriv, midFactor, chemPotL, chemPotR, choice)
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
        factor = FactorChoice(EigVal, EigValDagger, chemPotL, chemPotR, choice);
        
        Result = Result + Product*factor;
    end
    %disp('Finished calculation of the torque element.')
end

function [Result] = TorqueMatrixEV(Eigenvals, leftEVs, rightEVs, totalSysDeriv, midFactor, chemPotL, chemPotR, choice)
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
        factor = FactorChoice(EigVal, EigValDagger, chemPotL, chemPotR, choice);
        
        Result = Result + Product*factor;
    end
    %disp('Finished calculation of the torque element.')
end

%% calculate the factor
function [Factor] = FactorChoice(EigVal, EigValDagger, chemPotL, chemPotR, choice)
    if choice.conservative == true
        Factor = FactorElement(EigVal, EigValDagger, chemPotL) + FactorElement(EigVal, EigValDagger, chemPotR);
    elseif choice.nonconservative == true
        Factor = FactorElement(EigVal, EigValDagger, chemPotL) - FactorElement(EigVal, EigValDagger, chemPotR);
    elseif choice.left == true
        Factor = FactorElement(EigVal, EigValDagger, chemPotL);
    elseif choice.right == true
        Factor = FactorElement(EigVal, EigValDagger, chemPotR);
    end
end

function [result] = FactorElement(eig1, eig2, chemPot)
    %pot = chemPot(1);
    factor = 1/(eig1 -eig2);
    element1 = log(chemPot - eig1);
    element2 = log(chemPot - eig2);
    result = factor*(element1 - element2);
end

%% total torque in the linear transport approximation
function [TotalResult] = TorqueLin(Energy, totalSystem, totalSysDeriv, gammaL, gammaR, choice)
%calculates the transport through a molecule in the linear transport approximation
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
        TotalResult = TorqueZeroTemp(Energy, totalSystem, totalSysDeriv, midFactor);
    else
        ResultL = TorqueZeroTemp(Energy, totalSystem, totalSysDeriv, gammaL);
        ResultR = TorqueZeroTemp(Energy, totalSystem, totalSysDeriv, gammaR);
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

%% helping functions
function [chemPots] = setupPots(voltages)
    chemPots = struct('left', [], 'right', []);
    for j = 1:length(voltages)
        chemPotL = voltages(j)/2;
        chemPotR = -1*voltages(j)/2;
        chemPots(j) = struct('left', chemPotL, 'right', chemPotR);
    end
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