function [Results] = TransInt(totalSystem, gammaL, gammaR, voltages, options)
% calculate the transmission through a molecule for zero temperature
arguments
    totalSystem
    gammaL
    gammaR
    voltages
    options.Schur = true
    options.linearResponse = true
    options.print = false
end
    %disp('Starting calculation of the current.')
    if options.linearResponse == false
        chemPots = setupPots(voltages);
        if options.Schur == false
            % compute the Eigenvectors and the Eigenvalues of the system
            [Eigenvals, leftEVs, rightEVs] = getEigenvectors(totalSystem);%, checkMore=true);
        elseif options.Schur == true
            % compute the Schur decomposition of the System's pseudo Hamiltonian
            [Diag, upperTriag, SchurVec] = getSchur(totalSystem);
        end
        Results = zeros(1, length(chemPots));
        for i = 1:length(chemPots)
            chemPotL = chemPots(i).left;
            chemPotR = chemPots(i).right;
            if options.Schur == false
                Matrix = TransmissionMatrixEV(Eigenvals, leftEVs, rightEVs, gammaL, gammaR, chemPotL, chemPotR);
            elseif options.Schur == true
                Matrix = TransmissionMatrixSchur(Diag, upperTriag, SchurVec, gammaL, gammaR, chemPotL, chemPotR);
            end
            Results(i) = real(trace(Matrix));
            if options.print == true
                disp(['Voltage: ', num2str(chemPotL - chemPotR), ', j=', num2str(i)])
            end
        end
    elseif options.linearResponse == true
        Energies = voltages;
        Results = zeros(1, length(Energies));
        for i = 1:length(Energies)
            Matrix = TransmissionLin(Energies(i), totalSystem, gammaL, gammaR);
            Results(i) = trace(real(Matrix));
            if options.print == true
                disp(['Energy: ', num2str(Energies(i)), ', j=', num2str(i)])
            end
        end
    end
    %disp('Finished calculation of the current.')
end

%% total transmission for finite voltages
function [Result] = TransmissionMatrixSchur(Diag, upperTriag, SchurVec, gammaL, gammaR, chemPotL, chemPotR)
    index = struct('i', [], 'j', [], ...
                    'Eigenval', [], 'EigenvalD', [], ...
                    'leftEV', [], 'leftEVD', [], ...
                    'rightEV', [], 'rightEVD', []);
    for i = 1:length(Eigenvals)
        for j = 1:length(Eigenvals)
            idx = (i-1)*length(Eigenvals) + j;
            index(idx).i = i;
            index(idx).j = j;
            % set the normal variables
            index(idx).Eigenval = Eigenvals(i,i);
            index(idx).leftEV = leftEVs(:,i)';
            index(idx).rightEV = rightEVs(:,i);
            % set the daggered variables
            index(idx).EigenvalD = Eigenvals(j,j)';
            index(idx).leftEVD = leftEVs(:,j);
            index(idx).rightEVD = rightEVs(:,j)';
        end
    end

    %disp('Starting calculation of the transmission element.')
    Result = 0;
    parfor idx = 1:length(index)
        % get the normal left and right Eigenvectors
        leftEV = index(idx).leftEV;
        rightEV = index(idx).rightEV;
        
        % get the daggered Eigenvectors
        leftEVdagger = index(idx).leftEVD;
        rightEVdagger = index(idx).rightEVD;
            
        % compute the matrix element for chosen i and j
        ProductLeft = rightEV;
        ProductMid = leftEV * gammaR * leftEVdagger;
        ProductRight = rightEVdagger * gammaL;
            
        Product = ProductLeft * ProductMid * ProductRight;
        
        % compute the additional matrix element
        EigVal = index(idx).Eigenval;
        EigValDagger = index(idx).EigenvalD;
        factor = FactorElement(EigVal, EigValDagger, chemPotL) - FactorElement(EigVal, EigValDagger, chemPotR);
        
        Result = Result + Product*factor;
    end
    %disp('Finished calculation of the transmission element.')
end

function [Result] = TransmissionMatrixEV(Eigenvals, leftEVs, rightEVs, gammaL, gammaR, chemPotL, chemPotR)
    index = struct('i', [], 'j', [], ...
                    'Eigenval', [], 'EigenvalD', [], ...
                    'leftEV', [], 'leftEVD', [], ...
                    'rightEV', [], 'rightEVD', []);
    for i = 1:length(Eigenvals)
        for j = 1:length(Eigenvals)
            idx = (i-1)*length(Eigenvals) + j;
            index(idx).i = i;
            index(idx).j = j;
            % set the normal variables
            index(idx).Eigenval = Eigenvals(i,i);
            index(idx).leftEV = leftEVs(:,i)';
            index(idx).rightEV = rightEVs(:,i);
            % set the daggered variables
            index(idx).EigenvalD = Eigenvals(j,j)';
            index(idx).leftEVD = leftEVs(:,j);
            index(idx).rightEVD = rightEVs(:,j)';
        end
    end

    %disp('Starting calculation of the transmission element.')
    Result = 0;
    parfor idx = 1:length(index)
        % get the normal left and right Eigenvectors
        leftEV = index(idx).leftEV;
        rightEV = index(idx).rightEV;
        
        % get the daggered Eigenvectors
        leftEVdagger = index(idx).leftEVD;
        rightEVdagger = index(idx).rightEVD;
            
        % compute the matrix element for chosen i and j
        ProductLeft = rightEV;
        ProductMid = leftEV * gammaR * leftEVdagger;
        ProductRight = rightEVdagger * gammaL;
            
        Product = ProductLeft * ProductMid * ProductRight;
        
        % compute the additional matrix element
        EigVal = index(idx).Eigenval;
        EigValDagger = index(idx).EigenvalD;
        factor = FactorElement(EigVal, EigValDagger, chemPotL) - FactorElement(EigVal, EigValDagger, chemPotR);
        
        Result = Result + Product*factor;
    end
    %disp('Finished calculation of the transmission element.')
end

%% calculate the factor
function [result] = FactorElement(eig1, eig2, chemPot)
    if eig1 ~= eig2
        factor = 1/(eig1 - eig2);
        element1 = log(chemPot - eig1);
        element2 = log(chemPot - eig2);
        result = factor*(element1 - element2);
    else
        result = -1/(chemPot - eig1);
    end
end

%% total transmission in the linear transport approximation
function [Result] = TransmissionLin(Energy, totalSystem, gammaL, gammaR)
    % calculate the Greens Function
    GreensFuncInv = Energy*eye(length(totalSystem)) - totalSystem;
    GreensFunc = inv(GreensFuncInv);
    
    % calculate the matrix product
    Result = GreensFunc * gammaL * GreensFunc' * gammaR;
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