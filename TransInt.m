function [Results] = TransInt(totalSystem, gammaL, gammaR, voltages, options)
% calculate the transmission through a molecule for zero temperature
arguments
    totalSystem
    gammaL
    gammaR
    voltages
    options.linearResponse = true
    options.print = false
end
    %disp('Starting calculation of the current.')
    if options.linearResponse == false
        chemPots = setupPots(voltages);
        % compute the Schur decomposition of the System's pseudo Hamiltonian
        [Diag, upperTriag, SchurVec] = getSchur(totalSystem);
        Results = zeros(1, length(chemPots));
        for i = 1:length(chemPots)
            chemPotL = chemPots(i).left;
            chemPotR = chemPots(i).right;
            Matrix = TransmissionMatrix(Diag, upperTriag, SchurVec, gammaL, gammaR, chemPotL, chemPotR);
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

%% total transmission for finite voltages using Schur
function [Result] = TransmissionMatrix(Diag, upperTriag, SchurVec, gammaL, gammaR, chemPotL, chemPotR)
    index = struct('index1', [], 'index2', [], ...
                    'Eigenval', [], 'EigenvalD', [], ...
                    'leftEV', [], 'leftEVD', [], ...
                    'rightEV', [], 'rightEVD', []);
    idx = 0;
    for i = 1:length(Diag)
    for j = i:length(Diag)
        for k = 1:length(Diag)
        for l = k:length(Diag)
            idx = idx + 1;
            index(idx).index1 = [i, j];
            index(idx).index1 = [k, l];
            % set the normal variables
            for m = 1:length(Diag)-1
                index(idx).Eigenval = Eigenvals(i,i);
                index(idx).leftEV = leftEVs(:,i)';
                index(idx).rightEV = rightEVs(:,i);
            end
            % set the daggered variables
            for m = 1:length(Diag)-1
                index(idx).EigenvalD = Eigenvals(j,j)';
                index(idx).leftEVD = leftEVs(:,j);
                index(idx).rightEVD = rightEVs(:,j)';
            end
        end
        end
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

function [] = getElement(Diag, upperTriag, SchurVec)
arguments
    Diag
    upperTriag
    SchurVec
end
    index = struct('row', [], 'column', [], 'power', [], ...
                    'Eigenval', [], 'EigenvalD', [], ...
                    'leftEV', [], 'leftEVD', [], ...
                    'rightEV', [], 'rightEVD', []);
    idx = 0;
    for row = 1:length(Diag)
        for column = row:length(Diag)
            for power = 0:length(Diag)-1
                idx = idx + 1;
                % set the indices
                index(idx).row = row;
                index(idx).column = column;
                index(idx).power = power;
                % set the variables
                index(idx).path = [];
                index(idx).leftEV = leftEVs(:,i)';
                index(idx).rightEV = rightEVs(:,i);
            end
        end
    end
end

function [] = getPath(row, column, power)
    for row = 1:length(size)
            
        end
    end
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
%calculates the transport through a molecule in the linear transport approximation
arguments
    Energy
    totalSystem
    gammaL
    gammaR
end
    % G*B*Gt*C
    % GreensFunc * gammaL * GreensFunc' * gammaR
    GreensInv = Energy*eye(length(totalSystem)) - totalSystem;
    % F = decomposition(GreensInv,'lu');
    F = decomposition(GreensInv,'lu');  % create reusable LU object (works for sparse/dense)
    % Y = F \ B;
    Y = F \ gammaL;
    % W = C * Y;
    W = gammaR * Y;
    % Z = F' \ W;
    Z = F' \ W;                         % uses transpose of factorization
    % t = trace(Z);
    Result = Z;
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