function [Results] = TransInt(totalSystem, gammaL, gammaR, Eigenvals, leftEVs, rightEVs, chemPots, options)
% calculate the transmission through a molecule for zero temperature
arguments
    totalSystem
    gammaL
    gammaR
    Eigenvals
    leftEVs
    rightEVs
    chemPots
    options.linearResponse = false
end
    %disp('Starting calculation of the current.')
    Results = zeros(1, length(chemPots));
    for i = 1:length(chemPots)
        chemPotL = chemPots(i).left;
        chemPotR = chemPots(i).right;
        TotalResult = TransmissionMatrix(Eigenvals, leftEVs, rightEVs, gammaL, gammaR, chemPotL, chemPotR);
        Results(i) = real(trace(TotalResult));
        
        voltage = chemPotL - chemPotR;
        disp(['Voltage: ', num2str(voltage), ', j=', num2str(i)])
    end
    %disp('Finished calculation of the current.')
end

%% total transmission for finite voltages
function [Result] = TransmissionMatrix(Eigenvals, leftEVs, rightEVs, gammaL, gammaR, chemPotL, chemPotR)
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
        factor = factorElement(EigVal, EigValDagger, chemPotL) - factorElement(EigVal, EigValDagger, chemPotR);
        
        Result = Result + Product*factor;
    end
    %disp('Finished calculation of the transmission element.')
end

%% calculate the factor
function [result] = factorElement(eig1, eig2, chemPot)
    if eig1 ~= eig2
        factor = 1/(eig1 - eig2);
        element1 = log(chemPot - eig1);
        element2 = log(chemPot - eig2);
        result = factor*(element1 - element2);
    else
        result = -1/(chemPot - eig1);
    end
end