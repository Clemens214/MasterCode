function [Result] = transmission(totalSystem, gammaL, gammaR, Eigenvals, leftEVs, rightEVs, chemPot, options)
% calculate the transmission through a molecule for zero temperature
arguments
    totalSystem
    gammaL
    gammaR
    Eigenvals
    leftEVs
    rightEVs
    chemPot
    options.linearResponse = false
end
    %disp('Starting calculation of the current.')
    if options.linearResponse == true
        Energy = chemPot;
        TotalResult = Transmission(Energy, totalSystem, gammaL, gammaR);
    elseif options.linearResponse == false
        chemPotL = chemPot;
        chemPotR = -1*chemPot;
        TotalResult = currentElement(Eigenvals, leftEVs, rightEVs, gammaL, gammaR, chemPotL, chemPotR);
    end
    Result = real(trace(TotalResult));
    %disp('Finished calculation of the current.')
end

%% calculate the Transmission for zero temperature
function [Result] = Transmission(Energy, totalSystem, gammaL, gammaR)
    % calculate the Transport through the molecule

    % calculate the Greens Function
    GreensFuncInv = Energy*eye(length(totalSystem)) - totalSystem;
    GreensFunc = inv(GreensFuncInv);
    
    % calculate the matrix product
    Result = GreensFunc * gammaL * GreensFunc' * gammaR;
end

%% 
function [Result] = currentElement(Eigenvals, leftEVs, rightEVs, gammaL, gammaR, chemPotL, chemPotR)
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

    %disp('Starting calculation of the transmission element.')
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
        ProductLeft = rightEV;
        ProductMid = leftEV * gammaR * leftEVdagger;
        ProductRight = rightEVdagger * gammaL;
            
        Product = ProductLeft * ProductMid * ProductRight;
        
        % compute the additional matrix element
        factor = factorElement(EigVal, EigValDagger, chemPotL) - factorElement(EigVal, EigValDagger, chemPotR);
        
        Result = Result + Product*factor;
    end
    %disp('Finished calculation of the transmission element.')
end

%% calculate the factor for temperatures equal to zero
function [result] = factorElement(eig1, eig2, chemPot)
    %pot = chemPot(1);
    factor = 1/(eig1 -eig2);
    element1 = log(chemPot - eig1);
    element2 = log(chemPot - eig2);
    result = factor*(element1 - element2);
end