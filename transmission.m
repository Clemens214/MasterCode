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
    %disp('Starting calculation of the current element.')
    Result = 0;
    for i = 1:length(leftEVs)
        % get the normal left and right Eigenvectors
        EigVal = Eigenvals(i,i);
        leftEV = leftEVs(:,i)';
        rightEV = rightEVs(:,i);
        
        ProductLeft = rightEV;
        ProductMidLeft = leftEV * gammaR;
        for j = 1:length(leftEVs)
            % get the daggered Eigenvectors
            EigValDagger = Eigenvals(j,j)';
            leftEVdagger = leftEVs(:,j);
            rightEVdagger = rightEVs(:,j)';
            
            % compute the matrix element for chosen i and j
            ProductMid = ProductMidLeft * leftEVdagger;
            ProductRight = rightEVdagger * gammaL;
            
            Product = ProductLeft * ProductMid * ProductRight;
            
            % compute the additional matrix element
            factor = factorElement(EigVal, EigValDagger, chemPotL) - factorElement(EigVal, EigValDagger, chemPotR);
            
            Result = Result + Product*factor;
        end
    end
    %disp('Finished calculation of the current element.')
end

%% calculate the factor for temperatures equal to zero
function [result] = factorElement(eig1, eig2, chemPot)
    %pot = chemPot(1);
    factor = 1/(eig1 -eig2);
    element1 = log(chemPot - eig1);
    element2 = log(chemPot - eig2);
    result = factor*(element1 - element2);
end