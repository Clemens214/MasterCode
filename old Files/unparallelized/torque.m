function [Result, varargout] = torque(totalSystem, Eigenvals, leftEVs, rightEVs, totalSysDeriv, gammaL, gammaR, chemPots)
    %calculates the torque experienced by a molecule as calculated in the
    %PhD-Thesis

    TorqueTotal = 0;
    TotalMatrix = zeros(length(Eigenvals), length(Eigenvals));
    %disp('Starting calculation of the torque.')
    gammas = {gammaL, gammaR};
    for i = 1:length(gammas)
        %disp(['Index, i=',num2str(i)])
        gamma = gammas{i};
        chemPot = chemPots{i};
        
        % determine whether the calculation is to be done using linear response
        linearResponse = false;
        if linearResponse == false
            TorqueMatrix = Torque(Eigenvals, leftEVs, rightEVs, totalSysDeriv, gamma, chemPot);
        else
            TorqueMatrix = TorqueZeroTemp(totalSystem, totalSysDeriv, gamma);
        end
        TorqueElement = trace(TorqueMatrix);
        
        TotalMatrix = TotalMatrix + TorqueMatrix;
        TorqueTotal = TorqueTotal + TorqueElement;
    end
    %disp('Finished calculation of the torque.')
    Result = real(TorqueTotal);
    varargout{1} = real(TotalMatrix);
end

%% transmission at finite temperature
function [Result] = Torque(Eigenvals, leftEVs, rightEVs, totalSysDeriv, gamma, chemPot)
    %disp('Starting calculation of the current element.')
    Result = 0;
    for i = 1:length(leftEVs)
        % get the normal left and right Eigenvectors
        EigVal = Eigenvals(i,i);
        leftEV = leftEVs(:,i)';
        rightEV = rightEVs(:,i);
        
        ProductLeft = totalSysDeriv * rightEV;
        ProductMidLeft = leftEV * gamma;
        for j = 1:length(leftEVs)
            % get the daggered Eigenvectors
            EigValDagger = Eigenvals(j,j)';
            leftEVdagger = leftEVs(:,j);
            rightEVdagger = rightEVs(:,j)';
            
            % compute the matrix element for chosen i and j
            ProductMid = ProductMidLeft * leftEVdagger;
            ProductRight = rightEVdagger;
            
            Product = ProductLeft * ProductMid * ProductRight;
            
            % compute the additional matrix element
            factor = matrixElement(EigVal, EigValDagger, chemPot);
            
            Result = Result + Product*factor;
        end
    end
    %disp('Finished calculation of the current element.')
end

function [result] = matrixElement(eig1, eig2, chemPot)
    factor = 1/(eig1 -eig2);
    element1 = log(chemPot - eig1);
    element2 = log(chemPot - eig2);
    result = factor*(element1 - element2);
end

%% transmission at zero temperature
function [Transmission] = TorqueZeroTemp(totalSystem, totalSysDeriv, gamma)
    % calculate the Greens Function
    GreensFuncInv = eye(length(totalSystem)) - totalSystem;
    GreensFunc = inv(GreensFuncInv);
    
    % calculate the matrix product
    Transmission = totalSysDeriv * GreensFunc * gamma * GreensFunc';
end