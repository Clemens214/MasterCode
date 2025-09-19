function [ResultTotal] = torque(totalSystem, totalSysDeriv, gammaL, gammaR, Eigenvals, leftEVs, rightEVs, Temp)
    %calculates the torque experienced by a molecule
    
    TempVal = false;
    %disp('Starting calculation of the torque.')
    ResultTotal = 0;
    gammas = [gammaL, gammaR];
    for i = 1:length(gammas)
        gamma = gammas(i);
        if Temp == 0
            TransmissionMatrix = zeroTempTransmission(totalSystem, totalSysDeriv, gamma);
        else
            TransmissionMatrix = Transmission(Eigenvals, leftEVs, rightEVs, totalSysDeriv, gamma, Temp);
        end
        Result = trace(TransmissionMatrix);
        
        ResultTotal = ResultTotal + Result;
    end

    if Temp ~= 0 && TempVal == true
        TempResult = ResultTotal/Temp;
    end
    %disp('Finished calculation of the torque.')
end

%% transmission at finite temperature
function [Result] = Transmission(Eigenvals, leftEVs, rightEVs, totalSysDeriv, gamma, Temp)
    %disp('Starting calculation of the Hurwitz Zeta values.')
    Hurwitz = zeros(length(Eigenvals));
    HurwitzDagger = zeros(length(Eigenvals));
    for k = 1:length(Eigenvals)
        EigVal = Eigenvals(k,k);
        EigValDagger = Eigenvals(k,k)';
        % compute the Hurwitz Zeta functions
        beta = 1/Temp;
        Hurwitz(k) = HurwitzZeta(EigVal, 1/beta);
        HurwitzDagger(k) = HurwitzZeta(EigValDagger, 1/beta);
    end
    %disp('Finished calculation of the Hurwitz Zeta values.')
    
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
            factor = matrixElement(EigVal, EigValDagger, Hurwitz(i), HurwitzDagger(j));
            
            Result = Result + Product*factor;
        end
    end
    %disp('Finished calculation of the current element.')
end

function [result] = matrixElement(eig1, eig2, HurwitzEig1, HurwitzEig2)
    factor = 1/(eig1 -eig2);
    result = factor*HurwitzEig1 - factor*HurwitzEig2;
end

function [result] = HurwitzZeta (x, beta)
    factor = sign(imag(x))/(2*pi*1j*beta);
    value = 0.5+factor*x;
    result = -1*factor * hurwitzZeta(1, value);
end

%% transmission at zero temperature
function [Transmission] = zeroTempTransmission(totalSystem, totalSysDeriv, gamma)
    % calculate the Greens Function
    GreensFuncInv = eye(length(totalSystem)) - totalSystem;
    GreensFunc = inv(GreensFuncInv);
    
    % calculate the matrix product
    Transmission = totalSysDeriv * GreensFunc * gamma * GreensFunc';
end