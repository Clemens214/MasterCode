function [Result, varargout] = transport(totalSystem, gammaL, gammaR, Eigenvals, leftEVs, rightEVs, Temp, value, step)
    %calculates the transport through a molecule as calculated in the
    %Project-Praktikum

    % define the states to be used for calculating the current
    stateLeft = zeros(1,length(Eigenvals));
    stateLeft(value) = 1;
    stateRight = zeros(1,length(Eigenvals))';
    stateRight(value+step) = 1;
    
    %compute the current
    TempVal = false;
    midFactor = -1*(gammaL - gammaR);
    %disp('Starting calculation of the current.')
    if Temp == 0
        TransmissionMatrix = zeroTempTransmission(totalSystem, midFactor);
    else
        TransmissionMatrix = Transmission(Eigenvals, leftEVs, rightEVs, midFactor, Temp);
    end
    TransmissionElement = stateLeft * TransmissionMatrix * stateRight;

    Result = imag(TransmissionElement);
    if Temp ~= 0 && TempVal == true
        TempResult = Result/Temp;
    end
    varargout{1} = imag(TransmissionMatrix);
    %disp('Finished calculation of the current.')
end

%% transmission at finite temperature

%function [Result] = Transmission(stateLeft, stateRight, Eigenvals, leftEVs, rightEVs, midFactor, Temp)
        %ProductLeft = stateLeft * rightEV;
        %ProductRight = rightEVdagger * stateRight;
function [Result] = Transmission(Eigenvals, leftEVs, rightEVs, midFactor, Temp)
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
        
        ProductLeft = rightEV;
        ProductMidLeft = leftEV * midFactor;
        for j = 1:length(leftEVs)
            % get the daggered Eigenvectors
            EigValDagger = Eigenvals(j,j)';
            leftEVdagger = leftEVs(:,j);
            rightEVdagger = rightEVs(:,j)';
            
            % compute the matrix element for chosen i and j
            ProductMid = ProductMidLeft * leftEVdagger;
            ProductRight = rightEVdagger;
            
            %compute the total product of the matrices
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
    result = -1*factor * hurwitzZeta(2, value);
end

%% transmission at zero temperature
function [Transmission] = zeroTempTransmission(totalSystem, midFactor)
    % calculate the Greens Function
    GreensFuncInv = eye(length(totalSystem)) - totalSystem;
    GreensFunc = inv(GreensFuncInv);
    
    % calculate the matrix product
    Transmission = GreensFunc * midFactor * GreensFunc';
end