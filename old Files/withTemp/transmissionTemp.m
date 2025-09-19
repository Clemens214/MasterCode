function [Result, varargout] = transmission(totalSystem, gammaL, gammaR, Eigenvals, leftEVs, rightEVs, Temp)
    %calculates the transport through by a molecule as calculated in the
    %Phd Thesis
    
    TempVal = false;
    %disp('Starting calculation of the current.')
    if Temp == 0
        TransmissionMatrix = zeroTempTransmission(totalSystem, gammaL, gammaR);
    else
        TransmissionMatrix = Transmission(Eigenvals, leftEVs, rightEVs, gammaL, gammaR, Temp);
    end
    TransmissionElement = trace(TransmissionMatrix);

    Result = real(TransmissionElement);
    if Temp ~= 0 && TempVal == true
        TempResult = Result/Temp;
    end
    varargout{1} = real(TransmissionMatrix);
    %disp('Finished calculation of the current.')
end

%% transmission at finite temperature
function [Result] = Transmission(Eigenvals, leftEVs, rightEVs, gammaL, gammaR, Temp)
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
        %disp(['Hurwitz = ',num2str(Hurwitz(k)),' Hurwitz+ = ',num2str(HurwitzDagger(k))])
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
            factor = matrixElement(EigVal, EigValDagger, Hurwitz(i), HurwitzDagger(j));
            
            Result = Result + Product*factor;
            %disp(['i = ',num2str(i),', j = ',num2str(j),', Result = ',num2str(real(trace(Result)))])
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
    %result = -1*factor * hurwitzZeta(1, value);
    result = -1*factor * hurwitzZeta(2, value);
end

%% transmission at zero temperature
function [Transmission] = zeroTempTransmission(totalSystem, gammaL, gammaR)
    % calculate the Greens Function
    GreensFuncInv = eye(length(totalSystem)) - totalSystem;
    GreensFunc = inv(GreensFuncInv);
    
    % calculate the matrix product
    Transmission = GreensFunc * gammaR * GreensFunc' * gammaL;
end