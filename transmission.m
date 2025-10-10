function [Result] = transmission(totalSystem, gammaL, gammaR, Eigenvals, leftEVs, rightEVs, Temp, chemPot, options)
% calculate the transmission through a molecule
arguments
    totalSystem
    gammaL
    gammaR
    Eigenvals
    leftEVs
    rightEVs
    Temp
    chemPot
    options.value = 0
end
    %disp('Starting calculation of the current.')
    if Temp == 0
        TotalResult = Transmission(chemPot, totalSystem, gammaL, gammaR);
    else
        TotalResult = currentElement(Eigenvals, leftEVs, rightEVs, gammaL, gammaR, Temp, chemPot);
    end
    Result = trace(TotalResult);
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
function [Result] = currentElement(Eigenvals, leftEVs, rightEVs, gammaL, gammaR, Temp, chemPot)
    %disp('Starting calculation of the Hurwitz Zeta values.')
    Hurwitz = zeros(length(Eigenvals));
    HurwitzDagger = zeros(length(Eigenvals));
    for k = 1:length(Eigenvals)
        EigVal = Eigenvals(k,k) - chemPot;
        EigValDagger = Eigenvals(k,k)' - chemPot;
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
        ProductMidLeft = leftEV * gammaL;
        for j = 1:length(leftEVs)
            % get the daggered Eigenvectors
            EigValDagger = Eigenvals(j,j)';
            leftEVdagger = leftEVs(:,j);
            rightEVdagger = rightEVs(:,j)';
            
            % compute the matrix element for chosen i and j
            ProductMid = ProductMidLeft * leftEVdagger;
            ProductRight = rightEVdagger * gammaR;
            
            Product = ProductLeft * ProductMid * ProductRight;
            
            % compute the additional matrix element
            factor = factorElement(EigVal, EigValDagger, Hurwitz(i), HurwitzDagger(j));
            
            Result = Result + Product*factor;
        end
    end
    %disp('Finished calculation of the current element.')
end

%% calculate the factor for temperatures other than zero
function [result] = factorElement(eig1, eig2, HurwitzEig1, HurwitzEig2)
    factor = 1/(eig1 -eig2);
    result = factor*HurwitzEig1 - factor*HurwitzEig2;
end

function [result] = HurwitzZeta(x, beta)
    factor = sign(imag(x))/(2*pi*1j*beta);
    value = 0.5+factor*x;
    result = -1*factor * hurwitzZeta(2, value);
end

