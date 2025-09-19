function [entropyResult, particleResult, energyResult, ProductResult] = currentQuick(totalSystem, gammaL, gammaR, Eigenvals, leftEVs, rightEVs, Temp, chemPot, value)
    % define the states to be used for calculating the current
    stateLeft = zeros(1,length(Eigenvals));
    stateLeft(value+1) = 1;
    stateRight = zeros(1,length(Eigenvals))';
    stateRight(value+2) = 1;
    
    %compute the entropy current
    midFactor = gammaL - gammaR;
    %disp('Starting calculation of the current.')
    if Temp == 0
        ProductResult = imag(Transmission(chemPot, totalSystem, midFactor, value));
        particleResult = ProductResult;
        energyResult = chemPot*ProductResult;
        entropyResult = energyResult - chemPot*particleResult;
    else
        [entropyResult, particleResult, energyResult, ProductResult] = currentElement(stateLeft, stateRight, Eigenvals, leftEVs, rightEVs, midFactor, Temp, chemPot);
        entropyResult = imag(entropyResult);
        particleResult = imag(particleResult);
        energyResult = imag(energyResult);
        ProductResult = imag(ProductResult);
    end
    %disp('Finished calculation of the current.')
end

%% 
function [TransmissionEM] = Transmission (chemPot, totalSystemEM, midFactor, value)
    % calculate the Transmissions
    GreensFuncEM = GreensFunc(chemPot, totalSystemEM);
    TransmissionEM = Transport(GreensFuncEM, midFactor, value);
end

function [arrayInv] = GreensFunc (chemPot, totalSystem)
    array = eye(length(totalSystem))*chemPot - totalSystem;
    arrayInv = inv(array);
end

function [T] = Transport (GreensFunc, midFactor, value)
    % calculate the matrix product
    transport = GreensFunc * midFactor * GreensFunc';
    % calculate the current
    T = transport(value+1, value+2);
end

%% 
function [entropyResult, particleResult, energyResult, ProductResult] = currentElement (stateLeft, stateRight, Eigenvals, leftEVs, rightEVs, midFactor, Temp, chemPot)
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
    particleResult = 0;
    energyResult = 0;
    entropyResult = 0;
    ProductResult = 0;
    for i = 1:length(leftEVs)
        % get the normal left and right Eigenvectors
        EigVal = Eigenvals(i,i);
        leftEV = leftEVs(:,i)';
        rightEV = rightEVs(:,i);
        
        ProductLeft = stateLeft * rightEV;
        ProductMidLeft = leftEV * midFactor;
        for j = 1:length(leftEVs)
            % get the daggered Eigenvectors
            EigValDagger = Eigenvals(j,j)';
            leftEVdagger = leftEVs(:,j);
            rightEVdagger = rightEVs(:,j)';
            
            % compute the matrix element for chosen i and j
            ProductMid = ProductMidLeft * leftEVdagger;
            ProductRight = rightEVdagger * stateRight;
            
            Product = ProductLeft * ProductMid * ProductRight;
            
            % compute the additional matrix element
            factorEnergy = energyElementQuick(EigVal, EigValDagger, Hurwitz(i), HurwitzDagger(j));
            factorParticle = particleElementQuick(EigVal, EigValDagger, Hurwitz(i), HurwitzDagger(j));
            factor = factorEnergy - chemPot*factorParticle;
            
            particleResult = particleResult + Product*factorParticle;
            energyResult = energyResult + Product*factorEnergy;
            entropyResult = entropyResult + Product*factor*1/Temp;

            ProductResult = ProductResult + Product;
        end
    end
    %disp('Finished calculation of the current element.')
end

%%
function [result] = particleElementQuick (eig1, eig2, HurwitzEig1, HurwitzEig2)
    factor = 1/(eig1 -eig2);
    result = factor*HurwitzEig1 - factor*HurwitzEig2;
end

function [result] = energyElementQuick (eig1, eig2, HurwitzEig1, HurwitzEig2)
    factor = 1/(eig1 -eig2);
    result = eig1*factor*HurwitzEig1 - eig2*factor*HurwitzEig2;
end

function [result] = HurwitzZeta (x, beta)
    factor = sign(imag(x))/(2*pi*1j*beta);
    value = 0.5+factor*x;
    result = -1*factor * hurwitzZeta(2, value);
end