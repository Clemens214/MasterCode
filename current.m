function [entropyResultImag, particleResultImag, energyResultImag] = current(totalSystem, gammaL, gammaR, Eigenvals, leftEVs, rightEVs, Temp, chemPot, value)
    % define the states to be used for calculating the current
    stateLeft = totalSystem(:,value)';
    stateRight = totalSystem(:,value+1);
    
    %compute the entropy current
    midFactor = gammaL - gammaR;
    %disp('Starting calculation of the current.')
    [entropyResult, particleResult, energyResult, entropyAll, particleAll, energyAll] = currentElement(stateLeft, stateRight, Eigenvals, leftEVs, rightEVs, midFactor, Temp, chemPot);
    
    entropyResultImag = imag(entropyResult);
    particleResultImag = imag(particleResult);
    energyResultImag = imag(energyResult);
    
    entropyImag = imag(entropyAll);
    particleImag = imag(particleAll);
    energyImag = imag(energyAll);
    %disp('Finished calculation of the current.')
end

%% 
function [entropyResult, particleResult, energyResult, entropyAll, particleAll, energyAll] = currentElement (stateLeft, stateRight, Eigenvals, leftEVs, rightEVs, midFactor, Temp, chemPot)
    %disp('Starting calculation of the Hurwitz Zeta values.')
    Hurwitz = zeros(length(Eigenvals));
    HurwitzDagger = zeros(length(Eigenvals));
    for k = 1:length(Eigenvals)
        EigVal = Eigenvals(k,k) - chemPot;
        EigValDagger = Eigenvals(k,k)' - chemPot;
        % compute the Hurwitz Zeta functions
        beta = 1/Temp;
        Hurwitz(k) = HurwitzZeta(EigVal, beta);
        HurwitzDagger(k) = HurwitzZeta(EigValDagger, beta);
    end
    %disp('Finished calculation of the Hurwitz Zeta values.')
    
    %disp('Starting calculation of the current element.')
    particleResult = 0;
    energyResult = 0;
    entropyResult = 0;
    
    particleAll = zeros(length(leftEVs), length(leftEVs));
    energyAll = zeros(length(leftEVs), length(leftEVs));
    entropyAll = zeros(length(leftEVs), length(leftEVs));
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
            entropyResult = entropyResult + Product*factor;

            particleAll(i,j) = Product*factorParticle;
            energyAll(i,j) = Product*factorEnergy;
            entropyAll(i,j) = Product*factor;
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

%% 
function [result] = particleElement(eig1, eig2, T)
    beta = 1/T;
    factor = 1/(eig1 -eig2);
    result = factor*HurwitzZeta(eig1, beta) - factor*HurwitzZeta(eig2, beta);
end

function [result] = energyElement (eig1, eig2, T)
    beta = 1/T;
    factor = 1/(eig1 -eig2);
    result = eig1*factor*HurwitzZeta(eig1, beta) - eig2*factor*HurwitzZeta(eig2, beta);
end

%% 
function [result] = HurwitzZeta(x, beta)
    factor = sign(imag(x))/(2*pi*1j*beta);
    value = 0.5+factor*x;
    result = -1*factor * hurwitzZeta(2, value);
end