function [Result] = transport(totalSystem, gammaL, gammaR, Eigenvals, leftEVs, rightEVs, Temp, chemPot, options)
% calculates the transport through a molecule
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
    % define the states to be used for calculating the current
    stateLeft = zeros(1,length(Eigenvals));
    stateLeft(options.value+1) = 1;
    stateRight = zeros(1,length(Eigenvals))';
    stateRight(options.value+2) = 1;
    
    %compute the current
    midFactor = gammaL - gammaR;
    %disp('Starting calculation of the current.')
    if Temp == 0
        TotalResult = Transport(chemPot, totalSystem, midFactor, value);
    else
        TotalResult = currentElement(stateLeft, stateRight, Eigenvals, leftEVs, rightEVs, midFactor, Temp, chemPot);
    end
    Result = imag(TotalResult);
    %disp('Finished calculation of the current.')
end

%% 
function [T] = Transport (Energy, totalSystem, midFactor, value)
    % calculate the Transport through the molecule

    % calculate the Greens Function
    GreensFuncInv = Energy*eye(length(totalSystem)) - totalSystem;
    GreensFunc = inv(GreensFuncInv);
    
    % calculate the matrix product
    Transport = GreensFunc * midFactor * GreensFunc';
    % calculate the current
    T = Transport(value+1, value+2);
end

%% 
function [Result] = currentElement (stateLeft, stateRight, Eigenvals, leftEVs, rightEVs, midFactor, Temp, chemPot)
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
            factorParticle = particleElementQuick(EigVal, EigValDagger, Hurwitz(i), HurwitzDagger(j));
            
            Result = Result + Product*factorParticle;
        end
    end
    %disp('Finished calculation of the current element.')
end

%%
function [result] = particleElementQuick (eig1, eig2, HurwitzEig1, HurwitzEig2)
    factor = 1/(eig1 -eig2);
    result = factor*HurwitzEig1 - factor*HurwitzEig2;
end

function [result] = HurwitzZeta (x, beta)
    factor = sign(imag(x))/(2*pi*1j*beta);
    value = 0.5+factor*x;
    result = -1*factor * hurwitzZeta(2, value);
end