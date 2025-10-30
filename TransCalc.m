function [Results] = TransCalc(totalSystem, gammaL, gammaR, Eigenvals, leftEVs, rightEVs, chemPots, options)
% calculate the transmission through a molecule for zero temperature
arguments
    totalSystem
    gammaL
    gammaR
    Eigenvals
    leftEVs
    rightEVs
    chemPots
    options.linearResponse = false
end
    %disp('Starting calculation of the current.')
    if options.linearResponse == true
        Energies = getEnergies(chemPots);
        Results = Transmission(Energies, totalSystem, gammaL, gammaR);
    elseif options.linearResponse == false && options.integrate == true
        Results = integrate(chemPots, totalSystem, gammaL, gammaR);
    end
    %disp('Finished calculation of the current.')
end

%% integrate the transmission
function [Results] = integrate(chemPots, totalSystem, gammaL, gammaR, options)
arguments
    chemPots
    totalSystem
    gammaL
    gammaR
    options.stepMult = 10
end
    Energies = getEnergies(chemPots);
    Diffs = zeros(1, length(Energies)-1);
    for i = 2:length(Energies)
        Diffs(i) = Energies(i) - Energies(i-1); 
    end
    lcd = LowestCommonDenominator(Diffs);
    stepSize = (1/lcd) / options.stepMult;

    % calculate the transmissions
    evalPoints = makeList(max(Energies), min(Energies), stepSize);
    values = Transmission(evalPoints, totalSystem, gammaL, gammaR);

    % calculate the integrals
    Results = zeros(1, length(chemPots));
    for i = 1:length(chemPots)
        fermiFunc = getFermiFunc(evalPoints, chemPots(i).left) - getFermiFunc(evalPoints, chemPots(i).right);
        % calculate the Result
        yData = fermiFunc .* values;
        if length(evalPoints) > 1
            Results(i) = trapz(evalPoints, yData);
        elseif isscalar(evalPoints)
            Results(i) = 0;
        end
        disp(['Voltage: ', num2str(chemPots(i).left-chemPots(i).right), ', j=', num2str(i)])
    end
end

function [fermiFunc] = getFermiFunc(evalPoints, chemPot, Temp)
arguments
    evalPoints
    chemPot
    Temp = 0
end
    fermiFunc = zeros(size(evalPoints));
    for i = 1:length(evalPoints)
        E = evalPoints(i);
        if Temp ~= 0
            fermiFunc(i) = 1/(exp((E-chemPot)/Temp)+1);
        elseif Temp == 0
            if E <= chemPot
                fermiFunc(i) = 1;
            else
                fermiFunc(i) = 0;
            end
        end
    end
end

%% total transmission in the linear transport approximation
function [Results] = Transmission(Energies, totalSystem, gammaL, gammaR)
    %calculates the transport through a molecule in the linear transport approximation
    arguments
        Energies
        totalSystem
        gammaL
        gammaR
    end
    % calculate the transport matrix and the trace
    Traces = zeros(1, length(Energies));
    for i = 1:length(Energies)
        %Matrix = TransmissionZeroTemp(Energies(i), totalSystem, gammaL, gammaR);
        Matrix = TransmissionAlt(Energies(i), totalSystem, gammaL, gammaR);
        Traces(i) = trace(real(Matrix));
    end
    % return the results
    Results = Traces;
end

function [Result] = TransmissionAlt(Energy, totalSystem, gammaL, gammaR)
    % G*B*Gt*C
    % F = decomposition(GreensInv,'lu');
    % Y = F \ B;
    % W = C * Y;
    % Z = F' \ W;
    % t = trace(Z);
    
    % GreensFunc * gammaL * GreensFunc' * gammaR
    GreensInv = Energy*eye(length(totalSystem)) - totalSystem;
    F = decomposition(GreensInv,'lu');    % create reusable LU object (works for sparse/dense)
    
    Y = F \ gammaL;
    W = gammaR * Y;
    Z = F' \ W;                   % uses transpose of factorization
    % t = trace(Z);
    Result = Z;
end

%% helping functions
function [Filtered] = getEnergies(chemPots)
    Energies = zeros(1, length(chemPots)*2);
    for i = 1:length(chemPots)
        Energies(2*i-1) = chemPots(i).left;
        Energies(2*i) = chemPots(i).right;
    end
    Sorted = sort(Energies);
    Filtered = unique(Sorted);
end

function [values] = makeList(maxVal, minVal, stepVal)
    arguments
        maxVal
        minVal
        stepVal
    end
    numVal = (maxVal-minVal)/stepVal+1;
    values = linspace(minVal, maxVal, numVal);
end

%% ------------------------------ deprecated functions ------------------------------
%% total transmission for finite voltages
function [Result] = TransmissionMatrix(Eigenvals, leftEVs, rightEVs, gammaL, gammaR, chemPotL, chemPotR)
    index = struct('i', [], 'j', [], ...
                    'Eigenval', [], 'EigenvalD', [], ...
                    'leftEV', [], 'leftEVD', [], ...
                    'rightEV', [], 'rightEVD', []);
    for i = 1:length(Eigenvals)
        for j = 1:length(Eigenvals)
            idx = (i-1)*length(Eigenvals) + j;
            index(idx).i = i;
            index(idx).j = j;
            % set the normal variables
            index(idx).Eigenval = Eigenvals(i,i);
            index(idx).leftEV = leftEVs(:,i)';
            index(idx).rightEV = rightEVs(:,i);
            % set the daggered variables
            index(idx).EigenvalD = Eigenvals(j,j)';
            index(idx).leftEVD = leftEVs(:,j);
            index(idx).rightEVD = rightEVs(:,j)';
        end
    end

    %disp('Starting calculation of the transmission element.')
    Result = 0;
    parfor idx = 1:length(index)
        % get the normal left and right Eigenvectors
        leftEV = index(idx).leftEV;
        rightEV = index(idx).rightEV;
        
        % get the daggered Eigenvectors
        leftEVdagger = index(idx).leftEVD;
        rightEVdagger = index(idx).rightEVD;
            
        % compute the matrix element for chosen i and j
        ProductLeft = rightEV;
        ProductMid = leftEV * gammaR * leftEVdagger;
        ProductRight = rightEVdagger * gammaL;
            
        Product = ProductLeft * ProductMid * ProductRight;
        
        % compute the additional matrix element
        EigVal = index(idx).Eigenval;
        EigValDagger = index(idx).EigenvalD;
        factor = factorElement(EigVal, EigValDagger, chemPotL) - factorElement(EigVal, EigValDagger, chemPotR);
        
        Result = Result + Product*factor;
    end
    %disp('Finished calculation of the transmission element.')
end

function [result] = factorElement(eig1, eig2, chemPot)
    if eig1 ~= eig2
        factor = 1/(eig1 - eig2);
        element1 = log(chemPot - eig1);
        element2 = log(chemPot - eig2);
        result = factor*(element1 - element2);
    else
        result = -1/(chemPot - eig1);
    end
end

%% total transmission in the linear transport approximation
function [Result] = TransmissionZeroTemp(Energy, totalSystem, gammaL, gammaR)
    % calculate the Greens Function
    GreensFuncInv = Energy*eye(length(totalSystem)) - totalSystem;
    GreensFunc = inv(GreensFuncInv);
    
    % calculate the matrix product
    Result = GreensFunc * gammaL * GreensFunc' * gammaR;
end