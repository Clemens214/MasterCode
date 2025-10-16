function [Result] = torque(totalSystem, totalSysDeriv, gammaL, gammaR, Eigenvals, leftEVs, rightEVs, chemPot, options)
% calculate the torque through a molecule for zero temperature
arguments
    totalSystem
    totalSysDeriv
    gammaL
    gammaR
    Eigenvals
    leftEVs
    rightEVs
    chemPot
    options.linearResponse = false
    options.conservative = false
    options.nonconservative = false
    options.left = false
    options.right = false
end
    choice = struct('conservative', options.conservative, ...
                    'nonconservative', options.nonconservative, ...
                    'left', options.left, 'right', options.right);
    %disp('Starting calculation of the torque.')
    if options.linearResponse == true
        Energy = chemPot;
        TotalResult = choiceLin(Energy, totalSystem, totalSysDeriv, gammaL, gammaR, choice);
    elseif options.linearResponse == false
        chemPotL = chemPot;
        chemPotR = -1*chemPot;
        TotalResult = choiceCalc(Eigenvals, leftEVs, rightEVs, totalSysDeriv, gammaL, gammaR, chemPotL, chemPotR, choice);
    end
    Result = real(trace(TotalResult));
    %disp('Finished calculation of the torque.')
end

%% calculate the Torque for zero temperature
function [TotalResult] = choiceLin(Energy, totalSystem, totalSysDeriv, gammaL, gammaR, choice)
    if choice.conservative == true
        midFactor = gammaL + gammaR;
        TotalResult = Transmission(Energy, totalSystem, totalSysDeriv, midFactor);
    elseif choice.nonconservative == true
        midFactor = gammaL - gammaR;
        TotalResult = Transmission(Energy, totalSystem, totalSysDeriv, midFactor);
    elseif choice.left == true
        midFactor = gammaL;
        TotalResult = Transmission(Energy, totalSystem, totalSysDeriv, midFactor);
    elseif choice.right == true
        midFactor = gammaR;
        TotalResult = Transmission(Energy, totalSystem, totalSysDeriv, midFactor);
    else
        ResultL = Transmission(Energy, totalSystem, totalSysDeriv, gammaL);
        ResultR = Transmission(Energy, totalSystem, totalSysDeriv, gammaR);
        TotalResult = ResultL + ResultR;
    end
end

function [Result] = Transmission(Energy, totalSystem, totalSysDeriv, midFactor)
    % calculate the Transport through the molecule

    % calculate the Greens Function
    GreensFuncInv = Energy*eye(length(totalSystem)) - totalSystem;
    GreensFunc = inv(GreensFuncInv);
    
    % calculate the matrix product
    Result = totalSysDeriv * GreensFunc * midFactor * GreensFunc';
end

%% 
function [TotalResult] = choiceCalc(Eigenvals, leftEVs, rightEVs, totalSysDeriv, gammaL, gammaR, chemPotL, chemPotR, choice)
    if choice.conservative == true
        midFactor = gammaL + gammaR;
        TotalResult = currentElement(Eigenvals, leftEVs, rightEVs, totalSysDeriv, midFactor, chemPotL, chemPotR, choice);
    elseif choice.nonconservative == true
        midFactor = gammaL + gammaR;
        TotalResult = currentElement(Eigenvals, leftEVs, rightEVs, totalSysDeriv, midFactor, chemPotL, chemPotR, choice);
    elseif choice.left == true
        midFactor = gammaL;
        TotalResult = currentElement(Eigenvals, leftEVs, rightEVs, totalSysDeriv, midFactor, chemPotL, chemPotR, choice);
    elseif choice.right == true
        midFactor = gammaR;
        TotalResult = currentElement(Eigenvals, leftEVs, rightEVs, totalSysDeriv, midFactor, chemPotL, chemPotR, choice);
    else
        choiceL = choice;
        choiceL.left = true;
        ResultL = currentElement(Eigenvals, leftEVs, rightEVs, totalSysDeriv, gammaL, chemPotL, chemPotR, choiceL);
        choiceR = choice;
        choiceR.left = true;
        ResultR = currentElement(Eigenvals, leftEVs, rightEVs, totalSysDeriv, gammaR, chemPotL, chemPotR, choiceR);
        TotalResult = ResultL + ResultR;
    end
end

function [Result] = currentElement(Eigenvals, leftEVs, rightEVs, totalSysDeriv, midFactor, chemPotL, chemPotR, choice)
    index = struct('i', [], 'j', [], ...
                    'Eigenval', [], 'EigenvalD', [], ...
                    'leftEV', [], 'leftEVD', [], ...
                    'rightEV', [], 'rightEVD', []);
    for i = 1:length(Eigenvals)
        for j = 1:length(Eigenvals)
            idx = (i-1)*length(Eigenvals) + j;
            index(idx).i = i;
            index(idx).j = j;
            index(idx).Eigenval = Eigenvals(i,i);
            index(idx).leftEV = leftEVs(:,i)';
            index(idx).rightEV = rightEVs(:,i);
            index(idx).EigenvalD = Eigenvals(j,j)';
            index(idx).leftEVD = leftEVs(:,j);
            index(idx).rightEVD = rightEVs(:,j)';
        end
    end

    %disp('Starting calculation of the torque element.')
    Result = 0;
    parfor idx = 1:length(index)
        % get the normal left and right Eigenvectors
        EigVal = index(idx).Eigenval;
        leftEV = index(idx).leftEV;
        rightEV = index(idx).rightEV;
        
        % get the daggered Eigenvectors
        EigValDagger = index(idx).EigenvalD;
        leftEVdagger = index(idx).leftEVD;
        rightEVdagger = index(idx).rightEVD;
            
        % compute the matrix element for chosen i and j
        ProductLeft = totalSysDeriv * rightEV;
        ProductMid = leftEV * midFactor * leftEVdagger;
        ProductRight = rightEVdagger;
        
        Product = ProductLeft * ProductMid * ProductRight;
        
        % compute the additional matrix element
        factor = choiceFactor(EigVal, EigValDagger, chemPotL, chemPotR, choice);
        
        Result = Result + Product*factor;
    end
    %disp('Finished calculation of the torque element.')
end

%% calculate the factor
function [Factor] = choiceFactor(EigVal, EigValDagger, chemPotL, chemPotR, choice)
    if choice.conservative == true
        Factor = factorElement(EigVal, EigValDagger, chemPotL) + factorElement(EigVal, EigValDagger, chemPotR);
    elseif choice.nonconservative == true
        Factor = factorElement(EigVal, EigValDagger, chemPotL) - factorElement(EigVal, EigValDagger, chemPotR);
    elseif choice.left == true
        Factor = factorElement(EigVal, EigValDagger, chemPotL);
    elseif choice.right == true
        Factor = factorElement(EigVal, EigValDagger, chemPotR);
    end
end

function [result] = factorElement(eig1, eig2, chemPot)
    %pot = chemPot(1);
    factor = 1/(eig1 -eig2);
    element1 = log(chemPot - eig1);
    element2 = log(chemPot - eig2);
    result = factor*(element1 - element2);
end