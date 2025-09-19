function [Factors] = FactorCalc(Eigenvals, chemPots, options)
    %calculates the transport through by a molecule as calculated in the PhD-Thesis
    arguments
        Eigenvals
        chemPots
        options.linearApprox = false
    end
    
    index = struct('EigenVal', [], 'EigenValD', []);
    for i = 1:length(Eigenvals)
        for j = 1:length(Eigenvals)
            idx = (i-1)*length(Eigenvals) + j;
            % set all the different variables
            EigenVal = Eigenvals(i,i);
            EigenValD = Eigenvals(j,j)';
            % populate the struct
            index(idx) = struct('EigenVal', EigenVal, 'EigenValD', EigenValD);
        end
    end
    
    Factors = Factor(index, chemPots);
end

%% calculate the factors
function [Factors] = Factor(index, chemPots)
    [chemPotL, chemPotR] = chemPots{:};
    %disp('Starting calculation of the factors.')
    Factors = zeros(1, numel(index));
    parfor idx = 1:numel(index)
        % get the Eigenvalues
        EigVal = index(idx).EigenVal;
        EigValDagger = index(idx).EigenValD;
        
        % compute the additional matrix element
        Factors(idx) = factorElement(EigVal, EigValDagger, chemPotL) - factorElement(EigVal, EigValDagger, chemPotR);
    end
    %disp('Finished calculation of the factors.')
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