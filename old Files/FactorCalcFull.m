function [Result, varargout] = FactorCalcFull(Eigenvals, chemPots, options)
    %calculates the transport through by a molecule as calculated in the PhD-Thesis
    arguments
        Eigenvals
        chemPots
        options.Transmission = false
        options.Torque = false
        options.conservative = false
        options.nonconservative = false
    end
    % generate and populate the struct used in the parfor loop
    index = struct('i', [], 'j', [], ...
                    'EigenVal', [], 'EigenValD', []);
    for i = 1:length(Eigenvals)
        for j = 1:length(Eigenvals)
            idx = (i-1)*length(Eigenvals) + j;
            % set the normal variables
            EigenVal = Eigenvals(i,i);
            % set the daggered variables
            EigenValD = Eigenvals(j,j)';
            % populate the struct
            index(idx) = struct('i', i, 'j', j, ...
                                'EigenVal', EigenVal, 'EigenValD', EigenValD);
        end
    end
    [chemPotL, chemPotR] = chemPots{:};
    
    if options.Transmission == true && options.Torque == false
        MatricesTrans = matrixTrans(index, gammaL, gammaR);
        Result = MatricesTrans;
    elseif options.Transmission == false && options.Torque == true
        MatricesTorque = matrixTorque(index, gammaL, gammaR);
        Result = MatricesTorque;
    end
end

function [Factors, varargout] = FactorCalc(Eigenvals, chemPots, options)


    % calculate the factors
    Factors = zeros(1, numel(index));
    FactorsAdd = zeros(1, numel(index));
    Factors = struct('i', [], 'j', [], 'factor', []);
    FactorsAdd = struct('i', [], 'j', [], 'factor', []);
    parfor idx = 1:numel(index)
        % get the Eigenvalues
        EigVal = index(idx).EigenVal;
        EigValDagger = index(idx).EigenValD;

        % compute the factor
        [Factors(idx).factor, FactorsAdd(idx).factor] = FactorChoice(EigVal, EigValDagger, chemPotL, chemPotR, options)

        % return the indices
        Factors(idx).i = index(idx).i;
        Factors(idx).j = index(idx).j;
        FactorsAdd(idx).i = index(idx).i;
        FactorsAdd(idx).j = index(idx).j;
    end
    varargout{1} = FactorsAdd;
end

function [Matrices] = matrixTrans(index, gammaL, gammaR)
    Matrices = matrixTransFull(index, gammaL, gammaR);
end

function [Matrices] = matrixTorque(index, gammaL, gammaR)
    Matrices = matrixTorqueFull(index, totalSysDeriv, gammaL, gammaR);
end

%% calculate the matrices and return the full matrices
function [Matrices] = matrixTransFull(index, gammaL, gammaR)
    %disp('Starting calculation of the current matrices.')
    Matrices = cell(1, numel(index));
    parfor idx = 1:numel(index)
        % get the normal Eigenvectors
        leftEV = index(idx).leftEV;
        rightEV = index(idx).rightEV;
        
        % get the daggered Eigenvectors
        leftEVdagger = index(idx).leftEVD;
        rightEVdagger = index(idx).rightEVD;
        
        % compute the matrix element for chosen i and j
        ProductLeft = rightEV;
        ProductMid = leftEV * gammaR * leftEVdagger;
        ProductRight = rightEVdagger * gammaL;
        
        Matrices{idx} = ProductLeft * ProductMid * ProductRight;
    end
    %disp('Finished calculation of the current matrices.')
end

function [Matrices] = matrixTorqueFull(index, totalSysDeriv, gamma)
    %disp('Starting calculation of the torque matrices.')
    Matrices = cell(1, numel(index));
    parfor idx = 1:numel(index)
        % get the normal Eigenvectors
        leftEV = index(idx).leftEV;
        rightEV = index(idx).rightEV;
        
        % get the daggered Eigenvectors
        leftEVdagger = index(idx).leftEVD;
        rightEVdagger = index(idx).rightEVD;
        
        % compute the matrix element for chosen i and j
        ProductLeft = totalSysDeriv * rightEV;
        ProductMid = leftEV * gamma * leftEVdagger;
        ProductRight = rightEVdagger;
        
        Matrices{idx} = ProductLeft * ProductMid * ProductRight;
    end
    %disp('Finished calculation of the torque matrices.')
end