function [Transmissions, varargout] = TransCalc(SystemVals, chemPots, angles, voltages)
    %calculates the transport through by a molecule as in the PhD thesis
    arguments
        SystemVals
        chemPots
        angles
        voltages
    end
    % calculate the matrices in the presence of an applied voltage
    Transmissions = cell(1, length(angles));
    Traces = cell(1, length(angles));
    for i = 1:length(angles)
        gammaL = SystemVals(i).gammaL;
        gammaR = SystemVals(i).gammaR;
        leftEVs = SystemVals(i).leftEVs;
        rightEVs = SystemVals(i).rightEVs;
        
        % compute the matrices
        disp(['Start calculation of the transmission matrices. Angle = ',num2str(angles(i)),', i=',num2str(i)])
        TransMatrices = MatrixCalc(leftEVs, rightEVs, gammaL, gammaR);
        disp(['Finished calculation of the transmission matrices. Angle = ',num2str(angles(i)),', i=',num2str(i)])
        
        Transmissions{i} = zeros(1, length(voltages));
        Traces{i} = cell(1, length(voltages));
        for j = 1:length(voltages)
            % compute the factors
            EigenVals = SystemVals(i).EigenVals;
            chemPot = {chemPots(j).left, chemPots(j).right};
            Factors = FactorCalc(EigenVals, chemPot);
            
            % compute the transmission
            disp('Start calculation of the transmission.')
            [Transmission, TransTraces] = ResultCalc(TransMatrices, Factors);
            Transmissions{i}(j) = Transmission;
            Traces{i}{j} = TransTraces;
            disp(['Finished calculation of the transmission. ',num2str(Transmission), ...
                ' Angle = ',num2str(angles(i)),', i=',num2str(i), ...
                ', Voltage = ',num2str(voltages(j)),', j=',num2str(j)])
        end
    end
    varargout{1} = Traces;
end

function [Result, varargout] = ResultCalc(Matrices, Factors, options)
    %calculates the transport through a molecule
    arguments
        Matrices
        Factors
        options.evalIndex = 1
    end
    % calculate the transport matrix
    Size = size(Matrices);
    Matrix = zeros(Size(1), 1);
    for idx = 1:length(Matrices)
        Matrix = Matrix + Matrices(:,idx)*Factors(idx);
    end

    % calculate the traces
    Traces = allTrace(real(Matrix));
    
    % calculate the final results
    Result = sum(real(Matrix));
    varargout{1} = Traces;
    varargout{2} = real(Matrix);
end

%% calculate the matrices
function [Matrices] = MatrixCalc(leftEVs, rightEVs, gammaL, gammaR)
    %calculates all the different matrix elements
    arguments
        leftEVs
        rightEVs
        gammaL
        gammaR
    end
    % generate and populate the struct used in the parfor loop
    index = struct('leftEV', [], 'leftEVD', [], ...
                    'rightEV', [], 'rightEVD', []);
    for i = 1:length(leftEVs)
        for j = 1:length(leftEVs)
            idx = (i-1)*length(leftEVs) + j;
            % set all the different variables
            leftEV = leftEVs(:,i)';
            leftEVD = leftEVs(:,j);
            rightEV = rightEVs(:,i);
            rightEVD = rightEVs(:,j)';
            % populate the struct
            index(idx) = struct('leftEV', leftEV, 'leftEVD', leftEVD, ...
                                'rightEV', rightEV, 'rightEVD', rightEVD);
        end
    end
    % calculate the matrices
    Matrices = MatrixTrans(index, gammaL, gammaR);
end

function [Matrices] = MatrixTrans(index, gammaL, gammaR)
    Matrices = zeros(length(gammaL), length(index));
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
        
        Product = ProductLeft * ProductMid * ProductRight;
        Matrices(:, idx) = diag(Product);
    end
end

%% calculate the facctors
function [Factors] = FactorCalc(Eigenvals, chemPots)
    %calculates the transport through by a molecule as calculated in the PhD-Thesis
    arguments
        Eigenvals
        chemPots
    end
    % generate and populate the struct used in the parfor loop
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
    [chemPotL, chemPotR] = chemPots{:};
    % calculate the factors
    Factors = zeros(1, numel(index));
    parfor idx = 1:numel(index)
        % get the Eigenvalues
        EigVal = index(idx).EigenVal;
        EigValDagger = index(idx).EigenValD;
        
        % compute the factor
        Factors(idx) = factorElement(EigVal, EigValDagger, chemPotL) - factorElement(EigVal, EigValDagger, chemPotR);
    end
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

%% helping functions
function [Transmissions] = allTrace(Matrix)
    Transmissions = zeros(1, floor(length(Matrix)/2));
    for i = 1:floor(length(Matrix)/2)
        startIndex = i;
        endIndex = length(Matrix) - (i-1);
        matrix = Matrix(startIndex:endIndex);
        Transmissions(i) = sum(matrix);
    end
end