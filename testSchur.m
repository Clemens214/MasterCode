function [] = testSchur(totalSystem, omega)
%UNTITLED Summary of this function goes here
arguments
    totalSystem
    omega
end 
    [Diag, upperTriag, SchurVec] = getSchur(totalSystem);
    
    GreensInv = omega*eye(length(totalSystem)) - totalSystem;
    GreensFunc = inv(GreensInv);
    
    gammaL = eye(size(Diag));
    gammaR = eye(size(Diag));
    Result = MatrixCalc(omega, Diag, upperTriag, SchurVec, gammaL, gammaR);

    GreensCheck = GreensCalc(omega, Diag, upperTriag, SchurVec);
    
    [Test, ~, maxDiff] = checkResult(GreensFunc, GreensCheck);
    if all(Test)
        disp(['The Greens functions DO match! Maximum Difference: ', num2str(maxDiff)])
    else
        disp(['The Greens functions do NOT match! Maximum Difference: ', num2str(maxDiff)])
    end
end

function [Result] = GreensCalc(omega, Diag, upperTriag, SchurVec, options)
arguments
    omega
    Diag
    upperTriag
    SchurVec
    options.check = true
end
    indices = getIndices(length(Diag));
    %disp('Starting calculation of the Green's Function.')
    Matrix = zeros(size(Diag));
    Matrices = cell(1, length(Diag));
    for i = 1:numel(Matrices)
        Matrices{i} = zeros(length(Diag), length(Diag));
    end
    for idx = 1:numel(indices)
        factor = indices(idx).factor;
        values = zeros(1, numel(factor));
        for i = 1:height(factor)
            [factors, vals] = getVals(factor(i).power, factor(i).paths, Diag, upperTriag);
            denominators = omega - vals;
            elements = factors ./ denominators;
            values(i) = sum(elements);
        end
        row = indices(idx).row;
        column = indices(idx).column;
        Matrix(row, column) = sum(values);
        for i = 1:numel(values)
            Matrices{i}(row, column) = values(i);
        end
    end
    invDiag = (omega*eye(length(Diag)) - Diag);
    ResultSchur = invDiag \ Matrix;
    Result = SchurVec * ResultSchur * SchurVec';
    %disp('Finished calculation of the Green's Function.')
    if options.check == true
        for i = 1:numel(Matrices)
            Test = ((omega*eye(length(Diag)) - Diag) \ upperTriag)^(idx-1);
            if Matrices{i} ~= Test
                disp(['The matrices do NOT match!', ' power: ', num2str(power)])
            end
        end
    end
end

function [Result] = MatrixCalc(omega, Diag, upperTriag, SchurVec, gammaL, gammaR, options)
arguments
    omega
    Diag
    upperTriag
    SchurVec
    gammaL
    gammaR
    options.check = true
end
    indices = getIndices(length(Diag));
    combis = combinations(indices, indices);
    combis1 = combis{:, 1};
    combis2 = combis{:, 2};
    Matrix = zeros(size(Diag));
    for idx = 1:height(combis)
    %parfor idx = 1:height(combis)
        % get the normal matrix
        index1 = combis1(idx);
        matrix1 = zeros(size(Diag));
        matrix1(index1.row, index1.column) = 1;
        
        % get the daggered matrix
        index2 = combis2(idx);
        matrix2 = zeros(size(Diag));
        matrix2(index2.row, index2.column) = 1;
            
        % compute the matrix element for chosen matrices
        Product = gammaL * matrix1 * gammaR * matrix2';
        
        % compute the factor for the chosen matrices
        Factors = FactorCalc(omega, Diag, upperTriag, SchurVec, index1, index2);
        
        Result = Result + Product*factor;
    end
end

function [Result] = FactorCalc(omega, Diag, upperTriag, SchurVec, index1, index2, options)
arguments
    omega
    Diag
    upperTriag
    SchurVec
    index1
    index2
    options.check = true
end
    % get the combinations of the paths
    combis = combinations(index1.factor, index2.factor);
    combis1 = combis{:, 1};
    combis2 = combis{:, 2};
    parfor idx = 1:height(combis)
        [factors1, vals1] = getVals(combis1(idx).power, combis1(idx).paths, Diag, upperTriag);
        [factors2, vals2] = getVals(combis2(idx).power, combis2(idx).paths, Diag, upperTriag);
        factors = [factors1, factors2];
        vals = [vals1, vals2];
        for i = 1:length(vals)
            
        end
    end
end

function [factors, vals] = getVals(power, path, Diag, upperTriag)
    % define the values from the diagonal matrix
    vals = zeros(1, length(path));
    vals(1) = Diag(path(1), path(1));
    % define the factors from the upper triangular matrix
    factors = zeros(1, length(path));
    factors(1) = 1;
    % populate the lists
    for i = 2:length(path)
        vals(i) = Diag(path(i-1), path(i-1));
        if power == 0
            I = eye(path(k), path(k+1));
            factors(i) =  I(path(i-1), path(i));
        else
            factors(i) = upperTriag(path(i-1), path(i));
        end
    end
end

function [values] = getIndices(size)
    index = struct('index', [], 'row', [], 'column', [], 'power', [], 'paths', []);
    idx = 0;
    jdx = 0;
    for row = 1:size
        for column = row:size
            idx = idx+1;
            for power = 0:abs(row-column)
                jdx = jdx+1;
                index(jdx).index = idx;
                index(jdx).row = row;
                index(jdx).column = column;
                index(jdx).power = power;
            end
        end
    end
    parfor idx = 1:numel(index)
        row = index(idx).row;
        column = index(idx).column;
        power = index(idx).power;
        index(idx).paths = getPaths(row, column, power);
    end
    values = struct('row', [], 'column', [], 'factor', []);
    for i = 1:numel(index)
        idx = index(i).index;
        row = index(i).row;
        column = index(i).column;
        filt = index([index(:).row] == row);
        filt = filt([filt(:).column] == column);
        % set the variables
        values(idx).row = row;
        values(idx).column = column;
        factor = struct('power', [], 'paths', []);
        jdx = 0;
        for j = 1:length(filt)
            paths = filt(j).paths;
            for k = 1:height(filt(j).paths)
                jdx = jdx+1;
                factor(jdx).power = filt(j).power;
                factor(jdx).paths = paths(k, :);
            end
        end
        values(idx).factor = factor;
    end
end

function [paths] = getPaths(row, column, power)
    range = row+1 : column-1;
    if power > 0
        if isscalar(range)
            centers = range;
        else
            centers = nchoosek(range, power-1);
        end
        paths = zeros(size(centers, 1), size(centers, 2)+2);
        for idx = 1:size(centers, 1)
            paths(idx, :) = [row, centers(idx, :), column];
        end
    else
        paths = [row, column];
    end
    %disp(['row: ', num2str(row), ', column: ', num2str(column), ', power: ', num2str(power)])
end

%% helping functions

%% checking functions
function [Test, Diff, maxDiff] = checkResult(GreensFunc, GreensCheck, options)
arguments
    GreensFunc
    GreensCheck
    options.returnAbs = true
    options.Tolerance = 1e-10;
    %Tolerance = 1e-10;
end
    Test = isapprox(GreensFunc, GreensCheck, AbsoluteTolerance=options.Tolerance);
    Diff = GreensFunc - GreensCheck;
    if options.returnAbs == true
        Diff = abs(Diff);
    end
    maxDiff = max(max(abs(Diff)));
end