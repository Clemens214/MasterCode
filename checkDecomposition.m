function [] = checkDecomposition(totalSystem, omega)
%UNTITLED Summary of this function goes here
arguments
    totalSystem
    omega
end
    GreensInv = omega*eye(length(totalSystem)) - totalSystem;
    GreensFunc = inv(GreensInv);
    TestFunc = GreensFunc * GreensFunc;
    gammaL = eye(size(totalSystem));
    gammaR = eye(size(totalSystem));
    
    chemPotL = 1;
    chemPotR = -1;
    
    % compute the Eigenvectors and the Eigenvalues of the system
    [Eigenvals, leftEVs, rightEVs] = getEigenvectors(totalSystem);
    EigFunc = EigCalc(omega, Eigenvals, leftEVs, rightEVs, gammaL, gammaR, chemPotL, chemPotR);
    [TestEig, ~, maxDiffEig] = checkResult(TestFunc, EigFunc);
    if all(TestEig)
        disp(['The Eigenvalue decomposition DOES match the Greens function! Maximum Difference: ', num2str(maxDiffEig)])
    else
        disp(['The Eigenvalue decomposition does NOT match the Greens function! Maximum Difference: ', num2str(maxDiffEig)])
    end
    
    % compute the Schur decomposition of the system
    [Diag, upperTriag, SchurVec] = getSchur(totalSystem);
    SchurFunc = SchurCalc(omega, Diag, upperTriag, SchurVec, gammaL, gammaR, chemPotL, chemPotR);
    [TestSchur, ~, maxDiffSchur] = checkResult(TestFunc, SchurFunc);
    if all(TestSchur)
        disp(['The Schur decomposition DOES match the Greens function! Maximum Difference: ', num2str(maxDiffSchur)])
    else
        disp(['The Schur decomposition does NOT match the Greens function! Maximum Difference: ', num2str(maxDiffSchur)])
    end
    
    % check the results
    [Test, ~, maxDiff] = checkResult(EigFunc, SchurFunc);
    if all(Test)
        disp(['The multiplied Greens functions DO match! Maximum Difference: ', num2str(maxDiff)])
    else
        disp(['The multiplied Greens functions do NOT match! Maximum Difference: ', num2str(maxDiff)])
    end
end

%% calculate the Eigenvalue decomposition
function [Result] = EigCalc(omega, Eigenvals, leftEVs, rightEVs, gammaL, gammaR, chemPotL, chemPotR, options)
arguments
    omega
    Eigenvals
    leftEVs
    rightEVs
    gammaL
    gammaR
    chemPotL
    chemPotR
    options.check = true
end
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
        Factor = (omega - EigVal) * (omega - EigValDagger);
        %Factor = FactorElement(EigVal, EigValDagger, chemPotL) - FactorElement(EigVal, EigValDagger, chemPotR);
        
        Result = Result + Product*Factor;
    end
    %disp('Finished calculation of the transmission element.')
end

%% calculate the Schur decomposition
function [Result] = SchurCalc(omega, Diag, upperTriag, SchurVec, gammaL, gammaR, chemPotL, chemPotR, options)
arguments
    omega
    Diag
    upperTriag
    SchurVec
    gammaL
    gammaR
    chemPotL
    chemPotR
    options.check = true
end
    indices = getIndices(length(Diag));
    combis = combinations(indices, indices);
    combis1 = combis{:, 1};
    combis2 = combis{:, 2};
    Result = zeros(size(Diag));
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
        disp(['Matrix: idx=', num2str(idx)])
        Factor = FactorCalc(omega, Diag, upperTriag, index1, index2, chemPotL, chemPotR);
        
        Result = Result + Product*Factor;
    end
    Result = SchurVec * Result * SchurVec';
end

function [Result] = FactorCalc(omega, Diag, upperTriag, index1, index2, chemPotL, chemPotR, options)
arguments
    omega
    Diag
    upperTriag
    index1
    index2
    chemPotL
    chemPotR
    options.check = true
end
    % get the combinations of the paths
    combis = combinations(index1.factor, index2.factor);
    combis1 = combis{:, 1};
    combis2 = combis{:, 2};
    values = zeros(1, height(combis));
    for idx = 1:height(combis)
        [factors1, vals1] = getVals(combis1(idx).power, combis1(idx).paths, Diag, upperTriag);
        [factors2, vals2] = getVals(combis2(idx).power, combis2(idx).paths, Diag, upperTriag);
        vals = [vals1, conj(vals2)];
        factors = [factors1, conj(factors2)];
        % perform the partial fraction decomposition
        disp(['Factor: idx=', num2str(idx)])
        [constants, values, orders] = partialFraction(factors, vals);
        % calculate the result
        denominators = (omega - values).^orders;
        elements = constants ./ denominators;
        element = sum(elements);
        % save the result
        values(idx) = element;
    end
    Result = sum(values);
end

function [Constants, Denominators, orders] = partialFraction(factors, vals)
arguments
    factors
    vals
end
    % identify the unique values
    [values, ~, idx] = unique(vals);
    % get the orders of the values
    orders = groupcounts(idx);
    if numel(vals) ~= sum(orders)
        error('The orders and the number of elements in the partial fraction do NOT match!')
    end
    % calculate the coefficients for the different values
    PolyVals = zeros(numel(vals), numel(vals));
    Denominators = zeros(1, numel(vals));
    Orders = zeros(1, numel(vals));
    column = 0;
    for i = 1:numel(values)
        value = values(i);
        order = orders(i);
        % calculate the coefficents for the value
        for j = 1:order
            % get the values in the polynome
            removed = 0;
            mask = true(1, numel(vals));
            for k = 1:numel(vals)
                if value == vals(k) && removed < j
                    mask(k) = false;
                    removed = removed + 1;
                end
            end
            % calculate the coefficients
            Constants = poly(vals(mask));
            disp(Constants)
            % return the coefficients
            column = column + 1;
            PolyVals(:, column) = [zeros(1, numel(values)-j), Constants].';
            test = [zeros(1, numel(values)-j), Constants].';
            disp(['Size column=',num2str(numel(vals)),', Size Result=',num2str(numel(test))])
            Denominators(column) = value;
            Orders(column) = order;
        end
    end
    % calculate the numerators, includung the factor
    Factor = prod(factors);
    Numerators = [zeros(1, numel(vals)-1), Factor];
    % calculate the coefficients
    Constants = PolyVals \ Numerators.';
    Constants = Constants.';
end

%% helping functions
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

function [factors, vals] = getVals(power, path, Diag, upperTriag)
    if power == 0
        vals = Diag(path(1), path(1));
        factors = 1;
    else
        % define the values from the diagonal matrix
        vals = zeros(1, length(path)-1);
        vals(1) = Diag(path(1), path(1));
        % define the factors from the upper triangular matrix
        factors = zeros(1, length(path)-1);
        factors(1) = 1;
        % populate the lists
        for i = 2:length(path)
            vals(i) = Diag(path(i-1), path(i-1));
            factors(i) = upperTriag(path(i-1), path(i));
        end
    end
end

function [Result] = FactorElement(eig1, eig2, chemPot)
    if eig1 ~= eig2
        Factor = 1/(eig1 - eig2);
        Element1 = log(chemPot - eig1);
        Element2 = log(chemPot - eig2);
        Result = Factor*(Element1 - Element2);
    else
        Result = -1/(chemPot - eig1);
    end
end

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