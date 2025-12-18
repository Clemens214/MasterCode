function [] = testSchur(totalSystem, omega)
%UNTITLED Summary of this function goes here
arguments
    totalSystem
    omega
end 
    [Diag, upperTriag, SchurVec] = getSchur(totalSystem);
    
    GreensInv = omega*eye(length(totalSystem)) - totalSystem;
    GreensFunc = inv(GreensInv);
    
    GreensCheck = GreensCalc(omega, Diag, upperTriag, SchurVec);
    
    [Test, ~, maxDiff] = checkResult(GreensFunc, GreensCheck);
    if all(Test)
        disp("The Green's functions DO match! Maximum Difference: ", num2str(maxDiff))
    else
        disp("The Green's functions do NOT match! Maximum Difference: ", num2str(maxDiff))
    end
end

function [Result] = GreensCalc(omega, Diag, upperTriag, SchurVec)
    index = struct('row', [], 'column', [], 'power', [], 'paths', []);
    idx = 0;
    for row = 1:length(Diag)
        for column = row:length(Diag)
            for power = 0:abs(row-column)
                idx = idx + 1;
                % set the indices
                index(idx).row = row;
                index(idx).column = column;
                index(idx).power = power;
                % set the variables
                matrix = zeros(length(Diag), length(Diag));
                matrix(row, column) = 1;
                index(idx).matrix = matrix;
                index(idx).paths = getPaths(row, column, power);
            end
        end
    end
    Matrix = zeros(size(Diag));
    Matrices = cell(1, length(Diag));
    for i = 1:numel(Matrices)
        Matrices{i} = zeros(length(Diag), length(Diag));
    end
    for i = 1:numel(index)
        row = index(i).row;
        column = index(i).column;
        power = index(i).power;
        paths = index(i).paths;

        Factor = 0;
        for j = 1:height(paths)
            path = paths(j,:);
            Element = 1;
            for k = 1:length(path)-1
                if power == 0
                    I = eye(path(k), path(k+1));
                    val =  I(path(k), path(k+1));
                else
                    fac = omega - Diag(path(k), path(k));
                    val = 1/fac * upperTriag(path(k), path(k+1));
                end
                Element = Element * val;
            end
            Factor = Factor + Element;
        end
        Matrix(row, column) = Matrix(row, column) + Factor;
        MatricesMatrix = Matrices{power+1};
        MatricesMatrix(row, column) = Factor;
        Matrices{power+1} = MatricesMatrix;

        Test = ((omega*eye(length(Diag)) - Diag) \ upperTriag)^power;
        if Factor ~= Test(row, column)
            disp(['The factors do NOT match!', ' row: ', num2str(row), ', column: ', num2str(column), ', power: ', num2str(power)])
        end
    end
    for i = 1:numel(Matrices)
        Test = ((omega*eye(length(Diag)) - Diag) \ upperTriag)^(i-1);
        if Matrices{i} ~= Test
            disp(['The matrices do NOT match!', ' power: ', num2str(power)])
        end
    end
    invDiag = (omega*eye(length(Diag)) - Diag);
    ResultSchur = invDiag \ Matrix;
    Result = SchurVec * ResultSchur * SchurVec';
end

function [paths] = getPaths(row, column, power)
    disp(['row: ', num2str(row), ', column: ', num2str(column), ', power: ', num2str(power)])
    range = row+1 : column-1;
    if power > 0
        if length(range) == 1
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