function [Result, varargout] = MatrixCalcFull(leftEVs, rightEVs, gammaL, gammaR, options)
    %calculates all the different matrix elements
    arguments
        leftEVs
        rightEVs
        gammaL
        gammaR
        options.Transmission = false
        options.Torque = false
        options.conservative = false
        options.nonconservative
        options.saveMemory = true
    end
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
    
    if options.Transmission == true && options.Torque == false
        MatricesTrans = matrixTrans(index, gammaL, gammaR);
        Result = MatricesTrans;
    elseif options.Transmission == false && options.Torque == true
        MatricesTorque = matrixTorque(index, gammaL, gammaR);
        Result = MatricesTorque;
    end
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