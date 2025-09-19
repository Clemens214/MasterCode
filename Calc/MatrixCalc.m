function [Result, varargout] = MatrixCalc(leftEVs, rightEVs, gammaL, gammaR, totalSysDeriv, midFactor, options)
    %calculates all the different matrix elements
    arguments
        leftEVs
        rightEVs
        gammaL
        gammaR
        totalSysDeriv
        midFactor
        options.Transmission = false
        options.Torque = false
        options.left = false
        options.right = false
        options.conservative = false
        options.nonconservative = false
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
    if options.Transmission == true && options.Torque == false
        MatricesTrans = matrixTrans(index, gammaL, gammaR);
        Result = MatricesTrans;
    elseif options.Transmission == false && options.Torque == true
        MatricesTorque = matrixTorque(index, totalSysDeriv, gammaL, gammaR, ...
                                    left=options.left, right=options.right, ...
                                    conservative=options.conservative, nonconservative=options.nonconservative);
        Result = MatricesTorque;
    end
end

function [Matrices, varargout] = matrixTorque(index, totalSysDeriv, gammaL, gammaR, options)
    arguments
        index
        totalSysDeriv
        gammaL
        gammaR
        options.left = false
        options.right = false
        options.conservative = false
        options.nonconservative = false
    end
    if options.conservative == true && options.nonconservative == false
        midFactor = 
        
    elseif options.conservative == false && options.nonconservative == true

    else
        
    end
    Matrices = matrixTorqueCalc(index, totalSysDeriv, gammaL, gammaR);
end

%% calculate the matrices and only return the diagonals
function [Matrices] = matrixTrans(index, gammaL, gammaR)
    %disp('Starting calculation of the current matrices.')
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
        %disp(['Finished calculating matrix ', num2str(idx)])
    end
    %disp('Finished calculation of the current matrices.')
end

function [Matrices] = matrixTorqueCalc(index, totalSysDeriv, midFactor)
    %disp('Starting calculation of the torque matrices.')
    Matrices = zeros(length(gammaL), length(index));
    parfor idx = 1:numel(index)
        % get the normal Eigenvectors
        leftEV = index(idx).leftEV;
        rightEV = index(idx).rightEV;
        
        % get the daggered Eigenvectors
        leftEVdagger = index(idx).leftEVD;
        rightEVdagger = index(idx).rightEVD;
        
        % compute the matrix element for chosen i and j
        ProductLeft = totalSysDeriv * rightEV;
        ProductMid = leftEV * midFactor * leftEVdagger;
        ProductRight = rightEVdagger;
        
        Product = ProductLeft * ProductMid * ProductRight;
        Matrices(:, idx) = diag(Product);
    end
    %disp('Finished calculation of the torque matrices.')
end