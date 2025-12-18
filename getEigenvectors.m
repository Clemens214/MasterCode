function [Eigenvals, leftEVs, rightEVs, varargout] = getEigenvectors(totalSystem, options)
arguments
    totalSystem
    options.check = true
    options.checkMore = false
    options.returnAbs = true
end
    % [V,D,W] = eig(A)
    % returns a diagonal matrix D of eigenvalues 
    % and matrix V whose columns are the corresponding right eigenvectors, so that A*V = V*D.
    % also returns full matrix W whose columns are the corresponding left eigenvectors, so that W'*A = D*W'.
    [rightEVs, Eigenvals, leftEVs] = eig(totalSystem);
    
    [leftEVs, rightEVs, valuesOld, valuesNew] = normalize(leftEVs, rightEVs);
    
    %disp('Start checking the Eigenvectors.')
    if options.check == true
        % check the product
        [Product, maxOffDiag, minOnDiag] = TestProduct(Eigenvals, leftEVs, rightEVs);
        if options.checkMore == true
            disp(['maxOffDiag = ', num2str(maxOffDiag)])
            disp(['minOnDiag = ', num2str(minOnDiag)])
            absValuesOld = sort(abs(valuesOld));
            absValuesNew = sort(abs(valuesNew));
            if maxOffDiag > 10
                disp(['The values before normalization: \n', absValuesOld])
                disp(['The values after normalization: \n', absValuesNew])
            end
        end
        % return the values
        if options.returnAbs == true
            varargout{1} = abs(Product);
        else
            varargout{1} = Product;
        end
        if options.checkMore == true
            varargout{2} = maxOffDiag;
            varargout{3} = minOnDiag;
        end
    end
    if options.checkMore == true
        CorrVal = 2;
    else
        CorrVal = 0;
    end
    if options.check == true
        % check the eigenvectors
        [MatchLeft, MatchRight, DiffLeft, DiffRight] = TestEV(totalSystem, Eigenvals, leftEVs, rightEVs);
        if options.checkMore == true
            if all(MatchLeft) == true
                disp(['All the left eigenvectors match! Maximum Difference: ', num2str(max(max(DiffLeft)))]);
            end
            if all(MatchRight) == true
                disp(['All the right eigenvectors match! Maximum Difference: ', num2str(max(max(DiffRight)))]);
            end
        end
        % return the values
        varargout{2+CorrVal} = MatchLeft;
        varargout{3+CorrVal} = MatchRight;
        varargout{4+CorrVal} = DiffLeft;
        varargout{5+CorrVal} = DiffRight;
    end
    if options.checkMore == true
        checkResult(totalSystem)
    end
    %disp('Finished checking the Eigenvectors.')
end

function [leftEVs, rightEVs, varargout] = normalize (leftEVs, rightEVs)
    valuesOld = zeros(1, length(leftEVs));
    valuesNew = zeros(1, length(leftEVs));
    for i = 1:length(leftEVs)
        leftEV = leftEVs(:,i)';
        rightEV = rightEVs(:,i);
        valueOld = leftEV * rightEV;
        
        squareRoot = sqrt(valueOld);
        leftEVs(:,i) = leftEVs(:,i)/squareRoot';
        rightEVs(:,i) = rightEVs(:,i)/squareRoot;
        valueNew = leftEV * rightEV;

        valuesOld(i) = valueOld;
        valuesNew(i) = valueNew;
    end
    varargout{1} = valuesOld;
    varargout{2} = valuesNew;
end

%% checking functions
function [Product, maxOffDiag, minOnDiag] = TestProduct (Eigenvals, leftEVs, rightEVs)
    Product = zeros(length(Eigenvals), length(Eigenvals));
    for i = 1:length(Eigenvals)
        leftEV = leftEVs(:,i)';
        rightEV = rightEVs(:,i);
        Product = Product + (rightEV * leftEV);
    end
    maxOffDiag = 0;
    minOnDiag = 1;
    for i = 1:length(Eigenvals)
        for j = 1:length(Eigenvals)
            if i ==j && abs(Product(i,j)) < minOnDiag
                minOnDiag = abs(Product(i,i));
            elseif i ~= j && abs(Product(i,j)) > maxOffDiag
                maxOffDiag = abs(Product(i,j));
            end
        end
    end
end

function [TestLeft, TestRight, DiffLeft, DiffRight] = TestEV (Matrix, Eigenvals, leftEVs, rightEVs, options)
arguments
    Matrix 
    Eigenvals 
    leftEVs 
    rightEVs 
    options.returnAbs = true
    options.Tolerance = 1e-10;
end
    % check the left Eigenvectors
    % W = full matrix W whose columns are the corresponding left eigenvectors
    % W'*A = D*W'
    MatrixMultLeft = leftEVs' * Matrix;
    EigMultLeft = Eigenvals * leftEVs';
    TestLeft = isapprox(MatrixMultLeft, EigMultLeft, AbsoluteTolerance=options.Tolerance);
    DiffLeft = MatrixMultLeft-EigMultLeft;

    % check the right Eigenvectors
    % V = matrix V whose columns are the corresponding right eigenvectors
    % A*V = V*D
    MatrixMultRight = Matrix * rightEVs;
    EigMultRight = rightEVs * Eigenvals;
    TestRight = isapprox(MatrixMultRight, EigMultRight, AbsoluteTolerance=options.Tolerance);
    DiffRight = MatrixMultRight-EigMultRight;

    if options.returnAbs == true
        DiffLeft = abs(DiffLeft);
        DiffRight = abs(DiffRight);
    end
end

%% more checking functions
function [] = checkResult(totalSystem)
    % Example matrix (replace with your A)
    A = totalSystem;

    % 1) Compute right eigenvectors and eigenvalues
    [V, D] = eig(A);      % A * V = V * D
    
    % 2) Compute left eigenvectors (right eigenvectors of A.')
    [W, ~] = eig(A.');    % A.' * W = W * D   (columns of W are left eigenvectors)
    
    % 3) Enforce biorthonormality: W' * V should be the identity.
    %    Use matrix right-division to avoid explicit inv when possible:
    W = W / (W' * V);     % MATLAB: W * inv(W'*V). After this, W' * V ~ I
    
    % 4) Reconstruct A from the decomposition
    A_rec = V * D * W';   % uses transpose-conjugate of W (W' is conjugate transpose)
    
    % 5) Diagnostics: check biorthogonality and reconstruction error
    biorth_error = norm(W' * V - eye(size(A)), 'fro');     % should be ~0
    recon_error   = norm(A - A_rec, 'fro');                % should be small if diagonalizable
    
    fprintf('||W''*V - I||_F = %g\n', biorth_error);
    fprintf('||A - A_rec||_F = %g\n', recon_error);
    
    % Optional: also check condition number of V (large -> unstable)
    fprintf('cond(V) = %g\n', cond(V));

    % Conditioning of eigenvector matrix
    condV = cond(V);
    
    % Rank of eigenvector matrix
    rankV = rank(V);
    
    % Check if matrix is defective
    if rankV < size(A,1)
        disp('Matrix appears defective (not diagonalizable).');
    else
        disp('Matrix is diagonalizable but ill-conditioned.');
    end
    fprintf('Condition number of V: %g\n', condV);
end