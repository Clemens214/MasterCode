function [Eigenvals, leftEVs, rightEVs, varargout] = eigenvectors (totalSystem, options)
arguments
    totalSystem
    options.check = true
    options.checkMore = false
end
    %disp('Starting calculation of the Eigenvectors.')
    [rightEVs, Eigenvals, leftEVs] = eig(totalSystem);
    %disp('Finished calculation of the Eigenvectors.')
    
    %disp('Starting normalization of the Eigenvectors.')
    [leftEVs, rightEVs] = normalize(leftEVs, rightEVs);
    %disp('Finished normalization of the Eigenvectors.')
    
    %disp('Start checking the Eigenvectors.')
    if options.check == true
        [Product, maxOffDiag, minOnDiag] = TestProduct (Eigenvals, leftEVs, rightEVs);
        varargout{1} = Product;
        if options.checkMore == true
            CorrVal = 2;
            varargout{2} = maxOffDiag;
            disp(['maxOffDiag = ', num2str(maxOffDiag)])
            varargout{3} = minOnDiag;
            disp(['minOnDiag = ', num2str(minOnDiag)])
        else
            CorrVal = 0;
        end
        [MatchLeft, MatchRight, DiffLeft, DiffRight] = TestEV(totalSystem, Eigenvals, leftEVs, rightEVs);
        varargout{2+CorrVal} = MatchLeft;
        varargout{3+CorrVal} = MatchRight;
        varargout{4+CorrVal} = DiffLeft;
        varargout{5+CorrVal} = DiffRight;
    end
    %disp('Finished checking the Eigenvectors.')
end

function [leftEVs, rightEVs] = normalize (leftEVs, rightEVs)
    for i = 1:length(leftEVs)
        leftEV = leftEVs(:,i)';
        rightEV = rightEVs(:,i);
        value = leftEV * rightEV;
        
        squareRoot = sqrt(value);
        leftEVs(:,i) = leftEVs(:,i)/squareRoot';
        rightEVs(:,i) = rightEVs(:,i)/squareRoot;
    end
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

function [TestLeft, TestRight, DiffLeft, DiffRight] = TestEV (Matrix, Eigenvals, leftEVs, rightEVs)
    Tolerance = 1e-10;

    % check the left Eigenvectors
    % W = full matrix W whose columns are the corresponding left eigenvectors
    % W'*A = D*W'
    MatrixMultLeft = leftEVs' * Matrix;
    EigMultLeft = Eigenvals * leftEVs';
    TestLeft = isapprox(MatrixMultLeft, EigMultLeft, AbsoluteTolerance=Tolerance);
    DiffLeft = MatrixMultLeft-EigMultLeft;

    % check the right Eigenvectors
    % V = matrix V whose columns are the corresponding right eigenvectors
    % A*V = V*D
    MatrixMultRight = Matrix * rightEVs;
    EigMultRight = rightEVs * Eigenvals;
    TestRight = isapprox(MatrixMultRight, EigMultRight, AbsoluteTolerance=Tolerance);
    DiffRight = MatrixMultRight-EigMultRight;
end