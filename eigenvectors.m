function [Eigenvals, leftEVs, rightEVs, Product, MatchLeft, MatchRight, DiffLeft, DiffRight] = eigenvectors (totalSystem)
    %disp('Starting calculation of the Eigenvectors.')
    [rightEVs, Eigenvals, leftEVs] = eig(totalSystem);
    %disp('Finished calculation of the Eigenvectors.')
    
    %disp('Starting normalization of the Eigenvectors.')
    [leftEVs, rightEVs] = normalize(leftEVs, rightEVs);
    %disp('Finished normalization of the Eigenvectors.')
    
    %disp('Start checking the Eigenvectors.')
    Product = zeros(length(Eigenvals), length(Eigenvals));
    for i = 1:length(Eigenvals)
        leftEV = leftEVs(:,i)';
        rightEV = rightEVs(:,i);
        Product = Product + (rightEV * leftEV);
    end
    
    maxOffDiag = 0;
    minOnDiag = 1;
    for i = 1:length(Eigenvals)
        if abs(Product(i,i)) < minOnDiag
            minOnDiag = abs(Product(i,i));
        end
        for j = 1:length(Eigenvals)
            if i ~= j && abs(Product(i,j)) > maxOffDiag
                maxOffDiag = abs(Product(i,j));
            end
        end
    end
    %disp(['maxOffDiag = ', num2str(maxOffDiag)])
    %disp(['minOnDiag = ', num2str(minOnDiag)])

    [MatchLeft, MatchRight, DiffLeft, DiffRight] = TestEVquick(totalSystem, Eigenvals, leftEVs, rightEVs);
    %disp('Finished checking the Eigenvectors.')
end

%% 
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

%% 
function [TestLeft, TestRight, DiffLeft, DiffRight] = TestEVquick (Matrix, Eigenvals, leftEVs, rightEVs)
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

%% 
function [TestLeft, TestRight, DiffLeft, DiffRight] = TestEV (Matrix, Eigenvals, leftEVs, rightEVs)
    TestLeft = zeros(length(Eigenvals), length(Eigenvals));
    DiffLeft = zeros(length(Eigenvals), length(Eigenvals));
    TestRight = zeros(length(Eigenvals), length(Eigenvals));
    DiffRight = zeros(length(Eigenvals), length(Eigenvals));
    
    Tolerance = 1e-15;
    for i = 1:length(Eigenvals)
        % check the left Eigenvectors
        leftEV = leftEVs(i,:);
        MatrixMultLeft = leftEV * Matrix;
        EigMultLeft = Eigenvals(i) * leftEV;
        TestLeft(:,i) = isapprox(MatrixMultLeft, EigMultLeft, AbsoluteTolerance=Tolerance);
        DiffLeft(:,i) = MatrixMultLeft-EigMultLeft;
        % check the right Eigenvectors
        rightEV = rightEVs(:,i);
        MatrixMultRight = Matrix * rightEV;
        EigMultRight = Eigenvals(i) * rightEV;
        TestRight(:,i) = isapprox(MatrixMultRight, EigMultRight, AbsoluteTolerance=Tolerance);
        DiffRight(:,i) = MatrixMultRight-EigMultRight;
    end
end