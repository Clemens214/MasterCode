function [Eigenvals, leftEVs, rightEVs, varargout] = getEigenvectors(totalSystem, options)
% [V,D,W] = eig(A)
% returns a diagonal matrix D of eigenvalues 
% and matrix V whose columns are the corresponding right eigenvectors, so that A*V = V*D.
% also returns full matrix W whose columns are the corresponding left eigenvectors, so that W'*A = D*W'.
arguments
    totalSystem
    options.check = true
    options.checkMore = false
end
    % calculate the (left and right) Eigenvectors and Eigenvalues
    % normalize the Eigenvectors to enable reconstruction of the original matrix
    [rightEVs, Eigenvals, leftEVs] = eig(totalSystem);
    [leftEVs, rightEVs] = normalize(leftEVs, rightEVs);
    
    if options.check== true
        [Test, Diff, maxDiff] = checkResult(totalSystem, Eigenvals, leftEVs, rightEVs);
        varargout{1} = Test;
        varargout{2} = Diff;
        varargout{3} = maxDiff;
        % return the checking results
        allTest = all(all(Test));
        if options.checkMore == true && allTest == true
            disp(['The Eigenvalue decomposition DOES reproduce the Hamiltonian. Maximum Difference: ', num2str(maxDiff)])
        elseif options.checkMore == true && allTest == false
            disp(['The Eigenvalue decomposition does NOT reproduce the Hamiltonian. Maximum Difference: ', num2str(maxDiff)])
        elseif options.checkMore == false && allTest == false
            error(['The Eigenvalue decomposition does NOT reproduce the Hamiltonian. Maximum Difference: ', num2str(maxDiff)])
        end
    end
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
function [Test, Diff, maxDiff] = checkResult(totalSystem, Eigenvals, leftEVs, rightEVs, options)
% [V,D,W] = eig(A)
% returns a diagonal matrix D of eigenvalues 
% and matrix V whose columns are the corresponding right eigenvectors, so that A*V = V*D.
% also returns full matrix W whose columns are the corresponding left eigenvectors, so that W'*A = D*W'.
arguments
    totalSystem 
    Eigenvals 
    leftEVs 
    rightEVs
    options.returnAbs = true
    options.Tolerance = 1e-10;
    %Tolerance = 1e-10;
end
    Matrix = leftEVs' * Eigenvals * rightEVs;
    
    Test = isapprox(totalSystem, Matrix, AbsoluteTolerance=options.Tolerance);
    Diff = totalSystem - Matrix;
    if options.returnAbs == true
        Diff = abs(Diff);
    end
    maxDiff = max(max(abs(Diff)));
end