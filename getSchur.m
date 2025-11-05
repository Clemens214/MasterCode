function [Diag, upperTriag, SchurVec, varargout] = getSchur(totalSystem, options)
arguments
    totalSystem
    options.check = true
    options.checkMore = true
end
    % calculate the Schur-vectors
    % calculate the diagonal (Eigenvalues) and upper-triangular matrices
    [SchurVec, Triag] = schur(totalSystem, "complex");
    Diag = diag(diag(Triag));
    upperTriag = Traig - Diag;
    
    if options.check== true
        [Test, Diff, maxDiff] = checkResult(totalSystem, Diag, upperTriag, SchurVec, options);
        varargout{1} = Test;
        varargout{2} = Diff;
        varargout{3} = maxDiff;
        if options.checkMore == true
            allTest = all(Test);
            if allTest == true
                disp(['The Schur decomposition reproduces the Matrix. Maximum Difference: ', num2str(maxDiff)])
            else
                disp(['The Schur decomposition does NOT reproduce the Matrix. Maximum Difference: ', num2str(maxDiff)])
            end
        end
    end
end

%% checking functions
function [Test, Diff, maxDiff] = checkResult(totalSystem, Diag, upperTriag, SchurVec, options)
arguments
    totalSystem 
    Diag 
    upperTriag 
    SchurVec
    options.returnAbs = true
    options.Tolerance = 1e-10;
    %Tolerance = 1e-10;
end
    Matrix = SchurVec * (Diag + upperTriag) * SchurVec';
    
    Test = isapprox(totalSystem, Matrix, AbsoluteTolerance=options.Tolerance);
    Diff = MatrixMultLeft-EigMultLeft;
    if options.returnAbs == true
        Diff = abs(Diff);
    end
    maxDiff = max(max(abs(Diff)));
end