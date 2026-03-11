function [Diag, upperTriag, SchurVec, varargout] = getSchur(totalSystem, options)
arguments
    totalSystem
    options.check = true
    options.checkMore = false
end
    % calculate the Schur-vectors
    % calculate the diagonal (Eigenvalues) and upper-triangular matrices
    [SchurVec, Triag] = schur(totalSystem, "complex");
    Diag = diag(diag(Triag));
    upperTriag = Triag - Diag;
    
    if options.check== true
        [Test, Diff, maxDiff] = checkResult(totalSystem, Diag, upperTriag, SchurVec);
        varargout{1} = Test;
        varargout{2} = Diff;
        varargout{3} = maxDiff;
        % return the checking results
        allTest = all(all(Test));
        if options.checkMore == true && allTest == true
            disp(['The Schur decomposition DOES reproduce the Hamiltonian. Maximum Difference: ', num2str(maxDiff)])
        elseif options.checkMore == true && allTest == false
            disp(['The Schur decomposition does NOT reproduce the Hamiltonian. Maximum Difference: ', num2str(maxDiff)])
        elseif options.checkMore == false && allTest == false
            error(['The Schur decomposition does NOT reproduce the Hamiltonian. Maximum Difference: ', num2str(maxDiff)])
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
    Diff = totalSystem - Matrix;
    if options.returnAbs == true
        Diff = abs(Diff);
    end
    maxDiff = max(max(abs(Diff)));
end