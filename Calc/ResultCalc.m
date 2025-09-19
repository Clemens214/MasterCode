function [Result, varargout] = ResultCalc(Matrices, Factors, options)
    %calculates the transport/torque through by a molecule
    arguments
        Matrices
        Factors
        options.evalIndex = 1
        options.Transmission = false
        options.Torque = false
    end
    % calculate the transport/torque matrix and the traces
    Matrix = mult(Matrices, Factors);
    %Traces = allTrace(real(Matrix));
    Traces = allTrace(real(Matrix));
    
    % calculate the final results
    Result = Traces(options.evalIndex);
    varargout{1} = Traces;
    varargout{2} = real(Matrix);
end

%% transmission
function [Result] = mult(Matrices, Factors)
    Size = size(Matrices);
    Result = zeros(Size(1), 1);
    for idx = 1:length(Matrices)
        Result = Result + Matrices(:,idx)*Factors(idx);
        %disp(['Finished calculation of the Result for idx=',num2str(idx)])
    end
end

%% helping functions
function [Transmissions] = allTrace(Matrix)
    Transmissions = zeros(1, floor(length(Matrix)/2));
    for i = 1:floor(length(Matrix)/2)
        startIndex = i;
        endIndex = length(Matrix) - (i-1);
        matrix = Matrix(startIndex:endIndex);
        Transmissions(i) = sum(matrix);
    end
end