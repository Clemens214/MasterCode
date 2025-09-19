function [Result, varargout] = linearApprox(totalSystem, totalSysDeriv, gammaL, gammaR, options)
    %calculates the transport/torque through by a molecule in the linear transport approximation
    arguments
        totalSystem
        totalSysDeriv
        gammaL
        gammaR
        options.evalIndex = 1
        options.Transmission = false
        options.Torque = false
    end
    % calculate the transport/torque matrix and the traces
    Matrix = calc(totalSystem, totalSysDeriv, gammaL, gammaR, ...
                    Transmission=options.Transmission, Torque=options.Torque);
    Traces = allTrace(real(Matrix));
    
    % calculate the final results
    Result = Traces(options.evalIndex);
    varargout{1} = Traces;
    varargout{2} = real(Matrix);
end

%% linear transport approximation
function [Result, varargout] = calc(totalSystem, totalSysDeriv, gammaL, gammaR, options)
    % calculate the transmission and torque in the linear transport approximation
    arguments
        totalSystem
        totalSysDeriv
        gammaL
        gammaR
        options.Transmission = false
        options.Torque = false
    end
        
    if options.Transmission == true && options.Torque == false
        TransMatrix = TransmissionZeroTemp(totalSystem, gammaL, gammaR);
        
        Result = TransMatrix;
    elseif options.Transmission == false && options.Torque == true
        TorqueMatrixL = TorqueZeroTemp(totalSystem, totalSysDeriv, gammaL);
        TorqueMatrixR = TorqueZeroTemp(totalSystem, totalSysDeriv, gammaR);
        
        Result = TorqueMatrixL + TorqueMatrixR;
    else
        TransMatrix = TransmissionZeroTemp(totalSystem, gammaL, gammaR);
        TorqueMatrixL = TorqueZeroTemp(totalSystem, totalSysDeriv, gammaL);
        TorqueMatrixR = TorqueZeroTemp(totalSystem, totalSysDeriv, gammaR);
        
        Result = TransMatrix;
        varargout{1} = TorqueMatrixL + TorqueMatrixR;
    end
end

function [Transmission] = TransmissionZeroTemp(totalSystem, gammaL, gammaR)
    % calculate the Greens Function
    GreensFuncInv = eye(length(totalSystem)) - totalSystem;
    GreensFunc = inv(GreensFuncInv);
    
    % calculate the matrix product
    Transmission = GreensFunc * gammaR * GreensFunc' * gammaL;
end

function [Transmission] = TorqueZeroTemp(totalSystem, totalSysDeriv, gamma)
    % calculate the Greens Function
    GreensFuncInv = eye(length(totalSystem)) - totalSystem;
    GreensFunc = inv(GreensFuncInv);
    
    % calculate the matrix product
    Transmission = totalSysDeriv * GreensFunc * gamma * GreensFunc';
end

%% helping functions
function [Transmissions] = allTrace(Matrix)
    Transmissions = zeros(1, floor(length(Matrix)/2));
    for i = 1:floor(length(Matrix)/2)
        startIndex = i;
        endIndex = length(Matrix) - (i-1);
        matrix = Matrix(startIndex:endIndex, startIndex:endIndex);
        Transmissions(i) = trace(matrix);
    end
end