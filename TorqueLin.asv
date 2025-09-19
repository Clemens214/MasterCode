function [Torques, varargout] = TorqueLin(SystemVals, chemPots, angles, Energies, opt)
    %calculates the transport/torque through by a molecule in the linear transport approximation
    arguments
        SystemVals
        chemPots
        angles
        Energies
        opt.allVals = true
        opt.conservative = false
        opt.nonconservative = false
        opt.left = false
        opt.right = false
    end
    options = struct('conservative', opt.conservative, ...
                        'nonconservative', opt.nonconservative, ...
                        'left', opt.left, 'right', opt.right);
    % calculate the matrices in the linear transport approximation
    if opt.allVals == true
        Torques = cell(1, length(angles));
    else
        Torques = zeros(1, length(angles));
    end
    Traces = cell(1, length(angles));
    for i = 1:length(angles)
        totalSystem = SystemVals(i).totalSystem;
        totalSysDeriv = SystemVals(i).totalSysDeriv;
        gammaL = SystemVals(i).gammaL;
        gammaR = SystemVals(i).gammaR;
        % compute the transmission
        disp('Start calculation of the torque.')
        if opt.allVals == true
            [Torque, TorqueTraces] = AllCalc(Energies, totalSystem, totalSysDeriv, gammaL, gammaR, options);
            Torques{i} = Torque;
        else
            [Torque, TorqueTraces] = ValCalc(Energies, totalSystem, totalSysDeriv, gammaL, gammaR, options);
            Torques(i) = Torque;
        end
        Traces{i} = TorqueTraces;
        disp(['Finished calculation of the torque.', ...
            ' Angle = ',num2str(angles(i)),', i=',num2str(i)])
    end
    varargout{1} = Traces;
end

%% total transmission in the linear transport approximation
function [Result, varargout] = ValCalc(Energies, totalSystem, totalSysDeriv, gammaL, gammaR, options, opt)
    %calculates the torque experienced by a molecule in the linear transport approximation
    arguments
        Energies
        totalSystem
        totalSysDeriv
        gammaL
        gammaR
        options
        opt.EnergyVal = 0
        opt.evalIndex = 1
    end
    % calculate the torque matrix and the traces
    Matrix = ValMatrix(Energies, totalSystem, totalSysDeriv, gammaL, gammaR, options, ...
                        EnergyVal=opt.EnergyVal);
    Traces = allTrace(real(Matrix));
    
    % calculate the final results
    Result = Traces(opt.evalIndex);
    varargout{1} = Traces;
    varargout{2} = real(Matrix);
end

function [Result] = ValMatrix(Energies, totalSystem, totalSysDeriv, gammaL, gammaR, options, opt)
    % calculate the torque matrix in the linear transport approximation
    arguments
        Energies
        totalSystem
        totalSysDeriv
        gammaL
        gammaR
        options
        opt.EnergyVal = 0
    end
    EnergyVal = opt.EnergyVal;

    % calculate the torque
    Result = MatrixChoice(EnergyVal, totalSystem, totalSysDeriv, gammaL, gammaR, options);
end

%% torquances in the linear transport approximation
function [Results, varargout] = AllCalc(Energies, totalSystem, totalSysDeriv, gammaL, gammaR, options, opt)
    %calculates the torquances of a molecule in the linear transport approximation
    arguments
        Energies
        totalSystem
        totalSysDeriv
        gammaL
        gammaR
        options
        opt.evalIndex = 1
    end
    % calculate the torque matrix and the traces
    Matrices = AllMatrix(Energies, totalSystem, totalSysDeriv, gammaL, gammaR, options);
    
    % calculate the final results
    allResults = zeros(1, length(Matrices));
    allTraces = cell(1, length(Matrices));
    allMatrices = cell(1, length(Matrices));
    for i = 1:length(allMatrices)
        Traces = allTrace(real(Matrices{i}));
        allResults(i) = Traces(opt.evalIndex);
        allTraces{i} = Traces;
        allMatrices{i} = real(Matrices{i});
    end

    % return the results
    Results = allResults;
    varargout{2} = allTraces;
    varargout{3} = allMatrices;
end

function [Results] = AllMatrix(Energies, totalSystem, totalSysDeriv, gammaL, gammaR, options)
    % calculate the torquance matrices in the linear transport approximation
    arguments
        Energies
        totalSystem
        totalSysDeriv
        gammaL
        gammaR
        options
    end
    % calculate all the different Torque matrix values
    Results = cell(1, length(Energies));
    for i = 1:length(Energies)
        TorqueMatrix = MatrixChoice(Energies(i), totalSystem, totalSysDeriv, gammaL, gammaR, options);
        Results{i} = TorqueMatrix;
    end
end

%% calculate of the torque matrices
function [TorqueMatrix] = MatrixChoice(Energy, totalSystem, totalSysDeriv, gammaL, gammaR, options)
    arguments
        Energy
        totalSystem
        totalSysDeriv
        gammaL
        gammaR
        options
    end
    if options.conservative == true || options.nonconservative == true || options.left == true || options.right == true
        if options.conservative == true
            midFactor = gammaL + gammaR;
        elseif options.nonconservative == true
            midFactor = gammaL - gammaR;
        elseif options.left == true
            midFactor = gammaL;
        elseif options.right == true
            midFactor = gammaR;
        end
        TorqueMatrix = TorqueZeroTemp(Energy, totalSystem, totalSysDeriv, midFactor);
    else
        TorqueMatrixL = TorqueZeroTemp(Energy, totalSystem, totalSysDeriv, gammaL);
        TorqueMatrixR = TorqueZeroTemp(Energy, totalSystem, totalSysDeriv, gammaR);
        TorqueMatrix = TorqueMatrixL + TorqueMatrixR;
    end
end

function [Torque] = TorqueZeroTemp(Energy, totalSystem, totalSysDeriv, midFactor)
    % calculate the Greens Function
    GreensFuncInv = Energy*eye(length(totalSystem)) - totalSystem;
    GreensFunc = inv(GreensFuncInv);
    
    % calculate the matrix product
    Torque = totalSysDeriv * GreensFunc * midFactor * GreensFunc';
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