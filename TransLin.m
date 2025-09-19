function [Transmissions, varargout] = TransLin(SystemVals, chemPots, angles, Energies, option)
    %calculates the transport through by a molecule in the linear transport approximation
    arguments
        SystemVals
        chemPots
        angles
        Energies
        option.allVals = true
    end
    % calculate the matrices in the linear transport approximation
    if option.allVals == true
        Transmissions = cell(1, length(angles));
    else
        Transmissions = zeros(1, length(angles));
    end
    Traces = cell(1, length(angles));
    for i = 1:length(angles)
        totalSystem = SystemVals(i).totalSystem;
        gammaL = SystemVals(i).gammaL;
        gammaR = SystemVals(i).gammaR;
        % compute the transmission
        disp('Start calculation of the transmission.')
        if option.allVals == true
            [Transmission, TransTraces] = AllCalc(Energies, totalSystem, gammaL, gammaR);
            Transmissions{i} = Transmission;
        else
            [Transmission, TransTraces] = ValCalc(Energies, totalSystem, gammaL, gammaR);
            Transmissions(i) = Transmission;
        end
        Traces{i} = TransTraces;
        disp(['Finished calculation of the transmission.', ...
            ' Angle = ',num2str(angles(i)),', i=',num2str(i)])
    end
    varargout{1} = Traces;
end

%% total transmission in the linear transport approximation
function [Result, varargout] = ValCalc(Energies, totalSystem, gammaL, gammaR, options)
    %calculates the transport through a molecule in the linear transport approximation
    arguments
        Energies
        totalSystem
        gammaL
        gammaR
        options.EnergyVal = 0
        options.evalIndex = 1
    end
    % calculate the transport matrix and the traces
    Matrix = ValMatrix(Energies, totalSystem, gammaL, gammaR, ...
                        EnergyVal=options.EnergyVal);
    Traces = allTrace(real(Matrix));
    
    % calculate the final results
    Result = Traces(options.evalIndex);
    varargout{1} = Traces;
    varargout{2} = real(Matrix);
end

function [Result] = ValMatrix(Energies, totalSystem, gammaL, gammaR, options)
    % calculate the transmission in the linear transport approximation
    arguments
        Energies
        totalSystem
        gammaL
        gammaR
        options.EnergyVal = 0
    end
    EnergyVal = options.EnergyVal;
    
    % calculate the transmission
    Result = TransmissionZeroTemp(EnergyVal, totalSystem, gammaL, gammaR);
end

%% conductances in the linear transport approximation
function [Results, varargout] = AllCalc(Energies, totalSystem, gammaL, gammaR, options)
    %calculates the transport through a molecule in the linear transport approximation
    arguments
        Energies
        totalSystem
        gammaL
        gammaR
        options.evalIndex = 1
    end
    % calculate the transport matrix and the traces
    Matrices = AllMatrix(Energies, totalSystem, gammaL, gammaR);

    % calculate the final results
    allResults = zeros(1, length(Matrices));
    allTraces = cell(1, length(Matrices));
    allMatrices = cell(1, length(Matrices));
    for i = 1:length(allMatrices)
        Traces = allTrace(real(Matrices{i}));
        allResults(i) = Traces(options.evalIndex);
        allTraces{i} = Traces;
        allMatrices{i} = real(Matrices{i});
    end
    
    % return the results
    Results = allResults;
    varargout{2} = allTraces;
    varargout{3} = allMatrices;
end

function [Results] = AllMatrix(Energies, totalSystem, gammaL, gammaR)
    % calculate the conductance matrices in the linear transport approximation
    arguments
        Energies
        totalSystem
        gammaL
        gammaR
    end    
    % calculate all the different Transmission matrix values
    Results = cell(1, length(Energies));
    for i = 1:length(Energies)
        TransMatrix = TransmissionZeroTemp(Energies(i), totalSystem, gammaL, gammaR);
        Results{i} = TransMatrix;
    end
end

%% calculation of the transport matrices
function [Transmission] = TransmissionZeroTemp(Energy, totalSystem, gammaL, gammaR)
    % calculate the Greens Function
    GreensFuncInv = Energy*eye(length(totalSystem)) - totalSystem;
    GreensFunc = inv(GreensFuncInv);
    
    % calculate the matrix product
    Transmission = GreensFunc * gammaR * GreensFunc' * gammaL;
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