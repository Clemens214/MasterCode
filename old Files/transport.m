function [Result, varargout] = transport(totalSystem, gammaL, gammaR, value, step)
    %calculates the transport through a molecule as calculated in the
    %Project-Praktikum

    % define the states to be used for calculating the current
    stateLeft = zeros(1,length(totalSystem));
    stateLeft(value) = 1;
    stateRight = zeros(1,length(totalSystem))';
    stateRight(value+step) = 1;
    
    %compute the current
    midFactor = -1*(gammaL - gammaR);
    %disp('Starting calculation of the current.')
    TransmissionMatrix = Transmission(totalSystem, midFactor);
    %disp('Finished calculation of the current.')

    TransmissionElement = stateLeft * TransmissionMatrix * stateRight;
    Result = imag(TransmissionElement);

    varargout{1} = imag(TransmissionMatrix);
end

%% transmission at zero temperature
function [Transmission] = Transmission(totalSystem, midFactor)
    % calculate the Greens Function
    GreensFuncInv = eye(length(totalSystem)) - totalSystem;
    GreensFunc = inv(GreensFuncInv);
    
    % calculate the matrix product
    Transmission = GreensFunc * midFactor * GreensFunc';
end