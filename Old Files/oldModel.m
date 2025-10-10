function [] = oldModel ()
    [TransmissionsEM] = Transmission(omegas, sample, totalSystemEM, gammaL_EM, gammaR_EM, hoppingInter, hoppingLead, lengthSample, lengthLead);

end

function [TransmissionsEM, TransmissionsSI] = Transmission (omegas, totalSystemEM, gammaL_EM, gammaR_EM, lengthLead)
    %disp('Starting calculation of the Transmisson values.')
    % calculate the Transmissions
    TransmissionsEM = zeros(length(omegas));
    for i = 1:length(omegas)
        valueEM = lengthLead + 1;
        GreensFuncEM = GreensFunc(omegas(i), totalSystemEM);
        TransmissionsEM(i) = Transport(GreensFuncEM, gammaL_EM, gammaR_EM, valueEM);
    end
    %disp('Finished calculation of the Transmisson values.')
end

%%
function [arrayInv] = GreensFunc (omega, totalSystem)
    array = eye(length(totalSystem))*omega - totalSystem;
    arrayInv = inv(array);
end

%%
function [T] = Transport(GreensFunc, gammaL, gammaR, value)
    % calculate the matrix product
    transport = GreensFunc * (gammaL - gammaR) * GreensFunc';
    % calculate the current
    T = imag(transport(value, value+1)); 
end