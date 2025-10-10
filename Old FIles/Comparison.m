function [TransmissionsEM, TransmissionsSI] = Transmission (omegas, Hamiltonian, totalSystemEM, gammaL_EM, gammaR_EM, hoppingInter, hoppingBath, lengthSample, lengthLead)
    disp('Starting calculation of the Transmisson values.')
    % calculate the Transmissions
    TransmissionsSI = zeros(length(omegas));
    TransmissionsEM = zeros(length(omegas));
    for i = 1:length(omegas)
        valueSI = 1;
        [GreensFuncSI, totalSystemSI, gammaL_SI, gammaR_SI] = prepareSI(omegas(i), Hamiltonian, hoppingInter, hoppingBath, lengthSample);
        %GreensFuncSI = GreensFunc(omegas(i), totalSystemSI);
        TransmissionsSI(i) = Transport(GreensFuncSI, gammaL_SI, gammaR_SI, valueSI);
        
        valueEM = lengthLead + 1;
        GreensFuncEM = GreensFunc(omegas(i), totalSystemEM);
        TransmissionsEM(i) = Transport(GreensFuncEM, gammaL_EM, gammaR_EM, valueEM);
    end
    disp('Finished calculation of the Transmisson values.')
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