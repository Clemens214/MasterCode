function [Results] = TransCheck(Energies, totalSystem, GammaL, GammaR, options)
% calculate the transmission through a molecule for zero temperature
arguments
    Energies
    totalSystem
    GammaL
    GammaR
    options.alternative = true
end
    %disp('Starting calculation of the current.')
    Results = zeros(1, length(Energies));
    for i = 1:length(Energies)
        if options.alternative == true
            Matrix = TransmissionAlt(Energies(i), totalSystem, GammaL, GammaR);
        elseif options.alternative == false
            Matrix = TransmissionLin(Energies(i), totalSystem, GammaL, GammaR);
        end
        % return the Result
        Results(i) = trace(real(Matrix));
    end
    %disp('Finished calculation of the current.')

    plotConductance(Energies, Results)
end

function [Result] = TransmissionLin(Energy, totalSystem, gammaL, gammaR)
    % calculate the Greens Function
    GreensFuncInv = Energy*eye(size(totalSystem)) - totalSystem;
    GreensFunc = inv(GreensFuncInv);
    
    % calculate the matrix product
    Result = GreensFunc * gammaL * GreensFunc' * gammaR;
end

function [Result] = TransmissionAlt(Energy, totalSystem, gammaL, gammaR)
    % GreensFunc * gammaL * GreensFunc' * gammaR
    GreensInv = Energy*eye(size(totalSystem)) - totalSystem;
    F = decomposition(GreensInv,'lu');    % create reusable LU object (works for sparse/dense)
    % GreensFunc^{-1} * gammaL
    Y = F \ gammaL;
    % gammaR * GreensFunc^{-1} * gammaL
    W = gammaR * Y;
    % GreensFunc'^{-1} * gammaR * GreensFunc^{-1} * gammaL
    Z = F' \ W;
    % => GreensFunc^{-1} * gammaL * GreensFunc'^{-1} * gammaR
    Result = Z;
end

%% plotting function
function [] = plotConductance(Energies, Vals)
    figure('Name','Conductance', 'NumberTitle','off');
    % plot
    plot(Energies, Vals)
    % labels
    xlabel('Energy (units of t)');
    ylabel('T(E)');
    title('Conductance');
end