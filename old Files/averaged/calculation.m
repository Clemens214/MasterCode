function [Results, varargout] = calculation(sizeSample, orderSample, eigenenergy, hoppingsSample, sizeLead, leadVals, hoppingLead, hopping, voltages, angleDiff, options)
    % do the calculations for the transmission or torque
    arguments
        sizeSample
        orderSample
        eigenenergy
        hoppingsSample
        sizeLead
        leadVals
        hoppingLead
        hopping
        voltages
        angleDiff
        options.transmission = true
    end

    %preparing the sample's Hamiltonian
    sample = makeSample(eigenenergy, hoppingsSample, sizeSample, orderSample);
    
    Results = cell(1, length(angleDiff));
    ResultTraces = cell(1, length(angleDiff));
    ResultMatrix = cell(1, length(angleDiff));
    for i = 1:length(angleDiff)
        % define the hopping terms
        hoppingsInter = zeros(2, orderSample);
        if orderSample == 2
            hoppingsInter = hopping*[cos(angleDiff(i)), sin(angleDiff(i)); 1, 0];
            %hoppingsInter = hopping*[1, 0; 1, 0];
            hoppingsDeriv = hopping*[-1*sin(angleDiff(i)), cos(angleDiff(i)); 0, 0];
        elseif orderSample == 1
            hoppingsInter = hopping*[1; 1];
            hoppingsDeriv = hopping*[0; 0];
        end
        %disp(hoppingsInter(1,:))
    
        % preparing the Extended Molecule Hamiltonian
        [totalSystem, gammaL, gammaR] = makeSystem(sample, sizeSample, sizeLead, hoppingLead, hoppingsInter, leadVals);

        % prepare the derived Hamiltonian
        if options.transmission == false
            sampleDeriv = zeros(length(sample), length(sample));
            derivVals = {0, decay, offset};
            [totalSysDeriv, ~, ~] = makeSystem(sampleDeriv, sizeSample, sizeLead, 0, hoppingsDeriv, derivVals, check=false);
        end

        %compute the Eigenvectors and the Eigenvalues of the system
        [Eigenvals, leftEVs, rightEVs] = eigenvectors(totalSystem);
    
        Results{i} = zeros(1, length(voltages));
        ResultTraces{i} = cell(1, length(voltages));
        ResultMatrix{i} = cell(1, length(voltages));
        for j = 1:length(voltages)
            chemPotL = voltages(j)/2;
            chemPotR = -1*voltages(j)/2;
            chemPots = {chemPotL, chemPotR};
    
            evalIndex = floor(sizeLead/2);
            if options.transmission == true
                % compute the transmission
                disp('Starting calculation of the Transmission.')
                [Transmission, TransTrace, TransFull] = transmission(totalSystem, Eigenvals, leftEVs, rightEVs, gammaL, gammaR, chemPots);
                disp(['Finished calculation of the Transmission. Angle1 = ',num2str(angleDiff(i)),' i=',num2str(i), ...
                                                            ', Voltage = ',num2str(voltages(j)),' j=',num2str(j) ...
                                                            ', Transmission = ',num2str(Transmission)])
                Results{i}(j) = Transmission;
                ResultTraces{i}{j} = flip(TransTrace(1:sizeLead+1));
                ResultMatrix{i}{j} = cutMatrix(TransFull, sizeLead, length(sample));
            else
                % compute the torque
                disp('Starting calculation of the Torque.')
                [Torque, TorqueTrace, TorqueFull] = torque(totalSystem, Eigenvals, leftEVs, rightEVs, totalSysDeriv, gammaL, gammaR, chemPots);
                disp(['Finished calculation of the Torque. Angle1 = ',num2str(angleDiff(i)),' i=',num2str(i), ...
                                                            ', Voltage = ',num2str(voltages(j)),' j=',num2str(j) ...
                                                            ', Torque = ',num2str(Torque)])
                Results{i}(j) = Torque;
                ResultTraces{i}{j} = flip(TorqueTrace(1:sizeLead+1));
                ResultMatrix = cutMatrix(TorqueFull, sizeLead, length(sample));
            end
        end
    end
    varargout{1} = ResultTraces;
    varargout{2} = ResultMatrix;
end

%% helping function
function [matrix] = cutMatrix(Matrix, sizeLead, sizeSample)
    startIndex = sizeLead;
    endIndex = sizeLead+sizeSample + 1;
    matrix = Matrix(startIndex:endIndex, startIndex:endIndex);
    %disp(['Indices: from ',num2str(startIndex),' to ',num2str(endIndex)])
end