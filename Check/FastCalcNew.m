%% Set the variables
% variables for the sample
sizeSample = 48; %48
orderSample = 2;
eigenenergy = 0;
hopping = 1;
hoppingsSample = hopping*eye(orderSample);

% variables for the leads
sizeLead = 104; %128
maxVal = 1;
decay = 0.2;
offset = 32;
[leadVals, derivVals] = calcVals(maxVal, decay, offset);
hoppingLead = hopping;

% angles
angleMax = 2*pi;
angleStep = pi/8;
angles = makeList(angleMax, angleStep);

% voltages
voltMax = 4;
voltStep = 0.5;
voltages = makeList(voltMax, voltStep);

% energies
EnergyMax = 2.5;
EnergyStep = 0.1;
Energies = makeList(EnergyMax, EnergyStep, full=true);

%% Calculate all the inputs
disp('Start preparing all the variables.')
SystemVals = setupSystem(sizeSample, orderSample, sizeLead, eigenenergy, hopping, hoppingsSample, hoppingLead, leadVals, derivVals, angles);
chemPots = setupPots(voltages);
disp('Finished preparing all the variables.')

%% Compute the results
disp('Starting the calculation.')
Transmission = TransCalc(SystemVals, chemPots, angles, voltages);
%Torque = TorqueCalc(SystemVals, chemPots, angles, voltages);
%TorqueC = TorqueCalc(SystemVals, chemPots, angles, voltages, conservative=true);
%TorqueNC = TorqueCalc(SystemVals, chemPots, angles, voltages, nonconservative=true);
%TorqueL = TorqueCalc(SystemVals, chemPots, angles, voltages, conservative=true);
%TorqueR = TorqueCalc(SystemVals, chemPots, angles, voltages, conservative=true);
disp('Finished the calculation.')

%% plot the data
Plot(angles, voltages, Transmission, threeD=true, Num=1, Title='Transmission', Names={'Transmission'})
%Plot(angles, voltages, Torque, threeD=true, Num=1, Title='total Torque', Names={'Torque'})
%Plot(angles, voltages, TorqueC, threeD=true, Num=1, Title='conservative Torque', Names={'Torque'})
%Plot(angles, voltages, TorqueNC, threeD=true, Num=1, Title='nonconservative Torque', Names={'Torque'})
%Plot(angles, voltages, TorqueL, threeD=true, Num=1, Title='left Torque', Names={'Torque'})
%Plot(angles, voltages, TorqueR, threeD=true, Num=1, Title='right Torque', Names={'Torque'})
savefig("TransmissionFast.fig")

%% setup functions
function [SystemVals] = setupSystem(sizeSample, orderSample, sizeLead, eigenenergy, hopping, hoppingsSample, hoppingLead, leadVals, derivVals, angles)
    % prepare the sample's Hamiltonian
    sample = makeSample(eigenenergy, hoppingsSample, sizeSample, orderSample);
    % calculate the total Hamiltonian
    SystemVals = struct('totalSystem', [], 'totalSysDeriv', [], 'gammaL', [], 'gammaR', [], ...
                        'EigenVals', [], 'leftEVs', [], 'rightEVs', []);
    for i = 1:length(angles)
        % define the hopping terms betwenn the leads and the sample
        hoppingsInter = zeros(2, orderSample);
        if orderSample == 2
            hoppingsInter = hopping*[cos(angles(i)), sin(angles(i)); 1, 0];
            hoppingsDeriv = hopping*[-1*sin(angles(i)), cos(angles(i)); 0, 0];
        elseif orderSample == 1
            hoppingsInter = hopping*[1; 1];
            hoppingsDeriv = hopping*[0; 0];
        end
        % preparing the Extended Molecule Hamiltonian
        [totalSystem, gammaL, gammaR] = makeSystem(sample, sizeSample, orderSample, sizeLead, ...
                                                    hoppingLead, hoppingsInter, leadVals);
        % prepare the dreivative of the Hamiltonian
        sampleDeriv = zeros(length(sample), length(sample));
        hoppingDeriv = 0;
        [totalSysDeriv, ~, ~] = makeSystem(sampleDeriv, sizeSample, orderSample, sizeLead, ...
                                            hoppingDeriv, hoppingsDeriv, derivVals, check=false);
        % compute the Eigenvectors and the Eigenvalues of the system
        [EigenVals, leftEVs, rightEVs] = eigenvectors(totalSystem);
        % return the calculated values
        SystemVals(i) = struct('totalSystem', totalSystem, 'totalSysDeriv', totalSysDeriv, ...
                            'gammaL', gammaL, 'gammaR', gammaR, ...
                            'EigenVals', EigenVals, 'leftEVs', leftEVs, 'rightEVs', rightEVs);
    end
end

function [chemPots] = setupPots(voltages)
    chemPots = struct('left', [], 'right', []);
    for j = 1:length(voltages)
        chemPotL = voltages(j)/2;
        chemPotR = -1*voltages(j)/2;
        chemPots(j) = struct('left', chemPotL, 'right', chemPotR);
    end
end

%% helping functions
function [values] = makeList(maxVal, stepVal, options)
    arguments
        maxVal 
        stepVal 
        options.full = false
    end
    if options.full == false
        minVal = 0;
    else
        minVal = -1*maxVal;
    end
    numVal = (maxVal-minVal)/stepVal + 1;
    values = linspace(minVal, maxVal, numVal);
end

function [leadVals, derivVals] = calcVals(maxVal, decay, offset)
    arguments
        maxVal = 1
        decay = 0.3
        % offset should be at most half the length of the leads
        % normal: 32
        offset = 32
    end
    leadVals = {maxVal, decay, offset};
    derivVals = {0, decay, offset};
end

function [matrix] = cutMatrix(Matrix, sizeLead, sizeSample)
    startIndex = sizeLead;
    endIndex = sizeLead+sizeSample + 1;
    matrix = Matrix(startIndex:endIndex, startIndex:endIndex);
    %disp(['Indices: from ',num2str(startIndex),' to ',num2str(endIndex)])
end