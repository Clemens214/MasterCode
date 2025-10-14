%% Variables

% variables for the sample
sizeSample = 48;
orderSample = 2;
eigenenergy = 0;
hopping = 1;
hoppingsSample = hopping*eye(orderSample);

% variables for the leads
sizeLead = 104;
[leadVals, derivVals] = calcVals(maxVal = 1, decay = 0.2, offset = 32);
hoppingLead = hopping;

% variables for the hopping
angleMax = 2*pi;
angleStep = pi/8;
angles = makeList(angleMax, angleStep);

%variables for the calculation of the current
chemPotMax = 0;%1;
chemPotStep = 1;
chemPots = makeList(chemPotMax, chemPotStep, full=true);

%variables for the calculation of the transmission
omegaVal = 2.5;
omegaMax = omegaVal*hopping;
omegaStep = 0.05;
omegas = makeList(omegaMax, omegaStep, full=true);

%% Calculation

Transmission = zeros(length(angles), length(angles));
Torque = zeros(length(angles), length(angles));
for i = 1:length(angles)
for j = 1:length(angles)
    if orderSample == 1
        hoppingsInter = [hopping; hopping];
    elseif orderSample == 2
        hoppingsInter = [cos(angles(i)), sin(angles(i)); cos(angles(j)), sin(angles(j))];
    end
    
    % compute the Hamiltonian of the Sample
    sample = makeSample(eigenenergy, hoppingsSample, sizeSample,  orderSample);
    
    % preparing the Extended Molecule Hamiltonian
    [totalSystem, gammaL, gammaR] = makeSystemEM(sample, sizeSample, orderSample, sizeLead, hoppingLead, hoppingsInter, leadVals);
    [totalSysDeriv, ~, ~] = makeSystemEM(sample, sizeSample, orderSample, sizeLead, hoppingLead, hoppingsInter, derivVals, check=false);
    
    %compute the Eigenvectors and the Eigenvalues of the system
    disp('Starting calculation of the Eigenvectors.')
    [Eigenvals, leftEVs, rightEVs] = eigenvectors(totalSystem);
    disp('Finished calculation of the Eigenvectors.')
    
    omegas = [1];
    for idx = 1:length(omegas)
        %TransmissionResult = transmission(totalSystem, gammaL, gammaR, Eigenvals, leftEVs, rightEVs, omegas(idx), linearResponse=true);
        %Transmission(i, j) = TransmissionResult;

        TorqueResult = torque(totalSystem, totalSysDeriv, gammaL, gammaR, Eigenvals, leftEVs, rightEVs, omegas(idx), linearResponse=true, conservative=true);
        Torque(i, j) = TorqueResult;

        disp([  'Angle1: ', num2str(angles(i)), ', i=', num2str(i), ...
                ', Angle2: ', num2str(angles(j)), ', j=', num2str(j)])
    end
end
end

%% plot
%plotLin3D (angles, Transmission)
%indices = plotLin2D (2, 'Transmission', angles, Transmission);

plotLin3D (angles, Torque)
indices = plotLin2D (2, 'Torque', angles, Torque);

%% plotting functions
function [varargout] = plotLin2D (value, Title, angles, Vals)
    TransPlot = zeros(length(angles)*length(angles));
    angleDiff = zeros(length(angles)*length(angles));
    indices = zeros(length(angles)*length(angles));
    for i = 1:length(angles)
        for j = 1:length(angles)
            idx = (i-1)*length(angles) + j;
            TransPlot(idx) = Vals(i, j);
            angleDiff(idx) = angles(i) - angles(j);
            indices(idx) = idx;
        end
    end
    varargout{1} = indices;
    [angleSort, indices] = sort(angleDiff);
    TransSort = TransPlot(indices);
    % plot the data
    figure(value);
    plot(angleSort, TransSort)
    title(Title);
end

function [] = plotLin3D (angles, Vals)
    surf(angles, angles, Vals)
end

function [] = plotGraph (value, Title, omegas, Vals, chemPots)
    figure(value);
    title(Title);
    for i = 1:length(chemPots)
        plot(omegas, Vals{i})
        hold on
    end
    labels = strcat('chemPot = ',cellstr(num2str(chemPots.')));
    legend(labels)
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
    numVal = (maxVal-minVal)/stepVal+1;
    values = linspace(minVal, maxVal, numVal);
end

function [leadVals, derivVals] = calcVals(opt)%(maxVal, decay, offset)
    arguments
        opt.maxVal = 1
        opt.decay = 0.3
        % offset should be at most half the length of the leads
        opt.offset = 32
        % normal: 32
    end
    maxVal = opt.maxVal;
    decay = opt.decay;
    offset = opt.offset;
    leadVals = {maxVal, decay, offset};
    derivVals = {0, decay, offset};
end

function [] = saveVar(var, order)
    filename = append('Indices', int2str(order), '.mat');
    save(filename, "var")
end