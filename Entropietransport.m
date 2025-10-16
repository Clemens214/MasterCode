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
angleStep = pi/32;
angles = makeList(angleMax, angleStep);

%variables for the calculation of the current
voltageMax = 5;
voltageStep = 0.05;
voltages = makeList(voltageMax, voltageStep);

%% Calculation

Transmission = zeros(1, length(angles));
Torque = zeros(1, length(angles));
for i = 1:length(angles)
    if orderSample == 1
        hoppingsInter = [hopping; hopping];
        hoppingsDeriv = [0; 0];
    elseif orderSample == 2
        %hoppingsInter = [cos(angles(i)), sin(angles(i)); cos(angles(j)), sin(angles(j))];
        %hoppingsDeriv = [-1*sin(angles(i)), cos(angles(i)); -1*sin(angles(j)), cos(angles(j))];
        hoppingsInter = [cos(angles(i)), sin(angles(i)); 1, 0];
        hoppingsDeriv = [-1*sin(angles(i)), cos(angles(i)); 0, 0];
    end
    
    % compute the Hamiltonian of the Sample
    sample = makeSample(eigenenergy, hoppingsSample, sizeSample,  orderSample);
    
    % preparing the Extended Molecule Hamiltonian
    [totalSystem, gammaL, gammaR] = makeSystemEM(sample, sizeSample, orderSample, sizeLead, hoppingLead, hoppingsInter, leadVals);
    [totalSysDeriv, ~, ~] = makeSystemEM(zeros(size(sample)), sizeSample, orderSample, sizeLead, 0, hoppingsDeriv, derivVals, check=false);
    
    %compute the Eigenvectors and the Eigenvalues of the system
    disp('Starting calculation of the Eigenvectors.')
    [Eigenvals, leftEVs, rightEVs, Product] = eigenvectors(totalSystem, checkMore=true);
    disp('Finished calculation of the Eigenvectors.')
    
    voltages = [0];
    chemPots = setupPots(voltages);
    Transmission(i) = TransCalc(totalSystem, gammaL, gammaR, Eigenvals, leftEVs, rightEVs, chemPots, linearResponse=true);
    voltages = [2];
    chemPots = setupPots(voltages);
    chemPots(1).right = chemPots(1).left;
    Torque(i) = TorqueCalc(totalSystem, totalSysDeriv, gammaL, gammaR, Eigenvals, leftEVs, rightEVs, chemPots, linearResponse=true);
end

%% plot
plotAngle (1, 'Transmission', angles, Transmission);
plotAngle (2, 'Torque', angles, Torque);

%plotBoth (1, 'Transmission', angles, Transmission);
%plotBoth (3, 'Torque', angles, Torque);

%% chemPots
function [chemPots] = setupPots(voltages)
    chemPots = struct('left', [], 'right', []);
    for j = 1:length(voltages)
        chemPotL = voltages(j)/2;
        chemPotR = -1*voltages(j)/2;
        chemPots(j) = struct('left', chemPotL, 'right', chemPotR);
    end
end

%% plotting functions
function [] = plotAngle (value, Title, angles, Vals)
    figure(value);
    plot(angles, Vals)
    title(Title);
end

function [] = plotBoth(value, Title, angles, Vals)
    plotLin2D (value, Title, angles, Vals)
    plotLin3D (value, angles, Vals)
end

function [varargout] = plotLin2D (value, Title, angles, Vals)
    TransPlot = zeros(1, length(angles)*length(angles));
    angleDiff = zeros(1, length(angles)*length(angles));
    indices = zeros(1, length(angles)*length(angles));
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

function [] = plotLin3D (value, angles, Vals)
figure(value)
    surf(angles, angles, Vals)
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