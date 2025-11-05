%% Variables

% variables for the sample
sizeSample = 48;%48;
orderSample = 2;
eigenenergy = 0;
hopping = 1;
hoppingsSample = hopping*eye(orderSample);

% variables for the leads
sizeLead = 104;%104;
[leadVals, derivVals] = calcVals(maxVal = 1, decay = 0.2, offset = 32);
hoppingLead = hopping;

% variables for the hopping
angleMax = 2*pi;
angleStep = pi/16;
angles = makeList(angleMax, angleStep);

%variables for the calculation of the current
voltageMax = 5;
voltageStep = 0.5;
voltages = makeList(voltageMax, voltageStep);
chemPots = setupPots(voltages);
Energies = getEnergies(chemPots);

%% Calculation
Transmission = cell(1, length(angles));
Torque = cell(1, length(angles));
for i = 1:length(angles)
    if orderSample == 1
        hoppingsInter = [hopping; hopping];
        hoppingsDeriv = [0; 0];
    elseif orderSample == 2
        hoppingsInter = [cos(angles(i)), sin(angles(i)); 1, 0];
                        % cos(angles(j)), sin(angles(j))];
        hoppingsDeriv = [-1*sin(angles(i)), cos(angles(i)); 0, 0];
                        % -1*sin(angles(j)), cos(angles(j))];
    end
    
    % compute the Hamiltonian of the Sample
    sample = makeSample(eigenenergy, hoppingsSample, sizeSample,  orderSample);
    
    % preparing the Extended Molecule Hamiltonian
    [totalSystem, gammaL, gammaR, sigmaL, sigmaR] = makeSystemEM(sample, sizeSample, orderSample, sizeLead, hoppingLead, hoppingsInter, leadVals);
    totalSysDeriv = makeDeriv(sizeSample, orderSample, sizeLead, hoppingsDeriv, derivVals);

    checkMatrix(totalSystem);
    
    % calculating the values
    %Transmission{i} = TransCalc(totalSystem, gammaL, gammaR, chemPots);
    Torque{i} = TorqueInt(totalSystem, totalSysDeriv, gammaL, gammaR, chemPots);
    disp(['Angle: ', num2str(angles(i)), ', i=', num2str(i)])
end

%% plot
plot3D(1, 'Transmission', angles, voltages, Transmission)
%plot2D(2, 'Transmission', angles, voltages, Transmission)

plot3D(2, 'Torque', angles, voltages, Torque)

%plotAngle (1, 'Transmission', angles, Transmission);
%plotAngle (2, 'Torque', angles, Torque);

%% chemPots
function [chemPots] = setupPots(voltages)
    chemPots = struct('left', [], 'right', []);
    for j = 1:length(voltages)
        chemPotL = voltages(j)/2;
        chemPotR = -1*voltages(j)/2;
        chemPots(j) = struct('left', chemPotL, 'right', chemPotR);
    end
end

function [totalSysDeriv] = makeDeriv(sizeSample, orderSample, sizeLead, hoppingsDeriv, derivVals)
    sampleDeriv = zeros(sizeSample*orderSample, sizeSample*orderSample);
    hoppingDeriv = 0;
    [totalSysDeriv, ~, ~] = makeSystemEM(sampleDeriv, sizeSample, orderSample, sizeLead, hoppingDeriv, hoppingsDeriv, derivVals, check=false);
end

%% plotting functions
function [] = plotAngle (value, Title, angles, Vals)
    plotVals = cell2mat(Vals);
    figure(value);
    plot(angles, plotVals)
    title(Title);
end

function [] = plot2D (value, Title, angles, voltages, Data)
    TransPlot = cell(1, length(voltages));
    for i = 1:length(voltages)
        TransPlot{i} = zeros(1, length(angles));
        for j = 1:length(angles)
            TransPlot{i}(j) = Data{j}(i);
        end
    end
    % plot the data
    figure(value)
    hold on
    for i = 1:length(voltages)
        plot(angles, TransPlot{i});
    end
    hold off
    title(Title);
    labels = strcat('chemPot = ',cellstr(num2str(angles.')));
    legend(labels)
end

function [] = plot3D (value, Title, angles, voltages, Data)
    TransPlot = zeros(length(voltages), length(angles));
    for i = 1:length(Data)
        TransPlot(:, i) = Data{i}.';
    end
    figure(value)
    surf(angles, voltages, TransPlot)
    title(Title);
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

function [] = plotLin3D (value, Title, angles, Vals)
    figure(value)
    surf(angles, angles, Vals)
    title(Title);
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

function [Filtered] = getEnergies(chemPots)
    Energies = zeros(1, length(chemPots)*2);
    for i = 1:length(chemPots)
        Energies(2*i-1) = chemPots(i).left;
        Energies(2*i) = chemPots(i).right;
    end
    Sorted = sort(Energies);
    Filtered = unique(Sorted);
end

function [] = saveVar(var, order)
    filename = append('Indices', int2str(order), '.mat');
    save(filename, "var")
end

function [] = checkMatrix(totalSystem)
    % Example matrix (replace with your A)
    A = totalSystem;
    % 1) Compute right eigenvectors and eigenvalues
    [V, ~] = eig(A);      % A * V = V * D
    % Conditioning of eigenvector matrix
    condV = cond(V);
    % Rank of eigenvector matrix
    rankV = rank(V);
    % Check if matrix is defective
    if rankV < size(A,1)
        disp('Matrix appears defective (not diagonalizable).');
    else
        disp('Matrix is diagonalizable but ill-conditioned.');
    end
    fprintf('Condition number of V: %g\n', condV);
end