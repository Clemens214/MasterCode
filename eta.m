%% Variables

% variables for the sample
sizeSample = 5;
orderSample = 2;
energySample = 0;
hopping = 1;
hoppingsSample = hopping*eye(orderSample);
sampleVals = struct('size', sizeSample, 'order', orderSample, 'energy', energySample, 'hopping', hoppingsSample);

% variables for the leads
sizeLead = 1;
energyLead = energySample;
hoppingLead = hopping;
leadVals = struct('size', sizeLead, 'energy', energyLead, 'hopping', hoppingLead);

% variables for the hopping
angleMax = 2;
angleStep = 0.001;
anglesTick = makeList(angleMax, angleStep);
angles = pi*anglesTick;

%variables for the Energies
EnergyStep = 0.001;
EnergyMax = 2-EnergyStep;
Energies = makeList(EnergyMax, EnergyStep, full=true);
%Energies = makeList(EnergyMax, EnergyStep, full=false);

%variables for the voltages
voltageMax = 2*EnergyMax;
voltageStep = 0.01;
voltages = makeList(voltageMax, voltageStep);

etaMax = 1E-6;
etaMin = 0;
etas = [1E-16, 1E-14, 1E-12, 1E-10, 1E-8, 1E-6];

%% Calculation
Transmissions = cell(1, length(etas));
Torques = cell(1, length(etas));
Angulars = cell(1, length(etas));
for i = 1:length(etas)
    [Transmissions{i}, Torques{i}, Angulars{i}] = calc(angles, Energies, sampleVals, leadVals, etas(i));
end

%% Difference calculation
TransDiffs = zeros(1, length(Transmissions)); 
TorqueDiffs = zeros(1, length(Torques)); 
AngularDiffs = zeros(1, length(Angulars)); 
for i = 1:length(etas)
    TransDiffs(i) = calcDiff(Transmissions{i}, Transmissions{1}, absolute=true);
    TorqueDiffs(i) = calcDiff(Torques{i}, Torques{1}, absolute=true);
    AngularDiffs(i) = calcDiff(Angulars{i}, Angulars{1}, absolute=true);
end
Data = {TransDiffs, TorqueDiffs, AngularDiffs};
if true
    save('variables.mat')
end

%% plot
Plot('Difference', etas, Energies, Data, twoD=true, eta=true, Difference=true)

%% Calculating functions
function [Transmission, Torquance, Angular] = calc(angles, Energies, sampleVals, leadVals, eta)
    Transmission = cell(1, length(angles));
    Torquance = cell(1, length(angles));
    Angular = cell(1, length(angles));
    for i = 1:length(angles)
        if sampleVals.order == 1
            hoppingsInter = [1; 1];
            hoppingsDeriv = [0; 0];
        elseif sampleVals.order == 2
            hoppingsInter = [cos(angles(i)), sin(angles(i)); 1, 0];
                            %cos(angles(j)), sin(angles(j))];
            hoppingsDeriv = [-1*sin(angles(i)), cos(angles(i)); 0, 0];
                            %-1*sin(angles(j)), cos(angles(j))];
        end
        % compute the Hamiltonian of the Sample
        sample = makeSample(sampleVals.energy, sampleVals.hopping, sampleVals.size, sampleVals.order);
        % calculating the trace values
        Transmission{i} = TransCalc(sample, Energies, sampleVals, leadVals, hoppingsInter, eta=eta);
        Torquance{i} = TorqueCalc(sample, Energies, sampleVals, leadVals, hoppingsInter, hoppingsDeriv, eta=eta);
        Angular{i} = AngularCalc(sample, Energies, sampleVals, leadVals, hoppingsInter, eta=eta);
        disp(['Angle: ', num2str(angles(i)), ', i=', num2str(i)])
    end
end

function [Result, varargout] = calcDiff(Data, Comparison, options)
arguments
    Data 
    Comparison 
    options.absolute = true
    options.relative = false
end
    % lists of maxima
    maxAbs = zeros(1, length(Data));
    maxRel = zeros(1, length(Data));
    % lists of differences
    absDiffs = cell(1, length(Data));
    relDiffs = cell(1, length(Data));
    for i = 1:length(Data)
        absDiffs{i} = zeros(1, length(Data{i}));
        relDiffs{i} = zeros(1, length(Data{i}));
        for j = 1:length(Data{i})
            Diff = Data{i}(j) - Comparison{i}(j);
            %calculate the differences
            absDiffs{i}(j) = abs(Diff);
            relDiffs{i}(j) = abs(Diff)/abs(Comparison{i}(j));
        end
        maxAbs(i) = max(absDiffs{i});
        maxRel(i) = max(relDiffs{i});
    end
    if options.absolute == true && options.relative == false
        Result = max(maxAbs);
        varargout{1} = absDiffs;
    elseif options.relative == true
        Result = max(maxRel);
        varargout{1} = relDiffs;
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
    numVal = (maxVal-minVal)/stepVal+1;
    values = linspace(minVal, maxVal, numVal);
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