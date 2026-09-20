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
angleStep = 0.001;%1/8;
anglesTick = makeList(angleMax, angleStep);
angles = pi*anglesTick;

%variables for the Energies
EnergyMax = 1;
EnergyStep = 0.25;
%Energies = makeList(EnergyMax, EnergyStep, full=true);
Energies = makeList(EnergyMax, EnergyStep, full=false);

%variables for the voltages
voltageMax = 2*EnergyMax;
voltageStep = 0.01;
voltages = makeList(voltageMax, voltageStep);

%% Calculation
Transmission = cell(1, length(angles));
Torquance = cell(1, length(angles));
Angular = cell(1, length(angles));

for i = 1:length(angles)
    Sizes = [1, 2, 5, 10];
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
    sample = makeSample(energySample, hoppingsSample, sizeSample,  orderSample);
    
    % calculate the surface Green's function
    GreensL = zeros(1, length(Energies));
    GreensR = zeros(1, length(Energies));
    for j = 1:length(Energies)
        [~, ~, ~, ~, ~, ~, GreensL(j), GreensR(j)] = makeSystemSI (Energies(j), sample, 0, hoppingLead, hoppingsInter, hoppingsDeriv);
    end

    % calculating the trace values
    Transmission{i} = TransCalc(sample, Energies, sampleVals, leadVals, hoppingsInter);
    Torquance{i} = TorqueCalc(sample, Energies, sampleVals, leadVals, hoppingsInter, hoppingsDeriv);
    Angular{i} = AngularCalc(sample, Energies, sampleVals, leadVals, hoppingsInter);

    disp(['Angle: ', num2str(angles(i)), ', i=', num2str(i)])
end

%% plot
Plot('Transmission', anglesTick, Energies, Transmission, twoD=true, Value=true, Transmission=true)
Plot('Torquance', anglesTick, Energies, Torquance, twoD=true, Value=true, Torque=true)
Plot('Angular', anglesTick, Energies, Angular, twoD=true, Value=true, Angular=true)
disp('Test')

%% Greens
if false
    if orderSample == 1
        hoppingsInter = [hopping; hopping];
        hoppingsDeriv = [0; 0];
    elseif orderSample == 2
        hoppingsInter = [1, 0; 1, 0];
        hoppingsDeriv = [0, 1; 0, 0];
    end
    % compute the Hamiltonian of the Sample
    sample = makeSample(energySample, hoppingsSample, sizeSample,  orderSample);
    % calculate the surface Green's function
    GreensL = zeros(1, length(Energies));
    GreensR = zeros(1, length(Energies));
    for j = 1:length(Energies)
        [~, ~, ~, ~, ~, ~, GreensL(j), GreensR(j)] = makeSystemSI (Energies(j), sample, 0, hoppingLead, hoppingsInter, hoppingsDeriv);
    end
    % plot the surface Green's function
    plotGreens('GreensL', Energies, GreensL)
    plotGreens('GreensR', Energies, GreensR)
end

function [] = plotGreens(name, Energies, Greens)
    fig = figure(Name=name);
    hold on
    plot(Energies, real(Greens), linewidth=1);
    plot(Energies, imag(Greens), linewidth=1);
    hold off
end

%% chemPots
function [totalSysDeriv] = makeDeriv(sizeSample, orderSample, sizeLead, hoppingsDeriv)
    sampleDeriv = zeros(sizeSample*orderSample, sizeSample*orderSample);
    hoppingDeriv = 0;
    [totalSysDeriv, ~, ~] = makeSystemEM(sampleDeriv, sizeSample, orderSample, sizeLead, hoppingDeriv, hoppingsDeriv, maxVal=0, check=false);
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