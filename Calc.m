%% Variables

% variables for the sample
sizeSample = 48;
orderSample = 2;
energySample = 0;
hopping = 1;
hoppingsSample = hopping*eye(orderSample);
sampleVals = struct('size', sizeSample, 'order', orderSample, 'energy', energySample, 'hopping', hoppingsSample);

% variables for the leads
sizeLead = 104;
energyLead = energySample;
hoppingLead = hopping;
leadVals = struct('size', sizeLead, 'energy', energyLead, 'hopping', hoppingLead);

% variables for the hopping
angleMax = 2;
angleStep = 1/16;
anglesTick = makeList(angleMax, angleStep);
angles = makeList(pi*angleMax, pi*angleStep);

%variables for the voltages
voltageMax = 4.5;
voltageStep = 0.01;
voltages = makeList(voltageMax, voltageStep);

%variables for the Energies
EnergyMax = 2.5;
EnergyStep = 0.01;
Energies = makeList(EnergyMax, EnergyStep, full=true);

%sample = makeSample(energySample, hoppingsSample, sizeSample,  orderSample);
%checkDecomposition(sample, 0)

%% Calculation
Current = cell(1, length(angles));
Transmission = cell(1, length(angles));
Torque = cell(1, length(angles));
Torquance = cell(1, length(angles));
AngularInt = cell(1, length(angles));
Angular = cell(1, length(angles));
HelicityInt = cell(1, length(angles));
Helicity = cell(1, length(angles));
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
    
    %checkMatrix(totalSystem);
    
    % calculating the values
    %[Current{i}, Transmission{i}, EnergiesTrans] = TransCalc(sample, voltages, sampleVals, leadVals, hoppingsInter, linearResponse=false);
    [Torque{i}, Torquance{i}, EnergiesTorque] = TorqueCalc(sample, voltages, sampleVals, leadVals, hoppingsInter, hoppingsDeriv, EM=true, SI=false, linearResponse=false, nonconservative=true);
    %Angular{i} = AngularCalc(sample, Energies, sampleVals, leadVals, hoppingsInter);
    %Helicity{i} = HelicityCalc(sample, Energies, sampleVals, leadVals, hoppingsInter);
    disp(['Angle: ', num2str(angles(i)), ', i=', num2str(i)])
end

if true
    save('variables')
else
    load('variables')
end

%% plot
%Plot(1, anglesTick, [0], Transmission, twoD=true, Value=true, Title='Transmission')
%Plot(2, anglesTick, [1], Torque, twoD=true, Value=true, Title='Torque')

%Plot(1, angles, Energies, {Transmission, Torque}, Spectrum=true, Both=true, Title='Both')

%Plot(1, angles, voltages, Current, threeD=true, Title='Current')
%Plot(2, angles, EnergiesTrans, Transmission, threeD=true, Title='Transmission')
Plot(3, angles, voltages, Torque, threeD=true, Title='Torque')
Plot(4, angles, EnergiesTorque, Torquance, threeD=true, Title='Torquance')
%Plot(5, angles, voltages, AngularInt, threeD=true, Title='Angular Momentum, integrated')
%Plot(6, angles, EnergiesAngular, Angular, threeD=true, Title='Angular Momentum')
%Plot(7, angles, voltages, HelicityInt, threeD=true, Title='Helicity, integrated')
%Plot(8, angles, EnergiesHelicity, Helicity, threeD=true, Title='Helicity')

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