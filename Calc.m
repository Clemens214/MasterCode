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
angleMax = 2*pi;
angleStep = pi/64;
angles = makeList(angleMax, angleStep);
angles = pi/4;

%variables for the voltages
voltageMax = 2.5;
voltageStep = 0.01;
voltages = makeList(voltageMax, voltageStep);
voltages = 2;

%variables for the Energies
EnergyMax = 2.1;
EnergyStep = 0.01;
Energies = makeList(EnergyMax, EnergyStep, full=true);

sample = makeSample(energySample, hoppingsSample, sizeSample,  orderSample);

if false
    hoppingsSigma = [1, 0; 1, 0];
    SigmaL = zeros(1, length(Energies));
    SigmaR = zeros(1, length(Energies));
    for i = 1:length(Energies)
        [~, ~, ~, ~, sigmaL, sigmaR] = makeSystemSI(Energies(i), sample, energyLead, hoppingLead, hoppingsSigma);
        SigmaL(i) = sigmaL(1,1);
        SigmaR(i) = sigmaR(end,end);
        %disp(['Energy = ', num2str(Energies(i))])
    end
    figure(5)
    hold on
    plot(Energies, real(SigmaL))
    plot(Energies, imag(SigmaL))
    hold off
end

%checkDecomposition(sample, 0)

%% Calculation
Transmission = cell(1, length(angles));
Torque = cell(1, length(angles));
Angular = cell(1, length(angles));
Helicity = cell(1, length(angles));
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
    sample = makeSample(energySample, hoppingsSample, sizeSample,  orderSample);
    
    %checkMatrix(totalSystem);
    
    % calculating the values
    Transmission{i} = TransCalc(sample, Energies, sampleVals, leadVals, hoppingsInter, linearResponse=true ,EM=true);
    Torque{i} = TorqueCalc(sample, Energies, sampleVals, leadVals, hoppingsInter, hoppingsDeriv, linearResponse=true ,EM=true, conservative=true);
    Angular{i} = AngularCalc(sample, Energies, sampleVals, leadVals, hoppingsInter, linearResponse=true ,EM=true, nonconservative=true);
    Helicity{i} = HelicityCalc(sample, Energies, sampleVals, leadVals, hoppingsInter, linearResponse=true ,EM=true, nonconservative=true);
    disp(['Angle: ', num2str(angles(i)), ', i=', num2str(i)])
end

%% plot
Plot(1, angles, Energies, Transmission, twoD=true, Spectrum=true, Title='Transmission')
Plot(2, angles, Energies, Torque, twoD=true, Spectrum=true, Title='Torque')
Plot(3, angles, Energies, Angular, twoD=true, Spectrum=true, Title='Angular Momentum')
Plot(4, angles, Energies, Helicity, twoD=true, Spectrum=true, Title='Helicity')

%Plot(1, angles, Energies, {Transmission, Torque}, Spectrum=true, Both=true, Title='Both')

%Plot(2, angles, Energies, Transmission, Spectrum=true, twoD=true, Title='Transmission')
%Transmission = {TransmissionEM{1}, TransmissionSI{1}, TransmissionSI{1}-TransmissionEM{1}};
%plotSpectrum2D(1, 'Transmission', angles, Energies, Transmission)

%Plot(3, angles, Energies, Torque, Spectrum=true, twoD=true, Title='Torque')
%Torque = {TorqueEM{1}, TorqueSI{1}, TorqueSI{1}-TorqueEM{1}};
%plotSpectrum2D(2, 'Torque', angles, Energies, Torque)

%Plot(4, angles, voltages, Angular, threeD=true, Title='Angular Momentum')
%Angular = {AngularEM{1}, AngularSI{1}, AngularSI{1}-AngularEM{1}};
%plotSpectrum2D(3, 'Angular Momentum', angles, Energies, Angular)

%Plot(5, angles, voltages, Helicity, threeD=true, Title='Helicity')
%Helicity = {HelicityEM{1}, HelicitySI{1}, HelicitySI{1}-HelicityEM{1}};
%plotSpectrum2D(4, 'Helicity', angles, Energies, Helicity)

%Data = {Angular{1}, Helicity{1}, Torque{1}};
%Data = {Angular{1}, Helicity{1}};
%plotSpectrum2D(1, 'Angular Momentum vs Helicity vs Torquance', angles, Energies, Data)

function [] = plotSpectrum2D (value, Title, angles, voltages, Data)
    TransPlot = cell(1, length(angles));
    for i = 1:length(Data)
        TransPlot{i} = zeros(1, length(voltages));
        for j = 1:length(voltages)
            TransPlot{i}(j) = Data{i}(j);
        end
    end
    % plot the data
    figure(value)
    hold on
    %yyaxis left
    for i = 1:length(Data)%-1
        plot(voltages, TransPlot{i});
    end
    %ylabel('Nonconservative')
    if false
        yyaxis right
        plot(voltages, TransPlot{end});
        ylabel('Conservative')
    end
    hold off
    xlabel('Energy (units of t)');
    title(Title);
    %labels = strcat('Angle = ',cellstr(num2str(angles.')));
    labels = {'Extended Molecule'; 'Semi-infinite leads'; 'Difference'};
    legend(labels)
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