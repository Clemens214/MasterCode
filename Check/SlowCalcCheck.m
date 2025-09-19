%% Variables

% variables for the sample
sizeSample = 48; %48
orderSample = 2; %2;
eigenenergy = 0;
hopping = 1;
hoppingsSample = hopping*eye(orderSample);

% variables for the leads
sizeLead = 104; %128
maxVal = 1;
decay = 0.2; %0.3
offset = 32; %32
%offset should be at most half the length of the leads; normal: 32
leadVals = {maxVal, decay, offset};
hoppingLead = hopping;

% variables for the chemical potentials of the Leads
voltMax = 4;
voltStep = 1;
voltNum = voltMax/voltStep+1;
voltages = linspace(0, voltMax, voltNum);

% variables for the coupling between lead and sample
angleMax = 2*pi; %2*pi;
angleStep = pi/8;
angleNum = angleMax/angleStep+1;
angleDiff = linspace(0, angleMax, angleNum);

%% Calculation

%preparing the sample's Hamiltonian
sample = makeSample(eigenenergy, hoppingsSample, sizeSample, orderSample);

TotalTorques = cell(1, length(angleDiff));
LeftTorques = cell(1, length(angleDiff));
RightTorques = cell(1, length(angleDiff));
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
    sampleDeriv = zeros(length(sample), length(sample));
    derivVals = {0, decay, offset};
    [totalSysDeriv, ~, ~] = makeSystem(sampleDeriv, sizeSample, sizeLead, 0, hoppingsDeriv, derivVals, check=false);
    
    %compute the Eigenvectors and the Eigenvalues of the system
    [Eigenvals, leftEVs, rightEVs] = eigenvectors(totalSystem);
    
    TotalTorques{i} = zeros(1, length(voltages));
    LeftTorques{i} = zeros(1, length(voltages));
    RightTorques{i} = zeros(1, length(voltages));
    for j = 1:length(voltages)
        chemPotL = voltages(j)/2;
        chemPotR = -1*voltages(j)/2;
        chemPots = {chemPotL, chemPotR};
        
        evalIndex = floor(sizeLead/2);
        % compute the torque
        disp('Starting calculation of the Torque.')
        [TotalTorque, LeftTorque, RightTorque] = torqueCheck(totalSystem, Eigenvals, leftEVs, rightEVs, totalSysDeriv, gammaL, gammaR, chemPots);
        disp(['Finished calculation of the Torque. Angle1 = ',num2str(angleDiff(i)),' i=',num2str(i), ...
                                                    ', Voltage = ',num2str(voltages(j)),' j=',num2str(j) ...
                                                    ', Torque = ',num2str(TotalTorque)])
        TotalTorques{i}(j) = TotalTorque;
        LeftTorques{i}(j) = LeftTorque;
        RightTorques{i}(j) = RightTorque;
    end
end
save('Test.mat')
save("TotalTorques.mat","TotalTorques")
save("LeftTorques.mat","LeftTorques")
save("RightTorques.mat","RightTorques")

%% plot the data
figure(1)
totalPlot(angleDiff, voltages, TotalTorques)
savefig("TotalTorque.fig")

figure(2)
totalPlot(angleDiff, voltages, LeftTorques)
savefig("LeftTorque.fig")

figure(3)
totalPlot(angleDiff, voltages, RightTorques)
savefig("RightTorque.fig")

%% plotting functions
function [] = voltPlot(voltages, Transmissions)
    index1 = 1;
    index2 = 1;
    TransPlot =  Transmissions{index1, index2};
    plot(voltages, TransPlot)
end

function [] = diffPlot(angleDiff, Transmissions)
    TransPlot = zeros(1, length(angleDiff));
    for i = 1:length(Transmissions)
        TransPlot(i) = Transmissions{i}(1);
    end
    plot(angleDiff, TransPlot)
end

function [] = totalPlot(angleDiff, voltages, Transmissions)
    TransPlot = zeros(length(voltages), length(angleDiff));
    for i = 1:length(Transmissions)
        TransPlot(:, i) = Transmissions{i}.';
    end
    surf(angleDiff, voltages, TransPlot)
end

function [] = totalTracePlot(Traces, angleDiff, voltages)
    woZero = voltages;
    woZero(woZero==0) = []; % removes the zero entry from woZero
    
    plotFactor = 1;
    sizeMin = floor(sqrt(length(woZero)));
    sizeMax = ceil(sqrt(length(woZero)));
    indexPlot = 0;
    for i = 1:length(voltages)
        if voltages(i) ~= 0
            indexPlot = indexPlot + 1;
            subplot(sizeMin, sizeMax, indexPlot);
            % define the data used in the plot
            TransPlot = zeros(length(Traces{i}{1}), length(angleDiff));
            for j = 1:length(angleDiff)
                TransPlot(:, j) = Traces{j}{i}.';
            end
            index = 0:(length(Traces{1}{1}) -1);
            % plot the data
            endIndex = floor(length(index)*plotFactor);
            surf(angleDiff, index(1:endIndex), TransPlot(1:endIndex,:));
        end
    end
end

function [] = tracePlot(Traces, angleDiff, voltages)
    woZero = voltages;
    woZero(woZero==0) = []; % removes the zero entry from woZero
    
    plotFactor = 1;
    sizeMin = floor(sqrt(length(woZero)));
    sizeMax = ceil(sqrt(length(woZero)));
    indexPlot = 0;
    for i = 1:length(voltages)
        if voltages(i) ~= 0
            indexPlot = indexPlot + 1;
            subplot(sizeMin, sizeMax, indexPlot);
            % define the data used in the plot
            TransPlot = zeros(length(Traces{i}{1}), length(angleDiff));
            for j = 1:length(angleDiff)
                TransPlot(:, j) = Traces{j}{i}.';
            end
            index = 0:(length(Traces{1}{1}) -1);
            % plot the data
            endIndex = floor(length(index)*plotFactor);
            surf(angleDiff, index(1:endIndex), TransPlot(1:endIndex,:));
        end
    end
end

%% helping functions
function [matrix] = cutMatrix(Matrix, sizeLead, sizeSample)
    startIndex = sizeLead;
    endIndex = sizeLead+sizeSample + 1;
    matrix = Matrix(startIndex:endIndex, startIndex:endIndex);
    %disp(['Indices: from ',num2str(startIndex),' to ',num2str(endIndex)])
end