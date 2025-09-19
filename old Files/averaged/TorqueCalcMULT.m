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
voltMax = 4.5;
voltStep = 0.5;
voltNum = voltMax/voltStep+1;
voltages = linspace(0, voltMax, voltNum);

% variables for the coupling between lead and sample
angleMax = 2*pi; %2*pi;
angleStep = pi/10;
angleNum = angleMax/angleStep+1;
angleDiff = linspace(0, angleMax, angleNum);

avgNum = 1;

%% Calculation
AllResults = cell(1, avgNum);
AllTraces = cell(1, avgNum);
for i = 1:avgNum
    disp(['Average: ',num2str(i)])
    [Results, Traces] = calculation(sizeSample, orderSample, eigenenergy, hoppingsSample, sizeLead, leadVals, hoppingLead, hopping, voltages, angleDiff);
    AllResults{i} = Results;
    AllTraces{i} = Traces;
end

%% save the variables
save("variables.mat")
save("AllResults.mat", "AllResults")
save("AllTraces.mat", "AllTraces")

%% calculate the average
avgResults = average(AllResults);
avgTraces = average(AllTraces);

%% plot the data
%voltPlot(voltages, Transmissions)

%diffPlot(angleDiff, Transmissions)

figure(1)
totalPlot(angleDiff, voltages, avgResults)

figure(2)
totalTracePlot(angleDiff, voltages, avgTraces)

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

function [] = totalTracePlot(angleDiff, voltages, Traces)
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
function [Result] = average(Input)
    Result = cell(1, length(Input{1}));
    for i = 1:length(Input{1})
        Result{i} = cell(1, length(Input{1}{i}));
        for j = 1:length(Input{1}{i})
            if length(Input{1}{i}(j)) > 1
                Result{i}{j} = zeros(1, length(Input{1}{i}(j)));
                for k = 1:length(Input{1}{i}(j))
                    elements = zeros(1, length(Input));
                    for idx = 1:length(Input)
                        disp(['i=',num2str(i),' j=',num2str(j),' k=',num2str(idx)])
                        disp(['element = ', num2str(Input{idx}{i}(j))])
                        elements(idx) = Input{idx}{i}{j}(k);
                    end
                    Result{i}{j}(k) = mean(elements);
                end
            else
                elements = zeros(1, length(Input));
                for idx = 1:length(Input)
                    disp(['i=',num2str(i),' j=',num2str(j),' k=',num2str(idx)])
                    disp(['element = ', num2str(Input{idx}{i}(j))])
                    elements(idx) = Input{idx}{i}(j);
                end
                Result{i}{j} = mean(elements);
            end
        end
    end
end

function [matrix] = cutMatrix(Matrix, sizeLead, sizeSample)
    startIndex = sizeLead;
    endIndex = sizeLead+sizeSample + 1;
    matrix = Matrix(startIndex:endIndex, startIndex:endIndex);
    %disp(['Indices: from ',num2str(startIndex),' to ',num2str(endIndex)])
end