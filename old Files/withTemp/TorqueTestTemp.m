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
leadVals = [maxVal, decay, offset];
hoppingLead = hopping;

% variables for the coupling between lead and sample
angleMax = 2*pi;
angleStep = pi/10;
angleNum = angleMax/angleStep+1;
angles = linspace(0, angleMax, angleNum);
%angleDiff = linspace(0, angleMax, angleNum);

%variables for the calculation of the current
TempMax = 0; %.5; %2;
TempStep = 0.05; %0.05;
TempNum = TempMax/TempStep+1;
Temps = linspace(0, TempMax, TempNum);


%% Calculation

%preparing the sample's Hamiltonian
sample = makeSample(eigenenergy, hoppingsSample, sizeSample, orderSample);

angleDiff = zeros(length(angles), length(angles));
Transmissions = cell(length(angles), length(angles)); %= {cell([1, length(Temps)])};
for i = 1:length(angles)
    for j = 1:length(angles)
        indices = [i, j];
        hoppingsInter = zeros(2, orderSample);
        for k = 1:2
            sides = [i, j];
            if orderSample == 2
                hoppingsInter(k, :) = hopping*[cos(angles(sides(k))), sin(angles(sides(k)))];
            elseif orderSample == 1
                hoppingsInter(k) = hopping;
            end
        end
        %disp(hoppingsInter)

        % preparing the Extended Molecule Hamiltonian
        [totalSystem, gammaL, gammaR] = makeSystem(sample, sizeSample, sizeLead, hoppingLead, hoppingsInter, maxVal, decay, offset);
        
        %compute the Eigenvectors and the Eigenvalues of the system
        [Eigenvals, leftEVs, rightEVs, Product] = eigenvectors(totalSystem);
        
        angleDiff(i, j) = angles(i)-angles(j);
        Transmissions{i, j} = zeros(1, length(Temps));

        value = sizeLead + 1;
        step = orderSample;
        for k = 1:length(Temps)
            disp('Starting calculation of the Transmission.')
            [Result] = transport(totalSystem, gammaL, gammaR, Eigenvals, leftEVs, rightEVs, Temps(k), value, step);
            Transmissions{i, j}(k) = Result;
            disp(['Finished calculation of the Transmission. Angle1 = ',num2str(angles(i)),' i=',num2str(i), ...
                                                            ', Angle2 = ',num2str(angles(j)),' j=',num2str(j), ...
                                                            ', Temp = ',num2str(Temps(k)), ...
                                                            ', Transmission = ',num2str(Result)])
        end
    end
end

%anglePlot(Transmissions, angleMax, angleNum)

%TempPlot(Temps, Transmissions)

diffPlot(angleDiff, Transmissions)

%totalPlot(angleDiff, Temps, Transmissions)

%% plotting functions
function [] = totalPlot(angleDiff, Temps, Transmissions)
    [angleSort, TransSort] = plotSort(angleDiff, Transmissions);
    TransPlot = zeros(length(Temps), length(TransSort));
    for i = 1:length(TransSort)
        TransPlot(:, i) = TransSort{i}.';
    end
    surf(angleSort, Temps, TransPlot)
end

function [] = diffPlot(angleDiff, Transmissions)
    [angleSort, TransSort] = plotSort(angleDiff, Transmissions);
    TransPlot = zeros(1, length(angleSort));
    for i = 1:length(TransSort)
        TransPlot(i) = TransSort{i}(1);
    end
    plot(angleSort, TransPlot)
end

function [] = anglePlot(Transmissions, angleMax, angleNum)
    [angleX, angleY] = meshgrid(linspace(0, angleMax, angleNum), linspace(0, angleMax, angleNum));
    TransZ = zeros(length(Transmissions), length(Transmissions));
    for i = 1:length(Transmissions)
        for j = 1:length(Transmissions)
            TransZ(i, j) = Transmissions{i, j}(1);
        end
    end

    surf(angleX, angleY, TransZ)
end

function [] = TempPlot(Temps, Transmissions)
    angle1 = 0;
    angle2 = 0;
    TransPlot =  Transmissions{angle1, angle2};
    plot(Temps, TransPlot)
end

%% helper functions
function [angleSort, TransSort] = plotSort(angleDiff, Transmissions)
    angleSort = zeros(1, length(angleDiff)^2);
    TransSort = cell(1, length(angleDiff)^2);
    for i = 1:length(angleDiff)
        for j = 1:length(angleDiff)
            angleSort(1, (i-1)*length(angleDiff)+j) = angleDiff(i, j);
            TransSort{1, (i-1)*length(angleDiff)+j} = Transmissions{i, j};
        end
    end
    [angleSort, sortIndex] = sort(angleSort, 'ascend');
    TransSort = TransSort(sortIndex);
end