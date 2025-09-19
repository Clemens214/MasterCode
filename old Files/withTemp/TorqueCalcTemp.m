%% Variables

% variables for the sample
sizeSample = 4; %48
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
angleMax = 2*pi; %2*pi;
angleStep = pi/100;
angleNum = angleMax/angleStep+1;
angleDiff = linspace(0, angleMax, angleNum);

%variables for the calculation of the current
TempMax = 0; %.5; %2;
TempStep = 0.1; %0.05;
TempNum = TempMax/TempStep+1;
Temps = linspace(0, TempMax, TempNum);

%% Calculation

%preparing the sample's Hamiltonian
sample = makeSample(eigenenergy, hoppingsSample, sizeSample, orderSample);

TransmissionsOld = cell(1, length(angleDiff));
TransmissionsNew = cell(1, length(angleDiff));
for i = 1:length(angleDiff)
    hoppingsInter = zeros(2, orderSample);
    if orderSample == 2
        hoppingsInter = hopping*[cos(angleDiff(i)), sin(angleDiff(i)); 1, 0];
        %hoppingsInter = hopping*[1, 0; 1, 0];
    elseif orderSample == 1
        hoppingsInter = hopping*[1; 1];
    end
    %disp(hoppingsInter(1,:))
    
    % preparing the Extended Molecule Hamiltonian
    [totalSystem, gammaL, gammaR] = makeSystem(sample, sizeSample, sizeLead, hoppingLead, hoppingsInter, maxVal, decay, offset);
    
    checkHamiltonian(totalSystem)
    checkGamma(gammaL, 'gammaL')
    checkGamma(gammaR, 'gammaR')
    
    %compute the Eigenvectors and the Eigenvalues of the system
    [Eigenvals, leftEVs, rightEVs, Product] = eigenvectors(totalSystem);
    
    TransmissionsOld{i} = zeros(1, length(Temps));
    TransmissionsNew{i} = zeros(1, length(Temps));
    
    value = sizeLead + 1;
    step = orderSample;
    for k = 1:length(Temps)
        disp('Starting calculation of the Transmission.')
        [ResultOld, MatrixFullOld] = transport(totalSystem, gammaL, gammaR, Eigenvals, leftEVs, rightEVs, Temps(k), value ,step);
        disp('Finished old calculation, starting new calculation')
        [ResultNew, MatrixFullNew] = transmission(totalSystem, gammaL, gammaR, Eigenvals, leftEVs, rightEVs, Temps(k));
        %[Result] = torque(totalSystem, totalSysDeriv, gammaL, gammaR, Eigenvals, leftEVs, rightEVs, Temps(k));
        
        MatrixOld = cutMatrix(MatrixFullOld, sizeLead, length(sample));
        MatrixNew = cutMatrix(MatrixFullNew, sizeLead, length(sample));
        
        TransmissionsOld{i}(k) = ResultOld;
        TransmissionsNew{i}(k) = ResultNew;
        disp(['Finished calculation of the Transmission. Angle1 = ',num2str(angleDiff(i)),' i=',num2str(i), ...
                                                        ', Temp = ',num2str(Temps(k)), ...
                                                        ', Transmission = ',num2str(ResultOld),', ',num2str(ResultNew)])
    end
end

%TempPlot(Temps, TransmissionsOld)
%TempPlot(Temps, TransmissionsNew)

figure(1)
hold on
diffPlot(angleDiff, TransmissionsOld)
diffPlot(angleDiff, TransmissionsNew)
hold off

%totalPlot(angleDiff, Temps, Transmissions)

%% plotting functions
function [] = totalPlot(angleDiff, Temps, Transmissions)
    TransPlot = zeros(length(Temps), length(Transmissions));
    for i = 1:length(Transmissions)
        TransPlot(:, i) = Transmissions{i}.';
    end
    surf(angleDiff, Temps, TransPlot)
end

function [] = diffPlot(angleDiff, Transmissions)
    TransPlot = zeros(1, length(angleDiff));
    for i = 1:length(Transmissions)
        TransPlot(i) = Transmissions{i}(1);
    end
    plot(angleDiff, TransPlot)
end

function [] = TempPlot(Temps, Transmissions)
    index1 = 1;
    index2 = 1;
    TransPlot =  Transmissions{index1, index2};
    plot(Temps, TransPlot)
end

%% helping functions
function [] = checkHamiltonian(totalSystem)
    Diff = totalSystem - totalSystem';
    Diag = true;
    offDiag = true;
    for i = 1:length(Diff)
        for j = 1:length(Diff)
            if Diff(i,j) ~= 0 %not equal to zero
                if i==j
                    Diag = false;
                else
                    offDiag = false;
                end
            end
        end
    end

    if Diag == true && offDiag == true
        disp('The Hamiltonian is totally hermitian.')
    elseif Diag == true
        disp('The diagonal elements of the Hamiltonian are hermitian.')
    elseif offDiag == true
        disp('The offdiagonal elements of the Hamiltonian are hermitian.')
    else
        disp('No part of the Hamiltonian is hermitian.')
    end
end

function [] = checkGamma(gamma, name)
    [~, flag] = chol(gamma);
    if flag == 0
        disp([name, ' is symmetric positive definite. Flag = ', num2str(flag)])
    else
        disp([name, ' is not symmetric positive definite. Flag = ', num2str(flag)])
    end
end

function [matrix] = cutMatrix(Matrix, lengthLead, lengthSample)
    startIndex = lengthLead;
    endIndex = lengthLead+lengthSample + 1;
    matrix = Matrix(startIndex:endIndex, startIndex:endIndex);
    %disp(['Indices: from ',num2str(startIndex),' to ',num2str(endIndex)])
end