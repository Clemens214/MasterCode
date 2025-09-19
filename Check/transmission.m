function [Result, varargout] = transmission(totalSystem, Eigenvals, leftEVs, rightEVs, gammaL, gammaR, chemPots, options)
    %calculates the transport through by a molecule as calculated in the PhD-Thesis
    arguments
        totalSystem
        Eigenvals
        leftEVs
        rightEVs
        gammaL
        gammaR
        chemPots
        options.linearResponse = false
        options.evalIndex = 1
    end
    
    %disp('Starting calculation of the current.')
    % determine whether the calculation is to be done using linear response
    if options.linearResponse == false
        TransmissionMatrix = Transmission(Eigenvals, leftEVs, rightEVs, gammaL, gammaR, chemPots);
    else
        TransmissionMatrix = TransmissionZeroTemp(totalSystem, gammaL, gammaR);
    end
    %disp('Finished calculation of the current.')
    Traces = allTrace(real(TransmissionMatrix));
    Result = Traces(options.evalIndex);
    
    varargout{1} = Traces;
    varargout{2} = real(TransmissionMatrix);
end

%% transmission
function [Result] = Transmission(Eigenvals, leftEVs, rightEVs, gammaL, gammaR, chemPots)
    index = struct('i', [], 'j', []);
    for i = 1:length(leftEVs)
        for j = 1:length(leftEVs)
            idx = (i-1)*length(leftEVs) + j;
            index(idx) = struct('i', i, 'j', j);
        end
    end
    
    [chemPotL, chemPotR] = chemPots{:};
    %disp('Starting calculation of the current element.')
    Result = zeros(length(leftEVs), length(leftEVs));
    parfor idx = 1:numel(index)
        i = index(idx).i;
        j = index(idx).j;

        % get the normal Eigenvectors
        EigVal = Eigenvals(i,i);
        leftEV = leftEVs(:,i)';
        rightEV = rightEVs(:,i);
        %[EigVal, leftEV, rightEV] = makeReal(EigVal, leftEV, rightEV);
        
        % get the daggered Eigenvectors
        EigValDagger = Eigenvals(j,j)';
        leftEVdagger = leftEVs(:,j);
        rightEVdagger = rightEVs(:,j)';
        %[EigValDagger, leftEVdagger, rightEVdagger] = makeReal(EigValDagger, leftEVdagger, rightEVdagger);

        % compute the matrix element for chosen i and j
        ProductLeft = rightEV;
        ProductMid = leftEV * gammaR * leftEVdagger;
        ProductRight = rightEVdagger * gammaL;

        Product = ProductLeft * ProductMid * ProductRight;
        
        % compute the additional matrix element
        factor = matrixElement(EigVal, EigValDagger, chemPotL) - matrixElement(EigVal, EigValDagger, chemPotR);
        
        Result = Result + Product*factor;
        %disp(['i = ',num2str(i),', j = ',num2str(j),', Result = ',num2str(real(trace(Result)))])
    end
    %disp('Finished calculation of the current element.')
end

function [result] = matrixElement(eig1, eig2, chemPot)
    if eig1 ~= eig2
        factor = 1/(eig1 - eig2);
        element1 = log(chemPot - eig1);
        element2 = log(chemPot - eig2);
        result = factor*(element1 - element2);
    else
        result = -1/(chemPot - eig1);
    end
end

%% conductance
function [Transmission] = TransmissionZeroTemp(totalSystem, gammaL, gammaR)
    % calculate the Greens Function
    GreensFuncInv = eye(length(totalSystem)) - totalSystem;
    GreensFunc = inv(GreensFuncInv);
    
    % calculate the matrix product
    Transmission = GreensFunc * gammaR * GreensFunc' * gammaL;
end

%% helping functions
function [Transmissions] = allTrace(Matrix)
    Transmissions = zeros(1, floor(length(Matrix)/2));
    for i = 1:floor(length(Matrix)/2)
        startIndex = i;
        endIndex = length(Matrix) - (i-1);
        matrix = Matrix(startIndex:endIndex, startIndex:endIndex);
        Transmissions(i) = trace(matrix);
    end
end

function [Transmissions] = allTraceOld(Matrix, sizeLead, sizeSample)
    Transmissions = zeros(1, sizeLead+1);
    for i = 0:floor(length(Matrix)/2)
        startIndex = (sizeLead+1) - i;
        endIndex = sizeLead+sizeSample + i;
        matrix = Matrix(startIndex:endIndex, startIndex:endIndex);
        Transmissions(i+1) = trace(matrix);
    end
end