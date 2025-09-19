function [Result, varargout] = transmission(totalSystem, Eigenvals, leftEVs, rightEVs, gammaL, gammaR, chemPots)
    %calculates the transport through by a molecule as calculated in the
    %PhD-Thesis
    
    %disp('Starting calculation of the current.')
    % determine whether the calculation is to be done using linear response
    linearResponse = false;
    if linearResponse == false
        TransmissionMatrix = Transmission(Eigenvals, leftEVs, rightEVs, gammaL, gammaR, chemPots);
    else
        TransmissionMatrix = TransmissionZeroTemp(totalSystem, gammaL, gammaR);
    end
    %disp('Finished calculation of the current.')
    TransmissionElement = trace(TransmissionMatrix);
    Result = real(TransmissionElement);
    
    varargout{1} = real(TransmissionMatrix);
end

%% transmission
function [Result] = Transmission(Eigenvals, leftEVs, rightEVs, gammaL, gammaR, chemPots)
    %disp('Starting calculation of the current element.')
    Result = 0;
    for i = 1:length(leftEVs)
        % get the normal left and right Eigenvectors
        EigVal = Eigenvals(i,i);
        leftEV = leftEVs(:,i)';
        rightEV = rightEVs(:,i);
        
        ProductLeft = rightEV;
        ProductMidLeft = leftEV * gammaR;
        for j = 1:length(leftEVs)
            % get the daggered Eigenvectors
            EigValDagger = Eigenvals(j,j)';
            leftEVdagger = leftEVs(:,j);
            rightEVdagger = rightEVs(:,j)';
            
            % compute the matrix element for chosen i and j
            ProductMid = ProductMidLeft * leftEVdagger;
            ProductRight = rightEVdagger * gammaL;
            
            Product = ProductLeft * ProductMid * ProductRight;
            
            % compute the additional matrix element
            [chemPotL, chemPotR] = chemPots{:};
            factor = matrixElement(EigVal, EigValDagger, chemPotL) - matrixElement(EigVal, EigValDagger, chemPotR);
            
            Result = Result + Product*factor;
            %disp(['i = ',num2str(i),', j = ',num2str(j),', Result = ',num2str(real(trace(Result)))])
        end
    end
    %disp('Finished calculation of the current element.')
end

function [result] = matrixElement(eig1, eig2, chemPot)
    %pot = chemPot(1);
    factor = 1/(eig1 -eig2);
    element1 = log(chemPot - eig1);
    element2 = log(chemPot - eig2);
    result = factor*(element1 - element2);
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
function [matrix] = cutMatrix(Matrix, lengthLead, lengthSample)
    startIndex = lengthLead;
    endIndex = lengthLead+lengthSample + 1;
    matrix = Matrix(startIndex:endIndex, startIndex:endIndex);
    %disp(['Indices: from ',num2str(startIndex),' to ',num2str(endIndex)])
end