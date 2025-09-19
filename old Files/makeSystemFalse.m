function [totalSystem, gammaL, gammaR] = makeSystem(sample, sizeSystem, sizeLead, hoppingLead, hoppingsInter, maxVal, decay, offset)
    sizeSample = length(sample);
    sizeTotal = sizeSample+2*sizeLead;
    sizeSites = sizeSystem+2*sizeLead;
    
    % generate the Hamiltonians of the leads
    %[maxVal, decay, offset] = leadVals;
    left = false;
    right = true;

    sigmaL = makeLead(sizeTotal, sizeSites, sizeLead, maxVal, decay, offset, hoppingLead, left);
    sigmaR = makeLead(sizeTotal, sizeSites, sizeLead, maxVal, decay, offset, hoppingLead, right);

    leadLeft = sigmaL;
    leadRight = sigmaR;
    
    % compute the coupling strengths
    gammaL = 1j*(sigmaL - sigmaL');
    gammaR = 1j*(sigmaR - sigmaR');
    
    % generate the hopping matrices between the leads and the system
    interLeft = makeInter(sizeLead, sizeSample, hoppingsInter, left);
    interRight = makeInter(sizeLead, sizeSample, hoppingsInter, right);

    % generate the Hamiltonian of the total system
    totalSystem = combine(sizeSample, sizeLead, sample, leadLeft, leadRight, interLeft, interRight);
end

%% Functions
function [lead] = makeLead(sizeTotal, sizeSites, sizeLead, maxVal, decay, offset, hopping, right)
    orderSystem = sizeTotal/sizeSites;
    lead = zeros(sizeTotal, sizeTotal);

    for row = 1:sizeLead
        for column = 1:sizeLead
            rowSample = row;
            if right == true %'right'
                if row == column
                    lead(row, column) = 1j*maxVal/(1+exp(decay*(sizeTotal-rowSample-offset)));
                end
            else
                if row == column
                    lead(row, column) = 1j*maxVal/(1+exp(decay*(rowSample-offset)));
                elseif column == row+1 || column == row-1 %or
                    lead(row, column) = hopping;
                end
            end
        end
    end

    for row = sizeLead+1:sizeTotal-sizeLead
        for column = sizeLead+1:sizeTotal-sizeLead
            rowSample = floor(row/orderSystem)+1;
            if right == true %'right'
                if row == column
                    lead(row, column) = 1j*maxVal/(1+exp(decay*(sizeTotal-rowSample-offset)));
                end
            else
                if row == column
                    lead(row, column) = 1j*maxVal/(1+exp(decay*(rowSample-offset)));
                end
            end
        end
    end

    for row = sizeTotal-sizeLead+1:sizeTotal
        for column = sizeTotal-sizeLead+1:sizeTotal
            rowSample = row - (sizeTotal-sizeSites);
            if right == true %'right'
                if row == column
                    lead(row, column) = 1j*maxVal/(1+exp(decay*(sizeTotal-rowSample-offset)));
                elseif column == row+1 || column == row-1 %or
                    lead(row, column) = hopping;
                end
            else
                if row == column
                    lead(row, column) = 1j*maxVal/(1+exp(decay*(rowSample-offset)));
                end
            end
        end
    end
end

function [inter] = makeInter(sizeLead, sizeSystem, hoppingsInter, right)
    endVal = size(hoppingsInter);
    endVal = endVal(2);
    if right == true %'right'
        inter = zeros(sizeLead, sizeSystem);
        for i = 1:endVal
            inter(1, end-i+1) = hoppingsInter(2, i); %False formula
            %inter(1, end-endVal+i) = hoppingsInter(2, i); %Correct Formula
        end
    else
        inter = zeros(sizeSystem, sizeLead);
        for i = 1:endVal
            inter(i, end) = hoppingsInter(1, i);
        end
    end
end

function [totalSystem] = combine(sizeSample, sizeLead, sample, leadLeft, leadRight, interLeft, interRight)
    sizeTotal = length(leadLeft);
    totalSystem = leadLeft + leadRight;

    for row  = 1:sizeSample
        for column = 1:sizeSample
            rowTotal = row+sizeLead;
            columnTotal = column+sizeLead;
            if row == column
                totalSystem(rowTotal, columnTotal) = totalSystem(rowTotal, columnTotal) + sample(row, column);
            else
                totalSystem(rowTotal, columnTotal) = sample(row, column);
            end
        end
    end
    
    top = 1:sizeLead;
    mid = sizeLead+1:sizeTotal-sizeLead;
    bottom = sizeTotal-sizeLead+1:sizeTotal;

    totalSystem(mid, top) = interLeft;
    totalSystem(top, mid) = interLeft.'; %transpose
    
    totalSystem(bottom, mid) = interRight;
    totalSystem(mid, bottom) = interRight.'; %transpose
end