function [totalSystem, gammaL, gammaR] = prepareEM(sample, hopping, lengthSample, lengthTotal, maxVal, decay, offset)
    hoppingLead = hopping;
    hoppingInter = hopping;
    
    % compute the fermi functions
    sites = 1:lengthTotal;
    fermiFuncLeft = fermifuncLeft(sites, maxVal, decay, offset);
    fermiFuncRight = fermifuncRight(sites, maxVal, decay, offset, lengthTotal);
    
    % generate the Hamiltonians of the leads
    lengthLead = (lengthTotal-lengthSample)/2;
    sigmaL = makeHamiltonian(1j*fermiFuncLeft, hoppingLead, lengthLead, 'top left');
    sigmaR = makeHamiltonian(1j*fermiFuncRight, hoppingLead, lengthLead, 'bottom right');
    
    % compute the coupling strengths
    gammaL = 1j*(sigmaL - sigmaL');
    gammaR = 1j*(sigmaR - sigmaR');
    
    % generate the Hamiltonian of the total system
    totalSystem = combineH(sample, sigmaL, sigmaR, hoppingInter, lengthLead);
end

%% Functions

function [values] = fermifuncLeft(sites, maxVal, decay, offset)
    values = zeros(1, length(sites));
    for i = 1:length(sites)
        values(i) = maxVal/(1+exp(decay*(sites(i)+1-offset)));
    end
end

function [values] = fermifuncRight(sites, maxVal, decay, offset, lengthTotal)
    values = zeros(1, length(sites));
    for i = 1:length(sites)
        values(i) = maxVal/(1+exp(decay*(lengthTotal-sites(i)-offset)));
    end
end

%% 
function [arrayNew] = combineH(center, left, right, hoppingInter, size)
    sizeCenter = length(center);
    sizeLeft = length(left);
    sizeRight = length(right);
    sizeMax = max([sizeCenter, sizeLeft, sizeRight]);
    arrayNew = zeros(sizeMax, sizeMax);
    for row = size-1 : sizeMax-size
        for column = size-1 : sizeMax-size
           if row-size > 0 && column-size > 0
               arrayNew(row, column) = center(row-size, column-size);
           end
        end
    end
    values = [size, sizeMax-size];
    for row = 1:sizeMax
       for column = 1:sizeMax
            newFactor = 0;
            if row == values(1) && column == values(1)+1
                newFactor = hoppingInter;
            elseif row == values(1)+1 && column == values(1)
                newFactor = hoppingInter;
            end
            if row == values(2) && column == values(2)+1
                newFactor = hoppingInter;
            end
            if row == values(2)+1 && column == values(2)
                newFactor = hoppingInter;
            end
            arrayNew(row, column) = arrayNew(row, column) + left(row, column) + right(row, column) + newFactor;
        end
    end
end

