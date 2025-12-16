function [sample] = makeSample(eigenenergy, hoppings, size,  order)
% generate the sample's Hamiltonian
arguments
    eigenenergy
    hoppings
    size
    order = 1
end
    sizeSample = size*order;

    sample = zeros(size, size);
    for row = 1:sizeSample
        for column = 1:sizeSample
            if row == column
                sample(row, column) = eigenenergy;
            elseif abs(row-column) <= order
                rowMod = mod(row, order)+1;
                columnMod = mod(column, order)+1;
                sample(row, column) = hoppings(rowMod, columnMod);
            end
        end
    end
end

function [sample] = makeSampleNew(eigenenergy, hoppings, size, order)
% generate the sample's Hamiltonian
arguments
    eigenenergy
    hoppings
    size
    order = 1
end
    % include the eigenenergy
    smallEnergy = eye(size, size);
    % include the hopping
    smallHopping = zeros(size, size);
    for row = 1:size
        for column = 1:size
            if column == row+1 || row == column+1
                smallHopping(row, column) = 1;
            end
        end
    end
    % generate the sample's Hamiltonian
    sizeCenter = size*order;
    sample = zeros(sizeCenter, sizeCenter);
    for row = 1:size
        for column = 1:size
            if smallEnergy(row, column) ~= 0
                range = order*(row-1)+1 : order*row;
                sample(range, range) = eye(order, order)*eigenenergy;
            end
            if smallHopping(row, column) ~= 0
                range = order*(row-1)+1 : order*row;
                sample(range, range) = -1*hoppings;
            end
        end
    end
end