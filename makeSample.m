function [sample] = makeSample(eigenenergy, hoppings, size,  order)
% generate the sample's Hamiltonian
arguments
    eigenenergy
    hoppings
    size
    order = 2
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
            rangeRow = order*(row-1)+1 : order*row;
            rangeColumn = order*(column-1)+1 : order*column;
            if smallEnergy(row, column) ~= 0
                sample(rangeRow, rangeColumn) = eye(order, order)*eigenenergy;
            end
            if smallHopping(row, column) ~= 0
                sample(rangeRow, rangeColumn) = hoppings;
            end
        end
    end
end