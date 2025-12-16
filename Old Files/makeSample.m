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