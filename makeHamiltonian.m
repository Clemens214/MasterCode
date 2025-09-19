function [H] = makeHamiltonian(eigenenergies, hopping, size, alignment)
    if size == 0
        size = length(eigenenergies);
    end
    sizeLen = length(eigenenergies);
    maxSize = max([sizeLen, size]);
    
    if contains (alignment, 'top')
        startRow = 0;
    elseif contains (alignment, 'bottom')
        startRow = maxSize - size;
    end
    if contains(alignment, 'left')
        startColumn = 0;
    elseif contains(alignment, 'right')
        startColumn = maxSize - size;
    end
    
    H = zeros(maxSize, maxSize);
    for row = 1:length(eigenenergies)
        for column = 1:length(eigenenergies)
            if row == column
                H(row, column) = eigenenergies(row);
            end
        end
    end
    
    for row = 1:size
        for column = 1:size
            if column == row+1 || column == row-1
                H(startRow+row, startColumn+column) = hopping;
            end
        end
    end
end