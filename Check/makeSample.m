function [sample] = makeSample(eigenenergy, hoppingsSample, sizeSample, orderSample)
    maxSize = sizeSample*orderSample;

    sample = zeros(maxSize, maxSize);
    for row = 1:maxSize
        for column = 1:maxSize
            if row == column
                sample(row, column) = eigenenergy;
            elseif abs(row-column) <= orderSample
                k = mod(row, orderSample)+1;
                l = mod(column, orderSample)+1;
                sample(row, column) = hoppingsSample(k, l);
            end
        end
    end
end