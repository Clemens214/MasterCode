function [lcd, varargout] = LowestCommonDenominator(nums, options)
arguments
%   nums - Vector or matrix of non-integer numeric values
    nums
%   tol  - (optional) Tolerance for rational approximation (default: 1e-6)
    options.tol = 1e-6
end
    % Ensure nums is a vector
    nums = nums(:)';

    % Convert to rational approximations
    [numerator, denominator] = rat(nums);
    %[numerator, denominator] = rat(nums, tol);

    % Find the lowest common denominator
    lcd = denominator(1);
    for i = 2:length(denominator)
        lcd = lcm(lcd, denominator(i));
    end

    % Compute equivalent integer numerators
    num_ints = numerator .* (lcd ./ denominator);

    % Reshape numerators back to input shape
    varargout{1} = reshape(num_ints, size(nums));
end
