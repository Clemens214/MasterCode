function TF = isapprox(A,B,varargin)
% Syntax:
%     TF = isapprox(A,B)
%     TF = isapprox(A,B,toleranceLevel)
%     TF = isapprox(A,B,Name=Value)
%
%     Tolerance Levels:
%         "verytight"
%         "tight"
%         "loose"
%         "veryloose"
%
%     Name-Value Arguments:
%         AbsoluteTolerance
%         RelativeTolerance
%
% For more information, see documentation

%   Copyright 2024 The MathWorks, Inc.

[A,B,abstol,reltol] = parseInputs(A,B,varargin{:});
TF = A == B | 0 <= max(abstol,reltol.*max(abs(A),abs(B))) - abs(A-B);
end

% -------------------------------------------------------------------------

function [A,B,abstol,reltol] = parseInputs(A,B,varargin)

if ~isnumeric(A) || ~isnumeric(B)
    error(message("MATLAB:isapprox:NonNumericInput"));
end

checkForPrecisionLoss(A);
checkForPrecisionLoss(B);

inputIsSingle = isa(A,"single") || isa(B,"single");
A = double(A);
B = double(B);

validToleranceLevels = ["verytight","tight","loose","veryloose"];

if numel(varargin) < 2
    if inputIsSingle
                          % "verytight","tight","loose","veryloose"
        toleranceValues =      [5e-7      1e-6    1e-4     1e-2];
    else
                          % "verytight","tight","loose","veryloose"
        toleranceValues =      [1e-15     1e-12   1e-8     1e-3];
    end
    % Default of tolerance level syntax
    if isempty(varargin)
        levelMask = [true false false false]; % Default is the same as "verytight"
    else
        levelMask = matlab.internal.math.checkInputName(varargin{1},validToleranceLevels);
        if nnz(levelMask) ~= 1
            error(message("MATLAB:isapprox:InvalidToleranceLevel"));
        end
    end

    abstol = toleranceValues(levelMask);
    reltol = abstol;
else
    % Name-value argument syntax
    if rem(numel(varargin),2) == 1
        if any(matlab.internal.math.checkInputName(varargin{1},validToleranceLevels))
            error(message("MATLAB:isapprox:LevelWithOtherOption"));
        else
            error(message("MATLAB:isapprox:KeyWithoutValue"));
        end
    end

    abstol = 0;
    reltol = 0;
    validNames = ["AbsoluteTolerance","RelativeTolerance"];

    for ii = 1:2:numel(varargin)
        nameMask = matlab.internal.math.checkInputName(varargin{ii},validNames);
        if nameMask(1) % AbsoluteTolerance
            abstol = varargin{ii + 1};
            if isempty(abstol) || ~isfloat(abstol) || ~isreal(abstol) || ~allfinite(abstol) ...
                    || any(abstol < 0,"all")
                error(message("MATLAB:isapprox:InvalidAbsTol"));
            end
            abstol = double(abstol);
        elseif nameMask(2) % RelativeTolerance
            reltol = varargin{ii + 1};
            if isempty(reltol) || ~isfloat(reltol) || ~isreal(reltol) || ~allfinite(reltol) ...
                    || any(reltol < 0,"all") || any(reltol >= 1,"all")
                error(message("MATLAB:isapprox:InvalidRelTol"));
            end
            reltol = double(reltol);
        else
            error(message("MATLAB:isapprox:InvalidName"));
        end
    end
end
end

% -------------------------------------------------------------------------

function checkForPrecisionLoss(data)
if (isa(data,"int64") || isa(data,"uint64")) && any(abs(data) > flintmax,"all")
    error(message("MATLAB:isapprox:LossOfPrecision"));
end
end

