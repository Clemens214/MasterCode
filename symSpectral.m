%% Variables

% variables for the sample
sizeSample = 1;
orderSample = 2;
energySample = 0;
hopping = sym("t", "real");
assumptions(hopping)
hoppingsSample = hopping*eye(orderSample);

sample = makeSample(energySample, hoppingsSample, sizeSample,  orderSample);
disp(sample);

% variables for the leads
sizeLead = 0;
energyLead = energySample;
hoppingLead = hopping;

% variables for the hopping
angle = sym("a", "positive");
angleDiff = sym("d", "positive");

angleL = angle + angleDiff;
angleR = angle;

hoppingsInter = [cos(angleL), sin(angleL); cos(angleR), sin(angleR)];

%hoppingsDeriv = [-1*sin(angleDiff), cos(angleDiff); 0, 0];
hoppingsDeriv = [-1*sin(angle), cos(angle); -1*sin(angle), cos(angle)];

% variables for the Energy
Energy = sym("w", "real");
assumptions(Energy)

eta = sym("eta", "positive");

%% Calculate the total Greens Function of the Lead
GreensLead = CalcGreensElement(Energy, hoppingLead, energyLead, eta=eta);
GreensSimp = simplify(subs(GreensLead, {eta}, {0}));

%% Generate the total Hamiltonian of the System
GreensFunc = sym("G");
[totalSystem, gammaL, gammaR, totalSysDeriv, sigmaL, sigmaR] = makeSystem(Energy, sample, energySample, hoppingLead, hoppingsInter, hoppingsDeriv, Greens=GreensLead);
totalSystem = simplify(totalSystem);
gammaL = simplify(gammaL);
gammaR = simplify(gammaR);
gammaC = simplify(gammaL + gammaR);
gammaNC = simplify(gammaL - gammaR);

TransOp = gammaR;
TorqueOp = totalSysDeriv;
AngularOp = AngularOperator(length(totalSystem), sizeLead, orderSample);
HelicityOp = HelicityOperator(length(totalSystem), sizeLead, orderSample);

%% Simulate the spectral function
Greens1 = sym("A");
Greens2 = sym("C");
GreensCross = sym("B");

GreensSim = sym(zeros(2, 2));
GreensSim(1, 1) = Greens1;
GreensSim(2, 2) = Greens2;
GreensSim(1, 2) = GreensCross;
GreensSim(2, 1) = GreensCross;

Sigma = makeSigma(length(sample), sizeLead, hoppingLead, hoppingsInter, GreensFunc, right=true);
Gamma = -1j*(Sigma - Sigma');

Spectral = GreensSim * Sigma * GreensSim';
SpectralSimp = simplify(Spectral);

%% Calculate the observables
GreensFunc = CalcGreens (Energy, totalSystem);

disp('Starting calculation of spectral functions.')
spectralFuncL = CalcSpectral(GreensFunc, gammaL, gammaR, left=true);
disp('Finished calculation of left spectral function.')
spectralFuncR = CalcSpectral(GreensFunc, gammaL, gammaR, right=true);
disp('Finished calculation of right spectral function.')
spectralFuncC = CalcSpectral(GreensFunc, gammaL, gammaR, conservative=true);
disp('Finished calculation of conservative spectral function.')
spectralFuncNC = CalcSpectral(GreensFunc, gammaL, gammaR, nonconservative=true);
disp('Finished calculation of nonconservative spectral function.')

[Transmission, TransMatrix] = Calc(TransOp, spectralFuncL);
disp('Finished calculation of the Transmission.')
[TorqueC, TorqueMatrixC] = Calc(TorqueOp, spectralFuncC);
disp('Finished calculation of the conservative Torque.')
[TorqueNC, TorqueMatrixNC] = Calc(TorqueOp, spectralFuncNC);
disp('Finished calculation of the nonconservative Torque.')
[Angular, AngularMatrix] = Calc(AngularOp, spectralFuncC);
disp('Finished calculation of the Angular Momentum.')
[Helicity, HelicityMatrix] = Calc(HelicityOp, spectralFuncC);
disp('Finished calculation of the Helicity.')

%% substitute the values
TransSimp = simplify(subs(Transmission, {hopping, eta}, {1, 0}));
disp('Finished substitution of the Transmission.')
TorqueSimp = simplify(subs(Torque, {hopping, eta}, {1, 0}));
disp('Finished substitution of the Torque.')
AngularSimp = simplify(subs(Angular, {hopping, eta}, {1, 0}));
disp('Finished substitution of the Angular Momentum.')
HelicitySimp = simplify(subs(Helicity, {hopping, eta}, {1, 0}));
disp('Finished substitution of the Helicity.')

if true
    TransVal = simplify(subs(TransSimp, {angle}, {pi/4}));
    disp('Finished substitution of the Transmission.')
    TorqueVal = simplify(subs(TorqueSimp, {angle}, {pi/4}));
    disp('Finished substitution of the Torque.')
    AngularVal = simplify(subs(AngularSimp, {angle}, {pi/4}));
    disp('Finished substitution of the Angular Momentum.')
    HelicityVal = simplify(subs(HelicitySimp, {angle}, {pi/4}));
    disp('Finished substitution of the Helicity.')
    if true
        TransVal = simplify(subs(TransVal, {Energy}, {0}));
        disp('Finished substitution of the Transmission.')
        TorqueVal = simplify(subs(TorqueVal, {Energy}, {0}));
        disp('Finished substitution of the Torque.')
        AngularVal = simplify(subs(AngularVal, {Energy}, {0}));
        disp('Finished substitution of the Angular Momentum.')
        HelicityVal = simplify(subs(HelicityVal, {Energy}, {0}));
        disp('Finished substitution of the Helicity.')
    end
end

%% functions used in calculating the observables
function [Result, varargout] = Calc(operator, SpectralFunc)
arguments
    operator
    SpectralFunc
end
    % do the matrix multiplication
    Matrix = operator * SpectralFunc;
    MatrixSimp = simplify(Matrix);
    % calculate the trace
    Value = trace(MatrixSimp);
    ValueSimp = simplify(Value);
    % return the results
    Result = ValueSimp;
    varargout{1} = MatrixSimp;
end

function [Result, varargout] = CalcSpectral (GreensFunc, gammaL, gammaR, choice)
arguments
    GreensFunc
    gammaL
    gammaR
    choice.conservative = false
    choice.nonconservative = false
    choice.left = false
    choice.right = false
end
    % define the midfactor
    if choice.conservative == true
        midFactor = gammaL + gammaR;
    elseif choice.nonconservative == true
        midFactor = gammaL - gammaR;
    elseif choice.left == true
        midFactor = gammaL;
    elseif choice.right == true
        midFactor = gammaR;
    end
    % do the matrix multiplication
    Matrix = GreensFunc * midFactor * GreensFunc';
    MatrixSimp = simplify(Matrix);
    % return the results
    Result = MatrixSimp;
end

function [GreensSimp] = CalcGreens (Energy, totalSystem)
arguments
    Energy
    totalSystem
end
    % calculate the Greens Function
    GreensInv = Energy*eye(length(totalSystem)) - totalSystem;
    GreensFunc = inv(GreensInv);
    GreensSimp = simplify(GreensFunc); 
end

%% functions used in generating the operators
function [operator] = AngularOperator(sizeTotal, sizeLead, order)
arguments
    sizeTotal
    sizeLead
    order = 2
end
    Pauli = [0, -1j; 1j, 0];
    % generate the angular momentum operator
    sizeCenter = sizeTotal - 2*sizeLead;
    sizeSample = sizeCenter/order;
    small = eye(sizeSample, sizeSample);
    sample = zeros(sizeCenter, sizeCenter);
    for row = 1:sizeSample
        for column = 1:sizeSample
            if small(row, column) ~= 0
                range = order*(row-1)+1 : order*row;
                sample(range, range) = small(row, column) * Pauli;
            end
        end
    end
    % embed the operator in total Hamiltonian
    operator = zeros(sizeTotal, sizeTotal);
    mid = sizeLead+1 : sizeTotal-sizeLead;
    operator(mid, mid) = sample;
    operator = sym(operator);
end

function [operator] = HelicityOperator(sizeTotal, sizeLead, order)
arguments
    sizeTotal
    sizeLead
    order = 2
end
    Pauli = [0, -1j; 1j, 0];
    % generate the helicity operator
    sizeCenter = sizeTotal - 2*sizeLead;
    sizeSample = sizeCenter/order;
    small = zeros(sizeSample, sizeSample);
    for row = 1:sizeSample
        for column = 1:sizeSample
            if column == row+1
                small(row, column) = 1;
            elseif row == column+1
                small(row, column) = -1;
            end
        end
    end
    sample = zeros(sizeCenter, sizeCenter);
    for row = 1:sizeSample
        for column = 1:sizeSample
            rangeRow = order*(row-1)+1 : order*row;
            rangeColumn = order*(column-1)+1 : order*column;
            if small(row, column) ~= 0
                sample(rangeRow, rangeColumn) = 1j/2 * small(row, column) * Pauli;
            end
        end
    end
    % embed the operator in total Hamiltonian
    operator = zeros(sizeTotal, sizeTotal);
    mid = sizeLead+1 : sizeTotal-sizeLead;
    operator(mid, mid) = sample;
    operator = sym(operator);
end

%% functions used in generating the sample
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
    Energies = eye(order, order)*eigenenergy;
    % include the hopping
    smallHopping = zeros(size, size);
    for row = 1:size
        for column = 1:size
            if column == row+1 || row == column+1
                smallHopping(row, column) = 1;
            end
        end
    end
    % generate the symbolic Hamiltonian
    sizeCenter = size*order;
    sample = sym(zeros(sizeCenter, sizeCenter));
    for row = 1:size
        for column = 1:size
            rangeRow = order*(row-1)+1 : order*row;
            rangeColumn = order*(column-1)+1 : order*column;
            if smallEnergy(row, column) ~= 0
                sample(rangeRow, rangeColumn) = sym(Energies);
            end
            if smallHopping(row, column) ~= 0
                sample(rangeRow, rangeColumn) = sym(hoppings);
            end
        end
    end
    sample = simplify(sample);
end

%% functions used in generating the total Hamiltonian
function [totalSystem, GammaL, GammaR, varargout] = makeSystem (Energy, sample, eigenenergy, hoppingLead, hoppingsInter, hoppingsDeriv, options)
arguments
    Energy
    sample
    eigenenergy
    hoppingLead
    hoppingsInter
    hoppingsDeriv = zeros(size(hoppingsInter))
    options.Greens = sym("G")
    options.size = 0
end
    sizeSample = length(sample);
    sizeExtra = options.size;
    
    % compute the self-energies of the leads
    GreensFunc = options.Greens;
    [SigmaL, GreensL] = makeSigma (sizeSample, sizeExtra, hoppingLead, hoppingsInter, GreensFunc, left=true);
    [SigmaR, GreensR] = makeSigma (sizeSample, sizeExtra, hoppingLead, hoppingsInter, GreensFunc, right=true);

    % generate the Hamiltonian of the total system
    totalHamiltonian = combineSystem(sample, sizeExtra, eigenenergy, hoppingsInter, hoppingLead);
    totalSystem = totalHamiltonian + SigmaL + SigmaR;
    
    % compute the coupling strengths of the leads
    GammaL = -1j*(SigmaL - SigmaL');
    GammaR = -1j*(SigmaR - SigmaR');
    
    % generate the derivative of the Hamiltonian
    energyDeriv = 0;
    hoppingDeriv = 0;
    sampleDeriv = zeros(size(sample));
    SigmaDerivL = makeSigma (sizeSample, sizeExtra, hoppingDeriv, hoppingsDeriv, GreensFunc, left=true);
    SigmaDerivR = makeSigma (sizeSample, sizeExtra, hoppingDeriv, hoppingsDeriv, GreensFunc, right=true);
    totalHamDeriv = combineSystem(sampleDeriv, sizeExtra, energyDeriv, hoppingsDeriv, hoppingDeriv);
    totalSysDeriv = totalHamDeriv + SigmaDerivL + SigmaDerivR;

    % return the results
    varargout{1} = simplify(totalSysDeriv);
    varargout{2} = simplify(SigmaL);
    varargout{3} = simplify(SigmaR);
    varargout{4} = simplify(GreensL);
    varargout{5} = simplify(GreensR);
end

function [totalSystem] = combineSystem (sample, sizeExtra, eigenenergy, hoppingsInter, hoppingLead)
arguments
    sample
    sizeExtra
    eigenenergy
    hoppingsInter
    hoppingLead
end
    sizeSample = length(sample);
    sizeTotal = sizeSample + 2*sizeExtra;
    totalSystem = sym(zeros(sizeTotal, sizeTotal));
    % generate the sections
    top = 1 : sizeExtra;
    mid = sizeExtra+1 : sizeTotal-sizeExtra;
    bottom = sizeTotal-sizeExtra+1 : sizeTotal;
    % include the sample
    totalSystem(mid, mid) = sample;
    % include the environment
    if sizeExtra > 0
        % include the leads
        totalSystem(top, top) = makeLead(eigenenergy, hoppingLead, sizeExtra, left=true);
        totalSystem(bottom, bottom) = makeLead(eigenenergy, hoppingLead, sizeExtra, right=true);
        % include the hopping matrices
        % include the left side
        interLeft = makeInter(sizeSample, sizeExtra, hoppingsInter, left=true);
        totalSystem(mid, top) = interLeft;
        totalSystem(top, mid) = interLeft.';
        % include the right side
        interRight = makeInter(sizeSample, sizeExtra, hoppingsInter, right=true);
        totalSystem(bottom, mid) = interRight;
        totalSystem(mid, bottom) = interRight.';
    end
end

function [lead] = makeLead(eigenenergy, hopping, size, options)
arguments
    eigenenergy
    hopping
    size
    options.left = false
    options.right = false
end
    lead = sym(zeros(size, size));
    for row = 1 : size
        for column = 1 : size
            if row == column
                lead(row, column) = sym(eigenenergy);
            elseif abs(row-column) <= 1 && options.left == true && options.right == false
                lead(row, column) = sym(hopping);
            end
        end
    end
end

function [inter] = makeInter(sizeSample, sizeExtra, hoppingsInter, options)
% generate the hopping matrix between the lead and the sample
arguments
    sizeSample
    sizeExtra
    hoppingsInter
    options.left = false
    options.right = false
end
    endVal = size(hoppingsInter);
    endVal = endVal(2);
    if options.left == true && options.right == false
        inter = sym(zeros(sizeSample, sizeExtra));
        for i = 1:endVal
            inter(i, end) = sym(hoppingsInter(1, i));
        end
    elseif options.left == false && options.right == true
        inter = sym(zeros(sizeExtra, sizeSample));
        for i = 1:endVal
            inter(1, sizeSample-endVal+i) = sym(hoppingsInter(2, i));
        end
    end
end

function [SigmaTotal, varargout] = makeSigma (sizeSample, sizeExtra, hoppingLead, hoppingsInter, GreensFunc, options)
arguments
    sizeSample
    sizeExtra
    hoppingLead
    hoppingsInter
    GreensFunc = sym("G")
    options.left = false
    options.right = false
end
    sizeTotal = sizeSample + 2*sizeExtra;
    % calculate the self- energy
    if sizeExtra >= 1
        hoppings = hoppingLead;
    elseif sizeExtra == 0 && options.left == true
        hoppings = hoppingsInter(1, :);
    elseif sizeExtra == 0 && options.right == true
        hoppings = hoppingsInter(2, :);
    end
    Sigma = sym(zeros(length(hoppings), length(hoppings)));
    for i = 1:length(hoppings)
        for j = 1:length(hoppings)
            Sigma(i, j) = hoppings(i) * GreensFunc * hoppings(j);
        end
    end
    % place the self- energy matrix
    sizeSigma = length(Sigma);
    if options.left == true
        range = 1 : sizeSigma;
    elseif options.right == true
        range = sizeTotal-sizeSigma+1 : sizeTotal;
    end
    SigmaTotal = sym(zeros(sizeTotal, sizeTotal));
    SigmaTotal(range, range) = Sigma;
    varargout{1} = GreensFunc;
end

function [Result, varargout] = CalcGreensElement (w, t, eig, options)
arguments
    w
    t
    eig
    %options.eta = 1E-12
    options.eta = sym("eta")
end
    eta = options.eta;
    % calculate parts of the Greens function
    val = w - eig + 1j*eta;
    factor = val/(2*t^2);
    root = sqrt(1-( 4*t^2/val^2 ));
    % calculate the final Greens function elements
    G_plus = factor*(1 + root);
    G_minus = factor*(1 - root);
    % return the results
    Result = conj(G_minus);
    varargout{1} = G_plus;
end