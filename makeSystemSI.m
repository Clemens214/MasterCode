function [totalSystem, GammaL, GammaR, varargout] = makeSystemSI (Energy, sample, eigenenergy, hoppingLead, hoppingsInter, hoppingsDeriv, options)
arguments
    Energy
    sample
    eigenenergy
    hoppingLead
    hoppingsInter
    hoppingsDeriv = zeros(size(hoppingsInter))
    options.size = 0
    options.check = false
    options.checkMore = false
end
    sizeSample = length(sample);
    sizeExtra = options.size;

    hoppingLead = -1*hoppingLead;
    chiral_term = 0.55;
    epsilon = 100000;
    theta = pi/4;

    rangeL = 1 : 2;
    rangeM = sizeSample-3 : sizeSample-2;
    rangeR = sizeSample-1 : sizeSample;
    
    % compute the self-energies of the leads
    [SigmaL, GreensL] = makeSigma (sizeSample, sizeExtra, Energy, hoppingLead, hoppingsInter, eigenenergy, left=true);
    [SigmaR, GreensR] = makeSigma (sizeSample, sizeExtra, Energy, hoppingLead, hoppingsInter, eigenenergy, right=true);
    
    eta = 0.1;
    SigmaL(rangeL, rangeL) = [-1j*eta, 0; 0, -1j*eta];
    SigmaR(rangeR, rangeR) = [-1j*eta, 0; 0, -1j*eta];

    % generate the Hamiltonian of the total system
<<<<<<< Updated upstream
    totalHamiltonian = combineSystem(sample, sizeExtra, eigenenergy, hoppingsInter, hoppingLead);
=======
    totalHamiltonian = makeSystem(sample, sizeExtra, eigenenergy, hoppingsInter, hoppingLead);

    % site 1: end group term
    h1 = [0, 0; 0, 1]*epsilon;
    % site N-1: chiral term
    rc = rotation_matrix(theta - pi/6);
    chiral_group = chiral_term *( rc * h1 * rc.' ) /epsilon;
    % site N: end group term
    r = rotation_matrix(theta);
    rot_h1 = ( r * h1 * r.' );

    totalHamiltonian = sample;
    totalHamiltonian(rangeL, rangeL) = h1;
    totalHamiltonian(rangeM, rangeM) = chiral_group;
    totalHamiltonian(rangeR, rangeR) = rot_h1;

>>>>>>> Stashed changes
    totalSystem = totalHamiltonian + SigmaL + SigmaR;
    
    % compute the coupling strengths of the leads
    GammaL = -1j*(SigmaL - SigmaL');
    GammaR = -1j*(SigmaR - SigmaR');

    GammaL = 1j*(SigmaL - SigmaL');
    GammaR =-1j*(SigmaR - SigmaR');
    
    % check the results
    if options.check == true
        checkGamma(GammaL, 'GammaL')
        checkGamma(GammaR, 'GammaR')
        checkHamiltonian(totalSystem)
    end
    
    % generate the derivative of the Hamiltonian
    sampleDeriv = zeros(size(sample));
<<<<<<< Updated upstream
    totalSysDeriv = combineSystem(sampleDeriv, sizeExtra, 0, hoppingsDeriv, 0);
=======
    totalSysDeriv = makeSystem(sampleDeriv, sizeExtra, 0, hoppingsDeriv, 0);
    
    % site N-1: chiral term
    rc = rotation_matrix(theta - pi/6);
    drc = drotation_matrix(theta - pi/6);
    chiral_group = chiral_term *( drc * h1 * rc.' ) /epsilon + chiral_term *( rc * h1 * drc.' ) /epsilon;
    % site N: end group term
    dr = drotation_matrix(theta);
    r = rotation_matrix(theta);
    rot_h1 = ( dr * h1 * r.') + ( r * h1 * dr.' );

    totalSysDeriv = sampleDeriv;
    totalSysDeriv(rangeM, rangeM) = chiral_group;
    totalSysDeriv(rangeR, rangeR) = rot_h1;
>>>>>>> Stashed changes

    % return the results
    varargout{1} = totalSysDeriv;
    varargout{2} = SigmaL;
    varargout{3} = SigmaR;
    varargout{4} = GreensL;
    varargout{5} = GreensR;
end

% rotation matrix
function [matrix] = rotation_matrix(angle)
    c = cos(angle);
    s = sin(angle);
    matrix = [c, -s; s, c];
end

% d(rotation matrix) / d theta
function [matrix] = drotation_matrix(angle)
    dc = -sin(angle);
    ds = cos(angle);
    matrix = [dc, -ds; ds, dc];
end

%% functions used in generating the total Hamiltonian
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
    totalSystem = zeros(sizeTotal, sizeTotal);
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
    lead = zeros(size, size);
    for row = 1 : size
        for column = 1 : size
            if row == column
                lead(row, column) = eigenenergy;
            elseif abs(row-column) <= 1 && options.left == true && options.right == false
                lead(row, column) = hopping;
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
        inter = zeros(sizeSample, sizeExtra);
        for i = 1:endVal
            inter(i, end) = hoppingsInter(1, i);
        end
    elseif options.left == false && options.right == true
        inter = zeros(sizeExtra, sizeSample);
        for i = 1:endVal
            inter(1, sizeSample-endVal+i) = hoppingsInter(2, i);
        end
    end
end

%% functions used to calculate the self- energies
function [SigmaTotal, varargout] = makeSigma (sizeSample, sizeExtra, Energy, hoppingLead, hoppingsInter, eigenenergy, options)
arguments
    sizeSample
    sizeExtra
    Energy
    hoppingLead
    hoppingsInter
    eigenenergy
    options.left = false
    options.right = false
end
    sizeTotal = sizeSample + 2*sizeExtra;
    % return the self-energy matrix
    [Sigma, GreensFunc] = CalcSigma(Energy, sizeExtra, hoppingLead, hoppingsInter, eigenenergy, left=options.left, right=options.right);
    % place the self- energy matrix
    sizeSigma = length(Sigma);
    if options.left == true
        range = 1 : sizeSigma;
    elseif options.right == true
        range = sizeTotal-sizeSigma+1 : sizeTotal;
    end
    SigmaTotal = zeros(sizeTotal);
    SigmaTotal(range, range) = Sigma;
    varargout{1} = GreensFunc;
end

function [Sigma, varargout] = CalcSigma (Energy, sizeExtra, hoppingLead, hoppingsInter, eigenenergy, options)
arguments
    Energy
    sizeExtra
    hoppingLead
    hoppingsInter
    eigenenergy
    options.left = false
    options.right = false
end
    % calculate the self- energy
    if sizeExtra >= 1
        hoppings = hoppingLead;
    elseif sizeExtra == 0 && options.left == true
        hoppings = hoppingsInter(1, :);
    elseif sizeExtra == 0 && options.right == true
        hoppings = hoppingsInter(2, :);
    end
    Sigma = zeros(length(hoppings), length(hoppings));
    GreensFunc = CalcGreens(Energy, hoppingLead, eigenenergy);
    for i = 1:length(hoppings)
        for j = 1:length(hoppings)
            Sigma(i, j) = hoppings(i) * GreensFunc * hoppings(j);
        end
    end
    varargout{1} = GreensFunc;
end

function [Result, varargout] = CalcGreens (w, t, eig, options)
arguments
    w
    t
    eig
    options.eta = 1E-12
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
    Result = real(G_minus) + -1j*imag(G_minus);
    varargout{1} = G_plus;
end

%% checking functions
function [] = checkGamma(gamma, name)
    [~, flag] = chol(gamma);
    if flag == 0
        disp([name, ' is symmetric positive definite. Flag = ', num2str(flag)])
    else
        disp([name, ' is not symmetric positive definite. Flag = ', num2str(flag)])
    end
end

function [] = checkHamiltonian(totalSystem)
    Diff = totalSystem - totalSystem';
    Diag = true;
    offDiag = true;
    for i = 1:length(Diff)
        for j = 1:length(Diff)
            if Diff(i,j) ~= 0 %not equal to zero
                if i==j
                    Diag = false;
                else
                    offDiag = false;
                end
            end
        end
    end
    % return the result
    if Diag == true && offDiag == true
        disp('The Hamiltonian is totally hermitian.')
    elseif Diag == true
        disp('The diagonal elements of the Hamiltonian are hermitian.')
    elseif offDiag == true
        disp('The offdiagonal elements of the Hamiltonian are hermitian.')
    else
        disp('No part of the Hamiltonian is hermitian.')
    end
end