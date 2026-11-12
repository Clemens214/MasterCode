%% Variables

% variables for the sample (the analytical model requires a single site)
sizeSample = 1;
orderSample = 2;
energySample = 0;
hopping = 1;
hoppingsSample = hopping*eye(orderSample);
sampleVals = struct('size', sizeSample, 'order', orderSample, 'energy', energySample, 'hopping', hoppingsSample);

% variables for the leads
sizeLead = 1;
energyLead = energySample;
hoppingLead = hopping;
leadVals = struct('size', sizeLead, 'energy', energyLead, 'hopping', hoppingLead);

% variables for the coupling between the leads and the sample
% the simplified expressions of main.tex require hoppingInter = hoppingLead
hoppingInter = hopping;

% variables for the hopping
angleMax = 0.5;
angleStep = 0.1;
anglesTick = makeList(angleMax, angleStep);
angles = pi*anglesTick;

%variables for the Energies
EnergyMax = 2;
EnergyStep = 0.001;
Energies = makeList(EnergyMax, EnergyStep, full=true);

%% Calculation
% compare the energy spectra for the indicated relative twists
if true
    [TransNum, TorqueNum, AngularNum] = calc(angles, Energies, sampleVals, leadVals, hoppingInter);
    [TransAna, TorqueAna, AngularAna] = calcAnal(angles, Energies, sampleVals, leadVals, hoppingInter);
    checkDiff('Transmission', anglesTick, TransNum, TransAna)
    checkDiff('Torque', anglesTick, TorqueNum, TorqueAna)
    checkDiff('Angular', anglesTick, AngularNum, AngularAna)
    PlotCompare('TransmissionCompare', anglesTick, Energies, TransNum, TransAna, Transmission=true)
    PlotCompare('TorquanceCompare', anglesTick, Energies, TorqueNum, TorqueAna, Torque=true)
    PlotCompare('AngularCompare', anglesTick, Energies, AngularNum, AngularAna, Angular=true)
end

% compare the closed-form expressions with the block-matrix expressions
if false
    [TransBlock, TorqueBlock, AngularBlock] = calcAnal(angles, Energies, sampleVals, leadVals, hoppingInter, simplified=false);
    [TransSimp, TorqueSimp, AngularSimp] = calcAnal(angles, Energies, sampleVals, leadVals, hoppingInter, simplified=true);
    checkDiff('Transmission (block vs. simplified)', anglesTick, TransBlock, TransSimp)
    checkDiff('Torque (block vs. simplified)', anglesTick, TorqueBlock, TorqueSimp)
    checkDiff('Angular (block vs. simplified)', anglesTick, AngularBlock, AngularSimp)
end
disp('Test')

%% Calculating functions
function [Transmission, Torquance, Angular] = calc(angles, Energies, sampleVals, leadVals, hoppingInter)
% calculate the numerical traces using the semi-infinite leads
    Transmission = cell(1, length(angles));
    Torquance = cell(1, length(angles));
    Angular = cell(1, length(angles));
    for i = 1:length(angles)
        if sampleVals.order == 1
            hoppingsInter = hoppingInter*[1; 1];
            hoppingsDeriv = [0; 0];
        elseif sampleVals.order == 2
            hoppingsInter = hoppingInter*[cos(angles(i)), sin(angles(i)); 1, 0];
            % the torque is derived with respect to the relative twist,
            % such that only the left coupling contributes
            hoppingsDeriv = hoppingInter*[-1*sin(angles(i)), cos(angles(i)); 0, 0];
        end
        % compute the Hamiltonian of the Sample
        sample = makeSample(sampleVals.energy, sampleVals.hopping, sampleVals.size, sampleVals.order);
        % calculating the trace values
        Transmission{i} = TransCalc(sample, Energies, sampleVals, leadVals, hoppingsInter);
        Torquance{i} = TorqueCalc(sample, Energies, sampleVals, leadVals, hoppingsInter, hoppingsDeriv, conservative=true);
        Angular{i} = AngularCalc(sample, Energies, sampleVals, leadVals, hoppingsInter, nonconservative=true);
        disp(['Numerical; Angle: ', num2str(angles(i)), ', i=', num2str(i)])
    end
end

function [Transmission, Torquance, Angular] = calcAnal(angles, Energies, sampleVals, leadVals, hoppingInter, options)
% calculate the analytical traces of the four-state model
arguments
    angles
    Energies
    sampleVals
    leadVals
    hoppingInter
    options.simplified = true
    options.eta = 0
end
    if sampleVals.size ~= 1 || sampleVals.order ~= 2
        error('The analytical model requires sizeSample = 1 and orderSample = 2.')
    end
    if options.simplified == true && hoppingInter ~= leadVals.hopping
        disp('The simplified expressions require hoppingInter = hoppingLead, using the block expressions instead.')
        options.simplified = false;
    end
    Transmission = cell(1, length(angles));
    Torquance = cell(1, length(angles));
    Angular = cell(1, length(angles));
    for i = 1:length(angles)
        Transmission{i} = TransAnal(angles(i), Energies, sampleVals, leadVals, hoppingInter, simplified=options.simplified, eta=options.eta);
        Torquance{i} = TorqueAnal(angles(i), Energies, sampleVals, leadVals, hoppingInter, simplified=options.simplified, eta=options.eta);
        Angular{i} = AngularAnal(angles(i), Energies, sampleVals, leadVals, hoppingInter, simplified=options.simplified, eta=options.eta);
        disp(['Analytical; Angle: ', num2str(angles(i)), ', i=', num2str(i)])
    end
end

%% analytical functions
function [Results] = TransAnal(angle, Energies, sampleVals, leadVals, hoppingInter, options)
% transmission T(w) = Gamma_L Gamma_R |G_41|^2, eq. (Transmission_analytical)
arguments
    angle
    Energies
    sampleVals
    leadVals
    hoppingInter
    options.simplified = true
    options.eta = 0
end
    Results = zeros(1, length(Energies));
    for i = 1:length(Energies)
        if options.simplified == true
            % closed form for hoppingInter = hoppingLead = t, eq. (Transmission_simplified)
            [GreensSurf, gamma, Denom, t] = SimpleParts(angle, Energies(i), sampleVals, leadVals, eta=options.eta);
            Results(i) = 4*t^8 * gamma^2 * abs(GreensSurf)^4 * Energies(i)^2 * cos(angle)^2 / abs(Denom)^2;
        else
            % block-matrix form, eq. (Greens_Analytical)
            Blocks = GreensBlocks(angle, Energies(i), sampleVals, leadVals, hoppingInter, eta=options.eta);
            Results(i) = Blocks.Gamma^2 * abs(Blocks.G14)^2;
        end
    end
end

function [Results] = TorqueAnal(angle, Energies, sampleVals, leadVals, hoppingInter, options)
% symmetric torque trace <tau>_s = <tau>_L + <tau>_R with tau = -dH/d(Delta theta)
arguments
    angle
    Energies
    sampleVals
    leadVals
    hoppingInter
    options.simplified = true
    options.eta = 0
end
    Results = zeros(1, length(Energies));
    for i = 1:length(Energies)
        if options.simplified == true
            % closed form for hoppingInter = hoppingLead = t, <tau>_L = <tau>_R
            [GreensSurf, gamma, Denom, t] = SimpleParts(angle, Energies(i), sampleVals, leadVals, eta=options.eta);
            TorqueL = -2*t^2*gamma/(2*pi) * t^4 * Energies(i) * sin(2*angle) * real(GreensSurf^2) / abs(Denom)^2;
            Results(i) = 2*TorqueL;
        else
            % block-matrix form: <tau>_l = -Gamma/pi * Re( conj(G_1n) * dt_L' * G_mn )
            Blocks = GreensBlocks(angle, Energies(i), sampleVals, leadVals, hoppingInter, eta=options.eta);
            hoppingDeriv = hoppingInter*[-1*sin(angle); cos(angle)];
            TorqueL = -1*Blocks.Gamma/pi * real( conj(Blocks.G11) * (hoppingDeriv.' * Blocks.Gm1) );
            TorqueR = -1*Blocks.Gamma/pi * real( conj(Blocks.G14) * (hoppingDeriv.' * Blocks.Gm4) );
            Results(i) = TorqueL + TorqueR;
        end
    end
end

function [Results] = AngularAnal(angle, Energies, sampleVals, leadVals, hoppingInter, options)
% antisymmetric angular momentum trace <L_z>_a = <L_z>_L - <L_z>_R
arguments
    angle
    Energies
    sampleVals
    leadVals
    hoppingInter
    options.simplified = true
    options.eta = 0
end
    Pauli = [0, -1j; 1j, 0];
    Results = zeros(1, length(Energies));
    for i = 1:length(Energies)
        if options.simplified == true
            % closed form for hoppingInter = hoppingLead = t, <L_z>_R = -<L_z>_L
            [GreensSurf, gamma, Denom, t] = SimpleParts(angle, Energies(i), sampleVals, leadVals, eta=options.eta);
            AngularL = 4*t^6 * gamma^2 / (2*pi) * real(GreensSurf) * sin(2*angle) / abs(Denom)^2;
            Results(i) = 2*AngularL;
        else
            % block-matrix form: <L_z>_l = Gamma/(2 pi) * G_mn' * sigma_2 * G_mn
            Blocks = GreensBlocks(angle, Energies(i), sampleVals, leadVals, hoppingInter, eta=options.eta);
            AngularL = Blocks.Gamma/(2*pi) * real(Blocks.Gm1' * Pauli * Blocks.Gm1);
            AngularR = Blocks.Gamma/(2*pi) * real(Blocks.Gm4' * Pauli * Blocks.Gm4);
            Results(i) = AngularL - AngularR;
        end
    end
end

%% analytical helping functions
function [Blocks] = GreensBlocks(angle, Energy, sampleVals, leadVals, hoppingInter, options)
% blocks of the embedded Green's function with theta_L = angle, theta_R = 0
arguments
    angle
    Energy
    sampleVals
    leadVals
    hoppingInter
    options.eta = 0
end
    eta = 1j*options.eta;
    t = leadVals.hopping;
    GreensSurf = CalcGreens(Energy, t, leadVals.energy, eta=options.eta);
    % coupling vectors, eq. (Coupling_vector_analytical)
    hoppingL = hoppingInter*[cos(angle); sin(angle)];
    hoppingR = hoppingInter*[1; 0];
    % inverse of the modified sub-matrices, eqs. (subH_left) and (subH_Right)
    subInvL = GreensSurf;
    subInvR = GreensSurf;
    % Schur complement and its inverse, eqs. (Schur_complement) and (Schur_inverse)
    Schur = (Energy + eta - sampleVals.energy)*eye(2) - subInvL*(hoppingL*hoppingL.') - subInvR*(hoppingR*hoppingR.');
    SchurInv = [Schur(2,2), -1*Schur(1,2); -1*Schur(2,1), Schur(1,1)] / (Schur(1,1)*Schur(2,2) - Schur(1,2)*Schur(2,1));
    % blocks of the embedded Green's function, eq. (Greens_Analytical)
    Blocks.G11 = subInvL + subInvL * (hoppingL.' * SchurInv * hoppingL) * subInvL;
    Blocks.G14 = subInvL * (hoppingL.' * SchurInv * hoppingR) * subInvR;
    Blocks.Gm1 = SchurInv * hoppingL * subInvL;
    Blocks.Gm4 = SchurInv * hoppingR * subInvR;
    % broadening, eq. (Gamma_analytical)
    Blocks.Gamma = -2*t^2 * imag(GreensSurf);
end

function [GreensSurf, gamma, Denom, t] = SimpleParts(angle, Energy, sampleVals, leadVals, options)
% common parts of the closed-form expressions for hoppingInter = hoppingLead = t
arguments
    angle
    Energy
    sampleVals
    leadVals
    options.eta = 0
end
    t = leadVals.hopping;
    GreensSurf = CalcGreens(Energy, t, leadVals.energy, eta=options.eta);
    % gamma(w) = Im(g), eq. (Gamma_analytical)
    gamma = imag(GreensSurf);
    % determinant of the Schur complement, eq. (Schur_inverse)
    Denom = (1/GreensSurf)^2 - (t^2*GreensSurf)^2 * cos(angle)^2;
end

%% plotting functions
function [] = PlotCompare(name, angles, Energies, DataNum, DataAna, options)
% plot the numerical (solid) and analytical (dashed) spectra
arguments
    name
    angles
    Energies
    DataNum
    DataAna
    options.Transmission = false
    options.Torque = false
    options.Angular = false
end
    colororder({'b', 'r', 'g', 'm', 'c', 'y', 'k'});
    Palette = colororder();
    fig = figure(Name=name);
    hold on
    for i = 1:length(angles)
        color = Palette(mod(i-1, size(Palette, 1))+1, :);
        plot(Energies, DataNum{i}, '-', Color=color, linewidth=1);
        plot(Energies, DataAna{i}, '--', Color=color, linewidth=1);
    end
    hold off
    grid off
    % set the labels of the plot
    xlabel('\omega [|t|]');
    if options.Transmission == true
        ylabel('T(\omega)');
    elseif options.Torque == true
        ylabel('\tau_s(\omega)');
    elseif options.Angular == true
        ylabel('L_{z,a}(\omega) [1/t]');
    end
    labels = cell(1, 2*length(angles));
    for i = 1:length(angles)
        labels{2*i-1} = ['\Delta\theta=', num2str(angles(i)), '\pi (num.)'];
        labels{2*i} = ['\Delta\theta=', num2str(angles(i)), '\pi (anal.)'];
    end
    legend(labels, 'Location','northoutside', 'NumColumns', 2);
    % set the limits of the plot
    xlim([min(Energies), max(Energies)])
    yMax = max(cellfun(@(x) max(x, [], 'omitnan'), [DataNum, DataAna]));
    yMin = min(cellfun(@(x) min(x, [], 'omitnan'), [DataNum, DataAna]));
    if options.Transmission == true
        ylim([yMin, 1.1*yMax])
    elseif yMax > yMin
        ylim([1.1*yMin, 1.1*yMax])
    end
    resizeFig(fig)
    % change font size
    fontsize(12,"points")
    fontname("Helvetica")
    % export the plot
    cleanfigure;
    filename = strcat(name, '.tex');
    matlab2tikz(filename)
end

function [] = resizeFig (figure)
arguments
    figure
end
    % get the conversion factor
    set(figure, 'PaperPositionMode','auto')
    Width = 14.53; %cm
    Convert = get(0, 'ScreenPixelsPerInch') / 2.54;
    % set the figure size
    PosOld = get(figure, 'Position');
    PosNew = [PosOld(1)/Convert, PosOld(2)/Convert, Width, Width];
    set(figure, 'Units','centimeters', 'Position', PosNew)
end

%% checking functions
function [] = checkDiff(name, angles, DataNum, DataAna)
% print the maximum difference between the numerical and analytical traces
    for i = 1:length(angles)
        Diff = max(abs(DataNum{i} - DataAna{i}), [], 'omitnan');
        Scale = max(abs(DataNum{i}), [], 'omitnan');
        disp([name, '; Angle: ', num2str(angles(i)), 'pi, max. difference = ', num2str(Diff), ', max. value = ', num2str(Scale)])
    end
end

%% helping functions
function [values] = makeList(maxVal, stepVal, options)
    arguments
        maxVal
        stepVal
        options.full = false
    end
    if options.full == false
        minVal = 0;
    else
        minVal = -1*maxVal;
    end
    numVal = (maxVal-minVal)/stepVal+1;
    values = linspace(minVal, maxVal, numVal);
end

function [Result, varargout] = CalcGreens (w, t, eig, options)
arguments
    w
    t
    eig
    options.eta = 0
end
    eta = options.eta;
    % calculate parts of the Greens function
    z = w - eig + 1j*eta;
    factor = z/(2*t^2);
    root = sqrt(1-( 4*t^2/z^2 ));
    % calculate the final Greens function elements
    G_minus = factor*(1 - root);
    G_plus = factor*(1 + root);
    % return the results and ensure a negative imaginary part
    Result = real(G_minus) - 1j*abs(imag(G_minus));
    varargout{1} = real(G_plus) - 1j*abs(imag(G_plus));
end