function [G_NEGF, G_BW, varargout] = ResonantCheck(Energies, EnergiesDot, totalSystem, GammaL, GammaR, SigmaL, SigmaR, options)
arguments
    Energies
    EnergiesDot
    totalSystem
    GammaL
    GammaR
    SigmaL
    SigmaR
    options.fermiEnergy = 0
    options.eigEnergy = 0
    options.idxSample = ceil(length(totalSystem)/2)
end
    Energy = options.fermiEnergy;
    EnergyDot = options.eigEnergy;
    idxSample = options.idxSample;
    totalHamiltonian = totalSystem - SigmaL - SigmaR;
    
    % calculate the Extended Molecule Transmission and Conductance
    [T_NEGF] = calcEM(Energies, totalSystem, GammaL, GammaR);
    [G_NEGF] = sweepEM(EnergiesDot, totalSystem, GammaL, GammaR, Energy, idxSample);
    
    % calculate the Breit-Wigner Transmission
    [SigmaL_eff, GammaL_eff] = calcBWside(Energies, totalHamiltonian, SigmaL, EnergyDot, idxSample);
    [SigmaR_eff, GammaR_eff] = calcBWside(Energies, totalHamiltonian, SigmaR, EnergyDot, idxSample);
    
    [T_BW, Delta] = calcBWboth(Energies, SigmaL_eff, SigmaR_eff, GammaL_eff, GammaR_eff, EnergyDot);

    % calculate the Breit-Wigner Conductance
    [SigmaL_tmp, GammaL_tmp] = sweepBWside(EnergiesDot, totalHamiltonian, SigmaL, Energy, idxSample);
    [SigmaR_tmp, GammaR_tmp] = sweepBWside(EnergiesDot, totalHamiltonian, SigmaR, Energy, idxSample);
    
    [G_BW] = sweepBWboth(EnergiesDot, SigmaL_tmp, SigmaR_tmp, GammaL_tmp, GammaR_tmp, Energy);

    % other outputs
    varargout{1} = T_NEGF;
    varargout{2} = T_BW;
    varargout{3} = GammaL_eff;
    varargout{4} = GammaR_eff;
    varargout{5} = Delta;
end

%% -------------------- Other functions: EM -------------------- 
function [T_NEGF] = calcEM(Energies, totalSystem, GammaL, GammaR)
    eta = 1j*1e-12;
    % Energy loop: compute GR (full) & projections
    T_NEGF = zeros(size(Energies));    % NEGF full transmission
    for ii = 1:length(Energies)
        Energy = Energies(ii);
        % Full Green's function with both absorbers
        GreensR = ((Energy+eta)*eye(size(totalSystem)) - totalSystem) \ eye(size(totalSystem));
        GreensA = GreensR';
        % NEGF transmission (matrix formula)
        T_NEGF(ii) = real(trace(GammaL * GreensR * GammaR * GreensA));
    end
end

function [G_NEGF] = sweepEM(EnergiesDot, totalSystem, GammaL, GammaR, Energy, idxSample)
    eta = 1j*1e-12;
    % We'll compute T(EF) for each gate (shift dot onsite)
    G_NEGF = zeros(size(EnergiesDot));
    for k = 1:length(EnergiesDot)
        EnergyDot = EnergiesDot(k);
        % Update dot onsite in H and recompute dot GR at EF for NEGF (full)
        H = totalSystem;
        H(idxSample, idxSample) = H(idxSample, idxSample) + EnergyDot;
        % Full Green's function with both absorbers at EF
        GreensR = ((Energy+eta)*eye(size(H)) - H) \ eye(size(H));
        GreensA = GreensR';
        % NEGF conductance (matrix formula)
        G_NEGF(k) = real(trace(GammaL * GreensR * GammaR * GreensA));
    end
end

%% -------------------- Other functions: BW -------------------- 
function [SigmaDot, GammaDot] = calcBWside(Energies, totalHamiltonian, Sigma, EnergyDot, idxSample)
    eta = 1j*1e-12;
    SigmaDot = zeros(size(Energies));
    GammaDot = zeros(size(Energies));
    for ii = 1:length(Energies)
        Energy = Energies(ii);
        %compute the dot Green's function for one-side-only absorber
        GreensR = ((Energy+eta)*eye(size(totalHamiltonian)) - totalHamiltonian - Sigma) \ eye(size(totalHamiltonian));
        GreensDot = GreensR(idxSample, idxSample);
        % effective self-energy due to lead as seen by dot
        SigmaDot(ii) = Energy - EnergyDot - 1/GreensDot;
        % Extract Gammas and shift (energy dependent)
        GammaDot(ii) = -2 * imag(SigmaDot(ii));
    end
end

function [T_BW, Delta] = calcBWboth(Energies, SigmaL, SigmaR, GammaL, GammaR, EnergyDot)
    T_BW   = zeros(size(Energies));
    Delta = zeros(size(Energies));
    for ii = 1:length(Energies)
        Energy = Energies(ii);
        % Effective total self-energy from both sides
        Sigma = SigmaL(ii) + SigmaR(ii);
        % Extract Gammas (energy dependent)
        Gamma = GammaL(ii) + GammaR(ii);
        % Extract shift (energy dependent)
        Delta(ii) = real(Sigma);  
        % level shift = Re Sigma_eff_total
        
        % Breit-Wigner transmission using extracted Gammas and shift
        T_BW(ii) = (GammaL(ii) * GammaR(ii)) / ((Energy - EnergyDot - Delta(ii))^2 + (Gamma/2)^2);
    end
end

function [SigmaDot, GammaDot] = sweepBWside(EnergiesDot, totalHamiltonian, Sigma, Energy, idxSample)
    eta = 1j*1e-12;
    SigmaDot = zeros(size(EnergiesDot));
    GammaDot = zeros(size(EnergiesDot));
    for k = 1:length(EnergiesDot)
        EnergyDot = EnergiesDot(k);
        % Update dot onsite in H and recompute dot GR at EF for NEGF (full)
        H = totalHamiltonian;
        H(idxSample, idxSample) = H(idxSample, idxSample) + EnergyDot;
        % compute the dot Green's function for one-side-only absorber
        GreensR = ((Energy+eta)*eye(size(H)) - H - Sigma) \ eye(size(H));
        GreensDot = GreensR(idxSample, idxSample);
        % effective self-energy due to lead as seen by dot
        SigmaDot(k) = Energy - EnergyDot - 1/GreensDot;
        % Extract Gammas and shift (energy dependent)
        GammaDot(k) = -2 * imag(SigmaDot(k));
    end
end

function [G_BW, varargout] = sweepBWboth(EnergiesDot, SigmaL, SigmaR, GammaL, GammaR, Energy)
    G_BW = zeros(size(EnergiesDot));
    Delta = zeros(size(EnergiesDot));
    for ii = 1:length(EnergiesDot)
        EnergyDot = EnergiesDot(ii);
        % Effective total self-energy from both sides
        Sigma = SigmaL(ii) + SigmaR(ii);
        % Extract Gammas (energy dependent)
        Gamma = GammaL(ii) + GammaR(ii);
        % Extract shift (energy dependent)
        Delta(ii) = real(Sigma);
        % level shift = Re Sigma_eff_total
    
        % Breit-Wigner conductance using extracted Gammas and shift
        G_BW(ii) = (GammaL(ii) * GammaR(ii)) / ((Energy - EnergyDot - Delta(ii))^2 + (Gamma/2)^2);
    end
    varargout{1} = Delta;
end