function [sample] = makeSample(eigenenergy, hoppings, size, order, options)
% generate the sample's Hamiltonian
%
% Optional axially chiral end group (order = 2 only):
%   The p_x/p_y degeneracy is lifted along a local "director" d(phi) on the
%   first chiralSites bonds/sites of one end. The director spirals with a
%   fixed angle per bond/site,
%       phi_b = phase + (b-1)*twist ,   b = 1, ..., chiralSites ,
%   and enters through the traceless quadrupole matrix
%       D(phi) = [cos(2phi), sin(2phi); sin(2phi), -cos(2phi)] = 2 d d^T - I .
%   endType = "hopping": bond b carries the hopping block
%       T_b = hoppings * (I + anisotropy * D(phi_b)),
%     i.e. t(1+kappa) along the director, t(1-kappa) perpendicular to it.
%     This keeps the sublattice (particle-hole) symmetry of the chain.
%   endType = "onsite": site s carries the on-site term anisotropy * D(phi_s)
%     (crystal-field splitting by the substituents).
%   shift: additional orbital-independent on-site energy on the first
%     chiralSites sites of the end group (breaks particle-hole symmetry).
%   side: "left" (default) or "right" end of the chain.
%
%   twist = 0           -> anisotropic but achiral end group
%   twist -> -twist     -> enantiomer (mirror image)
%   anisotropy = 0      -> original rotationally isotropic chain
arguments
    eigenenergy
    hoppings
    size
    order = 2
    options.chiralSites = 0
    options.anisotropy = 0
    options.twist = 0
    options.phase = 0
    options.endType {mustBeMember(options.endType, ["hopping", "onsite"])} = "hopping"
    options.shift = 0
    options.side {mustBeMember(options.side, ["left", "right"])} = "left"
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

    % include the axially chiral end group
    M = options.chiralSites;
    if M > 0
        if order ~= 2
            error('makeSample:order', 'The chiral end group requires order = 2 (p_x, p_y).')
        end
        kappa = options.anisotropy;
        quad = @(phi) [cos(2*phi), sin(2*phi); sin(2*phi), -cos(2*phi)];
        if options.endType == "hopping"
            if M > size-1
                error('makeSample:sites', 'A hopping end group needs chiralSites <= size-1.')
            end
            for b = 1:M
                phi = options.phase + (b-1)*options.twist;
                if options.side == "left"
                    n = b;              % bond between sites n and n+1
                else
                    n = size - b;
                end
                rows = 2*n-1 : 2*n;
                cols = 2*n+1 : 2*n+2;
                block = hoppings * (eye(2) + kappa*quad(phi));
                sample(rows, cols) = block;
                sample(cols, rows) = block.';
            end
        elseif options.endType == "onsite"
            if M > size
                error('makeSample:sites', 'An on-site end group needs chiralSites <= size.')
            end
            for s = 1:M
                phi = options.phase + (s-1)*options.twist;
                if options.side == "left"
                    n = s;
                else
                    n = size + 1 - s;
                end
                range = 2*n-1 : 2*n;
                sample(range, range) = sample(range, range) + kappa*quad(phi);
            end
        end
        % orbital-independent shift of the end-group sites
        if options.shift ~= 0
            for s = 1:min(M, size)
                if options.side == "left"
                    n = s;
                else
                    n = size + 1 - s;
                end
                range = 2*n-1 : 2*n;
                sample(range, range) = sample(range, range) + options.shift*eye(2);
            end
        end
    end
end
