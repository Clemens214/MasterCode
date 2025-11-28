function [] = makeSystemSI (sample, sizeSample, orderSample, sizeLead, hoppingLead, hoppingsInter, leadVals, options)
arguments
    sample
    hoppingsInter
    eta = 1E-12;
    options.check = true
    options.checkMore = false
end
    hoppingsLeft = hoppingsInter(1, :);
    hoppingsRight = hoppingsInter(2, :);

    sigmaL = zeros(size(sample));
    sigmaR = zeros(size(sample));

    for 1: length()
    end
end

function [G, totalSystem, gammaL, gammaR] = prepareSI (w, Hamiltonian, t_C, t_B, N)
    eta = 0.000000001; %0.000000001
    sigmaL = zeros(N, N);
    sigmaR = zeros(N, N);
    
    sigmaL(1,1) = t_C^2 *(w +eta*1i)*(1-sqrt(1-(4*t_B^2)/(w+eta *1i)^2))/(2*t_B^2);
    sigmaR(N,N) = t_C^2 *(w +eta*1i)*(1-sqrt(1-(4*t_B^2)/(w+eta *1i)^2))/(2*t_B^2);
    
    sigmaLdagger = sigmaL';
    sigmaRdagger = sigmaR';
    
    gammaR = 1j*(sigmaR - sigmaRdagger);
    gammaL = 1j*(sigmaL - sigmaLdagger);
    
    totalSystem = Hamiltonian + sigmaL + sigmaR;

    invG = (w + 1i*eta)*eye(N)- Hamiltonian - sigmaL -sigmaR;
    G = inv(invG);
end