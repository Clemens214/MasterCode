function I = integrate_AGBGt(A, H, B, f, wmin, wmax, nquad)
%INTEGRATE_AGBGT  Numerically integrate  f(w) * A*G(w)*B*G(w)'  over [wmin,wmax]
%
%   I = integrate_AGBGT(A, H, B, f, wmin, wmax)
%   I = integrate_AGBGT(A, H, B, f, wmin, wmax, nquad)
%
%   Inputs:
%     A, B : matrices (same size as H)
%     H    : non-Hermitian system matrix (n x n)
%     f    : function handle for scalar weighting f(w)
%     wmin, wmax : integration limits (real scalars)
%     nquad (optional) : number of quadrature points (default = 32)
%
%   Output:
%     I : matrix integral  ∫ f(w) A G(w) B G(w)' dw
%
%   Notes:
%     - Uses Schur decomposition for stability: H = Q*T*Q'
%     - Computes G(w) via triangular back-substitution
%     - Integration done via Gauss–Legendre quadrature
%
%   Example:
%     n = 4;
%     H = randn(n) + 1i*randn(n);  % non-Hermitian
%     A = eye(n); B = eye(n);
%     f = @(w) exp(-w.^2);  % Gaussian weight
%     I = integrate_AGBGT(A, H, B, f, -5, 5);

    if nargin < 7, nquad = 32; end

    % 1. Schur decomposition for numerical stability
    [Q, T] = schur(H, 'complex');
    Atil = Q' * A * Q;
    Btil = Q' * B * Q;
    n = size(H, 1);

    % 2. Gauss–Legendre quadrature nodes and weights
    [x, w] = gauss_legendre(nquad);
    % Map from [-1,1] to [wmin,wmax]
    wm = (wmax + wmin)/2;
    wd = (wmax - wmin)/2;
    wpts = wm + wd*x;
    weights = wd*w;

    % 3. Accumulate integral in transformed basis
    I_tilde = zeros(n, n);
    for k = 1:nquad
        wk = wpts(k);
        fw = f(wk);

        % Compute G(wk) = (wk I - H)^(-1) stably via Schur form
        % Solve (wk I - T) X = I by triangular back-substitution
        X = solve_triangular(wk, T);

        Gtil = X;   % = (wk I - T)^(-1)
        term = fw * (Atil * Gtil * Btil * Gtil');
        I_tilde = I_tilde + weights(k) * term;
    end

    % 4. Transform back to original basis
    I = Q * I_tilde * Q';
end


%% --- Helper: Solve (wI - T) X = I for upper-triangular T ---
function X = solve_triangular(w, T)
    n = size(T,1);
    X = zeros(n);
    for j = 1:n
        % Solve column j of X
        rhs = zeros(n,1); rhs(j) = 1;
        for i = n:-1:1
            s = 0;
            if i < n
                s = T(i,i+1:n) * X(i+1:n,j);
            end
            X(i,j) = (rhs(i) + s) / (w - T(i,i));
        end
    end
end

%% --- Helper: Gauss–Legendre nodes and weights ---
function [x, w] = gauss_legendre(n)
    % Compute nodes/weights on [-1,1]
    beta = 0.5 ./ sqrt(1 - (2*(1:n-1)).^(-2));
    T = diag(beta,1) + diag(beta,-1);
    [V, D] = eig(T);
    x = diag(D);
    [x, idx] = sort(x);
    V = V(:, idx);
    w = 2 * V(1,:).^2;
end