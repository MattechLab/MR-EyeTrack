% Joint gaze-bin reconstruction with spatial TV + temporal TV across bins.
%
% Extends bmSteva (MattechLab/monalisa) to reconstruct nGaze bins in a single
% joint ADMM solve.  The objective is:
%
%   min_X  sum_g (1/2)||F_g C x_g - y_g||^2_HY
%        + delta_s * sum_g ||grad_s x_g||_1     (spatial TV per bin)
%        + delta_t * ||D_t X||_1                (temporal TV across bins)
%
% where D_t applies cyclic finite differences along the gaze dimension.
%
% ADMM variable split:
%   z_s(:,:,g)  = grad_s x_g          (spatial auxiliary, per bin)
%   z_t(:,g)    = (D_t X)(:,g)        (temporal auxiliary, across bins)
%
% Usage:
%   X_cell = bmSteva_4fr(X_init, z_s, u_s, z_t, u_t, ...
%                        y, ve, C, Gu, Gut, frSize, ...
%                        delta_s, rho_s, delta_t, rho_t, ...
%                        nCGD, ve_max, nIter)
%
% Inputs:
%   X_init   - [nPt_u x nGaze] initial image (each col = one gaze bin);
%              pass [] to start from zeros.
%   z_s, u_s - [nPt_u x imDim x nGaze] spatial ADMM auxiliary / dual;
%              pass [] to auto-initialise.
%   z_t, u_t - [nPt_u x nGaze] temporal ADMM auxiliary / dual;
%              pass [] to auto-initialise.
%   y        - {nGaze x 1} cell; y{g} is the k-space data matrix for bin g.
%   ve       - {nGaze x 1} cell; ve{g} is the volume-element vector for bin g.
%   C        - coil sensitivity maps (shared across bins).
%   Gu, Gut  - {nGaze x 1} cell; Gu{g} / Gut{g} are the forward / adjoint
%              NUFFT structs for bin g (output of bmTraj2SparseMat, unwrapped).
%   frSize   - output image size (1 x imDim); pass [] to use N_u from Gu{1}.
%   delta_s  - spatial TV weight (scalar or [lower; upper] for schedule).
%   rho_s    - spatial ADMM step (scalar or [lower; upper]).
%   delta_t  - temporal TV weight.
%   rho_t    - temporal ADMM step.
%   nCGD     - conjugate-gradient iterations per ADMM outer iteration.
%   ve_max   - k-space volume-element clipping threshold; pass [] to auto.
%   nIter    - number of outer ADMM iterations.
%
% Output:
%   X_cell   - {nGaze x 1} cell of reconstructed 3-D images (single complex).

function X_cell = bmSteva_4fr(X, z_s, u_s, z_t, u_t, ...
                               y, ve, C, Gu, Gut, frSize, ...
                               delta_s, rho_s, delta_t, rho_t, ...
                               nCGD, ve_max, nIter)

myEps = 10*eps('single');
nGaze = numel(y);

% --- Common grid parameters (all bins share the same Cartesian grid) -------
N_u    = double(int32(Gu{1}.N_u(:)'));
frSize = double(frSize(:)');
if isempty(frSize)
    frSize = N_u;
end
nPt_u = prod(frSize(:));
imDim = numel(N_u);
dK_u  = double(single(Gu{1}.d_u(:)'));
dX_u  = single((1./single(dK_u))./single(N_u));

HX = single(prod(dX_u(:)));   % image-space integration weight
HZ = single(prod(dX_u(:)));   % gradient-space integration weight (= HX)

% --- Coil maps & NUFFT deapodisation kernel (shared across bins) -----------
nCh = size(y{1}, 2);
C   = single(bmColReshape(C, N_u));
KFC      = single(bmKF(          C,       N_u, frSize, dK_u, nCh, ...
                       Gu{1}.kernel_type, Gu{1}.nWin, Gu{1}.kernelParam));
KFC_conj = single(bmKF_conj(conj(C),     N_u, frSize, dK_u, nCh, ...
                       Gu{1}.kernel_type, Gu{1}.nWin, Gu{1}.kernelParam));

% --- Per-bin k-space metric HY ---------------------------------------------
HY = cell(nGaze, 1);
for g = 1:nGaze
    y{g} = single(y{g});
    if isempty(ve_max)
        vm = max(ve{g}(:));
    else
        vm = ve_max;
    end
    HY{g} = min(single(bmY_ve_reshape(ve{g}, size(y{g}))), single(vm));
end

% --- Initialise stacked image X as [nPt_u x nGaze] -------------------------
if isempty(X)
    X = complex(zeros(nPt_u, nGaze, 'single'));
else
    Xtmp = complex(zeros(nPt_u, nGaze, 'single'));
    for g = 1:nGaze
        Xtmp(:, g) = single(bmColReshape(X(:, g), frSize));
    end
    X = Xtmp;
end

% --- Initialise z_s, u_s [nPt_u x imDim x nGaze] --------------------------
if isempty(z_s)
    z_s = complex(zeros(nPt_u, imDim, nGaze, 'single'));
    for g = 1:nGaze
        z_s(:, :, g) = bmBackGradient(X(:, g), frSize, dX_u);
    end
end
if isempty(u_s)
    u_s = complex(zeros(nPt_u, imDim, nGaze, 'single'));
end

% --- Initialise z_t, u_t [nPt_u x nGaze] ----------------------------------
if isempty(z_t)
    z_t = private_Dt(X, nGaze);
end
if isempty(u_t)
    u_t = complex(zeros(nPt_u, nGaze, 'single'));
end

[delta_s, rho_s] = private_init_delta_rho(delta_s, rho_s, nIter);
[delta_t, rho_t] = private_init_delta_rho(delta_t, rho_t, nIter);

% ==========================================================================
% ADMM loop
% ==========================================================================
for c = 1:nIter

    % -----------------------------------------------------------------------
    % x-update: nCGD steps of conjugate-gradient descent
    % -----------------------------------------------------------------------

    % Compute residuals at current X
    res_y   = cell(nGaze, 1);
    res_z_s = cell(nGaze, 1);
    for g = 1:nGaze
        res_y{g}   = y{g} - bmShanna(X(:,g), Gu{g}, KFC, frSize, 'MATLAB');
        res_z_s{g} = (z_s(:,:,g) - u_s(:,:,g)) - bmBackGradient(X(:,g), frSize, dX_u);
    end
    D_tX    = private_Dt(X, nGaze);
    res_z_t = (z_t - u_t) - D_tX;

    % Initial CG search direction = negative gradient of augmented Lagrangian
    dagA_res_next = complex(zeros(nPt_u, nGaze, 'single'));
    for g = 1:nGaze
        dagM_g  = (1/HX) * bmNakatsha(HY{g}.*res_y{g}, Gut{g}, KFC_conj, true, frSize, 'MATLAB');
        dagFs_g = rho_s(1,c) * (1/HX) * bmBackGradientT(HZ*res_z_s{g}, frSize, dX_u);
        dagFt_g = rho_t(1,c) * (1/HX) * private_DtT_col(HZ*res_z_t, g, nGaze);
        dagA_res_next(:, g) = dagM_g + dagFs_g + dagFt_g;
    end

    p_next            = dagA_res_next;
    sqn_dagA_res_next = real(dagA_res_next(:)' * (HX * dagA_res_next(:)));

    for i = 1:nCGD

        res_y_curr        = res_y;
        res_z_s_curr      = res_z_s;
        res_z_t_curr      = res_z_t;
        sqn_dagA_res_curr = sqn_dagA_res_next;
        p_curr            = p_next;

        if sqn_dagA_res_curr < myEps
            break;
        end

        % Apply the joint normal operator to p_curr
        Mp_curr   = cell(nGaze, 1);
        Fp_s_curr = cell(nGaze, 1);
        sqn_Mp    = 0;
        sqn_Fp_s  = 0;
        for g = 1:nGaze
            Mp_curr{g}   = bmShanna(p_curr(:,g), Gu{g}, KFC, frSize, 'MATLAB');
            Fp_s_curr{g} = bmBackGradient(p_curr(:,g), frSize, dX_u);
            sqn_Mp   = sqn_Mp   + real(Mp_curr{g}(:)' * (HY{g}(:) .* Mp_curr{g}(:)));
            sqn_Fp_s = sqn_Fp_s + real(Fp_s_curr{g}(:)' * (rho_s(1,c)*HZ*Fp_s_curr{g}(:)));
        end
        Fp_t_curr   = private_Dt(p_curr, nGaze);
        sqn_Fp_t    = real(Fp_t_curr(:)' * (rho_t(1,c)*HZ*Fp_t_curr(:)));
        sqn_Ap_curr = sqn_Mp + sqn_Fp_s + sqn_Fp_t;

        a = sqn_dagA_res_curr / sqn_Ap_curr;
        X = X + a*p_curr;

        if i == nCGD
            break;
        end

        % Propagate residuals (avoid redundant NUFFT calls)
        for g = 1:nGaze
            res_y{g}   = res_y_curr{g}   - a*Mp_curr{g};
            res_z_s{g} = res_z_s_curr{g} - a*Fp_s_curr{g};
        end
        res_z_t = res_z_t_curr - a*Fp_t_curr;

        % Recompute gradient for next CG direction
        dagA_res_next = complex(zeros(nPt_u, nGaze, 'single'));
        for g = 1:nGaze
            dagM_g  = (1/HX) * bmNakatsha(HY{g}.*res_y{g}, Gut{g}, KFC_conj, true, frSize, 'MATLAB');
            dagFs_g = rho_s(1,c) * (1/HX) * bmBackGradientT(HZ*res_z_s{g}, frSize, dX_u);
            dagFt_g = rho_t(1,c) * (1/HX) * private_DtT_col(HZ*res_z_t, g, nGaze);
            dagA_res_next(:, g) = dagM_g + dagFs_g + dagFt_g;
        end
        sqn_dagA_res_next = real(dagA_res_next(:)' * (HX * dagA_res_next(:)));
        b      = sqn_dagA_res_next / sqn_dagA_res_curr;
        p_next = dagA_res_next + b*p_curr;

    end % CG loop

    % -----------------------------------------------------------------------
    % z_s-update: per-bin spatial soft-threshold
    % -----------------------------------------------------------------------
    if rho_s(1,c) > 0
        for g = 1:nGaze
            bGx_u      = bmBackGradient(X(:,g), frSize, dX_u) + u_s(:,:,g);
            z_s(:,:,g) = bmProx_oneNorm(bGx_u, delta_s(1,c)/rho_s(1,c));
            u_s(:,:,g) = bGx_u - z_s(:,:,g);
        end
    end

    % -----------------------------------------------------------------------
    % z_t-update: temporal soft-threshold on cyclic gaze differences
    % Skip when rho_t = 0 (temporal coupling fully relaxed this iteration).
    % -----------------------------------------------------------------------
    if rho_t(1,c) > 0
        D_tX_u = private_Dt(X, nGaze) + u_t;
        z_t     = bmProx_oneNorm(D_tX_u, delta_t(1,c)/rho_t(1,c));
        u_t     = D_tX_u - z_t;
    end

    % --- monitoring --------------------------------------------------------
    R    = 0;
    TV_s = 0;
    for g = 1:nGaze
        tmp = y{g} - bmShanna(X(:,g), Gu{g}, KFC, frSize, 'MATLAB');
        R   = R + real(tmp(:)' * (HY{g}(:) .* tmp(:)));
        tmp  = bmBackGradient(X(:,g), frSize, dX_u);
        TV_s = TV_s + HZ*sum(abs(real(tmp(:)))) + HZ*sum(abs(imag(tmp(:))));
    end
    tmp  = private_Dt(X, nGaze);
    TV_t = HZ*sum(abs(real(tmp(:)))) + HZ*sum(abs(imag(tmp(:))));
    fprintf('Iter %03d/%03d | R=%.4e | TV_s=%.4e | TV_t=%.4e\n', c, nIter, R, TV_s, TV_t);

end % ADMM loop

% --- Reshape each bin back to frSize and return as cell -------------------
X_cell = cell(nGaze, 1);
for g = 1:nGaze
    X_cell{g} = bmBlockReshape(X(:, g), frSize);
end

end % function


% ==========================================================================
% Private helpers
% ==========================================================================

% Cyclic forward differences along the gaze dimension.
% DtX(:, g) = X(:, mod(g, nGaze)+1) - X(:, g)
function DtX = private_Dt(X, nGaze)
DtX = complex(zeros(size(X), 'single'));
for g = 1:nGaze
    g_next    = mod(g, nGaze) + 1;
    DtX(:, g) = X(:, g_next) - X(:, g);
end
end

% Column g of D_t^T * V:  V(:, g_prev) - V(:, g)
% g_prev = cyclic predecessor of g  (mod(g-2, nGaze)+1 in 1-based indexing)
function out = private_DtT_col(V, g, nGaze)
g_prev = mod(g-2, nGaze) + 1;
out    = V(:, g_prev) - V(:, g);
end

% Schedule delta and rho as row vectors of length nIter.
%   scalar          → constant across all iterations
%   [lo; hi]        → linear ramp from lo to hi  (warmup / continuation)
%   [d1;d2;...;dN]  → arbitrary per-iteration schedule (N must equal nIter)
function [delta, rho] = private_init_delta_rho(delta, rho, nIter)
delta = single(abs(delta(:)));
rho   = single(abs(rho(:)));
if numel(delta) == 1
    delta = repmat(delta, 1, nIter);
elseif numel(delta) == 2
    delta = single(linspace(delta(1), delta(2), nIter));
% else: full nIter-length schedule passed through as-is
end
delta = delta(:)';
if numel(rho) == 1
    rho = repmat(rho, 1, nIter);
elseif numel(rho) == 2
    rho = single(linspace(rho(1), rho(2), nIter));
end
rho = rho(:)';
end
