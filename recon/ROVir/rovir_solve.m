function out = rovir_solve(Xmat, roiIdx, intIdx, opts)
% ROVIR_SOLVE  Region-optimised virtual coil transform.
%
%   out = rovir_solve(Xmat, roiIdx, intIdx, opts)
%
% Finds a linear recombination of the physical coils that concentrates signal
% in a region of interest and suppresses it elsewhere, by solving the
% generalised eigenvalue problem
%
%       A v = lambda B v,     A = sum_{v in ROI} x_v x_v^H
%                             B = sum_{v in INT} x_v x_v^H
%
% The leading eigenvector maximises the signal-to-interference ratio; the
% eigenvectors sorted by descending eigenvalue give a ranked basis.
%
% Reference: Kim et al., "Region-optimized virtual (ROVir) coils", MRM 2021.
% Follows the MR RawDeface implementation (github.com/chiew-group/MRI-defacer),
% with the ROI and interference roles swapped: there the brain is kept and the
% face suppressed, here the orbits are kept and the rest of the head suppressed.
%
% Inputs
%   Xmat   [nVox x nCh] complex per-coil image values, one row per voxel
%   roiIdx linear indices (into the voxel dimension) of the region to keep
%   intIdx linear indices of the region to suppress
%   opts   struct, optional fields:
%            .lambda  Tikhonov load on B, as a fraction of mean(diag(B))
%                     (default 1e-4). B is formed from tens of thousands of
%                     voxels so it is normally well conditioned, but coils
%                     that see almost no interference can still make it
%                     near-singular, and the Cholesky factor below needs it
%                     positive definite.
%
% Output struct
%   V          [nCh x nCh] eigenvectors, columns sorted by descending SIR,
%              each normalised to unit l2 norm
%   eigval     [nCh x 1] corresponding generalised eigenvalues
%   roiEnergy  [nCh x 1] ROI signal energy of each virtual coil
%   intEnergy  [nCh x 1] interference energy of each virtual coil
%   sir        [nCh x 1] roiEnergy ./ intEnergy
%   roiRetain  [nCh x 1] % of ROI signal retained by the leading n coils,
%              Frobenius metric ||P A P||_F / ||A||_F -- what MR RawDeface
%              reports, kept here for comparability with the abstract
%   intRetain  [nCh x 1] same, for the interference region
%   roiEnergy_retain [nCh x 1] % of ROI signal *energy*, trace(P A)/trace(A)
%   intEnergy_retain [nCh x 1] same, for the interference region
%   A, B       the covariance matrices
%
% Prefer the *_retain energy curves when judging the result. The Frobenius
% metric squares the eigenvalue contributions a second time, so it is
% dominated by the leading eigenvector and reads far more favourably than the
% actual fraction of signal energy retained: on sub-015 at n = 10 it reports
% 0.35% interference retained where the true energy figure is 0.29% -- close
% here, but the two diverge sharply at larger n.
%
% CONVENTION (easy to get wrong, and silent when wrong):
%   A = X^H X, which pairs with forming virtual coils by PLAIN transpose,
%       y_virt = V.' * y        and       x_virt = X * V
%   NOT V' * y / X * conj(V). Using the conjugate here still produces a
%   plausible-looking image, but the retained subspace no longer matches the
%   one the eigenproblem optimised, and the suppression silently vanishes.
%   The energy identity  sum_{v in ROI} ||V.' x_v||^2 == trace(V^H A V)
%   is the check that catches it; R2 asserts it.
%
% Note on the retention curves: they are computed on the *orthonormalised*
% span of the leading n eigenvectors, because that is the subspace the data is
% actually projected onto. The generalised eigenvectors are B-orthogonal, not
% orthonormal, so using them directly would misstate the retained energy.

if nargin < 4 || isempty(opts); opts = struct(); end
if ~isfield(opts, 'lambda'); opts.lambda = 1e-4; end

nCh = size(Xmat, 2);

%% Covariance matrices

Xroi = Xmat(roiIdx, :);
Xint = Xmat(intIdx, :);

A = double(Xroi' * Xroi);
B = double(Xint' * Xint);

% Force exact Hermitian symmetry (the products are Hermitian up to roundoff,
% and chol/eig below are picky about it)
A = (A + A') / 2;
B = (B + B') / 2;

%% Regularise and solve

Breg = B + opts.lambda * mean(real(diag(B))) * eye(nCh);

% Cholesky whitening rather than eig(A,B) directly: this is what scipy's
% eigh(A, B) does internally, and it keeps the problem Hermitian throughout
% instead of handing MATLAB a general (non-symmetric) pencil.
[L, p] = chol(Breg, 'lower');
if p ~= 0
    error(['Interference covariance is not positive definite even after ' ...
           'Tikhonov loading (lambda = %g). Increase opts.lambda, or check ' ...
           'that the interference mask is non-empty.'], opts.lambda);
end

M = L \ A / L';
M = (M + M') / 2;
[W, D] = eig(M);
eigval = real(diag(D));

V = L' \ W;

% Unit-norm columns, then sort by descending eigenvalue (descending SIR)
V = V ./ vecnorm(V, 2, 1);
[eigval, ord] = sort(eigval, 'descend');
V = V(:, ord);

%% Per-virtual-coil energies

roiEnergy = real(diag(V' * A * V));
intEnergy = real(diag(V' * B * V));
sir       = roiEnergy ./ (intEnergy + eps);

%% Cumulative retention curves

normA = norm(A, 'fro');
normB = norm(B, 'fro');
trA   = real(trace(A));
trB   = real(trace(B));
roiRetain = zeros(nCh, 1);
intRetain = zeros(nCh, 1);
roiEnergy_retain = zeros(nCh, 1);
intEnergy_retain = zeros(nCh, 1);
for n = 1:nCh
    Q = orth(V(:, 1:n));          % orthonormal basis for the retained subspace
    P = Q * Q';                   % orthogonal projector onto it
    roiRetain(n) = norm(P * A * P, 'fro') / normA * 100;
    intRetain(n) = norm(P * B * P, 'fro') / normB * 100;
    roiEnergy_retain(n) = real(trace(Q' * A * Q)) / trA * 100;
    intEnergy_retain(n) = real(trace(Q' * B * Q)) / trB * 100;
end

out = struct('V', V, 'eigval', eigval, 'roiEnergy', roiEnergy, ...
             'intEnergy', intEnergy, 'sir', sir, 'roiRetain', roiRetain, ...
             'intRetain', intRetain, 'roiEnergy_retain', roiEnergy_retain, ...
             'intEnergy_retain', intEnergy_retain, ...
             'A', A, 'B', B, 'lambda', opts.lambda);
end
