%% R2 - Compute the ROVir transform from the per-coil gridded images
%
% Consumes woBin/x0_noC_48.mat (per-coil images, produced by
% recon/4-Recon/S02_chuv_woC.m with matrix_size = 48) and ROVir/masks.mat
% (from R1), and produces the coil recombination matrix V.
%
% The transform is computed once per subject from the *unbinned* data, which
% has the full readout count and therefore the best SNR. It is then applied
% unchanged to every gaze bin in R3: coil geometry does not depend on where the
% subject is looking.
%
% Outputs -> data/study/sub-NNN/recon/ROVir/rovir_transform.mat  (+ QC pngs)

clc; close all;

if ~exist('subject_num',  'var'); subject_num  = 15;    end
% 'ROVir'  = SIR-ranked generalised eigenproblem (suppress the head)
% 'ROI-PCA'= ROI-energy-ranked (lambda -> inf); same solver, max compression.
% Outputs go to recon/<variant>/ so the two never mix.
if ~exist('variant',      'var'); variant     = 'ROVir'; end
if ~exist('matrix_size',  'var'); matrix_size  = 48;    end
if ~exist('lambda',       'var'); lambda       = 1e-4;  end
% Default criterion: suppress the rest of the head to <= selectThr % of its
% energy, and let the channel count fall out of that. This encodes the actual
% goal ("keep the eyes, drop everything else, shrink the data") directly.
% 'operating_knee' is available but sits at n~37 here, only a 1.4x reduction.
if ~exist('selectMethod', 'var'); selectMethod = 'int_budget'; end
if ~exist('selectThr',    'var'); selectThr    = 2;     end  % % interference energy budget
if ~exist('nvManual',     'var'); nvManual     = [];    end  % override, e.g. 12

addpath(genpath('/home/debi/jaime/repos/MR-EyeTrack/recon'));

baseDir    = '/home/debi/jaime/repos/MR-EyeTrack/data/study';
subjectStr = sprintf('sub-%03d', subject_num);
reconDir   = fullfile(baseDir, subjectStr, 'recon');
rovirDir   = fullfile(reconDir, variant);
qcDir      = fullfile(rovirDir, 'qc');
if ~exist(qcDir, 'dir'); mkdir(qcDir); end

%% Load masks and per-coil images

M = load(fullfile(rovirDir, 'masks.mat'));
roiMask = M.roiMask; intMask = M.intMask; headMask = M.headMask;
N = size(roiMask, 1);

x0Path = fullfile(reconDir, 'woBin', sprintf('x0_noC_%d.mat', matrix_size));
P = load(x0Path, 'x0');
x0 = P.x0;
nCh = numel(x0);
fprintf('Loaded %s\n  %d coils, %s per coil\n', x0Path, nCh, mat2str(size(x0{1})));

if ~isequal(size(x0{1}), size(roiMask))
    error('Per-coil images are %s but masks are %s - regenerate masks at the matching grid.', ...
          mat2str(size(x0{1})), mat2str(size(roiMask)));
end

% [nVox x nCh]
Xmat = zeros(numel(roiMask), nCh, 'single');
for c = 1:nCh
    Xmat(:, c) = x0{c}(:);
end
clear x0 P;

%% Solve

fprintf('\nSolving ROVir (%d ROI voxels, %d interference voxels, lambda=%g)...\n', ...
        nnz(roiMask), nnz(intMask), lambda);
R = rovir_solve(Xmat, find(roiMask), find(intMask), struct('lambda', lambda));
fprintf('  eigenvalue range: %.3e .. %.3e\n', min(R.eigval), max(R.eigval));
fprintf('  SIR of best virtual coil: %.2f   (worst: %.4f)\n', max(R.sir), min(R.sir));

%% Choose how many virtual coils to keep

switch selectMethod
    case 'operating_knee'
        % Knee of the (interference retained, ROI retained) operating curve:
        % the point furthest from the chord joining its two endpoints. This is
        % the best trade of eye signal kept against head signal suppressed,
        % without having to pick an arbitrary percentage.
        xq = R.intEnergy_retain(:) / 100; y = R.roiEnergy_retain(:) / 100;
        d = abs((y(end)-y(1))*xq - (xq(end)-xq(1))*y + xq(end)*y(1) - y(end)*xq(1));
        [~, nv] = max(d);
        crit = 'knee of the ROI-vs-interference operating curve';
    case 'roi_retained'
        nv = find(R.roiEnergy_retain >= selectThr, 1, 'first');
        if isempty(nv); nv = nCh; end
        crit = sprintf('first n retaining >= %g%% of ROI signal energy', selectThr);
    case 'int_budget'
        nv = find(R.intEnergy_retain > selectThr, 1, 'first');
        if isempty(nv); nv = nCh; else; nv = max(nv - 1, 1); end
        crit = sprintf('largest n keeping interference energy <= %g%%', selectThr);
    case 'sir'
        % elbow of the SIR curve: point of maximum distance to the chord
        % joining its two endpoints
        y = R.sir(:) / max(R.sir); xq = (1:nCh)' / nCh;
        d = abs((y(end)-y(1))*xq - (xq(end)-xq(1))*y + xq(end)*y(1) - y(end)*xq(1));
        [~, nv] = max(d);
        crit = 'elbow of the SIR curve';
    otherwise
        error('Unknown selectMethod "%s"', selectMethod);
end
if ~isempty(nvManual); nv = nvManual; crit = 'manual override'; end

fprintf('\nRetaining %d of %d virtual coils (%s)\n', nv, nCh, crit);
fprintf('  ROI energy retained          : %6.2f %%\n', R.roiEnergy_retain(nv));
fprintf('  Interference energy retained : %6.2f %%\n', R.intEnergy_retain(nv));
fprintf('  contrast gain (ROI/int)      : %6.2fx\n', ...
        R.roiEnergy_retain(nv)/max(R.intEnergy_retain(nv), eps));
fprintf('  [Frobenius metric, as the abstract reports it: ROI %.2f%%, int %.2f%%]\n', ...
        R.roiRetain(nv), R.intRetain(nv));
fprintf('  k-space channel reduction    : %d -> %d  (%.1fx smaller)\n', ...
        nCh, nv, nCh/nv);

% Full trade-off table, so the choice above can be overruled from the data
fprintf('\n  %4s %10s %10s %10s\n', 'n', 'ROI%', 'int%', 'gain');
for n = [4 6 8 10 12 15 20 25 30 35 40 nCh]
    if n > nCh; continue; end
    fprintf('  %4d %10.2f %10.2f %9.1fx\n', n, R.roiEnergy_retain(n), ...
            R.intEnergy_retain(n), R.roiEnergy_retain(n)/max(R.intEnergy_retain(n), eps));
end

% Orthonormalise the retained span. The generalised eigenvectors are
% B-orthogonal, so feeding them to the recon directly would correlate the
% noise between virtual channels; monalisa's least-squares/CS recon assumes
% roughly white noise across channels. orth() keeps the same subspace.
Vret = orth(R.V(:, 1:nv));

% Guard the transpose/conjugate convention. A = X^H X pairs with x_virt = X*V;
% using conj(V) here instead silently destroys the suppression while still
% producing a believable image, so assert the energy identity rather than
% trusting it.
Xchk   = Xmat * Vret;
eDirect = sum(abs(Xchk(roiMask(:), :)).^2, 'all');
eQuad   = real(trace(Vret' * R.A * Vret));
relErr  = abs(eDirect - eQuad) / eQuad;
fprintf('\nConvention check: direct %.6e vs quadratic %.6e (rel. err %.2e)\n', ...
        eDirect, eQuad, relErr);
if relErr > 1e-4
    error(['ROVir convention mismatch (rel. err %.2e). x_virt = Xmat*Vret must ' ...
           'agree with trace(Vret^H A Vret); check for a stray conj/ctranspose.'], relErr);
end

%% Save

V         = R.V;
transform = struct('Vret', Vret, 'nv', nv, 'nCh', nCh, ...
                   'selectMethod', selectMethod, 'selectThr', selectThr, ...
                   'roiRetain', R.roiRetain, 'intRetain', R.intRetain, ...
                   'roiEnergy_retain', R.roiEnergy_retain, ...
                   'intEnergy_retain', R.intEnergy_retain, ...
                   'sir', R.sir, 'eigval', R.eigval, ...
                   'roiEnergy', R.roiEnergy, 'intEnergy', R.intEnergy, ...
                   'lambda', lambda, 'matrix_size', matrix_size, ...
                   'subject_num', subject_num, 'created', datetime('now'));

outPath = fullfile(rovirDir, 'rovir_transform.mat');
save(outPath, 'Vret', 'V', 'transform', '-v7.3');
fprintf('\nTransform saved: %s\n', outPath);

%% QC 1 - metric curves

f = figure('Visible', 'off', 'Position', [100 100 1200 850]);
tiledlayout(2, 2, 'Padding', 'compact');

nexttile;
semilogy(1:nCh, R.sir, '.-', 'MarkerSize', 14); hold on;
xline(nv, '--r', 'LineWidth', 1.5);
xlabel('virtual coil'); ylabel('SIR (ROI / interference)');
title('Signal-to-interference ratio per virtual coil'); grid on;

nexttile;
plot(1:nCh, R.roiEnergy_retain, '.-', 'MarkerSize', 14, 'Color', [0.1 0.6 0.2]); hold on;
plot(1:nCh, R.intEnergy_retain, '.-', 'MarkerSize', 14, 'Color', [0.85 0.2 0.2]);
xline(nv, '--k', 'LineWidth', 1.5);
xlabel('virtual coils retained'); ylabel('% retained');
legend({'eye ROI', 'rest of head'}, 'Location', 'east');
title('Cumulative signal ENERGY retention'); grid on; ylim([0 105]);

nexttile;
semilogy(1:nCh, R.roiEnergy, '.-', 'MarkerSize', 14, 'Color', [0.1 0.6 0.2]); hold on;
semilogy(1:nCh, R.intEnergy, '.-', 'MarkerSize', 14, 'Color', [0.85 0.2 0.2]);
xline(nv, '--k', 'LineWidth', 1.5);
xlabel('virtual coil'); ylabel('energy');
legend({'eye ROI', 'rest of head'}, 'Location', 'northeast');
title('Per-virtual-coil energy'); grid on;

nexttile;
plot(R.intEnergy_retain, R.roiEnergy_retain, '.-', 'MarkerSize', 14); hold on;
plot(R.intEnergy_retain(nv), R.roiEnergy_retain(nv), 'rp', 'MarkerSize', 18, 'MarkerFaceColor', 'r');
xlabel('% interference energy retained'); ylabel('% ROI energy retained');
title(sprintf('Operating curve (chosen: n=%d)', nv)); grid on;

qc1 = fullfile(qcDir, sprintf('R2_metrics_%s.png', subjectStr));
exportgraphics(f, qc1, 'Resolution', 130); close(f);
fprintf('QC metrics : %s\n', qc1);

%% QC 2 - RSS image before vs after, same window

Xv   = Xmat * Vret;                             % virtual-coil images (plain transpose)
rssO = reshape(sqrt(sum(abs(Xmat).^2, 2)), [N N N]);
rssV = reshape(sqrt(sum(abs(Xv).^2,  2)), [N N N]);

% Match the two on the ROI so the comparison shows suppression outside it
% rather than an overall gain difference.
sc   = mean(rssO(roiMask)) / mean(rssV(roiMask));
rssV = rssV * sc;

fprintf('\nMean magnitude (arbitrary units, ROI-matched):\n');
fprintf('  ROI        : original %.4g -> ROVir %.4g  (%.2fx)\n', ...
        mean(rssO(roiMask)), mean(rssV(roiMask)), mean(rssV(roiMask))/mean(rssO(roiMask)));
fprintf('  rest of head: original %.4g -> ROVir %.4g  (%.2fx)\n', ...
        mean(rssO(intMask)), mean(rssV(intMask)), mean(rssV(intMask))/mean(rssO(intMask)));

hi = prctile(rssO(headMask), 99);
zc = round(median(find(squeeze(any(any(roiMask,1),2)))));
zs = zc-4 : 2 : zc+4; zs = zs(zs >= 1 & zs <= N);
up = 8;
tile = zeros(2*N*up, numel(zs)*N*up, 'single');
for i = 1:numel(zs)
    k = zs(i);
    tile(1:N*up,        (i-1)*N*up+(1:N*up)) = kron(min(rssO(:,:,k)/hi, 1), ones(up));
    tile(N*up+(1:N*up), (i-1)*N*up+(1:N*up)) = kron(min(rssV(:,:,k)/hi, 1), ones(up));
end
qc2 = fullfile(qcDir, sprintf('R2_rss_before_after_%s.png', subjectStr));
imwrite(tile, qc2);
fprintf('QC images  : %s   (top row = all %d coils, bottom = %d ROVir coils, slices %s)\n', ...
        qc2, nCh, nv, mat2str(zs));
