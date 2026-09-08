%% R7 - Energy-retention table across channel counts, and cheap nv derivation
%
% Two jobs:
%
%  1) Write the full "how many channels do I need" table (every n from 1 to
%     nCh) as CSV + figure, from the transform R2 already computed.
%
%  2) Optionally derive a SMALLER-nv mitosius from an existing larger-nv one,
%     without re-reading the 16 GB raw .dat.
%
%     This works because the retained subspaces are nested: for n' < n,
%     span(V(:,1:n')) is contained in span(V(:,1:n)). So with
%     W = Vret_n' * Vret_n2   ([n x n2], Vret_n orthonormal),
%
%         y_n2 = W.' * y_n          C_n2 = C_n * W
%
%     i.e. a small matrix multiply on the already-compressed data. Only works
%     DOWNWARD (n2 <= n); a larger nv needs the raw data again via R3.
%
% Outputs -> recon/ROVir/nv_energy_table.csv, qc/R7_nv_curve.png
%            mitosius/ROVir_<n2>/<bin>/ and C_rovir_<n2>.mat when deriving

clc; clearvars -except subject_num deriveTo binName srcNv variant; close all;

addpath(genpath('/home/debi/jaime/repos/MR-EyeTrack/recon'));
addpath(genpath('/home/debi/MatTechLab/monalisa'));

if ~exist('subject_num', 'var'); subject_num = 15;      end
% 'ROVir'  = SIR-ranked generalised eigenproblem (suppress the head)
% 'ROI-PCA'= ROI-energy-ranked (lambda -> inf); same solver, max compression.
% Outputs go to recon/<variant>/ so the two never mix.
if ~exist('variant',      'var'); variant     = 'ROVir'; end
if ~exist('srcNv',       'var'); srcNv       = 20;      end  % existing mitosius
if ~exist('deriveTo',    'var'); deriveTo    = [];      end  % e.g. [6 10] or [] for table only
if ~exist('binName',     'var'); binName     = 'woBin'; end

baseDir    = '/home/debi/jaime/repos/MR-EyeTrack/data/study';
subjectStr = sprintf('sub-%03d', subject_num);
reconDir   = fullfile(baseDir, subjectStr, 'recon');
rovirDir   = fullfile(reconDir, variant);
qcDir      = fullfile(rovirDir, 'qc');
if ~exist(qcDir, 'dir'); mkdir(qcDir); end

T   = load(fullfile(rovirDir, 'rovir_transform.mat'), 'V', 'transform');
tr  = T.transform;
nCh = tr.nCh;
n   = (1:nCh)';

%% 1. Table

roiE = tr.roiEnergy_retain(:);
intE = tr.intEnergy_retain(:);
gain = roiE ./ max(intE, eps);

tbl = table(n, ...
            roiE, intE, gain, ...
            tr.roiRetain(:), tr.intRetain(:), ...
            nCh ./ n, 100 * n / nCh, ...
    'VariableNames', {'n_channels', 'roi_energy_pct', 'int_energy_pct', ...
                      'contrast_gain', 'roi_frobenius_pct', 'int_frobenius_pct', ...
                      'compression_ratio', 'channels_kept_pct'});

csvPath = fullfile(rovirDir, 'nv_energy_table.csv');
writetable(tbl, csvPath);
fprintf('Full table (%d rows) written: %s\n\n', height(tbl), csvPath);

fprintf('%4s %12s %12s %10s %12s\n', 'n', 'ROI energy%', 'head energy%', 'gain', 'compression');
for k = [4 6 8 10 12 15 20 25 30 35 40 46 nCh]
    if k > nCh; continue; end
    fprintf('%4d %12.2f %12.2f %9.1fx %11.2fx\n', k, roiE(k), intE(k), gain(k), nCh/k);
end

%% Figure

f = figure('Visible', 'off', 'Position', [100 100 1150 460]);
tiledlayout(1, 2, 'Padding', 'compact');

nexttile;
plot(n, roiE, '.-', 'MarkerSize', 13, 'Color', [0.10 0.60 0.25], 'LineWidth', 1.3); hold on;
plot(n, intE, '.-', 'MarkerSize', 13, 'Color', [0.85 0.20 0.20], 'LineWidth', 1.3);
xline(20, '--k', 'n=20', 'LineWidth', 1.3, 'LabelVerticalAlignment', 'top');
xlabel('virtual channels kept (n)'); ylabel('% of signal energy retained');
legend({'eye ROI (keep)', 'rest of head (suppress)'}, 'Location', 'northwest');
title('How much signal survives at each channel count'); grid on; ylim([0 102]);

nexttile;
semilogy(n, gain, '.-', 'MarkerSize', 13, 'Color', [0.20 0.35 0.75], 'LineWidth', 1.3); hold on;
xline(20, '--k', 'LineWidth', 1.3);
yline(1, ':k');
xlabel('virtual channels kept (n)'); ylabel('ROI / head energy ratio vs original');
title('Contrast gain (higher = more eye-selective)'); grid on;

exportgraphics(f, fullfile(qcDir, 'R7_nv_curve.png'), 'Resolution', 130); close(f);
fprintf('\nFigure: %s\n', fullfile(qcDir, 'R7_nv_curve.png'));

%% 2. Derive smaller-nv datasets (optional)

if isempty(deriveTo); fprintf('\nNo deriveTo requested - table only.\n'); return; end

srcM = fullfile(reconDir, 'mitosius', sprintf('%s_%d', variant, srcNv), binName);
if ~exist(srcM, 'dir'); error('Source mitosius missing: %s', srcM); end

VretSrc = orth(T.V(:, 1:srcNv));
load(fullfile(rovirDir, sprintf('C_rovir_%d.mat', srcNv)), 'C_rovir');
Csrc = C_rovir;

ySrc  = bmMitosius_load(srcM, 'y');
tSrc  = bmMitosius_load(srcM, 't');
veSrc = bmMitosius_load(srcM, 've');
fprintf('\nSource mitosius nv=%d loaded: y %s\n', srcNv, mat2str(size(ySrc{1})));

for n2 = deriveTo(:)'
    if n2 > srcNv
        warning('nv=%d > source nv=%d - cannot derive upward, skipping. Use R3.', n2, srcNv);
        continue
    end
    Vret2 = orth(T.V(:, 1:n2));
    W     = VretSrc' * Vret2;                       % [srcNv x n2]

    % Guard the nesting assumption: Vret2 must lie in span(VretSrc)
    resid = norm(Vret2 - VretSrc * W, 'fro') / norm(Vret2, 'fro');
    fprintf('\nnv=%d: subspace nesting residual %.2e\n', n2, resid);
    if resid > 1e-8
        error('Subspace not nested (residual %.2e) - derive from raw with R3 instead.', resid);
    end

    % Mitosius stores y as [nPt x nCh], i.e. already transposed relative to
    % y_virt = Vret.' * y. In that layout the recombination is a PLAIN
    % right-multiply: y_stored_n2 = y_stored_src * W. (No conj - same trap as
    % everywhere else in this pipeline, see rovir_solve.m.)
    y2 = { ySrc{1} * W };
    e1 = sum(abs(y2{1}).^2, 'all');
    e2 = real(trace(W' * (ySrc{1}' * ySrc{1}) * W));
    assert(abs(e1-e2)/e2 < 1e-4, 'Derived-channel energy mismatch (%.2e)', abs(e1-e2)/e2);

    C_rovir = reshape(reshape(Csrc, [], srcNv) * W, ...
                      [size(Csrc,1) size(Csrc,2) size(Csrc,3) n2]);
    Cfile = fullfile(rovirDir, sprintf('C_rovir_%d.mat', n2));
    save(Cfile, 'C_rovir', '-v7.3');

    outM = fullfile(reconDir, 'mitosius', sprintf('%s_%d', variant, n2), binName);
    if ~exist(outM, 'dir'); mkdir(outM); end
    bmMitosius_create(outM, y2, tSrc, veSrc);
    fprintf('  derived mitosius: %s\n', outM);
    fprintf('  coil maps       : %s\n', Cfile);
end
