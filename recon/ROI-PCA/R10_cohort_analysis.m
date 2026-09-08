%% R10 - Cohort analysis: does ROI-PCA hold across all subjects?
%
% For every subject with a finished ROI-PCA reconstruction, measure the eye ROI
% against the 52-channel reference and summarise across the cohort.
%
% Both comparisons are reported, but **the gridded one is the one to quote**.
% STEVA's `delta` is scale-dependent and the 52-channel references were
% normalised with the S3 interactive roipoly while ROI-PCA uses the ROI-based
% constant, so `delta = 1` is not equal regularisation between them. The
% gridded recons have no regularisation at all, so that confound cannot apply.
%
% Outputs -> recon/ROI-PCA/cohort_results.csv
%            recon/ROI-PCA/figures/figE_cohort_*.png

clc; clearvars -except subs nv; close all;

addpath(genpath('/home/debi/jaime/repos/MR-EyeTrack/recon'));
addpath(genpath('/home/debi/MatTechLab/monalisa'));

if ~exist('subs', 'var'); subs = 1:15; end
if ~exist('nv',   'var'); nv   = 4;    end

repoRoot = '/home/debi/jaime/repos/MR-EyeTrack';
base     = fullfile(repoRoot, 'data/study');
outDir   = fullfile(repoRoot, 'recon/ROI-PCA');
figDir   = fullfile(outDir, 'figures');
if ~exist(figDir, 'dir'); mkdir(figDir); end

rows = {}; keep = []; orbG = {}; orbR = {};

for s = subs
    sd = fullfile(base, sprintf('sub-%03d', s), 'recon');
    pd = fullfile(sd, 'ROI-PCA');
    fS = fullfile(pd, sprintf('x_steva_ROI-PCA_%d_woBin_nIter_20_delta_1.000.mat', nv));
    fG = fullfile(pd, sprintf('x0_ROI-PCA_%d_woBin.mat', nv));
    rS = fullfile(sd, 'woBin', 'x_steva_nIter_20_delta_1.000.mat');
    rG = fullfile(sd, 'woBin', 'x0.mat');
    if ~all(cellfun(@(p) exist(p,'file')==2, {fS,fG,rS,rG}))
        fprintf('sub-%03d: incomplete, skipping\n', s); continue
    end

    xoG = abs(single(grab(rG,'x0'))); xoS = abs(single(grab(rS,'x')));
    xpG = abs(single(grab(fG,'x0'))); xpS = abs(single(grab(fS,'x')));
    Nb  = size(xoG,1);

    M   = load(fullfile(pd,'masks.mat'),'roiMask','headMask');
    N48 = size(M.roiMask,1);
    roi = bmImResize(single(M.roiMask ),[N48 N48 N48],[Nb Nb Nb])>0.5;
    hd  = bmImResize(single(M.headMask),[N48 N48 N48],[Nb Nb Nb])>0.5;
    air = ~imdilate(hd, strel('sphere',6));

    T  = load(fullfile(pd,'rovir_transform.mat'),'transform');
    roiE = T.transform.roiEnergy_retain(nv);

    [nrG, ccG, snrG, snrRefG] = metrics(xpG, xoG, roi, air);
    [nrS, ccS, snrS, snrRefS] = metrics(xpS, xoS, roi, air);

    rows(end+1,:) = {s, nv, 52/nv, roiE, nrG, ccG, snrG, snrRefG, ...
                     nrS, ccS, snrS, snrRefS, 100*(snrG/snrRefG-1)};  %#ok<SAGROW>
    keep(end+1) = s;                                                   %#ok<SAGROW>
    orbG{end+1} = xoG * 1;                                             %#ok<SAGROW>
    orbR{end+1} = xpG * (mean(xoG(roi))/mean(xpG(roi)));               %#ok<SAGROW>
    fprintf('sub-%03d  gridded: NRMSE %.4f corr %.4f SNR %.1f (ref %.1f)\n', ...
            s, nrG, ccG, snrG, snrRefG);
end

if isempty(rows); error('No finished reconstructions found.'); end

T = cell2table(rows, 'VariableNames', {'subject','nv','compression', ...
    'roi_energy_pct','grid_nrmse','grid_corr','grid_snr','grid_snr_ref', ...
    'steva_nrmse','steva_corr','steva_snr','steva_snr_ref','snr_change_pct'});
writetable(T, fullfile(outDir, sprintf('cohort_results_nv%d.csv', nv)));

fprintf('\n===== COHORT SUMMARY (n = %d subjects, nv = %d, %.1fx) =====\n', ...
        height(T), nv, 52/nv);
f = @(v) sprintf('min %.4f  median %.4f  max %.4f', min(v), median(v), max(v));
fprintf('GRIDDED (quote these):\n');
fprintf('  ROI NRMSE     : %s\n', f(T.grid_nrmse));
fprintf('  ROI corr      : %s\n', f(T.grid_corr));
fprintf('  ROI SNR change: min %+.1f%%  median %+.1f%%  max %+.1f%%\n', ...
        min(T.snr_change_pct), median(T.snr_change_pct), max(T.snr_change_pct));
fprintf('STEVA (confounded by delta scaling - context only):\n');
fprintf('  ROI NRMSE     : %s\n', f(T.steva_nrmse));
fprintf('  ROI corr      : %s\n', f(T.steva_corr));
fprintf('\nWorst subject by gridded corr: sub-%03d (%.4f)\n', ...
        T.subject(T.grid_corr==min(T.grid_corr)), min(T.grid_corr));
fprintf('Table: %s\n', fullfile(outDir, sprintf('cohort_results_nv%d.csv', nv)));

%% Figure 1 - per-subject scatter

fh = figure('Visible','off','Position',[100 100 1200 420]);
tiledlayout(1,3,'Padding','compact');
nexttile; bar(T.subject, T.grid_corr); ylim([0.98 1]); grid on;
xlabel('subject'); ylabel('ROI correlation'); title('Gridded: fidelity vs 52-ch');
nexttile; bar(T.subject, T.snr_change_pct); grid on;
xlabel('subject'); ylabel('ROI SNR change (%)'); title('Gridded: SNR cost');
nexttile; plot(T.roi_energy_pct, T.grid_corr, 'o', 'MarkerFaceColor',[.2 .4 .8]); grid on;
xlabel('ROI energy retained (%)'); ylabel('ROI correlation');
title('Does energy predict fidelity?');
exportgraphics(fh, fullfile(figDir, sprintf('figE_cohort_metrics_nv%d.png', nv)), ...
               'Resolution', 130); close(fh);

%% Figure 2 - orbit montage, reference above ROI-PCA, one column per subject

nS = numel(keep); Nb = size(orbG{1},1);
zc = round(24 * Nb/48); r1 = 1:round(Nb*0.40); r2 = round(Nb*0.24):round(Nb*0.76);
fh = figure('Visible','off','Position',[100 100 130*nS 460]);
tl = tiledlayout(2, nS, 'Padding','compact','TileSpacing','none');
for i = 1:nS
    hi = prctile(orbG{i}(:), 99.8);
    nexttile(i);      imshow(min(orbG{i}(r1,r2,zc)/hi,1), []);
    title(sprintf('%03d', keep(i)), 'FontSize', 8);
    if i==1, ylabel('52-ch'); end
    nexttile(nS+i);   imshow(min(orbR{i}(r1,r2,zc)/hi,1), []);
end
title(tl, {sprintf('Cohort, gridded recon at slice 24 — top: 52 channels, bottom: ROI-PCA %d channels (%.0fx)', nv, 52/nv), ...
    'No regularisation in either row.'}, 'FontWeight','bold');
exportgraphics(fh, fullfile(figDir, sprintf('figF_cohort_orbits_nv%d.png', nv)), ...
               'Resolution', 130); close(fh);
fprintf('Figures: %s\n', figDir);

%% ---------------------------------------------------------------------------

function [nrmse, cc, snr, snrRef] = metrics(x, ref, roi, air)
    x      = x * (mean(ref(roi)) / mean(x(roi)));      % match inside the ROI
    nrmse  = sqrt(mean((x(roi)-ref(roi)).^2)) / mean(ref(roi));
    cc     = corr(double(x(roi)), double(ref(roi)));
    snr    = mean(x(roi))   / std(double(x(air)));
    snrRef = mean(ref(roi)) / std(double(ref(air)));
end

function v = grab(p, f)
    S = load(p, f); v = S.(f);
    if iscell(v); v = v{1}; end
end
