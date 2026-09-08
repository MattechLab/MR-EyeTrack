%% R9 - Compare every finished reconstruction across both variants
%
% Scans recon/ROVir/ and recon/ROI-PCA/ for STEVA volumes, measures each in the
% eye ROI against the 52-channel reference, and writes a summary table plus a
% figure. Works with whatever exists - missing channel counts are skipped.
%
% Outputs -> recon/ROI-PCA/figures/  and  recon/ROI-PCA/variant_comparison.csv

clc; clearvars -except subject_num; close all;

addpath(genpath('/home/debi/jaime/repos/MR-EyeTrack/recon'));
addpath(genpath('/home/debi/MatTechLab/monalisa'));

if ~exist('subject_num', 'var'); subject_num = 15; end

repoRoot   = '/home/debi/jaime/repos/MR-EyeTrack';
subjectStr = sprintf('sub-%03d', subject_num);
reconDir   = fullfile(repoRoot, 'data/study', subjectStr, 'recon');
outDir     = fullfile(reconDir, 'ROI-PCA');
figDir     = fullfile(repoRoot, 'recon/ROI-PCA/figures');
if ~exist(figDir, 'dir'); mkdir(figDir); end

%% Reference and masks

xo = grab(fullfile(reconDir, 'woBin', 'x_steva_nIter_20_delta_1.000.mat'), 'x');
xo = abs(single(xo));
Nb = size(xo, 1);

M    = load(fullfile(reconDir, 'ROVir', 'masks.mat'));
N48  = size(M.roiMask, 1);
roi  = bmImResize(single(M.roiMask ), [N48 N48 N48], [Nb Nb Nb]) > 0.5;
head = bmImResize(single(M.headMask), [N48 N48 N48], [Nb Nb Nb]) > 0.5;
air  = ~imdilate(head, strel('sphere', 6));
sigO = std(double(xo(air)));

%% Collect every STEVA volume from both variants

nCh   = 52;
rows  = {}; names = {}; vols = {};
for v = {'ROVir', 'ROI-PCA'}
    vd = fullfile(reconDir, v{1});
    if ~exist(vd, 'dir'); continue; end
    f = dir(fullfile(vd, sprintf('x_steva_%s_*_woBin_nIter_20_delta_1.000.mat', v{1})));
    for i = 1:numel(f)
        tok = regexp(f(i).name, sprintf('x_steva_%s_(\\d+)_', v{1}), 'tokens', 'once');
        if isempty(tok); continue; end
        nv = str2double(tok{1});
        x  = abs(single(grab(fullfile(vd, f(i).name), 'x')));
        x  = x * (mean(xo(roi)) / mean(x(roi)));       % match on the ROI
        sig = std(double(x(air)));
        rows(end+1, :) = {v{1}, nv, nCh/nv, ...
            sqrt(mean((x(roi)-xo(roi)).^2))/mean(xo(roi)), ...
            corr(double(x(roi)), double(xo(roi))), ...
            sig, mean(x(roi))/sig, ...
            mean(x(head & ~roi))/mean(x(roi))};        %#ok<SAGROW>
        names{end+1} = sprintf('%s n=%d', v{1}, nv);   %#ok<SAGROW>
        vols{end+1}  = x;                              %#ok<SAGROW>
    end
end

T = cell2table(rows, 'VariableNames', {'variant','nv','compression', ...
        'roi_nrmse','roi_corr','bg_sigma','roi_snr','head_over_roi'});
T = sortrows(T, {'variant','nv'});
writetable(T, fullfile(outDir, 'variant_comparison.csv'));

fprintf('\n%-16s %4s %7s %10s %9s %10s %9s %9s\n', 'recon','nv','compr', ...
        'ROI NRMSE','ROI corr','bg sigma','ROI SNR','head/ROI');
fprintf('%-16s %4d %6.1fx %10.4f %9.4f %10.3e %9.1f %9.3f\n', ...
        '52-ch reference', nCh, 1, 0, 1, sigO, mean(xo(roi))/sigO, ...
        mean(xo(head & ~roi))/mean(xo(roi)));
for i = 1:height(T)
    fprintf('%-16s %4d %6.1fx %10.4f %9.4f %10.3e %9.1f %9.3f\n', ...
        T.variant{i}, T.nv(i), T.compression(i), T.roi_nrmse(i), T.roi_corr(i), ...
        T.bg_sigma(i), T.roi_snr(i), T.head_over_roi(i));
end
fprintf('\nTable: %s\n', fullfile(outDir, 'variant_comparison.csv'));

%% Figure: orbits, reference on top then every variant

ord  = [find(strcmp(T.variant,'ROI-PCA'))', find(strcmp(T.variant,'ROVir'))'];
show = [{xo}, vols(ord)];
labs = [{'52-ch reference'}, names(ord)];
zs   = round((22:2:26) * Nb/N48);
r1   = 1:round(Nb*0.42); r2 = round(Nb*0.22):round(Nb*0.78);
hi   = prctile(xo(roi), 99.5);
nR   = numel(show);

f = figure('Visible','off','Position',[100 100 1050 260*nR]);
tl = tiledlayout(nR, numel(zs), 'Padding','compact','TileSpacing','compact');
for j = 1:nR
  for i = 1:numel(zs)
    nexttile((j-1)*numel(zs)+i);
    imshow(min(show{j}(r1,r2,zs(i))/hi, 1), []);
    if i==1
        text(-0.05, 0.5, labs{j}, 'Units','normalized','Rotation',90, ...
             'HorizontalAlignment','center','FontWeight','bold','FontSize',9);
    end
    if j==1, title(sprintf('slice %d', zs(i))); end
  end
end
title(tl, {'ROI-PCA vs ROVir vs the 52-channel reference, at the orbits', ...
    'All intensity-matched inside the eye ROI.'}, 'FontWeight','bold');
exportgraphics(f, fullfile(figDir,'figA_variant_orbits.png'), 'Resolution', 130); close(f);

% whole head, to show ROI-PCA does not deface
f = figure('Visible','off','Position',[100 100 1300 200*nR]);
zsW = round((20:3:32) * Nb/N48);
tl = tiledlayout(nR, numel(zsW), 'Padding','compact','TileSpacing','compact');
for j = 1:nR
  for i = 1:numel(zsW)
    nexttile((j-1)*numel(zsW)+i);
    imshow(min(show{j}(:,:,zsW(i))/hi, 1), []);
    if i==1
        text(-0.06, 0.5, labs{j}, 'Units','normalized','Rotation',90, ...
             'HorizontalAlignment','center','FontWeight','bold','FontSize',8);
    end
    if j==1, title(sprintf('slice %d', zsW(i))); end
  end
end
title(tl, {'Whole head - ROI-PCA keeps the head, ROVir discards it', ...
    'Same compressed data either way; only the ranking of the virtual channels differs.'}, ...
    'FontWeight','bold');
exportgraphics(f, fullfile(figDir,'figB_variant_wholehead.png'), 'Resolution', 110); close(f);
fprintf('Figures: %s\n', figDir);

function v = grab(p, f)
    S = load(p, f); v = S.(f);
    if iscell(v); v = v{1}; end
end
