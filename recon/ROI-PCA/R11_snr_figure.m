%% R11 - Raw ROI SNR per subject and channel count
%
% Absolute SNR, not percentage change: the percentage hides that subjects start
% from very different baselines (52-ch reference SNR ranges 4.3 to 8.7 here), so
% the same -25% means something different for each of them.
%
% Reads cohort_results_nv{4,6,8}.csv written by R10_cohort_analysis.m.
% Outputs -> recon/ROI-PCA/figures/figG_raw_snr.png  + snr_by_channels.csv

clc; clearvars; close all;

outDir = '/home/debi/jaime/repos/MR-EyeTrack/recon/ROI-PCA';
figDir = fullfile(outDir, 'figures');
if ~exist(figDir, 'dir'); mkdir(figDir); end

nvList = [4 6 8];
T = cell(1, numel(nvList));
for i = 1:numel(nvList)
    f = fullfile(outDir, sprintf('cohort_results_nv%d.csv', nvList(i)));
    if ~exist(f, 'file'); error('Missing %s - run R10 at nv=%d first.', f, nvList(i)); end
    T{i} = readtable(f);
end

subs   = T{1}.subject;
snrRef = T{1}.grid_snr_ref;                      % same 52-ch reference for all nv
snr    = [T{1}.grid_snr, T{2}.grid_snr, T{3}.grid_snr];

% The reference must be identical across the three tables; if not, the rows are
% not aligned and the comparison is meaningless.
for i = 2:numel(nvList)
    assert(isequal(T{i}.subject, subs), 'Subject order differs between tables.');
    assert(max(abs(T{i}.grid_snr_ref - snrRef)) < 1e-6, 'Reference SNR differs between tables.');
end

%% Table

out = table(subs, snrRef, snr(:,1), snr(:,2), snr(:,3), ...
    'VariableNames', {'subject','snr_52ch','snr_nv4','snr_nv6','snr_nv8'});
writetable(out, fullfile(outDir, 'snr_by_channels.csv'));

fprintf('%-9s %9s %8s %8s %8s\n','subject','52-ch','nv=8','nv=6','nv=4');
for k = 1:height(out)
    fprintf('sub-%03d   %9.2f %8.2f %8.2f %8.2f\n', ...
        out.subject(k), out.snr_52ch(k), out.snr_nv8(k), out.snr_nv6(k), out.snr_nv4(k));
end
fprintf('%-9s %9.2f %8.2f %8.2f %8.2f   <- median\n', 'MEDIAN', ...
    median(out.snr_52ch), median(out.snr_nv8), median(out.snr_nv6), median(out.snr_nv4));

%% Figure

f = figure('Visible','off','Position',[100 100 1250 520]);
tiledlayout(1, 2, 'Padding','compact');

nexttile;
b = bar(subs, [snrRef, snr(:,3), snr(:,2), snr(:,1)], 'grouped');
b(1).FaceColor = [0.25 0.25 0.25];
b(2).FaceColor = [0.15 0.55 0.80];
b(3).FaceColor = [0.45 0.72 0.35];
b(4).FaceColor = [0.90 0.55 0.20];
xlabel('subject'); ylabel('ROI SNR (gridded, absolute)');
legend({'52 channels','nv=8','nv=6','nv=4'}, 'Location','northwest');
title('Raw ROI SNR per subject and channel count'); grid on;
xticks(subs);

nexttile;
plot([52 8 6 4], [median(snrRef) median(snr(:,3)) median(snr(:,2)) median(snr(:,1))], ...
     'o-', 'LineWidth', 1.8, 'MarkerFaceColor','w', 'MarkerSize', 8); hold on;
for k = 1:numel(subs)
    plot([52 8 6 4], [snrRef(k) snr(k,3) snr(k,2) snr(k,1)], '-', ...
         'Color', [0.6 0.6 0.6 0.45]);
end
set(gca,'XDir','reverse','XScale','log');
xticks([4 6 8 52]); xticklabels({'4','6','8','52'});
xlabel('channels kept'); ylabel('ROI SNR (gridded, absolute)');
title('Each subject (grey) and the median (blue)'); grid on;

exportgraphics(f, fullfile(figDir,'figG_raw_snr.png'), 'Resolution', 130); close(f);
fprintf('\nFigure: %s\nTable : %s\n', fullfile(figDir,'figG_raw_snr.png'), ...
        fullfile(outDir,'snr_by_channels.csv'));
