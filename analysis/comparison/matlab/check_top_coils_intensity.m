%% Check Top Coils (Intensity-based) Across Subjects + Comparison with Energy-based

clc, clearvars, close all;

N       = 20;
baseDir = '/home/debi/jaime/repos/MR-EyeTrack/data/study';
outPath = fileparts(mfilename('fullpath'));

subDirs = dir(fullfile(baseDir, 'sub-*'));
subDirs = subDirs([subDirs.isdir]);
nSubs   = numel(subDirs);

%% Load intensity-based rankings
all_idx = nan(52, nSubs);
subNames = cell(1, nSubs);

for i = 1:nSubs
    subID      = subDirs(i).name;
    subNames{i} = subID;
    matFile = fullfile(baseDir, subID, 'recon', 'mitosius', 'woBin_comp', ...
                       'idx_coilSelection_intensity.mat');
    if exist(matFile, 'file')
        S = load(matFile, 'idx');
        if numel(S.idx) ~= 52
            warning('Unexpected idx size for %s', subID);
        end
        all_idx(:, i) = S.idx(:);
    else
        warning('File not found for %s', subID);
    end
end

%% Export CSV
idxTable = array2table(all_idx, 'VariableNames', subNames);
disp(idxTable(1:N, :));
outFile = fullfile(outPath, 'idx_coilSelection_intensity.csv');
writetable(idxTable, outFile);
fprintf('CSV written to: %s\n', fullfile(pwd, outFile));

%% Overlap heatmap
overlapMat = nan(nSubs, nSubs);
for i = 1:nSubs
    set_i = all_idx(1:N, i); set_i = set_i(~isnan(set_i));
    for j = 1:nSubs
        set_j = all_idx(1:N, j); set_j = set_j(~isnan(set_j));
        overlapMat(i, j) = numel(intersect(set_i, set_j));
    end
end

figure;
h = heatmap(subNames, subNames, overlapMat);
h.Title    = sprintf('Overlap of top %d coils — Intensity-based', N);
h.XLabel   = 'Subject'; h.YLabel = 'Subject';
h.Colormap = parula; h.ColorLimits = [0 N];

%% Mean rank per coil across subjects
coilRank = nan(52, nSubs);
for s = 1:nSubs
    idx_s = all_idx(:, s);
    for r = 1:numel(idx_s)
        if ~isnan(idx_s(r))
            coilRank(idx_s(r), s) = r;
        end
    end
end
meanRank = mean(coilRank, 2, 'omitnan');
[meanRankSorted, coilOrder] = sort(meanRank);
topN_coils_int    = coilOrder(1:N);
topN_meanRank_int = meanRankSorted(1:N);

topN_table = table((1:N)', topN_coils_int, topN_meanRank_int, ...
    'VariableNames', {'Index', 'Coil', 'MeanRank'});
disp(topN_table);

outMatFile = fullfile(outPath, sprintf('top%d_coils_intensity_N%d.mat', N, nSubs));
save(outMatFile, 'topN_coils_int', 'topN_meanRank_int');
fprintf('MAT file written to: %s\n', fullfile(pwd, outMatFile));

%% Compare intensity vs energy-based
energyCsv = fullfile(outPath, 'idx_coilSelection.csv');
if ~exist(energyCsv, 'file')
    warning('Energy-based CSV not found (%s). Skipping comparison.', energyCsv);
    return
end

all_idx_eng = table2array(readtable(energyCsv));

coilRank_eng = nan(52, nSubs);
for s = 1:nSubs
    idx_s = all_idx_eng(:, s);
    for r = 1:numel(idx_s)
        if ~isnan(idx_s(r))
            coilRank_eng(idx_s(r), s) = r;
        end
    end
end
meanRank_eng = mean(coilRank_eng, 2, 'omitnan');
[~, coilOrder_eng] = sort(meanRank_eng);
topN_coils_eng = coilOrder_eng(1:N);

% Overlap
shared = intersect(topN_coils_int, topN_coils_eng);
fprintf('\nTop-%d overlap between methods: %d / %d coils\n', N, numel(shared), N);
fprintf('Shared coils: '); fprintf('%d ', shared(:)'); fprintf('\n');
fprintf('Intensity-only: '); fprintf('%d ', setdiff(topN_coils_int, topN_coils_eng)'); fprintf('\n');
fprintf('Energy-only:    '); fprintf('%d ', setdiff(topN_coils_eng, topN_coils_int)'); fprintf('\n');

% Scatter plot: mean rank per coil, one method vs the other
figure; hold on;
scatter(meanRank_eng, meanRank, 30, [0.7 0.7 0.7], 'filled');
scatter(meanRank_eng(topN_coils_eng), meanRank(topN_coils_eng), 60, [0.2 0.4 0.8], 'filled', 'DisplayName', sprintf('Top-%d energy', N));
scatter(meanRank_eng(topN_coils_int), meanRank(topN_coils_int), 60, [0.9 0.3 0.2], 'filled', 'DisplayName', sprintf('Top-%d intensity', N));

% Coil labels for union of both top-N sets
for c = union(topN_coils_eng(:)', topN_coils_int(:)')
    text(meanRank_eng(c) + 0.3, meanRank(c) + 0.3, num2str(c), 'FontSize', 7);
end

% Identity line
ax = axis;
lims = [min(ax(1), ax(3)), max(ax(2), ax(4))];
plot(lims, lims, 'k--', 'HandleVisibility', 'off');

xlabel('Mean rank — energy-based');
ylabel('Mean rank — intensity-based');
title(sprintf('Coil mean rank: energy vs intensity (top-%d highlighted)', N));
legend('All coils', sprintf('Top-%d energy', N), sprintf('Top-%d intensity', N), 'Location', 'best');
grid on;
