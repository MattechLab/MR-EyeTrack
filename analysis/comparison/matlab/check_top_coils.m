%% Check Top Coils Across Subjects

clc, clearvars, close all;

N = 20;  % Number of top coils to consider

% Base directory
baseDir = '/home/debi/jaime/repos/MR-EyeTrack/data/study';

% Find subject folders (sub-XXX)
subDirs = dir(fullfile(baseDir, 'sub-*'));
subDirs = subDirs([subDirs.isdir]);

nSubs = numel(subDirs);

% Preallocate (52 coils x N subjects)
all_idx = nan(52, nSubs);
subNames = cell(1, nSubs);

for i = 1:nSubs
    subID = subDirs(i).name;
    subNames{i} = subID;

    matFile = fullfile( ...
        baseDir, subID, ...
        'recon', 'mitosius', 'woBin_comp', ...
        'idx_coilSelection.mat');

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

% Convert to table with subject-labeled columns
idxTable = array2table(all_idx, 'VariableNames', subNames);

% Display
disp(idxTable(1:20, :));

% Export it to csv
outPath = 'analysis/comparison/matlab';
outFile = fullfile(outPath, 'idx_coilSelection.csv');
writetable(idxTable, outFile);
fprintf('CSV written to: %s\n', fullfile(pwd, outFile));

%% Analyze Top Coils Frequency
% Number of top entries to compare

nSubs = size(all_idx, 2);

% Preallocate overlap matrix
overlapMat = nan(nSubs, nSubs);

for i = 1:nSubs
    set_i = all_idx(1:N, i);
    set_i = set_i(~isnan(set_i));   % remove NaNs

    for j = 1:nSubs
        set_j = all_idx(1:N, j);
        set_j = set_j(~isnan(set_j));   % remove NaNs

        % Number of shared values (order-independent)
        overlapMat(i, j) = numel(intersect(set_i, set_j));
    end
end

% Heatmap
figure;
h = heatmap(subNames, subNames, overlapMat);
h.Title = sprintf('Overlap of first %d idx values (order-independent)', N);
h.XLabel = 'Subject';
h.YLabel = 'Subject';
h.Colormap = parula;
h.ColorLimits = [0 N];


%% Coil ranking
[nRanks, nSubs] = size(all_idx);   % nRanks should be 52

% Initialize matrix: rows = coils, columns = subjects
coilRank = nan(52, nSubs);

% Build rank matrix
% coilRank(c, s) = rank position of coil c in subject s
for s = 1:nSubs
    idx = all_idx(:, s);
    idx = idx(:);

    for r = 1:numel(idx)
        if ~isnan(idx(r))
            coilRank(idx(r), s) = r;
        end
    end
end

% Average rank per coil across subjects (ignore NaNs)
meanRank = mean(coilRank, 2, 'omitnan');

% Sort coils by best (lowest) mean rank
[meanRankSorted, coilOrder] = sort(meanRank);

% Top N coils (subject-independent)
topN_coils = coilOrder(1:N);
topN_meanRank = meanRankSorted(1:N);

% Display nicely
topN_table = table((1:N)', topN_coils, topN_meanRank, ...
    'VariableNames', {'Index', 'Coil', 'MeanRank'});

disp(topN_table)

% Save topN_coils to .mat
outMatFile = fullfile(outPath, sprintf('top%d_coils_N15.mat', N));
save(outMatFile, 'topN_coils', 'topN_meanRank');
fprintf('MAT file written to: %s\n', fullfile(pwd, outMatFile));

%% Compare the order of N=5 with N=15
n5 = load('analysis/comparison/matlab/top20_coils_N5.mat');
n15 = load('analysis/comparison/matlab/top20_coils_N15.mat');

% Find the positions of the N=5 top coils in the N=15 list
[isInN15, locInN15] = ismember(n5.topN_coils, n15.topN_coils);

% Table showing N=5 coil, its position in N=15 (if present), and mean rank
compareTable = table((1:numel(n5.topN_coils))', n5.topN_coils(:), locInN15(:), ...
    'VariableNames', {'N5_Order', 'Coil', 'PositionInN15'});

disp('Comparison of N=5 top coils order in N=15 list:');
disp(compareTable);
