% Analyse overlap between motion-contaminated readouts and clean gaze bins.
%
% For each gaze direction (region 0-3), counts how many readouts that are
% assigned to that clean bin also fall in the motion period (0's in the
% no-mo mask).  Run from the MR-EyeTrack root directory.

subject_num = 2;
sub_str     = sprintf('sub-%03d', subject_num);

% Resolve repo root regardless of where the script is called from
script_dir = fileparts(mfilename('fullpath'));
repo_root  = fullfile(script_dir, '..', '..');
base_dir   = fullfile(repo_root, 'data', 'study');

nomo_path   = fullfile(base_dir, sub_str, 'recon', 'bins', 'no-mo', ...
                       'eMask_th0.75_winLen10.mat');
clean_dir   = fullfile(base_dir, sub_str, 'recon', 'bins', 'clean');
region_labels = {'up', 'down', 'left', 'right'};

% Load no-motion mask  (1 = no motion, 0 = motion)
tmp      = load(nomo_path);
fn       = fieldnames(tmp);
nomo     = logical(tmp.(fn{1}));   % 1×N
N        = numel(nomo);

n_motion = sum(~nomo);
fprintf('\n=== %s | N readouts = %d | motion readouts = %d (%.1f%%)\n\n', ...
        sub_str, N, n_motion, 100*n_motion/N);

fprintf('%-8s  %-10s  %-14s  %-14s  %-14s\n', ...
        'Region', 'Clean #', 'Clean+motion', 'Clean+no-mo', 'Motion frac (%)');
fprintf('%s\n', repmat('-', 1, 68));

results = struct();
for r = 0:3
    clean_path = fullfile(clean_dir, sprintf('eMask_th0.75_region%d.mat', r));
    tmp2   = load(clean_path);
    fn2    = fieldnames(tmp2);
    clean  = logical(tmp2.(fn2{1}));   % 1×N

    n_clean         = sum(clean);
    n_clean_motion  = sum(clean & ~nomo);   % assigned to bin BUT in motion
    n_clean_nomo    = sum(clean &  nomo);   % assigned AND motion-free

    frac = 100 * n_clean_motion / max(n_clean, 1);

    fprintf('%-8s  %-10d  %-14d  %-14d  %.2f\n', ...
            region_labels{r+1}, n_clean, n_clean_motion, n_clean_nomo, frac);

    results(r+1).region        = region_labels{r+1};
    results(r+1).n_clean       = n_clean;
    results(r+1).n_motion      = n_clean_motion;
    results(r+1).n_nomo        = n_clean_nomo;
    results(r+1).motion_frac   = frac;
end

% --- Summary figure ---
% figure('Name', sprintf('%s – motion/clean overlap', sub_str), ...
%        'NumberTitle', 'off', 'Color', 'w');

% fracs = [results.motion_frac];
% bar(fracs, 'FaceColor', [0.8 0.2 0.2]);
% set(gca, 'XTickLabel', region_labels, 'XTick', 1:4);
% ylabel('Motion-contaminated readouts in clean bin (%)');
% title(sprintf('%s: fraction of clean bin readouts that carry motion', sub_str));
% ylim([0 max(fracs)*1.2 + 1]);
% grid on;
