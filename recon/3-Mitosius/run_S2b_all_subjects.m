% Run S2b_coilSelection_intensity.m for all subjects interactively.
% ROI drawing (roipoly) is required for each subject — run from MATLAB desktop.
%
% Skip logic: if idx_coilSelection_intensity.mat already exists, subject is skipped.

clc; close all;

baseDir     = '/home/debi/jaime/repos/MR-EyeTrack/data/study';
matrix_size = 48;
nSubjects   = 15;

scriptDir = fileparts(mfilename('fullpath'));

for sub = 1:nSubjects

    clc, close all;

    subjectStr = sprintf('sub-%03d', sub);
    idxPath = fullfile(baseDir, subjectStr, 'recon', 'mitosius', 'woBin_comp', ...
                       'idx_coilSelection_intensity.mat');

    if exist(idxPath, 'file')
        fprintf('[skip] Subject %d: already processed.\n', sub);
        continue
    end

    x0Path = fullfile(baseDir, subjectStr, 'recon', 'woBin', ...
                      sprintf('x0_noC_%d.mat', matrix_size));
    if ~exist(x0Path, 'file')
        warning('x0_noC not found for subject %d — run S02 first. Skipping.', sub);
        continue
    end

    fprintf('\n=== Subject %d / %d ===\n', sub, nSubjects);
    subject_num = sub;
    run(fullfile(scriptDir, 'S2b_coilSelection_intensity.m'));

end

disp('All subjects processed.');
