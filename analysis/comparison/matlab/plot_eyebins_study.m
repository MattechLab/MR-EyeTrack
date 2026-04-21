%% Visual comparison with subplots
clc; clearvars; close all;
addpath(genpath('/home/debi/MatTechLab/monalisa'));

%% Config

% Variables
subject_num = 1;
mask_type = {'clean', 'clean_0.75'};
region_idx = 0;  % 0:up 1:down 2:left 3:right 4:center mask

% Base directory
baseDir = '/home/debi/jaime/repos/MR-EyeTrack/data/study';

% Construct subject folder name (zero-padded to 3 digits)
subjectStr = sprintf('sub-%03d', subject_num);

% Full path to dataset directories
datasetDir = fullfile(baseDir, subjectStr, 'rawdata');
reconDir   = fullfile(baseDir, subjectStr, 'recon');
binsDir1    = fullfile(reconDir, 'bins', mask_type{1});
binsDir2    = fullfile(reconDir, 'bins', mask_type{2});
ETDir      = fullfile(baseDir, subjectStr, 'eyemasks');

% Eye bins
mask_1_path = [binsDir1, '/eMask_th0.75_region', num2str(region_idx), '.mat'];
mask_2_path = [binsDir2, '/eMask_th0.75_region', num2str(region_idx), '.mat'];

disp(['mask_1_path: ', mask_1_path]);
disp(['mask_2_path: ', mask_2_path]);

% Load images
m1 = load(mask_1_path, 'eMaskN');
m2 = load(mask_2_path, 'eMaskN');

%% Plot
figure;
plot(m1.eMaskN, '.-');
title(['eMaskN 1D Plot - ', mask_type{1}]);
xlabel('Index');
ylabel('Mask Value');

figure;
plot(m2.eMaskN, '.-');
title(['eMaskN 1D Plot - ', mask_type{2}]);
xlabel('Index');
ylabel('Mask Value');
