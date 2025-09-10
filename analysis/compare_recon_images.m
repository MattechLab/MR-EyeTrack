%% Visual comparison with subplots
clc; clearvars;
addpath(genpath('/home/debi/jaime/repos/MR-EyeTrack/results'));
addpath(genpath('/home/debi/MatTechLab/monalisa'));  % for bmImage

%% Config

% Variables
subject_num = 2;
region_idx = 0; % 0:up 1:down 2:left 3:right 4:center mask
mask_type = {'clean', 'clean_0.95'};

% Paths
resultsDir = '/home/debi/jaime/repos/MR-EyeTrack/results';
x_all_path = [resultsDir, '/Sub00', num2str(subject_num), '/T1_LIBRE_woBinning/output/th8/x_steva_nIter_20_delta_1.000.mat'];
x_clean_path = [resultsDir, '/Sub00', num2str(subject_num), '/T1_LIBRE_Binning/output/', mask_type{1}, '/th8/x_steva_regionidx_', num2str(region_idx), '_nIter_20_delta_1.000.mat'];
x_clean_095_path = [resultsDir, '/Sub00', num2str(subject_num), '/T1_LIBRE_Binning/output/', mask_type{2}, '/th8/x_steva_regionidx_', num2str(region_idx), '_nIter_20_delta_1.000.mat'];
disp(['x_all_path: ', x_all_path]);
disp(['x_clean_path: ', x_clean_path]);
disp(['x_clean_095_path: ', x_clean_095_path]);

% Load images
x_all = load(x_all_path, 'x');
x_clean = load(x_clean_path, 'x');
x_clean_095 = load(x_clean_095_path, 'x');

%% Subplots
bmImage(cat(2, [x_all.x, x_clean.x, x_clean_095.x]), 'slice', 128, ... % Adjust slice index as needed
    'titles', {'Reconstruction without Mask', ['Reconstruction with Mask: ', mask_type{1}], ['Reconstruction with Mask: ', mask_type{2}]});
