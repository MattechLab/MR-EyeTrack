clc; close all;
if ~exist('subject_num', 'var'); subject_num = 3;  end
if ~exist('matrix_size', 'var'); matrix_size = 48; end

addpath(genpath('/home/debi/jaime/repos/MR-EyeTrack/recon'));
addpath(genpath('/home/debi/MatTechLab/monalisa'));

baseDir     = '/home/debi/jaime/repos/MR-EyeTrack/data/study';
subjectStr  = sprintf('sub-%03d', subject_num);
reconDir    = fullfile(baseDir, subjectStr, 'recon');

%% Load per-coil x0 (computed without C by S02_chuv_woC.m)
x0Path = fullfile(reconDir, 'woBin', sprintf('x0_noC_%d.mat', matrix_size));
disp(['Loading x0 from: ', x0Path]);

%% Run coil selection
weights_norm = coilSelectionEyesROI_intensity(reconDir, x0Path);
