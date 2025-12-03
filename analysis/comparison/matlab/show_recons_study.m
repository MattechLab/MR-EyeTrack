clc; clearvars; close all;

addpath(genpath('/home/debi/jaime/repos/MR-EyeTrack/recon'));
addpath(genpath('/home/debi/MatTechLab/monalisa'));
addpath(genpath('/home/debi/yiwei/forclone/pulseq'));

%% Config

% Variables
subject_num = 1;
region_idx = 0;
mask_type = 'clean';   % use char instead of string

Matrix_size = 240;
reconFov = 240;

% Pulseq
seqFolder = '/home/debi/jaime/repos/MR-EyeTrack/data/study/pulseq';
seqName_list = {
    'yj_seq2_t1w_libre_main_TR6.2ms_TE3.6ms_swap1_FA6_RF2_mreye_track_trajPTP_44_1872.seq', ...
    'yj0_seq8_t1w_libre_pre_TR6.2ms_TE3.6ms_swap1_FA4_RF2_rfmod2_trajPTP_nSeg88_nShot89.seq'};

% Base directory
baseDir = '/home/debi/jaime/repos/MR-EyeTrack/data/study/data';

% Construct subject folder name (zero-padded to 3 digits)
subjectStr = sprintf('sub-%03d', subject_num);

% Full path to dataset directories
subjectDir  = fullfile(baseDir, subjectStr);
rawDir      = fullfile(subjectDir, 'rawdata');
reconDir    = fullfile(subjectDir, 'recon');

%% Load and display images

% File paths
file1 = fullfile(reconDir, 'woBin', 'xrms.mat');
file2 = fullfile(reconDir, mask_type, 'x0', sprintf('x0_regionidx%d.mat', region_idx));
file3 = fullfile(reconDir, mask_type, 'x', sprintf('x_steva_regionidx_%d_nIter_20_delta_1.000.mat', region_idx));

% Load files
fprintf('Loading files for subject %d, region index %d...\n', subject_num, region_idx);

data1 = load(file1);
data2 = load(file2);
data3 = load(file3);

% Display images
bmImage(data1.xrms);
title(sprintf('xrms.mat\nSub-%03d', subject_num));

bmImage(data2.x0);
title(sprintf('x0 regionidx%d\nSub-%03d', region_idx, subject_num));

bmImage(data3.x);
title(sprintf('x steva regionidx %d\nSub-%03d', region_idx, subject_num));

fprintf('Display complete.\n');