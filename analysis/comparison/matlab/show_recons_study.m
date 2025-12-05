clc; clearvars; close all;

addpath(genpath('/home/debi/jaime/repos/MR-EyeTrack/recon'));
addpath(genpath('/home/debi/MatTechLab/monalisa'));
addpath(genpath('/home/debi/yiwei/forclone/pulseq'));

%% Config

% Variables
subject_num = 4;
region_idx = 0;
mask_type = 'clean';   % use char instead of string

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
% Load files
fprintf('Loading files for subject %d, region index %d...\n', subject_num, region_idx);

file1 = fullfile(reconDir, 'woBin', 'xrms.mat');
data1 = load(file1);
bmImage(data1.xrms);
title(sprintf('xrms.mat\nSub-%03d', subject_num));

file2 = fullfile(reconDir, 'woBin', 'xrms_HC.mat');
data2 = load(file2);
bmImage(data2.xrms);
title(sprintf('xrms HC\nSub-%03d', subject_num));

file3 = fullfile(reconDir, mask_type, 'x0', sprintf('x0_regionidx%d.mat', region_idx));
data3 = load(file3);
bmImage(data3.x0{1});
title(sprintf('x0 regionidx%d\nSub-%03d', region_idx, subject_num));

file4 = fullfile(reconDir, mask_type, 'x', sprintf('x_steva_regionidx_%d_nIter_20_delta_1.000.mat', region_idx));
data4 = load(file4);
bmImage(data4.x);
title(sprintf('x steva regionidx %d\nSub-%03d', region_idx, subject_num));

fprintf('Display complete.\n');