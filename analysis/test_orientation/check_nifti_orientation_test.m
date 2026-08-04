%% check_nifti_orientation_test.m -- MR-Eye Track
% Test script: compare the DICOM MPRAGE reference against the Twix-converted
% reconstruction for sub-XXX using check_nifti_orientation.
% Edit subjectNum or reconPath to test a different subject or NIfTI variant.

clc; clearvars; close all;

baseDir    = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'data', 'study');
subjectNum = 15;  % Change this to test a different subject
subjectStr = sprintf('sub-%03d', subjectNum);

refPath   = fullfile(baseDir, subjectStr, 'dicom', ...
    'csTFL_mp-rage_1mm-iso_CP_acc4.6_5_MR', [subjectStr '.nii.gz']);

% Switch reconPath to compare different NIfTI variants
reconPath = fullfile(baseDir, subjectStr, 'recon', 'woBin', ...
    'x_steva_nIter_20_delta_1.000_twix.nii.gz');

figTitle = sprintf('MPRAGE vs Twix recon — %s', subjectStr);

check_nifti_orientation(refPath, reconPath, figTitle);


%% --- AudioBOLD orientation check ---

clc; clearvars; close all;

audioBOLDRoot = '/mnt/filer01/MatTechLab/yannick.bovier/AudioBOLD';
subID         = 'sub-1';

refNifti    = fullfile(audioBOLDRoot, subID, 'anat', [subID '_T1w.nii.gz']);
yannickNii  = fullfile(audioBOLDRoot, subID, 'func', [subID '_task-audioPulseq_run-1_bold.nii.gz']);
myNii       = fullfile('/home/debi/Downloads', [subID '_BOLD_sense_recon.nii.gz']);

% Yannick's BOLD vs MPRAGE
check_nifti_orientation(refNifti, yannickNii, sprintf('MPRAGE vs Yannick BOLD — %s', subID));

% My conversion vs MPRAGE
check_nifti_orientation(refNifti, myNii, sprintf('MPRAGE vs my BOLD recon — %s', subID));


%% --- Yiwei 2.0 MR-Eye orientation check ---

clc; clearvars; close all;

filerRoot  = '/mnt/filer01/MatTechLab/yiwei.jia';
datasetDir = fullfile(filerRoot, 'datasets', '0005');
outDir     = '/home/debi/Downloads';

refNifti = fullfile(datasetDir, 'dicom_0005', 'csTFL_mp-rage_1mm-iso_CP_acc4.6_5_MR', 'mprage.nii.gz');
outNii1  = fullfile(outDir, 'yiwei_0005_MID00030_T1w_libre.nii.gz');
outNii2  = fullfile(outDir, 'yiwei_0005_MID00025_T2w_libre.nii.gz');

% Sample 1: T1w LIBRE vs MPRAGE
check_nifti_orientation(refNifti, outNii1, 'MPRAGE vs Yiwei T1w LIBRE (MID00030)');

% Sample 2: T2w LIBRE vs MPRAGE
check_nifti_orientation(refNifti, outNii2, 'MPRAGE vs Yiwei T2w LIBRE (MID00025)');
