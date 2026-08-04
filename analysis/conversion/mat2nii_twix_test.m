%% mat2nii_twix_test.m (MR-Eye Track)
% Test script: convert the sub-XXX woBin reconstruction to NIfTI using
% mat2nii_twix_function.  Edit subjectNum and the path variables below to
% run on a different subject or reconstruction file.

clc; clearvars; close all;

repoRoot = '/home/debi/jaime/repos/MR-EyeTrack';
addpath(genpath(fullfile(repoRoot, 'analysis')));
addpath(genpath(fullfile(repoRoot, 'recon')));
addpath(genpath('/home/debi/MatTechLab/monalisa'));
addpath(genpath('/home/debi/yiwei/forclone/pulseq'));

subjectNum = 4;  % Change this to test a different subject
subjectStr = sprintf('sub-%03d', subjectNum);
subDir     = fullfile(repoRoot, 'data/study', subjectStr);

% Input files
matFile  = fullfile(subDir, 'recon/woBin/x_steva_nIter_20_delta_1.000.mat');
twixFile = fullfile(subDir, 'rawdata', [subjectStr '_T1wLIBRE.dat']);
seqFile  = fullfile(repoRoot, 'data/study/pulseq', ...
    'yj_seq202_t1w_libre_main_TR8.0ms_TE3.6ms_swap1_FA6_RF2_mreye_track_trajPTP_44_1872_gdsp.seq');
refNifti = fullfile(subDir, 'dicom/csTFL_mp-rage_1mm-iso_CP_acc4.6_5_MR', ...
    [subjectStr '.nii.gz']);

% Output NIfTI
outputNiiGz = fullfile(subDir, 'recon/woBin/x_steva_nIter_20_delta_1.000_twix.nii.gz');

% Twix metadata cache — avoids re-reading the full raw .dat on repeated runs
twixMetaFile = fullfile(subDir, 'recon/woBin/twix_orientation_metadata.mat');

% raw=(−A,−R,+S) → permute dims 1↔2, then negate both
reorientFcn = @(v) flip(flip(permute(v, [2 1 3 (4:ndims(v))]), 1), 2);
mat2nii_twix(matFile, twixFile, seqFile, refNifti, outputNiiGz, twixMetaFile, reorientFcn);


%% --- AudioBOLD (Yannick) ---

clc; clearvars; close all;

repoRoot = '/home/debi/jaime/repos/MR-EyeTrack';
addpath(genpath(fullfile(repoRoot, 'analysis')));
addpath(genpath(fullfile(repoRoot, 'recon')));
addpath(genpath('/home/debi/MatTechLab/monalisa'));
addpath(genpath('/home/debi/yiwei/forclone/pulseq'));

% BIDS root: /mnt/filer01/MatTechLab/yannick.bovier/AudioBOLD/
%   MPRAGE NIfTI : <root>/sub-X/anat/sub-X_T1w.nii.gz
%   Reconstruction: sourcedata/sub-X/BOLD_sense_recon.mat
%                   field x = 1×40 cell of 240³ frames → stacked to 4D
%   Twix file    : sourcedata/sub-X/meas_*_MAIN.dat
%   Seq file     : sourcedata/sub-X/AudioBOLD*.seq
% Outputs saved to /home/debi/Downloads/

audioBOLDRoot = '/mnt/filer01/MatTechLab/yannick.bovier/AudioBOLD';
outDir        = '/home/debi/Downloads';
subjects      = [2,3];

for iSub = 1:numel(subjects)
    subID  = sprintf('sub-%d', subjects(iSub));
    subSrc = fullfile(audioBOLDRoot, 'sourcedata', subID);

    matFile  = fullfile(subSrc, 'BOLD_sense_recon.mat');
    refNifti = fullfile(audioBOLDRoot, subID, 'anat', [subID '_T1w.nii.gz']);

    % Find MAIN Twix file (meas_*_MAIN.dat)
    twixList = dir(fullfile(subSrc, '*_MAIN.dat'));
    assert(~isempty(twixList), 'No *_MAIN.dat found for %s', subID);
    twixFile = fullfile(subSrc, twixList(1).name);

    % Find AudioBOLD .seq file (skip the PreScan gre_rad one)
    seqList = dir(fullfile(subSrc, 'AudioBOLD*.seq'));
    assert(~isempty(seqList), 'No AudioBOLD*.seq found for %s', subID);
    seqFile = fullfile(subSrc, seqList(1).name);

    outputNiiGz  = fullfile(outDir, [subID '_BOLD_sense_recon.nii.gz']);
    twixMetaFile = fullfile(subSrc, 'twix_orientation_metadata.mat');

    % raw=(−R,−A,+S) → negate dim1 and dim2
    reorientFcn = @(v) flip(flip(v, 1), 2);
    fprintf('\n=== AudioBOLD %s ===\n', subID);
    mat2nii_twix(matFile, twixFile, seqFile, refNifti, outputNiiGz, twixMetaFile, reorientFcn);
end

%% --- Yiwei 2.0 MR-Eye ---

clc; clearvars; close all;

repoRoot = '/home/debi/jaime/repos/MR-EyeTrack';
addpath(genpath(fullfile(repoRoot, 'analysis')));
addpath(genpath(fullfile(repoRoot, 'recon')));
addpath(genpath('/home/debi/MatTechLab/monalisa'));
addpath(genpath('/home/debi/yiwei/forclone/pulseq'));

filerRoot  = '/mnt/filer01/MatTechLab/yiwei.jia';
datasetDir = fullfile(filerRoot, 'datasets', '0005');
outDir     = '/home/debi/Downloads';

% MPRAGE: shared by both samples
refNifti = fullfile(datasetDir, 'dicom_0005', 'csTFL_mp-rage_1mm-iso_CP_acc4.6_5_MR', 'mprage.nii.gz');

% --- Sample 1: T1w LIBRE (MID00030) ---
matFile1  = fullfile(filerRoot, 'recon_results', '0005', 'MID00030_recon', ...
    'T1_LIBRE_woBinning', 'rovir_ncoil25_ortho', 'output', ...
    'x_20260304_195803', 'x_Nx480_nIter15_delta_1.000.mat');
twixFile1 = fullfile(datasetDir, 'meas_MID00030_FID04261_t1w_libre.dat');
seqFile1  = fullfile(filerRoot, 'datasets', '0001', ...
    'yj_seq104_t1w_libre_main_TR8.0ms_TE3.6ms_swap1_FA6_Traj1_nSeg44_nShot1055_Fid0_mreye_2p0_gdsp.seq');
outNii1   = fullfile(outDir, 'yiwei_0005_MID00030_T1w_libre.nii.gz');
metaFile1 = fullfile(outDir, 'yiwei_0005_MID00030_twix_meta.mat');

% T1w raw=(−S,+R,−A) → cycle dims [2,3,1], negate dim2 and dim3
reorientFcn = @(v) flip(flip(permute(v, [2 3 1 (4:ndims(v))]), 2), 3);
fprintf('\n=== Yiwei Sample 1 — T1w LIBRE (MID00030) ===\n');
mat2nii_twix(matFile1, twixFile1, seqFile1, refNifti, outNii1, metaFile1, reorientFcn);

% --- Sample 2: T2w LIBRE (MID00025) ---
matFile2  = fullfile(filerRoot, 'recon_results', '0005', 'MID00025_recon', ...
    'T1_LIBRE_woBinning', 'rovir_ncoil25_ortho', 'output', ...
    'x_20260318_160957', 'x_Nx480_nIter15_delta_4.000.mat');
twixFile2 = fullfile(datasetDir, 'meas_MID00025_FID04256_t2w_libre.dat');
seqFile2  = fullfile(datasetDir, ...
    'yj_seq608_t2w_libre_main_TR8.0ms_TE3.6ms_swap1_FA6_Traj1_nSeg44_nShot1000_Fid0_t50_crusher_646p92.seq');
outNii2   = fullfile(outDir, 'yiwei_0005_MID00025_T2w_libre.nii.gz');
metaFile2 = fullfile(outDir, 'yiwei_0005_MID00025_twix_meta.mat');

% T2w raw=(−R,−S,−A) → swap dims 2↔3, negate dim 2
%
% The dim-3 flip that used to be here made the volume left-right mirrored.
% Note dim 3 is the left-right axis for this dataset, not dim 1: mat2nii_twix
% writes the reorientFcn output directly and takes the direction cosines from
% the reference, and Yiwei's MPRAGE is stored ('P','I','L'), so the written
% array runs A-P, S-I, L-R.  Removing the flip is what corrects the handedness.
%
% Found by analysis/test_orientation/tissue_check/lr_flip_test.py: registering
% the volume and a mirrored copy to the same-session T1w LIBRE and comparing
% MI/CC.  Visual checks in Mango cannot catch this -- a mirrored brain looks
% entirely plausible, and axis identity and rotation are unchanged by a mirror.
reorientFcn_t2w = @(v) flip(permute(v, [1 3 2 (4:ndims(v))]), 2);
fprintf('\n=== Yiwei Sample 2 — T2w LIBRE (MID00025) ===\n');
mat2nii_twix(matFile2, twixFile2, seqFile2, refNifti, outNii2, metaFile2, reorientFcn_t2w);

