%% mat2nii_twix_test.m
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

mat2nii_twix(matFile, twixFile, seqFile, refNifti, outputNiiGz, twixMetaFile);


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

    fprintf('\n=== AudioBOLD %s ===\n', subID);
    mat2nii_twix(matFile, twixFile, seqFile, refNifti, outputNiiGz, twixMetaFile);
end
