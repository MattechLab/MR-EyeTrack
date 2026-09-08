%% R12 - Convert the gridded (x0) volumes to NIfTI for the whole cohort
%
% Produces a matched set for visual comparison of the compression:
%
%   woBin/x0.nii.gz                       52-channel reference (gridded)
%   ROI-PCA/x0_ROI-PCA_<nv>_woBin.nii.gz  compressed, nv = 4, 6, 8
%
% All go through analysis/conversion/mat2nii_twix.m with the MR-EyeTrack
% reorientFcn, so every volume lands in the same world space and overlays
% directly in ITK-SNAP.
%
% Why x0 and not STEVA: the gridded recon has no regularisation, so what you see
% is what the channel reduction actually did. STEVA can hide channel loss as
% smoothing (see the README), which makes it the wrong thing to eyeball when
% judging compression.
%
% Existing .nii.gz are OVERWRITTEN - the ones in the repo predate the current
% converter and were made with different geometry handling, so they are not
% comparable with the ROI-PCA outputs.
%
% The Twix metadata cache (woBin/twix_orientation_metadata.mat) is written on
% first use per subject; only 3 of 15 had it, so the first conversion for each
% of the other 12 reads the raw .dat once (~2 min) and every later call is fast.

clc; clearvars -except subs nvList; close all;

addpath(genpath('/home/debi/jaime/repos/MR-EyeTrack/recon'));
addpath(genpath('/home/debi/jaime/repos/MR-EyeTrack/analysis'));
addpath(genpath('/home/debi/MatTechLab/monalisa'));
addpath(genpath('/home/debi/yiwei/forclone/pulseq'));

if ~exist('subs',   'var'); subs   = 1:15;   end
if ~exist('nvList', 'var'); nvList = [4 6 8]; end

repoRoot = '/home/debi/jaime/repos/MR-EyeTrack';
seqFile  = fullfile(repoRoot, 'data/study/pulseq', ...
    'yj_seq202_t1w_libre_main_TR8.0ms_TE3.6ms_swap1_FA6_RF2_mreye_track_trajPTP_44_1872_gdsp.seq');

% MR-EyeTrack LIBRE convention, raw axes (-A, -R, +S)
reorientFcn = @(v) flip(flip(permute(v, [2 1 3 (4:ndims(v))]), 1), 2);

nOK = 0; nFail = 0; failed = {};

for s = subs
    subjectStr = sprintf('sub-%03d', s);
    subjectDir = fullfile(repoRoot, 'data/study', subjectStr);
    reconDir   = fullfile(subjectDir, 'recon');
    pcaDir     = fullfile(reconDir, 'ROI-PCA');

    twixFile = fullfile(subjectDir, 'rawdata', ...
                        dir(fullfile(subjectDir,'rawdata','*_T1wLIBRE.dat')).name);
    refNifti = fullfile(subjectDir, 'dicom', ...
                        'csTFL_mp-rage_1mm-iso_CP_acc4.6_5_MR', [subjectStr '.nii.gz']);
    twixMeta = fullfile(reconDir, 'woBin', 'twix_orientation_metadata.mat');

    if exist(refNifti, 'file') ~= 2
        fprintf(2, '%s: reference NIfTI missing, skipping subject\n', subjectStr);
        nFail = nFail + 1; failed{end+1} = subjectStr; %#ok<SAGROW>
        continue
    end

    % 52-channel gridded reference, then each compressed channel count
    jobs = { fullfile(reconDir, 'woBin', 'x0.mat'), ...
             fullfile(reconDir, 'woBin', 'x0.nii.gz') };
    for nv = nvList
        jobs(end+1, :) = { ...
            fullfile(pcaDir, sprintf('x0_ROI-PCA_%d_woBin.mat', nv)), ...
            fullfile(pcaDir, sprintf('x0_ROI-PCA_%d_woBin.nii.gz', nv)) }; %#ok<SAGROW>
    end

    for j = 1:size(jobs, 1)
        inMat = jobs{j,1}; outNii = jobs{j,2};
        if exist(inMat, 'file') ~= 2
            fprintf('  %s: missing %s, skipping\n', subjectStr, inMat);
            continue
        end
        try
            if exist(outNii, 'file') == 2; delete(outNii); end   % force refresh
            mat2nii_twix(inMat, twixFile, seqFile, refNifti, outNii, twixMeta, reorientFcn);
            nOK = nOK + 1;
        catch ME
            fprintf(2, '  FAILED %s -> %s\n         %s\n', inMat, outNii, ME.message);
            nFail = nFail + 1; failed{end+1} = outNii; %#ok<SAGROW>
        end
    end
    fprintf('=== %s done ===\n', subjectStr);
end

fprintf('\nConverted %d volumes, %d failures\n', nOK, nFail);
if nFail > 0
    fprintf('Failures:\n'); fprintf('  %s\n', failed{:});
end
