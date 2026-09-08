%% R6 - Convert the ROVir reconstructions to NIfTI
%
% Uses analysis/conversion/mat2nii_twix.m, the same converter and the same
% MR-EyeTrack reorientFcn as the rest of the study, so the ROVir volumes land
% in the same world space as the existing woBin NIfTIs and can be overlaid on
% the MPRAGE directly.
%
%   reorientFcn (MR-EyeTrack LIBRE, raw axes (-A, -R, +S)):
%       @(v) flip(flip(permute(v, [2 1 3 (4:ndims(v))]), 1), 2)
%
% Reuses woBin/twix_orientation_metadata.mat as the Twix metadata cache so the
% 16 GB raw .dat is not re-read.
%
% Outputs -> data/study/sub-NNN/recon/ROVir/*.nii.gz

clc; clearvars -except subject_num nv binName variant only; close all;

addpath(genpath('/home/debi/jaime/repos/MR-EyeTrack/recon'));
addpath(genpath('/home/debi/jaime/repos/MR-EyeTrack/analysis'));
addpath(genpath('/home/debi/MatTechLab/monalisa'));
addpath(genpath('/home/debi/yiwei/forclone/pulseq'));

if ~exist('subject_num', 'var'); subject_num = 15;      end
% 'ROVir'  = SIR-ranked generalised eigenproblem (suppress the head)
% 'ROI-PCA'= ROI-energy-ranked (lambda -> inf); same solver, max compression.
% Outputs go to recon/<variant>/ so the two never mix.
if ~exist('variant',      'var'); variant     = 'ROVir'; end
if ~exist('nv',          'var'); nv          = 20;      end
if ~exist('binName',     'var'); binName     = 'woBin'; end

repoRoot   = '/home/debi/jaime/repos/MR-EyeTrack';
subjectStr = sprintf('sub-%03d', subject_num);
subjectDir = fullfile(repoRoot, 'data/study', subjectStr);
reconDir   = fullfile(subjectDir, 'recon');
rovirDir   = fullfile(reconDir, variant);

twixFile = fullfile(subjectDir, 'rawdata', dir(fullfile(subjectDir,'rawdata','*_T1wLIBRE.dat')).name);
seqFile  = fullfile(repoRoot, 'data/study/pulseq', ...
    'yj_seq202_t1w_libre_main_TR8.0ms_TE3.6ms_swap1_FA6_RF2_mreye_track_trajPTP_44_1872_gdsp.seq');
refNifti = fullfile(subjectDir, 'dicom/csTFL_mp-rage_1mm-iso_CP_acc4.6_5_MR', [subjectStr '.nii.gz']);
twixMeta = fullfile(reconDir, 'woBin', 'twix_orientation_metadata.mat');

assert(exist(refNifti,'file')==2, 'Reference NIfTI not found: %s', refNifti);

% MR-EyeTrack LIBRE convention (see analysis/conversion/mat2nii_twix.m header)
reorientFcn = @(v) flip(flip(permute(v, [2 1 3 (4:ndims(v))]), 1), 2);

% Convert EVERY reconstruction present, not just one nv. Any x_steva_* / x0_*
% in the ROVir folder is picked up, so a new channel count needs no edit here.
% Set `only` to a substring to restrict (e.g. only = 'ROVir_20').
if ~exist('only', 'var'); only = ''; end

cand = [dir(fullfile(rovirDir, 'x_steva_*.mat')); dir(fullfile(rovirDir, 'x0_*.mat'))];
jobs = cell(0, 2);
for i = 1:numel(cand)
    [~, base] = fileparts(cand(i).name);
    if ~isempty(only) && ~contains(base, only); continue; end
    % x0_rovir_<nv>.mat is R3's own gridded check recon. It duplicates R4's
    % x0_ROVir_<nv>_<bin>.mat and stores its volume under the field
    % "x0_rovir", which mat2nii_twix does not recognise. Skip it.
    if ~isempty(regexp(base, '^x0_rovir_\d+$', 'once')); continue; end
    jobs(end+1, :) = {fullfile(rovirDir, cand(i).name), ...
                      fullfile(rovirDir, [base '.nii.gz'])};   %#ok<SAGROW>
end

% The 52-channel STEVA reference, converted the same way, for side-by-side
refMat = fullfile(reconDir, 'woBin', 'x_steva_nIter_20_delta_1.000.mat');
if exist(refMat, 'file')
    jobs(end+1, :) = {refMat, fullfile(rovirDir, 'REF52_x_steva_nIter_20_delta_1.000.nii.gz')};
end

fprintf('%d volume(s) queued for conversion\n', size(jobs, 1));

for i = 1:size(jobs, 1)
    inMat = jobs{i,1}; outNii = jobs{i,2};
    if ~exist(inMat, 'file')
        fprintf('SKIP (missing): %s\n', inMat); continue
    end
    if exist(outNii, 'file')
        fprintf('SKIP (exists) : %s\n', outNii); continue
    end
    fprintf('\n--- %s\n  -> %s\n', inMat, outNii);
    % Don't let one unconvertible volume abort the whole batch
    try
        mat2nii_twix(inMat, twixFile, seqFile, refNifti, outNii, twixMeta, reorientFcn);
    catch ME
        fprintf(2, 'FAILED  %s\n        %s\n', inMat, ME.message);
    end
end

fprintf('\nDone. NIfTIs in %s\n', rovirDir);
