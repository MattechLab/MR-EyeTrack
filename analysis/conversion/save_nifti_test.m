%% Test save_nifti (Yannick) on sub-001 woBin reconstruction
clc; clearvars; close all;

repoRoot = '/home/debi/jaime/repos/MR-EyeTrack';
subID    = 'sub-015';

twixFile = fullfile(repoRoot, 'data/study', subID, 'rawdata', [subID '_T1wLIBRE.dat']);
matFile  = fullfile(repoRoot, 'data/study', subID, 'recon/woBin/x_steva_nIter_20_delta_1.000.mat');
outFile  = fullfile(repoRoot, 'data/study', subID, 'recon/woBin/x_steva_nIter_20_delta_1.000_yannick.nii.gz');

addpath(genpath('/home/debi/MatTechLab/monalisa'));
addpath(genpath(fullfile(repoRoot, 'analysis')));

assert(exist(twixFile, 'file') == 2, 'Twix file not found: %s', twixFile);
assert(exist(matFile,  'file') == 2, 'MAT file not found: %s',  matFile);

data = load(matFile);
vol  = single(abs(data.x));
fprintf('Loaded volume: %s\n', mat2str(size(vol)));

FoV  = 240;            % mm — verify against reader.acquisitionParams.FoV
TR   = 8e-3;           % s  — TR8.0ms from sequence name
dK_u = [1 1 1] / FoV;

save_nifti(twixFile, dK_u, vol, TR, outFile);
fprintf('Saved: %s\n', outFile);
