% =====================================================
% Author: Yiwei Jia
% Date: June 23
% ------------------------------------------------
% This script is used for generate binning mask 
% according to the ET mask, where
% sampling rate of ET mask: 1ms
% sampling rate of readouts: TR=8ms
% =====================================================
clearvars; clc;
addpath(genpath('/home/debi/jaime/repos/MR-EyeTrack/recon'));
addpath(genpath('/home/debi/jaime/repos/MR-EyeTrack/recon/Recon_scripts'));

%% Config

% Variables
subject_num = 2;
mask_type = "clean";

% Paths
datasetDir = ['/home/debi/jaime/repos/MR-EyeTrack/data/pilot/sub-0', num2str(subject_num), '/'];
reconDir = '/home/debi/jaime/repos/MR-EyeTrack/results';
otherDir = [reconDir, '/Sub00', num2str(subject_num),'/T1_LIBRE_Binning/other/', mask_type, '/'];
ETDir = '/home/debi/jaime/repos/MR-EyeTrack/data/pilot/masks_1206/';

% Check if the directory exists
if ~isfolder(otherDir)
    % If it doesn't exist, create it
    mkdir(otherDir);
    disp(['Directory created: ', otherDir]);
else
    disp(['Directory already exists: ', otherDir]);
end

%%
% nShotOff should be aligned with the case of woBinning 
nShotOff = 14; 
nSeg = 22;

% winLen: the length of the readout sliding window to determine the preservation.
% th_ratio: the ratio for thresholding the ET mask.
winLen = 10;
th_ratio = 3/4;

% To be specific:
% To determine if the current readout should be preserved or not,
% we set the sliding window of winLen=10, so that we can check the
% corresponding ET window of length winLen*int(TR) = 80 ET points
% if more than th_ratio*winLen*int(TR) = 60 ET points are true (i.e. located within the criterion region)
    % we will maintain the current readout, i.e. binningMask value = 1
% else
    % discard this readout, i.e. binningMask value = 0
%% Generate the full eMask
% This function will guide you manually select the raw data and ET masks
% to generate the ET-guided binning mask for monalisa recon
% If you'd like to generate 4 bins with 4 masks,
% please enter nBin=4, and select ET masks for 4 times

% Generate the full eMask (4 x N matrix)
eMask = eyeGenerateBinningWin(datasetDir, nShotOff, nSeg, th_ratio, ETDir, winLen, true);

% Saving data and Convert to Monalisa format
%--------------------------------------------------------------------------
region_idx_list = 0:3;  % 0:up 1:down 2:left 3:right 4:center mask

for region_idx = region_idx_list
    % Extract the row corresponding to this region
    single_eMask = eMask(region_idx + 1, :);  % +1 because MATLAB is 1-based indexing

    % Define the file path for this region
    eMaskFilePath = [otherDir, sprintf('eMask_th%.2f_region%i.mat', th_ratio, region_idx)];

    % Save this row into the .mat file (variable name is 'eMaskN')
    eMaskN = single_eMask;  % overwrite for saving clarity, or use different var name
    save(eMaskFilePath, 'eMaskN');

    % Display confirmation
    disp('eMask has been saved here:')
    disp(eMaskFilePath)
end



