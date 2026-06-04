% =====================================================
% Author: Yiwei Jia
% Date: June 23
% ------------------------------------------------
% This script is used for generate binning mask 
% according to the ET mask, where
% sampling rate of ET mask: 1ms
% sampling rate of readouts: TR=6.2ms
% =====================================================
clearvars; clc;
addpath(genpath('/home/debi/jaime/repos/MR-EyeTrack/recon'));
addpath(genpath('/home/debi/jaime/repos/MR-EyeTrack/recon/Recon_scripts'));
addpath(genpath('/home/debi/yiwei/forclone/pulseq'));
addpath(genpath('/home/debi/MatTechLab/monalisa'));

%% Config

% Variables
subject_num = 15;
mask_type = 'no-mo';

% Base paths
baseDir = '/home/debi/jaime/repos/MR-EyeTrack/data/study';

subjectStr = sprintf('sub-%03d', subject_num);

subjectDir  = fullfile(baseDir, subjectStr);
rawDir      = fullfile(subjectDir, 'rawdata');
ETDir       = fullfile(subjectDir, 'eyemasks');
reconDir    = fullfile(subjectDir, 'recon');
binsDir     = fullfile(reconDir, 'bins', mask_type, filesep);

% Directory existence check
if ~isfolder(binsDir)
    mkdir(binsDir);
    disp(['Directory created: ', binsDir]);
else
    disp(['Directory already exists: ', binsDir]);
end

%% Pulseq
% nShotOff should be aligned with the case of woBinning
seqFolder = fullfile(baseDir, 'pulseq');
seqName_list = {
    'yj_seq100_t1w_libre_pre_TR6.2ms_TE3.6ms_swap1_FA6_Traj1_nSeg44_nShot191_Fid0_mreye_2p0.seq', ...
    'yj_seq202_t1w_libre_main_TR8.0ms_TE3.6ms_swap1_FA6_RF2_mreye_track_trajPTP_44_1872_gdsp.seq'};

% Sequence parameters
seqFile = seqFolder + "/" + seqName_list{2};  % main sequence
seqParams = extract_seq_params(seqFile);

nShotOff = 14; 
nSeg = seqParams.nseg;
nShot = seqParams.nshot;

%% Generate the (single) full eMask
% This function will guide you manually select the raw data and ET masks
% to generate the ET-guided binning mask for monalisa recon
% If you'd like to generate 4 bins with 4 masks,
% please enter nBin=4, and select ET masks for 4 times

% winLen: the length of the readout sliding window to determine the preservation.
% th_ratio: the ratio for thresholding the ET mask.
winLen = 10;
th_ratio = 0.75;

% Generate the full eMask (4 x N matrix)
eMask = eyeGenerateBinningWin(rawDir, nShotOff, nSeg, th_ratio, ETDir, winLen, true, true, mask_type);

% Saving data and Convert to Monalisa format
%--------------------------------------------------------------------------
% Extract the row corresponding to this region
% single_eMask = eMask(region_idx + 1, :);  % +1 because MATLAB is 1-based indexing
single_eMask = eMask(1, :);  % for single eMask generation

% Define the file path for this region
eMaskFilePath = [binsDir, sprintf('eMask_th%.2f_winLen%i.mat', th_ratio, winLen)];

% Save this row into the .mat file (variable name is 'eMaskN')
eMaskN = single_eMask;  % overwrite for saving clarity, or use different var name
save(eMaskFilePath, 'eMaskN');

% Display confirmation
disp('eMask has been saved here:')
disp(eMaskFilePath)

% To be specific:
% To determine if the current readout should be preserved or not,
% we set the sliding window of winLen=10, so that we can check the
% corresponding ET window of length winLen*int(TR) = 10*6.2 = 62 ET points
% if more than th_ratio*winLen*int(TR) = 0.75*62 = 46.5 ET points are true (i.e. located within the criterion region)
    % we will maintain the current readout, i.e. binningMask value = 1
% else
    % discard this readout, i.e. binningMask value = 0
