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

%% Config

% Variables
subject_num = 4;
mask_type = 'no-mo';

% Base paths
baseDir = '/home/debi/jaime/repos/MR-EyeTrack/data/study/data';

subjectStr = sprintf('sub-%03d', subject_num);

subjectDir  = fullfile(baseDir, subjectStr);
rawDir      = fullfile(subjectDir, 'rawdata');
ETDir       = fullfile(subjectDir, 'EyeMasks');
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
seqFolder = '/home/debi/jaime/repos/MR-EyeTrack/data/study/pulseq';
seqName_list = {
    'yj_seq2_t1w_libre_main_TR6.2ms_TE3.6ms_swap1_FA6_RF2_mreye_track_trajPTP_44_1872.seq', ...
    'yj0_seq8_t1w_libre_pre_TR6.2ms_TE3.6ms_swap1_FA4_RF2_rfmod2_trajPTP_nSeg88_nShot89.seq'};

% Sequence parameters
seqFile = seqFolder + "/" + seqName_list{1};  % main sequence
seqParams = extract_seq_params(seqFile);

nShotOff = 14; 
nSeg = seqParams.nseg;

%% Generate the (single) full eMask
% This function will guide you manually select the raw data and ET masks
% to generate the ET-guided binning mask for monalisa recon
% If you'd like to generate 4 bins with 4 masks,
% please enter nBin=4, and select ET masks for 4 times

% winLen: the length of the readout sliding window to determine the preservation.
% th_ratio: the ratio for thresholding the ET mask.
winLen = 7;
th_ratio = 0.9;

% Generate the full eMask (4 x N matrix)
eMask = eyeGenerateBinningWin(rawDir, nShotOff, nSeg, th_ratio, ETDir, winLen, true);

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
