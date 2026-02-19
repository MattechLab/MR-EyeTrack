% =====================================================
% Author: Yiwei Jia
% Date: April 3
% ------------------------------------------------
% This script is used for generate mask 
% to eliminate the readouts in the non-steady state
% nShotOff * nSeg
% =====================================================
clearvars; clc; close all;
addpath(genpath('/home/debi/jaime/repos/MR-EyeTrack/recon'));
addpath(genpath('/home/debi/jaime/repos/MR-EyeTrack/recon/Recon_scripts'));
addpath(genpath('/home/debi/yiwei/forclone/pulseq'));

%% Config

% Variables
subject_num = 7;
mask_type = 'woBin';

% Base paths
baseDir = '/home/debi/jaime/repos/MR-EyeTrack/data/study';

% Pulseq
seqFolder = fullfile(baseDir, 'pulseq');
seqName_list = {
    'yj_seq100_t1w_libre_pre_TR6.2ms_TE3.6ms_swap1_FA6_Traj1_nSeg44_nShot191_Fid0_mreye_2p0.seq', ...
    'yj_seq202_t1w_libre_main_TR8.0ms_TE3.6ms_swap1_FA6_RF2_mreye_track_trajPTP_44_1872_gdsp.seq'};

% Construct subject folder name (zero-padded to 3 digits)
subjectStr = sprintf('sub-%03d', subject_num);

% Full path to dataset directories
subjectDir  = fullfile(baseDir, subjectStr);
rawDir      = fullfile(subjectDir, 'rawdata');
reconDir    = fullfile(subjectDir, 'recon');
binsDir     = fullfile(reconDir, 'bins', mask_type, filesep);

% Create binsDir if it doesn't exist
if ~isfolder(binsDir)
    mkdir(binsDir);
end

% Identify files by pattern
bodyCoilFile  = fullfile(rawDir, dir(fullfile(rawDir, '*_BC.dat')).name);
arrayCoilFile = fullfile(rawDir, dir(fullfile(rawDir, '*_HC.dat')).name);
measureFile   = fullfile(rawDir, dir(fullfile(rawDir, '*_T1wLIBRE.dat')).name);

% Display or use them
disp('Found files:');
disp(bodyCoilFile);
disp(arrayCoilFile);
disp(measureFile);

%% Step 1: Load the Raw Data

% Sequence parameters
seqFile = seqFolder + "/" + seqName_list{2};  % main sequence
seqParams = extract_seq_params(seqFile);

% Reader
autoFlag = true;  % Disable validation UI
reader = createRawDataReader(measureFile, autoFlag);
reader.acquisitionParams.nShot_off = 14;
nShotOff = reader.acquisitionParams.nShot_off;
reader.acquisitionParams.traj_type = 'pulseq';
reader.acquisitionParams.pulseqTrajFile_name = strcat(seqFile);
if isfield(seqParams, 'nshot')
    reader.acquisitionParams.nShot = seqParams.nshot;
    nShot = reader.acquisitionParams.nShot;
end
if isfield(seqParams, 'nseg')
    reader.acquisitionParams.nSeg = seqParams.nseg;
    nSeg = reader.acquisitionParams.nSeg;
end


%% Step 2: Create the eMask
eMask = ones(1, nShot*nSeg);
eMask = (eMask>0);
eMask(1:nShotOff*nSeg) = 0;
% Saving data and Convert to Monalisa format
%--------------------------------------------------------------------------    
eMaskFilePath = [binsDir,'eMask_woBin.mat'];

% Save the CMask to the .mat file
save(eMaskFilePath, 'eMask');
disp('eMask has been saved here:')
disp(eMaskFilePath)
