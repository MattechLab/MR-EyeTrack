% Example: rigid motion correction in a Mathilda recon pipeline
%
% Assumptions:
% - You already computed a 4D sequential reconstruction (x_cs)
% - Motion parameters are estimated using SPM
% - Motion is applied directly to k-space data and trajectory
% ============================================================
% The idea is to add two steps  
% 1. Estimate motion params from 4D image
% estimate_rigid_motion_params_spm() function
% 2. Correct t and y using apply_rigid_motion_to_kspace() function
% Here is my example:

%% !!! THIS SCRIPT DOES NOT RUN !!!
% You need to make your own 4D recon and have your own working
% reconstruction script, it's just an illustration of how / where to use
% the two provided functionser
addpath(genpath('/Users/mauroleidi/Desktop/MattechGit/monalisa/'));
addpath(genpath('/Users/mauroleidi/Desktop/MattechGit/pulseSeqYiwei/pulseq'));

pathTo4DImage = '/Users/mauroleidi/Desktop/MattechGit/flexyphy_paper_folder/HPC/motion_Estimation/motion_correction/data/x_cs_tres5s.mat';
PixelDims4DImage = [2 2 2 5];

motionParams = estimate_rigid_motion_params_spm( ...
    pathTo4DImage, ...
    'x_cs', ...
    PixelDims4DImage, ...
    'timeseries');

pathToSiemensRd = '/Volumes/sanDisk/5_december_2025/sub-01/raw/gre/original/meas_MID00330_FID19934_gre_std.dat';
pathToPulseqSeq = '/Volumes/sanDisk/5_december_2025/sub-01/raw/gre/original/gre_seq1_t1w_gre_main_TR4.5ms_TE2.0ms_swap0_FA12_RF0_nSeg34_nShot2440_gsm2_traj_original.seq';

% Coil sense and ROI mask for normalization
coilSenseDir = '/Volumes/sanDisk/5_december_2025/sub-01/derivatives/gre/coilSense';

%% ============================================================
% 2) Build raw data reader
% ============================================================

fprintf('\n=== Loading raw data ===\n');

reader = createRawDataReader(pathToSiemensRd, true);
reader.acquisitionParams.traj_type = 'pulseq';
reader.acquisitionParams.pulseqTrajFile_name = pathToPulseqSeq;

seqParams = extract_seq_params(pathToPulseqSeq);
if isfield(seqParams,'nshot')
    reader.acquisitionParams.nShot = seqParams.nshot;
end
if isfield(seqParams,'nseg')
    reader.acquisitionParams.nSeg = seqParams.nseg;
end

% Load raw data and trajectory
y = reader.readRawData(true, true);
t = bmTraj(reader.acquisitionParams);
ve = bmVolumeElement(t,'voronoi_full_radial3');

%% ============================================================
% 3) Load coil sensitivities and ROI
% ============================================================

load(fullfile(coilSenseDir,'coilSens_lowres.mat'), ...
     'coilSensMap','N_u_lowres');
load(fullfile(coilSenseDir,'roiMask_lowres.mat'), ...
     'roiMask','N_u_lowres');

N_u = [240 240 240];
FoV = 240;
dK_u = [1 1 1] / FoV;

C    = bmImResize(coilSensMap, N_u_lowres, N_u);
mask = bmImResize(single(roiMask), N_u_lowres, N_u) > 0.5;

%% ============================================================
% 4) Normalize raw data (same as standard recon)
% ============================================================

fprintf('=== Normalization ===\n');

x_norm = bmMathilda(y, t, ve, C, N_u, N_u, dK_u);
normalize_val = median(abs(x_norm(mask)));

y = y / normalize_val;

%% ============================================================
% 5) Apply rigid motion to k-space data and trajectory
% ============================================================

fprintf('=== Applying rigid motion to k-space ===\n');

% 5 seconds in between each frame in ms
motionTimeMs = ((0:size(motionParams,1)-1)*5 + 2.5) * 1000;

% get timesteamps of lines 
timestamps = reshape(reader.acquisitionParams.timestamp(:), [reader.acquisitionParams.nSeg, reader.acquisitionParams.nShot]);  % [seg, shot]
timestamps(1,:) = [];   % remove SI → [nSeg-1, nShot]
if reader.acquisitionParams.nShot_off > 0
    timestamps(:,1:reader.acquisitionParams.nShot_off) = [];
end
timestamps = timestamps(:);
% Convert them in milliseconds from the start of the acquiisitions
costTime = 2.5; % Siemens
timestampMs = double(timestamps) * costTime;
timestampMs = timestampMs - min(timestampMs);

[y_moco, t_moco] = apply_rigid_motion_to_kspace( ...
    y, t, timestampMs,motionParams, motionTimeMs );

ve_moco = bmVolumeElement(t_moco,'voronoi_full_radial3');

%% ============================================================
% 6) Final motion-corrected reconstruction (all lines)
% ============================================================

fprintf('=== Motion-corrected Mathilda recon ===\n');

x_moco = bmMathilda( ...
    y_moco, ...
    t_moco, ...
    ve_moco, ...
    C, ...
    N_u, ...
    N_u, ...
    dK_u );

%% ============================================================
% 7) Reference reconstruction (no motion correction)
% ============================================================

fprintf('=== Reference recon (no motion correction) ===\n');

x_nomoco = bmMathilda( ...
    y, ...
    t, ...
    ve, ...
    C, ...
    N_u, ...
    N_u, ...
    dK_u );

%% ============================================================
% 8) Quick visual comparison
% ============================================================
bmImage(cat(2,[x_nomoco,x_moco]))

%% ============================================================
% 9) 4-bin sequential reconstruction comparison
% ============================================================

fprintf('=== 4-bin reconstructions ===\n');

nLines = size(y,3);
nBins  = 4;
binSize = floor(nLines / nBins);

% Build bin masks
binMask = false(nBins, nLines);
for b = 1:nBins
    sIdx = (b-1)*binSize + 1;
    eIdx = min(b*binSize, nLines);
    binMask(b, sIdx:eIdx) = true;
end

%% --- Motion-corrected 4-bin ---
[y_bins_moco, t_bins_moco] = bmMitosis(y_moco, t_moco, binMask);
ve_bins_moco = bmVolumeElement(t_bins_moco,'voronoi_full_radial3');

x_bins_moco = cell(1,nBins);
for b = 1:nBins
    x_bins_moco{b} = bmMathilda( ...
        y_bins_moco{b}, t_bins_moco{b}, ve_bins_moco{b}, ...
        C, N_u, N_u, dK_u);
end

%% --- Non-motion-corrected 4-bin ---
[y_bins_nomoco, t_bins_nomoco] = bmMitosis(y, t, binMask);
ve_bins_nomoco = bmVolumeElement(t_bins_nomoco,'voronoi_full_radial3');

x_bins_nomoco = cell(1,nBins);
for b = 1:nBins
    x_bins_nomoco{b} = bmMathilda( ...
        y_bins_nomoco{b}, t_bins_nomoco{b}, ve_bins_nomoco{b}, ...
        C, N_u, N_u, dK_u);
end