% Motion correction for sub-005 (study_excluded)
% Adapted from example_script.m
%
% Uses xrms.mat (3D RMS image from woBin reconstruction) for SPM motion
% parameter estimation. Since xrms is a single frame, SPM returns zero
% motion parameters — the moco and no-moco reconstructions will match.
% Replace pathTo4DImage with a multi-frame cell-array .mat to get
% non-trivial correction.
% ============================================================

addpath(genpath('/home/debi/jaime/repos/MR-EyeTrack/recon'));
addpath(genpath('/home/debi/jaime/repos/MR-EyeTrack/analysis/motion_correction'));
addpath(genpath('/home/debi/MatTechLab/monalisa'));
addpath(genpath('/home/debi/yiwei/forclone/pulseq'));

%% === Config ===

subject_num = 5;
baseDir     = '/home/debi/jaime/repos/MR-EyeTrack/data/study_excluded';
studyDir    = '/home/debi/jaime/repos/MR-EyeTrack/data/study';

subjectStr = sprintf('sub-%03d', subject_num);
subjectDir = fullfile(baseDir, subjectStr);
rawDir     = fullfile(subjectDir, 'rawdata');
reconDir   = fullfile(subjectDir, 'recon');

seqFile = fullfile(studyDir, 'pulseq', ...
    'yj_seq202_t1w_libre_main_TR8.0ms_TE3.6ms_swap1_FA6_RF2_mreye_track_trajPTP_44_1872_gdsp.seq');

%% ============================================================
% 1) Estimate motion parameters from xrms
% ============================================================

% xrms is a 240x240x240 single real array (1 mm isotropic).
% estimate_rigid_motion_params_spm wraps it as a single-frame cell,
% so SPM realignment gives zero parameters (reference = itself).
pathTo4DImage    = fullfile(reconDir, 'woBin', 'xrms.mat');
PixelDims4DImage = [1 1 1 1];  % [dx dy dz dt] in mm and s

motionParams = estimate_rigid_motion_params_spm( ...
    pathTo4DImage, ...
    'xrms', ...
    PixelDims4DImage, ...
    'xrms_moco');

%% ============================================================
% 2) Build raw data reader
% ============================================================

fprintf('\n=== Loading raw data ===\n');

measureFile = fullfile(rawDir, dir(fullfile(rawDir, '*_T1wLIBRE.dat')).name);
seqParams   = extract_seq_params(seqFile);

reader = createRawDataReader(measureFile, true);
reader.acquisitionParams.nShot_off = 14;
reader.acquisitionParams.traj_type = 'pulseq';
reader.acquisitionParams.pulseqTrajFile_name = seqFile;
if isfield(seqParams, 'nshot')
    reader.acquisitionParams.nShot = seqParams.nshot;
end
if isfield(seqParams, 'nseg')
    reader.acquisitionParams.nSeg = seqParams.nseg;
end

y  = reader.readRawData(true, true);
t  = bmTraj(reader.acquisitionParams);
ve = bmVolumeElement(t, 'voronoi_full_radial3');

%% ============================================================
% 3) Load coil sensitivities
% ============================================================

fprintf('=== Loading coil maps ===\n');

FoV  = reader.acquisitionParams.FoV;
N_u  = [240 240 240];
dK_u = [1 1 1] / FoV;

load(fullfile(reconDir, 'C.mat'), 'C');
C = bmImResize(C, [48 48 48], N_u);

%% ============================================================
% 4) Normalize raw data
% ============================================================

fprintf('=== Normalization ===\n');

x_norm = bmMathilda(y, t, ve, C, N_u, N_u, dK_u);
centre = x_norm( ...
    round(N_u(1)/4):round(N_u(1)/4*3), ...
    round(N_u(2)/4):round(N_u(2)/4*3), ...
    round(N_u(3)/2));
normalize_val = mean(abs(centre(:)));

y = y / normalize_val;

%% ============================================================
% 5) Apply rigid motion to k-space data and trajectory
% ============================================================

fprintf('=== Applying rigid motion to k-space ===\n');

% Single xrms frame → motionTimeMs has one entry at t = 0 ms.
% apply_rigid_motion_to_kspace extrapolates flat → all lines get
% the same (zero) motion from the single reference frame.
motionTimeMs = zeros(size(motionParams, 1), 1);

% Timestamps per readout line
timestamps = reshape(reader.acquisitionParams.timestamp(:), ...
    [reader.acquisitionParams.nSeg, reader.acquisitionParams.nShot]);
timestamps(1, :) = [];  % remove SI row → [nSeg-1, nShot]
if reader.acquisitionParams.nShot_off > 0
    timestamps(:, 1:reader.acquisitionParams.nShot_off) = [];
end
timestamps = timestamps(:);

costTime    = 2.5;  % Siemens timestamp unit in ms
timestampMs = double(timestamps) * costTime;
timestampMs = timestampMs - min(timestampMs);

[y_moco, t_moco] = apply_rigid_motion_to_kspace( ...
    y, t, timestampMs, motionParams, motionTimeMs);

ve_moco = bmVolumeElement(t_moco, 'voronoi_full_radial3');

%% ============================================================
% 6) Motion-corrected reconstruction
% ============================================================

fprintf('=== Motion-corrected Mathilda recon ===\n');

x_moco = bmMathilda(y_moco, t_moco, ve_moco, C, N_u, N_u, dK_u);

%% ============================================================
% 7) Reference reconstruction (no motion correction)
% ============================================================

fprintf('=== Reference recon (no motion correction) ===\n');

x_nomoco = bmMathilda(y, t, ve, C, N_u, N_u, dK_u);

%% ============================================================
% 8) Visual comparison
% ============================================================

bmImage(cat(2, x_nomoco, x_moco));

%% ============================================================
% 9) Save results
% ============================================================

mocoDir = fullfile(reconDir, 'moco');
if ~exist(mocoDir, 'dir'); mkdir(mocoDir); end

save(fullfile(mocoDir, 'x_moco.mat'),   'x_moco',   '-v7.3');
save(fullfile(mocoDir, 'x_nomoco.mat'), 'x_nomoco', '-v7.3');
disp('Results saved to:');
disp(mocoDir);
