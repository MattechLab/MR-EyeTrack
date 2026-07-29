% Joint gaze-resolved reconstruction (spatial TV + temporal TV across bins).
%
% Runs a single ADMM solve over all 4 gaze bins simultaneously, coupling
% them via a finite-difference penalty in the gaze dimension (temporal TV).
% Static anatomy (brain, fat, sclera) benefits from 4x the data; only the
% eye region is allowed to differ across bins.
%
% Set subject_num and mask_type in the workspace before running, or accept
% the defaults below.
%
% Outputs (one file per gaze bin, variable name 'x' for mat2nii compatibility):
%   <reconDir>/<mask_type>/x_joint/
%     x_joint_regionidx_N_nIter_20_ds_1.000_dt_0.100.mat

clc; clearvars; close all;

addpath(genpath('/home/debi/jaime/repos/MR-EyeTrack/recon'));
addpath(genpath('/home/debi/MatTechLab/monalisa'));
addpath(genpath('/home/debi/yiwei/forclone/pulseq'));

%% Config

if ~exist('subject_num', 'var')
    subject_num = 3;
end
if ~exist('mask_type', 'var')
    mask_type = 'clean';
end

nIter   = 20;

% Delta schedules — one value per iteration; rho = 10*delta always.
% Spatial: ramp from 0.1 to 1 so data fidelity dominates early iterations.
% Temporal: ramp from 0 to 0.1 so bins stabilise independently before
%           cross-bin coupling is switched on.
delta_s = linspace(0.1, 1,   nIter)';
delta_t = linspace(0,   0.1, nIter)';
rho_s   = 10 * delta_s;
rho_t   = 10 * delta_t;
nCGD    = 4;

fprintf('Joint recon: subject=%d, mask=%s, ds=[%.3f→%.3f], dt=[%.3f→%.3f]\n', ...
        subject_num, mask_type, delta_s(1), delta_s(end), delta_t(1), delta_t(end));

baseDir    = '/home/debi/jaime/repos/MR-EyeTrack/data/study';
seqFolder  = fullfile(baseDir, 'pulseq');
seqName_list = { ...
    'yj_seq100_t1w_libre_pre_TR6.2ms_TE3.6ms_swap1_FA6_Traj1_nSeg44_nShot191_Fid0_mreye_2p0.seq', ...
    'yj_seq202_t1w_libre_main_TR8.0ms_TE3.6ms_swap1_FA6_RF2_mreye_track_trajPTP_44_1872_gdsp.seq'};

subjectStr = sprintf('sub-%03d', subject_num);
subjectDir = fullfile(baseDir, subjectStr);
rawDir     = fullfile(subjectDir, 'rawdata');
reconDir   = fullfile(subjectDir, 'recon');

xJointDir  = fullfile(reconDir, mask_type, 'x_joint');
if ~exist(xJointDir, 'dir')
    mkdir(xJointDir);
end

nRegions   = 4;      % 0=up, 1=down, 2=left, 3=right

% Check whether all 4 output files already exist
allExist = true;
for r = 0:nRegions-1
    xPath_r = fullfile(xJointDir, sprintf( ...
        'x_joint_regionidx_%i_nIter_%d_ds_%.3f_dt_%.3f.mat', ...
        r, nIter, delta_s(end), delta_t(end)));
    if ~exist(xPath_r, 'file')
        allExist = false;
        break;
    end
end
if allExist
    disp('All joint recon outputs already exist, skipping.')
    return
end

%% Load sequence parameters

seqFile   = seqFolder + "/" + seqName_list{2};
seqParams = extract_seq_params(seqFile);

autoFlag  = true;
arrayCoilFile = fullfile(rawDir, dir(fullfile(rawDir, '*_HC.dat')).name);
reader = createRawDataReader(arrayCoilFile, autoFlag);
reader.acquisitionParams.nShot_off  = 14;
reader.acquisitionParams.traj_type  = 'pulseq';
reader.acquisitionParams.pulseqTrajFile_name = strcat(seqFile);
if isfield(seqParams, 'nshot')
    reader.acquisitionParams.nShot = seqParams.nshot;
end
if isfield(seqParams, 'nseg')
    reader.acquisitionParams.nSeg = seqParams.nseg;
end

FoV         = reader.acquisitionParams.FoV;
voxel_size  = round(FoV/240);
matrix_size = 240;
N_u         = [matrix_size, matrix_size, matrix_size];
n_u         = N_u;
dK_u        = [1, 1, 1]./FoV;

ve_max = 10*prod(dK_u(:));

%% Load and rotate coil sensitivity maps

CfilePath = fullfile(reconDir, 'C.mat');
load(CfilePath, 'C');
disp(['C loaded from: ', CfilePath]);

rotationMap = containers.Map([2, -1], [3, -1]);
if isKey(rotationMap, subject_num)
    k = rotationMap(subject_num);
    C = rot90(C, k);
end
C = bmImResize(C, [48, 48, 48], N_u);
disp('C resized.')

%% Load all 4 gaze bins

x0Dir = fullfile(reconDir, mask_type, 'x0');
if ~exist(x0Dir, 'dir')
    mkdir(x0Dir);
end

y_all   = cell(nRegions, 1);
ve_all  = cell(nRegions, 1);
Gu_all  = cell(nRegions, 1);
Gut_all = cell(nRegions, 1);
X_init  = cell(nRegions, 1);  % initial estimates from gridding

for r = 0:nRegions-1
    mDir_r = fullfile(reconDir, 'mitosius', mask_type, ['mask_', num2str(r)]);
    fprintf('Loading bin %d from %s\n', r, mDir_r);

    y_r  = bmMitosius_load(mDir_r, 'y');
    t_r  = bmMitosius_load(mDir_r, 't');
    ve_r = bmMitosius_load(mDir_r, 've');

    % Unwrap the single-frame cell returned by bmMitosius_load / bmTraj2SparseMat
    y_all{r+1}  = y_r{1};
    ve_all{r+1} = ve_r{1};

    [Gu_r, Gut_r] = bmTraj2SparseMat(t_r, ve_r, N_u, dK_u);
    Gu_all{r+1}   = Gu_r{1};
    Gut_all{r+1}  = Gut_r{1};

    % Per-bin initial estimate (gridding recon, reused from S4 if available)
    x0Path_r = fullfile(x0Dir, ['x0_regionidx' num2str(r) '.mat']);
    if exist(x0Path_r, 'file')
        tmp = load(x0Path_r, 'x0');
        X_init{r+1} = tmp.x0{1};
        fprintf('  x0 loaded from %s\n', x0Path_r);
    else
        X_init{r+1} = bmMathilda(y_r{1}, t_r{1}, ve_r{1}, C, N_u, n_u, dK_u, [], [], [], []);
        x0 = X_init(r+1);
        save(x0Path_r, 'x0', '-v7.3');
        fprintf('  x0 computed and saved to %s\n', x0Path_r);
    end
end

disp('All bins loaded.')

%% Build stacked initial image X [nPt_u x nRegions]

nPt_u = prod(N_u);
X_stacked = complex(zeros(nPt_u, nRegions, 'single'));
for g = 1:nRegions
    X_stacked(:, g) = single(X_init{g}(:));
end

%% Joint ADMM reconstruction

disp('Starting joint reconstruction...')
X_joint = bmSteva_4fr(X_stacked, [], [], [], [], ...
                      y_all, ve_all, C, Gu_all, Gut_all, n_u, ...
                      delta_s, rho_s, delta_t, rho_t, ...
                      nCGD, ve_max, nIter);

% x = bmSteva(  x0{1}, [], [], y{1}, ve{1}, C, Gu{1}, Gut{1}, n_u, ...
%                                         delta, rho, nCGD, ve_max, ...
%                                         nIter, ...
%                                         bmWitnessInfo('steva_d0p1_r1_nCGD4', witness_ind));                      

%% Save one file per gaze bin (variable 'x' for mat2nii compatibility)

for r = 0:nRegions-1
    x = X_joint{r+1};
    xPath_r = fullfile(xJointDir, sprintf( ...
        'x_joint_regionidx_%i_nIter_%d_ds_%.3f_dt_%.3f.mat', ...
        r, nIter, delta_s(end), delta_t(end)));
    save(xPath_r, 'x', '-v7.3');
    fprintf('Saved bin %d -> %s\n', r, xPath_r);
end

disp('Joint reconstruction complete.')
close all;
