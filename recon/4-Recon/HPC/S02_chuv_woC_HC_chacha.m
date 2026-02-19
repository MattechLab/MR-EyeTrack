%% Init
clc; clearvars;
addpath(genpath('/home/jaime.barrancohernandez/shared_datasets/pulseq'))

%% Initialize the directories and acquire the Coil

% Parameters
subject_num = 10;
saveflag = 1;

% Base directory
baseDir = '/home/jaime.barrancohernandez/mnt/jaime.barranco/MR-EyeTrack/data/study';

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
reconDir = fullfile(subjectDir, 'recon');

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
seqFile = seqFolder + "/" + seqName_list{1};  % pre-scan
seqParams = extract_seq_params(seqFile);

% Reader
autoFlag = true;  % Disable validation UI
reader = createRawDataReader(arrayCoilFile, autoFlag);
reader.acquisitionParams.nShot_off = 14;
reader.acquisitionParams.traj_type = 'pulseq';
reader.acquisitionParams.pulseqTrajFile_name = strcat(seqFile);
if isfield(seqParams, 'nshot')
    reader.acquisitionParams.nShot = seqParams.nshot;
end
if isfield(seqParams, 'nseg')
    reader.acquisitionParams.nSeg = seqParams.nseg;
end
%
% Load the raw data and compute trajectory and volume elements
y_tot = reader.readRawData(true, true);  % Filter nShotOff and SI

%% === Load raw data and trajectory ===
t_tot = bmTraj(reader.acquisitionParams);
ve_tot = bmVolumeElement(t_tot, 'voronoi_full_radial3');

%% === Reconstruction configuration ===
matrix_size = 240;
N_u   = [matrix_size matrix_size matrix_size];
dK_u  = [1 1 1] / 240;

nCh = size(y_tot, 1);
disp(['Number of channels: ', num2str(nCh)]);

%% === Perform reconstruction per coil ===
x0 = cell(nCh, 1);

for iCh = 1:nCh
    x0{iCh} = bmMathilda(y_tot(iCh,:), t_tot, ve_tot, [], N_u, N_u, dK_u, [], [], [], []);
    disp(['Processing channel: ', num2str(iCh), '/', num2str(nCh)]);
end

bmImage(x0);

% x0Path = fullfile(reconDir, 'x0_noC.mat');
% if saveflag
%     save(x0Path, 'x0', '-v7.3');
%     disp('Saved:');
%     disp(x0Path);
% end

%% === Root-mean-square combination ===
% Root mean square across the channels
% Initialize an array to store sum of squared images
[nx, ny, nz] = size(x0{1}); % Get the dimensions (240,240,240)
numCoils = numel(x0); % Number of coils

sum_of_squares = zeros(nx, ny, nz, 'single');

% straightforward
% sum_of_squares = sum_of_squares + abs(x0{coil}).^2;
% eliminate extra square-root step
for coil = 1:numCoils
    sum_of_squares = sum_of_squares + real( x0{coil} .* conj(x0{coil}) );
end

% Compute the root mean square (RMS)
xrms = sqrt(sum_of_squares / numCoils);

bmImage(xrms);

%% Save RMS image

% Create the folder if it doesn't exist
if ~exist(reconDir, 'dir')
    mkdir(reconDir);
end

xrmsPath = fullfile(reconDir, 'woBin', 'xrms_HC.mat');
if ~exist(fullfile(reconDir, 'woBin'), 'dir')
    mkdir(fullfile(reconDir, 'woBin'));
end

if saveflag
    save(xrmsPath, 'xrms', '-v7.3');
    disp('Saved RMS image:');
    disp(xrmsPath)
end
