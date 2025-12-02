%% Init
clc; clearvars;

%% === Add paths ===
addpath(genpath('/Users/cag/Documents/forclone/Recon_scripts'));
addpath(genpath('/Users/cag/Documents/forclone/pulseq_v15'));
addpath(genpath('/Users/cag/Documents/forclone/monalisa'));

%% Initialize the directories and acquire the Coil

% Parameters
subject_num = 1;
saveflag = 1;

% Pulseq
seqFolder = '/home/debi/jaime/repos/MR-EyeTrack/data/study/pulseq';
seqName_list = {
    'yj_seq2_t1w_libre_main_TR6.2ms_TE3.6ms_swap1_FA6_RF2_mreye_track_trajPTP_44_1872.seq', ...
    'yj0_seq8_t1w_libre_pre_TR6.2ms_TE3.6ms_swap1_FA4_RF2_rfmod2_trajPTP_nSeg88_nShot89.seq'};

% Base directory
baseDir = '/home/debi/jaime/repos/MR-EyeTrack/data/study/data';

% Construct subject folder name (zero-padded to 3 digits)
subjectStr = sprintf('sub-%03d', subject_num);

% Full path to dataset directories
subjectDir  = fullfile(baseDir, subjectStr);
rawDir      = fullfile(subjectDir, 'rawdata');
reconDir = fullfile(subjectDir, 'recon');

% Get the list of meas_MID... files
files = dir(fullfile(rawDir, 'meas_MID*_FID*.dat'));

if numel(files) ~= 3
    warning('Expected 3 files, found %d', numel(files));
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
% =====================================================
% Helper function to extract sequence definitions
% =====================================================
function params = extract_seq_params(seqFile)
    params = struct();
    if ~isfile(seqFile), return; end
    try
        seq = mr.Sequence();
        seq.read(seqFile);
        defs = seq.definitions;
        keysList = keys(defs);
        for i = 1:numel(keysList)
            key = keysList{i};
            val = defs(key);
            cleanKey = regexprep(lower(key), '[^a-z0-9_]', '');
            params.(cleanKey) = val;
        end
    catch ME
        warning('Failed to read seq params from %s: %s', seqFile, ME.message);
    end
end

% Sequence parameters
seqFile = seqFolder + "/" + seqName_list{2};  % prescans
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

% bmImage(x0);

% x0Path = fullfile(reconDir, 'x0_noC.mat');
% if saveflag
%     save(x0Path, 'x0', '-v7.3');
%     disp('Saved:');
%     disp(x0Path);
% end

%% === Root-mean-square combination ===
[nx, ny, nz] = size(x0{1});
numCoils = numel(x0);

sum_of_squares = zeros(nx, ny, nz, 'single');

for coil = 1:numCoils
    sum_of_squares = sum_of_squares + real( x0{coil} .* conj(x0{coil}) );
end

xrms = sqrt(sum_of_squares / numCoils);

bmImage(xrms);

%% Save RMS image

% Create the folder if it doesn't exist
if ~exist(reconDir, 'dir')
    mkdir(reconDir);
end

xrmsPath = fullfile(reconDir, 'xrms_BC.mat');

if saveflag
    save(xrmsPath, 'xrms', '-v7.3');
    disp('Saved RMS image:');
    disp(xrmsPath)
end
