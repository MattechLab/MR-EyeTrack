%% Init
clc; clearvars;

%% === User Input Section ===
% Open file dialogs for sequence and raw data files

% Select sequence file
[seqName, seqFolder] = uigetfile('*.seq', 'Select sequence file');
if seqName == 0
    error('Sequence file selection was cancelled');
end

% Select raw data file
[measureFilename, rawDir] = uigetfile('*.dat', 'Select raw data file');
if measureFilename == 0
    error('Raw data file selection was cancelled');
end

% Create cell array with the selected sequence file
seqName_list = {seqName};

%% === Add paths ===
addpath(genpath('/home/debi/jaime/repos/MR-EyeTrack/recon'));
addpath(genpath('/home/debi/MatTechLab/monalisa'));
addpath(genpath('/home/debi/yiwei/forclone/pulseq'));

%% Initialize the directories and acquire the Coil

% Parameters
saveflag = 0;

% Full path to the measurement file
measureFile = fullfile(rawDir, measureFilename);

% Recon directory (same level as rawDir)
[parentDir, ~] = fileparts(rawDir);
reconDir = fullfile(parentDir, 'recon');

%% Step 1: Load the Raw Data

% Sequence file (already selected by user)
seqFile = fullfile(seqFolder, seqName);
seqParams = extract_seq_params(seqFile);

% Reader
autoFlag = true;  % Disable validation UI
reader = createRawDataReader(measureFile, autoFlag);
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
[nx, ny, nz] = size(x0{1});
numCoils = numel(x0);

sum_of_squares = zeros(nx, ny, nz, 'single');

for coil = 1:numCoils
    sum_of_squares = sum_of_squares + real( x0{coil} .* conj(x0{coil}) );
end

xrms = sqrt(sum_of_squares / numCoils);

bmImage(xrms);

%% Save RMS image

% Ask user if they want to save the image
saveChoice = questdlg('Do you want to save the RMS image?', 'Save Image', 'Yes', 'No', 'No');

if strcmp(saveChoice, 'Yes')
    % Open file dialog to select save location and filename
    [filename, pathname] = uiputfile('*.mat', 'Save RMS image as', 'xrms.mat');
    
    if filename ~= 0
        % Construct full path
        fullPath = fullfile(pathname, filename);
        
        % Create the folder if it doesn't exist
        if ~exist(pathname, 'dir')
            mkdir(pathname);
        end
        
        % Save the image
        save(fullPath, 'xrms', '-v7.3');
        disp('Saved RMS image:');
        disp(fullPath);
    else
        disp('Save cancelled by user');
    end
end
