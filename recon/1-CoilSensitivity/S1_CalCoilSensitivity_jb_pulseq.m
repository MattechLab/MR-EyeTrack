% =====================================================
% Author: Yiwei Jia
% Date: June 05
% ------------------------------------------------
% [Coil sensitivity] -> binning mask eMask -> Mitosius
% Update: this script is derived from Demo script
% by Mauro in Monalisa version Feb.5
% The old script has issue when running mask generation
% With readers, the param setting is more organized
% =====================================================
clc, clearvars;
addpath(genpath('/home/debi/jaime/repos/MR-EyeTrack/recon/Recon_scripts'));
addpath(genpath('/home/debi/MatTechLab/internal_monalisa'));
addpath(genpath('/home/debi/yiwei/forclone/pulseq'));

%% Initialize the directories and acquire the Coil

% Parameters
subject_num = 1;

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

%% Load and Configure Data
% Read data using the library's `createRawDataReader` function
% This readers makes the usage of Siemens and ISMRMRD files equivalent for
% the library

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

% Body Coil Reader
bodyCoilReader = createRawDataReader(bodyCoilFile, true);
bodyCoilReader.acquisitionParams.traj_type = 'pulseq';
bodyCoilReader.acquisitionParams.pulseqTrajFile_name = strcat(seqFile);
if isfield(seqParams, 'nshot')
    bodyCoilReader.acquisitionParams.nShot = seqParams.nshot;
end
if isfield(seqParams, 'nseg')
    bodyCoilReader.acquisitionParams.nSeg = seqParams.nseg;
end
bodyCoilReader.acquisitionParams.nShot_off = 14;

% Head Coil Reader
arrayCoilReader = createRawDataReader(arrayCoilFile, true);
arrayCoilReader.acquisitionParams.traj_type = 'pulseq';
arrayCoilReader.acquisitionParams.pulseqTrajFile_name = strcat(seqFile);
if isfield(seqParams, 'nshot')
    arrayCoilReader.acquisitionParams.nShot = seqParams.nshot;
end
if isfield(seqParams, 'nseg')
    arrayCoilReader.acquisitionParams.nSeg = seqParams.nseg;
end
arrayCoilReader.acquisitionParams.nShot_off = 14;

% Ensure consistency in number of shot-off points
nShotOff = arrayCoilReader.acquisitionParams.nShot_off;

%%
quickCalC = 1; 
%1: quickly calculate C by the function organized by Mauro Leidi
%0: explore the details of C calculation in monalisa, with comments

if quickCalC
    autoFlag=true;
    nIter = 5;
    C = mlComputeCoilSensitivity(bodyCoilReader, arrayCoilReader, [48,48,48], autoFlag, nIter);
else
    % Parameters
    dK_u = [1, 1, 1] ./ arrayCoilReader.acquisitionParams.FoV;   % Cartesian grid spacing
    N_u = [24, 24, 24];             % Adjust this value as needed, low resolution is sufficient
    % Compute Trajectory and Volume Elements
    [y_body, t, ve] = bmCoilSense_nonCart_data(bodyCoilReader, N_u);
    y_surface = bmCoilSense_nonCart_data(arrayCoilReader, N_u);
    
    % Compute the gridding matrices (subscript is a reminder of the result)
    % Gn is from uniform to Non-uniform
    % Gu is from non-uniform to Uniform
    % Gut is Gu transposed
    [Gn, Gu, Gut] = bmTraj2SparseMat(t, ve, N_u, dK_u);
    % Create Mask, we should select the box to mask out the background noise
    mask = bmCoilSense_nonCart_mask_automatic(y_body, Gn, false);
    
    close all;
    % Reference coil sensitivity using the body coils. This is used as 
    % % a reference to estiamte the sensitivity of each head coil
    [y_ref, C_ref] = bmCoilSense_nonCart_ref(y_body, Gn, mask, []); 
    
    % Estimate the coil sensitivity of each surface coil using one body coil
    % image as reference image C_c = (X_c./x_ref)
    C_array_prime = bmCoilSense_nonCart_primary(y_surface, y_ref, C_ref, Gn, ve, mask);
    
    % Do a recon, predending the selected body coil is one channel among the
    % others, and optimize the coil sensitivity estimate by alternating steps
    % Of gradient descent (X,C)
    nIter = 5; 
    [C, x] = bmCoilSense_nonCart_secondary(y_surface, C_array_prime, y_ref, C_ref, Gn, Gu, Gut, ve, nIter, false); 
    close all;
end
%

%% Save C into the folder

% bmImage(C)

saveCDir  = reconDir;
CfileName = 'C.mat';

% Create the folder if it doesn't exist
if ~exist(saveCDir, 'dir')
    mkdir(saveCDir);
end

% Full path to  C file
CfilePath = fullfile(saveCDir, CfileName);

% Save the matrix C to the .mat file
save(CfilePath, 'C');
disp('Coil sensitivity C has been saved here:')
disp(CfilePath)
