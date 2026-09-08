clc; clearvars; close all;

addpath(genpath('/home/debi/jaime/repos/MR-EyeTrack/recon'));
addpath(genpath('/home/debi/MatTechLab/monalisa'));
addpath(genpath('/home/debi/yiwei/forclone/pulseq'));

%% Config

% Variables
subject_num = 5;
mask_type = 'woBin';  % use char instead of string

matrix_size = 240;  % Max nominal spatial resolution

% Base directory
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

% Output paths (for x0 woBin)
x0Dir       = fullfile(reconDir, mask_type);
x0Path      = fullfile(x0Dir, 'x0.mat');
if ~exist(x0Dir, 'dir')
    mkdir(x0Dir);
end

% Mitosius output directory
mDir = fullfile(reconDir, 'mitosius', mask_type);
if ~exist(mDir, 'dir')
    mkdir(mDir);
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

% Directory existence check
if ~isfolder(binsDir)
    mkdir(binsDir);
    disp(['Directory created: ', binsDir]);
else
    disp(['Directory already exists: ', binsDir]);
end

%% Step 1: Load the Raw Data

% Sequence parameters
seqFile = seqFolder + "/" + seqName_list{2};  % main sequence
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

% Load the raw data and compute trajectory and volume elements
y_tot = reader.readRawData(true, true);  % Filter nShotOff and SI

%% Twix
% myTwix = bmTwix(measureFile);

%%
t_tot = bmTraj(reader.acquisitionParams);                 % Compute trajectory
ve_tot = bmVolumeElement(t_tot, 'voronoi_full_radial3');  % Volume elements

%% Step 2: Load Coil Sensitivity Maps

% Load the coil sensitivity previously measured
CfileName = 'C.mat';
CfilePath = fullfile(reconDir, CfileName);
load(CfilePath, 'C');  % Load sensitivity maps
disp(['C is loaded from:', CfilePath]);

%% Adjust grid size for coil sensitivity maps

FoV = reader.acquisitionParams.FoV;  % Field of View

% ==============================================
% Warning: due to the memory limit, all the voxel_size set on debi
% is always >= 1 to make sure the matrix size <=240
% voxel_size = round(FoV/240);
% So the mitosius saved on debi
% is the smaller than the full resolution.
% ===============================================
N_u = [matrix_size, matrix_size, matrix_size];
dK_u = [1, 1, 1]./FoV;

%% Rotate C
% Define rotation rules: subject_num → rotation k (rot90)
rotationMap = containers.Map( ...
    [2, -1], ...  % [subject_num, rotation] (-1 is clockwise, +1 is counter-clockwise)
    [3, -1] ...
);

% Apply rotation if subject exists in the map
if isKey(rotationMap, subject_num)
    k = rotationMap(subject_num);
    C = rot90(C, k);
    disp(size(C));
    bmImage(sqrt(sum(C.^2, 4)));
end

%% Resize C
C = bmImResize(C, [48, 48, 48], N_u);
% C = flip(flip(flip(C, 1), 2), 3);

%% Step 3: Normalize the Raw Data

if N_u > 240
    normalization = false;
else 
    normalization = true;
end
if normalization
    x_tot = bmMathilda(y_tot, t_tot, ve_tot, C, N_u, N_u, dK_u, [], [], [], []);
    x0=x_tot;
    bmImage(x0);
    
    % ROI
    % temp_im = getimage(gca);
    % temp_roi = roipoly;
    % normalize_val = mean(temp_im(temp_roi(:)));

    % Automatic ROI selection for normalization
    temp_im = x_tot(...
        round(matrix_size/4):round(matrix_size/4*3), ...
        round(matrix_size/4):round(matrix_size/4*3), ...
        round(matrix_size/2));
    normalize_val = mean(abs(temp_im(:)));
    
    bmImage(temp_im);

    % The normalize_val is super small, it is 5e-10, very small
    % again 3e-9
    % The value of one complex point is like: -0.0396 - 0.1162i
    disp(['normalize_val: ', num2str(normalize_val)])

    y_tot(1,1,123)
end
% only once !!!!
if real(y_tot)<1
    if normalization
        y_tot_norm = y_tot/normalize_val; 
        y_tot_norm(1,1,123)
    else
        y_tot_norm = y_tot/(2.5e-10);
        y_tot_norm(1,1,123)
    end
end

%% Save x0 recon woBin
save(x0Path, 'x0', '-v7.3');
disp('x0 has been saved here:')
disp(x0Path)

%% Prepare eye mask

eMaskFilePath = [binsDir, 'eMask_woBin'];

eyeMask = load(eMaskFilePath); 
fields = fieldnames(eyeMask);  % Get the field names
firstField = fields{1};  % Get the first field name
eyeMask = eyeMask.(firstField);  % Access the first field's value
disp(eMaskFilePath)
disp('is loaded!')
% Eliminate the first segment of all the spokes for accuracies

size_Mask = size(eyeMask);
nbins = size_Mask(1);
eyeMask = reshape(eyeMask, [nbins, reader.acquisitionParams.nSeg, reader.acquisitionParams.nShot]); 
eyeMask(:, 1, :) = [];

eyeMask(:, :, 1:reader.acquisitionParams.nShot_off) = []; 
eyeMask = bmPointReshape(eyeMask);

%% Run the mitosis function and compute volume elements

[y, t] = bmMitosis(y_tot_norm, t_tot, eyeMask); 
y = bmPermuteToCol(y); 
ve  = bmVolumeElement(t, 'voronoi_full_radial3'); 

% Save all the resulting datastructures on the disk. You are now ready
% to run your reconstruction

bmMitosius_create(mDir, y, t, ve); 
disp('Mitosius files are saved!')
disp(mDir)
