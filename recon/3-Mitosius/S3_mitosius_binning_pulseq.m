clc; clearvars; close all;

addpath(genpath('/home/debi/jaime/repos/MR-EyeTrack/recon'));
addpath(genpath('/home/debi/MatTechLab/monalisa'));
addpath(genpath('/home/debi/yiwei/forclone/pulseq'));

%% Config

% Variables
subject_num = 1;
mask_type = 'clean';   % use char instead of string

Matrix_size = 240;
reconFov = 240;

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
reconDir    = fullfile(subjectDir, 'recon');
binsDir     = fullfile(reconDir, 'bins', mask_type, filesep);

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

% Directory existence check
if ~isfolder(binsDir)
    mkdir(binsDir);
    disp(['Directory created: ', binsDir]);
else
    disp(['Directory already exists: ', binsDir]);
end

%% Step 1: Load the Raw Data

% Sequence parameters
seqFile = seqFolder + "/" + seqName_list{1};  % main sequence
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
voxel_size = round(FoV/240);
% So the mitosius saved on debi
% is the smaller than the full resolution.
% ===============================================
matrix_size = 240;  % Max nominal spatial resolution
N_u = [matrix_size, matrix_size, matrix_size];
dK_u = [1, 1, 1]./FoV;

%% Rotate C for sub-001
if subject_num == 1
    C_rot = rot90(C, 1);  % Rotate 90 degrees counter-clockwise
    disp(size(C_rot));
    % bmImage(C_rot);
    bmImage(sqrt(sum(C_rot.^2, 4)));
    C = C_rot;
end

%% Resize C
C = bmImResize(C, [48, 48, 48], N_u);
% C = flip(flip(flip(C, 1), 2), 3);

%% Step 3: Normalize the Raw Data

if N_u >240
    normalization = false;
else 
    normalization = true;
end
if normalization
    x_tot = bmMathilda(y_tot, t_tot, ve_tot, C, N_u, N_u, dK_u); 
    % x_perm = permute(x_tot, [2,3,1]);
    x0=x_tot;
    bmImage(x0);
    temp_im = getimage(gca);
    bmImage(temp_im); 
    temp_roi = roipoly; 
    normalize_val = mean(temp_im(temp_roi(:))); 
    % The normalize_val is super small, it is 5e-10, very small
    % again 3e-9
    % The value of one complex point is like: -0.0396 - 0.1162i
    disp('normalize_val')
    disp(normalize_val)
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

%% [OPTIONAL] Save x0 recon 

% x0Dir = fullfile(reconDir, ['Sub00', num2str(subject_num)], 'T1_LIBRE_Binning/output/');

% if ~isfolder(x0Dir)
%     % If it doesn't exist, create it
%     mkdir(x0Dir);
%     disp(['Directory created: ', x0Dir]);
% else
%     disp(['Directory already exists: ', x0Dir]);
% end
% x0Path = fullfile(x0Dir, 'x0.mat');

% % Save the x0 to the .mat file
% save(x0Path, 'x0', '-v7.3');
% disp('x0 has been saved here:')
% disp(x0Path)

%% Before running this cell, make sure the Mask is well-prepared.
% Load the masked coil sensitivity 

% if woBinning
% MaskFilePath = [binsDir, 'eMask_woBin.mat'];
% if withBinning

for region_idx = 0:3

    th_ratio = 0.75;
    mDir = fullfile(reconDir, 'mitosius', mask_type, ['mask_', num2str(region_idx)]);
    eMaskFilePath = fullfile(binsDir, sprintf('eMask_th%.2f_region%i.mat', th_ratio, region_idx));

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
    eyeMask(:, 1, :) = [];  % SI
    eyeMask(:, :, 1:reader.acquisitionParams.nShot_off) = [];  % SS
    eyeMask = bmPointReshape(eyeMask); 
        
    % Run the mitosis function and compute volume elements
    [y, t] = bmMitosis(y_tot_norm, t_tot, eyeMask); 
    y = bmPermuteToCol(y); 
    ve  = bmVolumeElement(t, 'voronoi_full_radial3' ); 

    % Save all the resulting datastructures on the disk. You are now ready
    % to run your reconstruction
    bmMitosius_create(mDir, y, t, ve); 
    disp('Mitosius files are saved!')
    disp(mDir)

end
