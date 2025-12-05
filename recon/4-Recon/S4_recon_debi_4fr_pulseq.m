clc; clearvars; close all;

addpath(genpath('/home/debi/jaime/repos/MR-EyeTrack/recon'));
addpath(genpath('/home/debi/MatTechLab/monalisa'));
addpath(genpath('/home/debi/yiwei/forclone/pulseq'));

%% Config

% Variables
subject_num = 4;
mask_type = 'clean';   % use char instead of string
region_idx = 1; % 0:up 1:down 2:left 3:right 4:center mask

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
mDir        = fullfile(reconDir, 'mitosius', mask_type, ['mask_', num2str(region_idx)]);

% Output paths (x path depends on nIter and delta, defined later)
x0Dir       = fullfile(reconDir, mask_type, 'x0');
x0Path      = fullfile(x0Dir, ['x0_regionidx' num2str(region_idx) '.mat']);
if ~exist(x0Dir, 'dir')
    mkdir(x0Dir);
end

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

%% Load mitosius
y   = bmMitosius_load(mDir, 'y'); 
t   = bmMitosius_load(mDir, 't'); 
ve  = bmMitosius_load(mDir, 've'); 

disp('Mitosius has been loaded!')

%% Load Coil Sensitivity Maps
CfileName = 'C.mat';
CfilePath = fullfile(reconDir, CfileName);
load(CfilePath, 'C');  % Load sensitivity maps
disp(['C is loaded from:', CfilePath]);

%% compileScript()
nFr     = 1; 
% best achivable resolution is 1/ N_u*dK_u If you have enough coverage
FoV = reader.acquisitionParams.FoV;  % Field of View

% ==============================================
% Warning: due to the memory limit, all the voxel_size set on debi
% is always >= 1 to make sure the matrix size <=240
voxel_size = round(FoV/240);
% So the mitosius saved on debi
% is the smaller than the full resolution.
% ===============================================
matrix_size = 240;  % Max nominal spatial resolution
N_u = [matrix_size, matrix_size, matrix_size]; % Matrix size: Size of the Virtual cartesian grid in the fourier space (regridding)
n_u = N_u; % Image size (output)
dK_u = [1, 1, 1]./FoV; % Spacing of the virtual cartesian grid

%% Rotate C for sub-001
if subject_num == 1
    C_rot = rot90(C, 1);  % Rotate 90 degrees counter-clockwise
    disp(size(C_rot));
    % bmImage(C_rot);
    bmImage(sqrt(sum(C_rot.^2, 4)));
    C = C_rot;
elseif subject_num == 4
    C_rot = rot90(C, -1);  % Rotate 90 degrees clockwise
    disp(size(C_rot));
    % bmImage(C_rot);
    bmImage(sqrt(sum(C_rot.^2, 4)));
    C = C_rot;
end

%% Resize C
C = bmImResize(C, [48, 48, 48], N_u);
disp('C resized!')

%%
x0 = cell(nFr, 1);
    for i = 1:nFr
        x0{i} = bmMathilda(y{i}, t{i}, ve{i}, C, N_u, n_u, dK_u, [], [], [], []);
    end
    % isequal(x0_p, x0)
    %
    bmImage(x0);

%% Save the x0 to the .mat file
save(x0Path, 'x0', '-v7.3');
disp('x0 has been saved here:')
disp(x0Path)

%%
[Gu, Gut] = bmTraj2SparseMat(t, ve, N_u, dK_u);

%% bmSteva
deltaArray = 1;

% nIter = 30; % iterations before stopping
nIter = 20; %20, 30
witness_ind = [];
delta = deltaArray(1);
% delta     = 0.1; %0.01, 0.1, 1
rho       = 10*delta;
nCGD      = 4;
ve_max    = 10*prod(dK_u(:));


x = bmSteva(  x0{1}, [], [], y{1}, ve{1}, C, Gu{1}, Gut{1}, n_u, ...
                                        delta, rho, nCGD, ve_max, ...
                                        nIter, ...
                                        bmWitnessInfo('steva_d0p1_r1_nCGD4', witness_ind));

bmImage(x)

% Save the x to the .mat file
xDir  = fullfile(reconDir, mask_type, 'x');
xPath = fullfile(xDir, sprintf('x_steva_regionidx_%i_nIter_%d_delta_%.3f.mat', region_idx, nIter, delta));
if ~exist(xDir, 'dir')
    mkdir(xDir);
end

save(xPath, 'x');
disp('x has been saved here:')
disp(xPath)

% %% .mat to .nii.gz
% image = load(xPath);

% % Define NIfTI metadata (optional but recommended for completeness) You can
% % adjust these properties according to your needs.
% nii_hdr = struct;  % Create default NIfTI header
% nii_hdr.ImageSize = size(image.x);
% nii_hdr.PixelDimensions = [0.5 0.5 0.5];  % Adjust these values if needed
% 
% % Write the NIfTI file niftiwrite(volume_data, nifti_file, nii_hdr);
% nifti_file = fullfile(xDir, sprintf('x_steva_regionidx%i_nIter%d_delta_%.3f.nii', region_idx, nIter, delta));
% niftiwrite(image.x, nifti_file);
% disp(['Data has been saved as a NIfTI file: ', nifti_file]);
% 
% %% Compress to .nii.gz
% gzip(nifti_file);
% 
% % (Optional) remove the uncompressed file delete(nifti_file);
