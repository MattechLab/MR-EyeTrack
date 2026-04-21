clc; clearvars; close all;

%% Config

% Variables
subject_num = 12;
mask_type = 'woBin';   % use char instead of string

% Base directory
baseDir = '/usr/src/app/data/study';

% Construct subject folder name (zero-padded to 3 digits)
subjectStr = sprintf('sub-%03d', subject_num);

% Full path to dataset directories
subjectDir  = fullfile(baseDir, subjectStr);
rawDir      = fullfile(subjectDir, 'rawdata');
reconDir    = fullfile(subjectDir, 'recon');
mDir        = fullfile(reconDir, 'mitosius', mask_type);

% Output paths (x path depends on nIter and delta, defined later)
x0Dir       = fullfile(reconDir, mask_type);
x0Path      = fullfile(x0Dir, 'x0.mat');
if ~exist(x0Dir, 'dir')
    mkdir(x0Dir);
end

% Get the list of meas_MID... files
files = dir(fullfile(rawDir, 'sub-*.dat'));

% Identify files by pattern
arrayCoilFile = fullfile(rawDir, dir(fullfile(rawDir, '*_HC.dat')).name);

% Display or use them
disp('Found files:');
disp(arrayCoilFile);

%% Step 1: Load the Raw Data
% Reader
autoFlag = true;  % Disable validation UI
reader = createRawDataReader(arrayCoilFile, autoFlag);

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
disp('C resized!')

%%
x0 = cell(nFr, 1);
    for i = 1:nFr
        x0{i} = bmMathilda(y{i}, t{i}, ve{i}, C, N_u, n_u, dK_u, [], [], [], []);
    end
    % isequal(x0_p, x0)
    % bmImage(x0);

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

% bmImage(x)

% Save the x to the .mat file
xDir  = fullfile(reconDir, mask_type);
xPath = fullfile(xDir, sprintf('x_steva_nIter_%d_delta_%.3f.mat', nIter, delta));
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
