%% The script is not changed yet... it cannot run on debi

clc; clearvars;
addpath(genpath('/home/debi/yiwei/forclone/Recon_scripts'));

%%
subject_num = 2;
region_idx = 0; % 0:up 1:down 2:left 3:right 4:center mask

datasetDir = ['/home/debi/jaime/repos/MR-EyeTrack/data/pilot/sub-0', num2str(subject_num), '/rawdata'];
reconDir = '/home/debi/jaime/tmp/250613_JB/';
otherDir = [reconDir, '/Sub00', num2str(subject_num),'/T1_LIBRE_Binning/other/'];
mDir = [reconDir, '/Sub00', num2str(subject_num),'/T1_LIBRE_Binning/mitosius/mask_', num2str(region_idx), '/'];
saveCDirList = {strcat('/Sub00',num2str(subject_num),'/T1_LIBRE_Binning/C/')...
    strcat('/Sub00',num2str(subject_num),'/T1_LIBRE_woBinning/C/')};

if subject_num == 1
    bodyCoilFile    = [datasetDir, '/meas_MID2400614_FID182868_BEAT_LIBREon_eye_BC_BC.dat'];
    arrayCoilFile   = [datasetDir, '/meas_MID00615_FID182869_BEAT_LIBREon_eye_HC_BC.dat'];
    measureFile     = [datasetDir, '/meas_MID00605_FID182859_BEAT_LIBREon_eye_(23_09_24).dat'];
elseif subject_num == 2
    bodyCoilFile    = [datasetDir, '/meas_MID00589_FID182843_BEAT_LIBREon_eye_BC_BC.dat'];
    arrayCoilFile   = [datasetDir, '/meas_MID00590_FID182844_BEAT_LIBREon_eye_HC_BC.dat'];
    measureFile     = [datasetDir, '/meas_MID00580_FID182834_BEAT_LIBREon_eye_(23_09_24).dat'];
elseif subject_num == 3
    bodyCoilFile    = [datasetDir, '/meas_MID00563_FID182817_BEAT_LIBREon_eye_BC_BC.dat'];
    arrayCoilFile   = [datasetDir, '/meas_MID00564_FID182818_BEAT_LIBREon_eye_HC_BC.dat'];
    measureFile     = [datasetDir, '/meas_MID00554_FID182808_BEAT_LIBREon_eye_(23_09_24).dat'];
end

%% reader
autoFlag = true;  % Disable validation UI
reader = createRawDataReader(bodyCoilFile, autoFlag);

%% Load mitosius
y   = bmMitosius_load(mDir, 'y'); 
t   = bmMitosius_load(mDir, 't'); 
ve  = bmMitosius_load(mDir, 've'); 

disp('Mitosius has been loaded!')

%% Load Coil Sensitivity Maps
saveCDir = [reconDir,saveCDirList{1}];
CfileName = 'C.mat';
CfilePath = fullfile(saveCDir, CfileName);
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

%%
C = bmImResize(C, [48, 48, 48], N_u);

%%
x0 = cell(nFr, 1);
    for i = 1:nFr
        x0{i} = bmMathilda(y{i}, t{i}, ve{i}, C, N_u, n_u, dK_u, [], [], [], []);
    end
    % isequal(x0_p, x0)
    %
    bmImage(x0);


%% save x0

x0Dir = [reconDir, '/Sub00',num2str(subject_num),'/T1_LIBRE_Binning/output/mask_', num2str(region_idx)];

if ~isfolder(x0Dir)
    % If it doesn't exist, create it
    mkdir(x0Dir);
    disp(['Directory created: ', x0Dir]);
else
    disp(['Directory already exists: ', x0Dir]);
end
x0Path = fullfile(x0Dir, 'x0.mat');
% Save the x0 to the .mat file
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

if isunix
    xDir = [reconDir, '/Sub00', num2str(subject_num), '/T1_LIBRE_Binning/output/th8/'];
else
    xDir = [reconDir, '\Sub00', num2str(subject_num), '\240821_recon\240']; 
end
if ~isfolder(xDir)
    % If it doesn't exist, create it
    mkdir(xDir);
    disp(['Directory created: ', xDir]);
else
    disp(['Directory already exists: ', xDir]);
end

xPath = fullfile(xDir, sprintf('x_steva_regionidx%i_nIter%d_delta_%.3f.mat', region_idx, nIter, delta));

% Save the x0 to the .mat file
save(xPath, 'x');
disp('x has been saved here:')
disp(xPath)

%% .mat to .nii.gz
image = load(xPath);

% Define NIfTI metadata (optional but recommended for completeness)
% You can adjust these properties according to your needs.
nii_hdr = struct;  % Create default NIfTI header
nii_hdr.ImageSize = size(image.x);
nii_hdr.PixelDimensions = [0.5 0.5 0.5];  % Adjust these values if needed

% Write the NIfTI file
% niftiwrite(volume_data, nifti_file, nii_hdr);
nifti_file = fullfile(xDir, sprintf('x_steva_regionidx%i_nIter%d_delta_%.3f', region_idx, nIter, delta));
niftiwrite(image.x, nifti_file);
disp(['Data has been saved as a NIfTI file: ', nifti_file]);

%% Compress to .nii.gz
gzip([nifti_file '.nii']);

% (Optional) remove the uncompressed file
delete([nifti_file '.nii']);
