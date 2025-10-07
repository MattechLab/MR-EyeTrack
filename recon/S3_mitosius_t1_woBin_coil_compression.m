%% Init
clc; clearvars;
addpath(genpath('/home/debi/jaime/repos/MR-EyeTrack/recon'));
addpath(genpath('/home/debi/MatTechLab/monalisa'));

%% Config

subject_num = 3;
mask_note = '_woFilt';  % '' or '_woFilt'
coilCompression = 0;  % 1: SVD-based compression; 0: coil selection based on energy
nChCompressed = 20;  % Number of virtual coils after compression

Matrix_size = 240;
reconFov = 240;

datasetDir = ['/home/debi/jaime/repos/MR-EyeTrack/data/pilot/sub-00', num2str(subject_num), '/rawdata/'];
reconDir = '/home/debi/jaime/repos/MR-EyeTrack/results';

x0Dir = [reconDir, '/Sub00', num2str(subject_num), '/T1_LIBRE_woBinning/output/x0', mask_note, '/'];

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

otherDir = [reconDir, strcat('/Sub00',num2str(subject_num),'/T1_LIBRE_woBinning/other/')];

% Check if the directory exists
if ~isfolder(otherDir)
    % If it doesn't exist, create it
    mkdir(otherDir);
    disp(['Directory created: ', otherDir]);
else
    disp(['Directory already exists: ', otherDir]);
end

if ~isfolder(x0Dir)
    % If it doesn't exist, create it
    mkdir(x0Dir);
    disp(['Directory created: ', x0Dir]);
else
    disp(['Directory already exists: ', x0Dir]);
end

saveCDir = [reconDir, strcat('/Sub00',num2str(subject_num),'/T1_LIBRE_woBinning/C/')];

%% Step 1: Load the Raw Data

autoFlag = true;  % Disable validation UI
reader = createRawDataReader(measureFile, autoFlag);
reader.acquisitionParams.traj_type = 'full_radial3_phylotaxis';
if strcmp(mask_note, '_woFilt')
    reader.acquisitionParams.nShot_off = 0;  % 0 for woFilt
    reader.acquisitionParams.selfNav_flag = 0;
    y_tot = reader.readRawData(false, false);  % Filter nShot_off (nShot_off*nSeg) and SI (1st nShot)
elseif strcmp(mask_note, '')
    reader.acquisitionParams.nShot_off = 14;  % 14 for with filtering
    reader.acquisitionParams.selfNav_flag = 1;
    y_tot = reader.readRawData(true, true);  % Filter nShot_off (nShot_off*nSeg) and SI (1st nShot)
end

%% Twix
myTwix = bmTwix(measureFile);

%%
t_tot = bmTraj(reader.acquisitionParams);                 % Compute trajectory
ve_tot = bmVolumeElement(t_tot, 'voronoi_full_radial3');  % Volume elements

%% Step 2: Load Coil Sensitivity Maps
% Load the coil sensitivity previously measured
CfileName = 'C.mat';
CfilePath = fullfile(saveCDir, CfileName);
load(CfilePath, 'C');  % Load sensitivity maps
disp(['C is loaded from:', CfilePath]);

%% Adjust grid size for coil sensitivity maps
FoV = reader.acquisitionParams.FoV;  % Field of View

% ==============================================
% Warning: due to the memory limit, all the voxel_size set on debi
% is always >= 1 to make sure the matrix size <=240
voxel_size = round(FoV/240);
% So the mitosius saved on debi
% is smaller than the full resolution.
% ===============================================
matrix_size = 240;  % Max nominal spatial resolution
N_u = [matrix_size, matrix_size, matrix_size];
n_u = [matrix_size, matrix_size, matrix_size];
dK_u = [1, 1, 1]./FoV;

%%
C_resized = bmImResize(C, [48, 48, 48], N_u);
% C = flip(flip(flip(C, 1), 2), 3);

%% Check Mathilda before compression
x0 = bmMathilda(y_tot, t_tot, ve_tot, C_resized, N_u, N_u, dK_u);
bmImage(x0)

% Save x0_woFilt
x0_Path = fullfile(x0Dir, 'x0.mat');
save(x0_Path, 'x0', '-v7.3');
disp('x0_Filt has been saved here:')
disp(x0_Path)

%% Coil compression
nCh = size(y_tot, 1);
nx = size(y_tot, 2);
ntviews = size(y_tot, 3);

if coilCompression == 1
    % SVD-based coil compression
    % --- Compress k-space data ---
    % Reshape to [nCh x (nx*ntviews)]
    D = reshape(y_tot, nx * ntviews, nCh);
    [U, S, V] = svd(D, 'econ');
    singular_values = diag(S);
    total_variance = sum(singular_values .^ 2);
    explained_variance = singular_values(1:nChCompressed) .^ 2;
    percentage_explained = (explained_variance / total_variance) * 100;
    y_tot_comp = reshape(D * V(:, 1:nChCompressed), nChCompressed, nx, ntviews);
    
    % Save the compressed k-space data
    y_outputDir = [reconDir, strcat('/Sub00',num2str(subject_num), '/T1_LIBRE_woBinning/y_tot_comp')];
    y_tot_comp_path = fullfile(y_outputDir, 'y_tot_comp.mat');
    if ~isfolder(y_outputDir)
        % If it doesn't exist, create it
        mkdir(y_outputDir);
        disp(['Directory created: ', y_outputDir]);
    else
        disp(['Directory already exists: ', outputDir]);
    end
    save(y_tot_comp_path, 'y_tot_comp', '-v7.3');
    disp(['y_tot_comp has been saved here: ', y_tot_comp_path]);

    % --- Compress coil sensitivity maps ---
    % Original C: [nx, ny, nz, nCh]
    [nxC, nyC, nzC, nChCheck] = size(C);
    if nChCheck ~= nCh
        error('Mismatch: C has %d channels but y_tot has %d', nChCheck, nCh);
    end
    
    % Reshape to [nVoxels x nCh]
    Cmat = reshape(C, nxC * nyC * nzC, nCh);     % [nVoxels x 64]
    % Apply compression
    Cmat_comp = Cmat * V(:, 1:nChCompressed);    % [nVoxels x 20] 
    % Reshape back to [nx x ny x nz x nChCompressed]
    C_comp = reshape(Cmat_comp, nxC, nyC, nzC, nChCompressed);

    % Save the compressed C
    C_outputDir = [reconDir, strcat('/Sub00',num2str(subject_num), '/T1_LIBRE_woBinning/C_comp')];
    if ~isfolder(C_outputDir)
        % If it doesn't exist, create it
        mkdir(C_outputDir);
        disp(['Directory created: ', C_outputDir]);
    else
        disp(['Directory already exists: ', C_outputDir]);
    end
    C_comp_path = fullfile(C_outputDir, 'C.mat');
    save(C_comp_path, 'C_comp', '-v7.3');
    disp(['C_comp has been saved here: ', C_comp_path]);

    % Plot the explained variance
    f=figure;
    f.Position = [100 100 500 800];
    plot(percentage_explained, '.-', 'LineWidth', 2, 'Color', 'r', 'MarkerSize', 20)
    xlabel('Virtual Coil nr.')
    ylabel('Explained [%]')
    title(sprintf('Coil Compression, Explanation for nChCompressed = %d is %.2f%%', nChCompressed, sum(percentage_explained(1:nChCompressed))))
    % Save the figure
    figPath = fullfile(outputDir, 'exp_var_20ch.png');
    saveas(f, figPath);
    disp(['Explained variance figure saved at: ', figPath]);

    % Set the folder for mitosius saving
    mDir = [reconDir, strcat('/Sub00', num2str(subject_num)),'/T1_LIBRE_woBinning/mitosius_comp/mask_', mask_note, '/'];
    
    % Set the x0Dir for saving x0
    x0Dir = [reconDir, '/Sub00', num2str(subject_num), '/T1_LIBRE_woBinning/output/x0_comp/'];

elseif coilCompression == 0
    % Coil selection based on energy
    % --- Compute coil energy ---
    % Energy per coil = sum over all samples of |signal|^2
    coilEnergy = squeeze(sum(abs(y_tot).^2, [2 3]));  % size: [nCh, 1]

    % --- Sort coils by energy ---
    [~, idxSorted] = sort(coilEnergy, 'descend');

    % Select indices of top coils
    idxKeep = idxSorted(1:nChCompressed);

    % --- Reduce y_tot and C accordingly ---
    y_tot_comp = y_tot(idxKeep, :, :);
    C_comp     = C(:,:,:, idxKeep);

    % Save the compressed k-space data
    y_outputDir = [reconDir, strcat('/Sub00',num2str(subject_num), '/T1_LIBRE_woBinning/y_tot_sel', mask_note)];
    y_tot_comp_path = fullfile(y_outputDir, 'y_tot_comp.mat');
    if ~isfolder(y_outputDir)
        % If it doesn't exist, create it
        mkdir(y_outputDir);
        disp(['Directory created: ', y_outputDir]);
    else
        disp(['Directory already exists: ', y_outputDir]);
    end
    save(y_tot_comp_path, 'y_tot_comp', '-v7.3');
    disp(['y_tot_comp has been saved here: ', y_tot_comp_path]);

    % Save the compressed C
    C_outputDir = [reconDir, strcat('/Sub00',num2str(subject_num), '/T1_LIBRE_woBinning/C_sel')];
    if ~isfolder(C_outputDir)
        % If it doesn't exist, create it
        mkdir(C_outputDir);
        disp(['Directory created: ', C_outputDir]);
    else
        disp(['Directory already exists: ', C_outputDir]);
    end
    C_comp_path = fullfile(C_outputDir, 'C_woFilt.mat');
    save(C_comp_path, 'C_comp', '-v7.3');
    disp(['C_comp has been saved here: ', C_comp_path]);

    % Set the folder for mitosius saving
    mDir = [reconDir, strcat('/Sub00', num2str(subject_num)),'/T1_LIBRE_woBinning/mitosius_sel/mask', mask_note, '/'];
    
    % Set the x0Dir for saving x0
    x0Dir = [reconDir, '/Sub00', num2str(subject_num), '/T1_LIBRE_woBinning/output/x0_sel/'];

end

%% [Optional] Read the new compressed raw data and coil sensitivity
% load(fullfile(y_outputDir, 'y_tot_comp.mat')); % y_tot_comp
% load(fullfile(C_outputDir, 'C.mat')); % C_comp

%%
C_comp_resized = bmImResize(C_comp, [48, 48, 48], N_u);
disp('C_comp resized!')

%% Check mathilda after compression
x0_comp = bmMathilda(y_tot_comp, t_tot, ve_tot, C_comp_resized, N_u, N_u, dK_u);
bmImage(x0_comp)

%% Step 3: Normalize the Raw Data
if N_u > 240
    normalization = false;
else 
    normalization = true;
end

if normalization
    x0_comp = bmMathilda(y_tot_comp, t_tot, ve_tot, C_comp_resized, N_u, N_u, dK_u); 
    x_perm = permute(x0_comp, [2,3,1]);
    bmImage(x_perm)
    temp_im = getimage(gca);  
    bmImage(temp_im); 
    temp_roi = roipoly; 
    normalize_val = mean(temp_im(temp_roi(:))); 
    % The normalize_val is super small, it is 5e-10, very small
    % again 3e-9
    % The value of one complex point is like: -0.0396 - 0.1162i
    disp('normalize_val')
    disp(normalize_val)
    y_tot_comp(1,1,123)
end

% only once !!!!
if real(y_tot_comp)<1
    if normalization
        y_tot_comp_norm = y_tot_comp/normalize_val; 
        y_tot_comp_norm(1,1,123)
    else
        y_tot_comp_norm = y_tot_comp/(2.5e-10); 
        y_tot_comp_norm(1,1,123)
    end
end

%% Save preliminar reconstructed image (mathilda)
if ~isfolder(x0Dir)
    % If it doesn't exist, create it
    mkdir(x0Dir);
    disp(['Directory created: ', x0Dir]);
else
    disp(['Directory already exists: ', x0Dir]);
end
x0Path = fullfile(x0Dir, ['x0_comp', mask_note, '.mat']);
% Save the x0 to the .mat file
save(x0Path, 'x0_comp', '-v7.3');
disp('x0 has been saved here:')
disp(x0Path)

%% Prepare eye mask

eMaskFilePath = [otherDir, 'eMask_woBin'];

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
% eyeMask(:, 1, :) = [];  % SI
% eyeMask(:, :, 1:reader.acquisitionParams.nShot_off) = [];  % SS
eyeMask = bmPointReshape(eyeMask);

%% Run the mitosis function and compute volume elements

[y, t] = bmMitosis(y_tot_comp_norm, t_tot, eyeMask);
y = bmPermuteToCol(y);
ve  = bmVolumeElement(t, 'voronoi_full_radial3');

% Save all the resulting datastructures on the disk. You are now ready
% to run your reconstruction

bmMitosius_create(mDir, y, t, ve);
disp('Mitosius files are saved!')
disp(mDir)
