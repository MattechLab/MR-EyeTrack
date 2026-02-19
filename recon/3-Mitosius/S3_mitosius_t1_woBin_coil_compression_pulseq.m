%% Init
clc; clearvars; close all;

addpath(genpath('/home/debi/jaime/repos/MR-EyeTrack/recon'));
addpath(genpath('/home/debi/MatTechLab/monalisa'));
addpath(genpath('/home/debi/yiwei/forclone/pulseq'));
addpath(genpath('/home/debi/jaime/repos/MR-EyeTrack/analysis/comparison/matlab'));

% Config

% Variables
subject_num = 5;
nChCompressed = 20;  % Number of virtual coils after compression
coilCompression = 3;  % 0: SVD-based compression; 1: coil selection based on energy; 2: coil selection with eye ROIs; 3: topN coils from analysis

mask_type = sprintf('woBin_comp/woBin_comp_%d', nChCompressed);   % use char instead of string

Matrix_size = 240;
reconFov = 240;

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
binsDir     = fullfile(reconDir, 'bins/woBin', filesep);

% Output paths (for x0 woBin)
x0Dir       = fullfile(reconDir, 'woBin', filesep);
if ~exist(x0Dir, 'dir')
    mkdir(x0Dir);
end
x0Path      = fullfile(x0Dir, 'x0.mat');

% Mitosius output directory
mDir = fullfile(reconDir, 'mitosius', mask_type);
if ~exist(mDir, 'dir')
    mkdir(mDir);
end

% Get the list of meas_MID... files
files = dir(fullfile(rawDir, 'sub-*.dat'));
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
myTwix = bmTwix(measureFile);

%% Trajectory and volume elements
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
C_resized = bmImResize(C, [48, 48, 48], N_u);
% C = flip(flip(flip(C, 1), 2), 3);
disp('C has been resized!')

%% Check Mathilda before compression
x0 = bmMathilda(y_tot, t_tot, ve_tot, C_resized, N_u, N_u, dK_u);
bmImage(x0)

% Save x0 recon woBin
save(x0Path, 'x0', '-v7.3');
disp('x0 has been saved here:')
disp(x0Path)


%% Coil compression

if coilCompression == 0

    % SVD-based coil compression
    nCh = size(y_tot, 1);
    nx = size(y_tot, 2);
    ntviews = size(y_tot, 3);

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
    y_outputDir = mDir;
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
    C_outputDir = mDir;
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
    figPath = fullfile(mDir, 'exp_var_20ch.png');
    saveas(f, figPath);
    disp(['Explained variance figure saved at: ', figPath]);

elseif coilCompression == 1

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

elseif coilCompression == 2

    % Coil selection with eye ROIs
    weights_norm = coilSelectionEyesROI(reconDir, C);

    % Sort coils by weights
    [sortedW, idx] = sort(weights_norm, 'descend');
   
    % Select indices of top coils
    idxKeep = idx(1:nChCompressed);

    % --- Reduce y_tot and C accordingly ---
    y_tot_comp = y_tot(idxKeep, :, :);
    C_comp     = C(:,:,:, idxKeep);

    % Save the compressed k-space data
    % y_outputDir = mDir;
    % y_tot_comp_path = fullfile(y_outputDir, sprintf('y_tot_comp_%d.mat', nChCompressed));
    % if ~isfolder(y_outputDir)
    %     % If it doesn't exist, create it
    %     mkdir(y_outputDir);
    %     disp(['Directory created: ', y_outputDir]);
    % else
    %     disp(['Directory already exists: ', y_outputDir]);
    % end
    % save(y_tot_comp_path, 'y_tot_comp', '-v7.3');
    % disp(['y_tot_comp has been saved here: ', y_tot_comp_path]);

    % Save the compressed C
    C_outputDir = reconDir;
    C_comp_path = fullfile(C_outputDir, sprintf('C_comp_%d.mat', nChCompressed));
    save(C_comp_path, 'C_comp', '-v7.3');
    disp(['C_comp has been saved here: ', C_comp_path]);

elseif coilCompression == 3

    % Load topN_coils
    topN_coils_path = fullfile('/home/debi/jaime/repos/MR-EyeTrack/analysis/comparison/matlab', sprintf('top%d_coils.mat', nChCompressed));
    load(topN_coils_path, 'topN_coils');
    idxKeep = topN_coils;

    % --- Reduce y_tot and C accordingly ---
    y_tot_comp = y_tot(idxKeep, :, :);
    C_comp     = C(:,:,:, idxKeep);

    % Save the compressed k-space data
    % y_outputDir = mDir;
    % y_tot_comp_path = fullfile(y_outputDir, sprintf('y_tot_comp_%d.mat', nChCompressed));
    % if ~isfolder(y_outputDir)
    %     % If it doesn't exist, create it
    %     mkdir(y_outputDir);
    %     disp(['Directory created: ', y_outputDir]);
    % else
    %     disp(['Directory already exists: ', y_outputDir]);
    % end
    % save(y_tot_comp_path, 'y_tot_comp', '-v7.3');
    % disp(['y_tot_comp has been saved here: ', y_tot_comp_path]);

    % Save the compressed C
    C_outputDir = reconDir;
    C_comp_path = fullfile(C_outputDir, sprintf('C_comp_%d.mat', nChCompressed));
    save(C_comp_path, 'C_comp', '-v7.3');
    disp(['C_comp has been saved here: ', C_comp_path]);

end

%% [Optional] Read the new compressed raw data and coil sensitivity
% load(fullfile(y_outputDir, sprintf('y_tot_comp_%d.mat', nChCompressed))); % y_tot_comp
% load(fullfile(C_outputDir, sprintf('C_comp_%d.mat', nChCompressed))); % C_comp

%% Resize C
C_comp_resized = bmImResize(C_comp, [48, 48, 48], N_u);
disp('C_comp has been resized!')

%% Step 3: Normalize the Raw Data
if N_u > 240
    normalization = false;
else 
    normalization = true;
end

if normalization
    x0_comp = bmMathilda(y_tot_comp, t_tot, ve_tot, C_comp_resized, N_u, N_u, dK_u);
    bmImage(x0_comp)
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
x0CompPath = fullfile(x0Dir, sprintf('x0_comp_%d.mat', nChCompressed));
% Save the x0 to the .mat file
save(x0CompPath, 'x0_comp', '-v7.3');
disp('x0 has been saved here:')
disp(x0CompPath)

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

[y, t] = bmMitosis(y_tot_comp_norm, t_tot, eyeMask);
y = bmPermuteToCol(y);
ve  = bmVolumeElement(t, 'voronoi_full_radial3');

% Save all the resulting datastructures on the disk. You are now ready
% to run your reconstruction

bmMitosius_create(mDir, y, t, ve);
disp('Mitosius files are saved!')
disp(mDir)

%% With binning
% binsDir = fullfile(reconDir, 'bins/clean', filesep);
% mask_type = sprintf('clean_comp/clean_comp_%d', nChCompressed);   % use char instead of string
% mDir = fullfile(reconDir, 'mitosius', mask_type);
% if ~exist(mDir, 'dir')
%     mkdir(mDir);
% end

% for region_idx = 0:3

%     th_ratio = 0.75;
%     mDir = fullfile(reconDir, 'mitosius', mask_type, ['mask_', num2str(region_idx)]);
%     eMaskFilePath = fullfile(binsDir, sprintf('eMask_th%.2f_region%i.mat', th_ratio, region_idx));

%     eyeMask = load(eMaskFilePath);
%     fields = fieldnames(eyeMask);  % Get the field names
%     firstField = fields{1};  % Get the first field name
%     eyeMask = eyeMask.(firstField);  % Access the first field's value
%     disp(eMaskFilePath)
%     disp('is loaded!')
    
%     % Eliminate the first segment of all the spokes for accuracies
%     size_Mask = size(eyeMask);
%     nbins = size_Mask(1);
%     eyeMask = reshape(eyeMask, [nbins, reader.acquisitionParams.nSeg, reader.acquisitionParams.nShot]); 
%     eyeMask(:, 1, :) = [];  % SI
%     eyeMask(:, :, 1:reader.acquisitionParams.nShot_off) = [];  % SS
%     eyeMask = bmPointReshape(eyeMask); 
        
%     % Run the mitosis function and compute volume elements
%     [y, t] = bmMitosis(y_tot_comp_norm, t_tot, eyeMask); 
%     y = bmPermuteToCol(y); 
%     ve  = bmVolumeElement(t, 'voronoi_full_radial3'); 

%     % Save all the resulting datastructures on the disk. You are now ready
%     % to run your reconstruction
%     bmMitosius_create(mDir, y, t, ve); 
%     disp('Mitosius files are saved!')
%     disp(mDir)

% end
