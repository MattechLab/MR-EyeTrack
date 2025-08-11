clear; clc; close all;

% =====================================================
% Author: Yiwei Jia
% Date: April 17
% ------------------------------------------------
% Coil sensitivity -> binning mask eMask -> [Mitosius]
% Update: the eyeGenerateBinning is replaced with eyeGenerateBinningWin
% Add a new txt file saved along side with the mask to log the details
% =====================================================

%%
subject_num = 1;

datasetDir = '/home/debi/jaime/repos/MR-EyeTrack/data/pilot/sub-01/rawdata/';
reconDir = '/home/debi/jaime/tmp/250613_JB/';

mask_note_list={'woBin','maskNote_for_sub2', 'maskNote_for_sub3'};

mask_note = mask_note_list{subject_num};
x0Dir = [reconDir, '/Sub00', num2str(subject_num), '/T1_LIBRE_woBinning/output/mask_', mask_note, '/'];

% common for all subjects
bodyCoilFile     = [datasetDir, '/meas_MID00614_FID182868_BEAT_LIBREon_eye_BC_BC.dat'];
arrayCoilFile    = [datasetDir, '/meas_MID00615_FID182869_BEAT_LIBREon_eye_HC_BC.dat'];
measureFile = [datasetDir, '/meas_MID00605_FID182859_BEAT_LIBREon_eye_(23_09_24).dat'];

otherDir = [reconDir, strcat('/Sub00',num2str(subject_num),'/T1_LIBRE_woBinning/other/')];

% Check if the directory exists
if ~isfolder(otherDir)
    % If it doesn't exist, create it
    mkdir(otherDir);
    disp(['Directory created: ', otherDir]);
else
    disp(['Directory already exists: ', otherDir]);
end

saveCDirList = {strcat('/Sub00',num2str(subject_num),'/T1_LIBRE_Binning/C/'),
    strcat('/Sub00',num2str(subject_num),'/T1_LIBRE_woBinning/C/')};

%% Step 1: Load the Raw Data

autoFlag = true;  % Disable validation UI
reader = createRawDataReader(measureFile, autoFlag);
reader.acquisitionParams.nShot_off = 14;
reader.acquisitionParams.traj_type = 'full_radial3_phylotaxis';
%
% Load the raw data and compute trajectory and volume elements
y_tot = reader.readRawData(true, true);  % Filter nshotoff and SI

%%
t_tot = bmTraj(reader.acquisitionParams);                 % Compute trajectory

ve_tot = bmVolumeElement(t_tot, 'voronoi_full_radial3');  % Volume elements

%% Step 2: Load Coil Sensitivity Maps

% Load the coil sensitivity previously measured
saveCDir = [reconDir,saveCDirList{2}];
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
% is the smaller than the full resolution.
% ===============================================
matrix_size = 240;  % Max nominal spatial resolution
N_u = [matrix_size, matrix_size, matrix_size];
dK_u = [1, 1, 1]./FoV;

%%
C = bmImResize(C, [48, 48, 48], N_u);
% C = flip(flip(flip(C, 1), 2), 3);

%% tmp: check HC

x0 = bmMathilda(y_tot, t_tot, ve_tot, C, N_u, N_u, dK_u);
bmImage(x0)

%% Coil compression (can be here or after normalization)
coilCompression = 1;
nChCompressed = 4;
nCh = size(y_tot, 1);
nx = size(y_tot, 2);
ntviews = size(y_tot, 3);
if coilCompression == 1
    D = reshape(y_tot, nx * ntviews, nCh);
    [U, S, V] = svd(D, 'econ');
    singular_values = diag(S);
    total_variance = sum(singular_values .^ 2);
    explained_variance = singular_values(1:nChCompressed) .^ 2;
    percentage_explained = (explained_variance / total_variance) * 100;
    y_tot = reshape(D * V(:, 1:nChCompressed), nChCompressed, nx, ntviews);
    % Plot the explained variance
    f=figure;
    f.Position = [100 100 500 800];
    plot(percentage_explained,'.-','LineWidth',2,'Color','r','MarkerSize',20)
    xlabel('Virtual Coil nr.')
    ylabel('Explained [%]')
    title(sprintf('Coil Compression, Explanation for nChCompressed = %d is %.2f%%', nChCompressed, sum(percentage_explained(1:nChCompressed))))
else
    % coil combination to save memory (for now just coil selection) EP test
    selected_coils = [1,2,4,5,6:2:size(kdata_raw,3)];
    kdata_raw=kdata_raw(:,:,selected_coils);
end

%% Step 3: Normalize the Raw Data

if N_u >240
    normalization = false;
else 
    normalization = true;
end
if normalization
    x_tot = bmMathilda(y_tot, t_tot, ve_tot, C, N_u, N_u, dK_u); 
    x_perm = permute(x_tot, [2,3,1]);
    x0=x_tot;
    %
    bmImage(x_perm)
    %
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
        y_tot = y_tot/normalize_val; 
        y_tot(1,1,123)
    else
        y_tot = y_tot/(2.5e-10); 
        y_tot(1,1,123)
    end
end

%%
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

%% Set the folder for mitosius saving

mDir = [reconDir, strcat('/Sub00', num2str(subject_num)),'/T1_LIBRE_woBinning/mitosius/mask_', mask_note, '/'];

%% Prepare eye mask

eMaskFilePath = [otherDir,'eMask_woBin'];

eyeMask = load(eMaskFilePath); 
fields = fieldnames(eyeMask);  % Get the field names
firstField = fields{1};  % Get the first field name
eyeMask = eyeMask.(firstField);  % Access the first field's value
disp(eMaskFilePath)
disp('is loaded!')
% Eleminate the first segment of all the spokes for accuracies

%
size_Mask = size(eyeMask);
nbins = size_Mask(1);
eyeMask = reshape(eyeMask, [nbins, reader.acquisitionParams.nSeg, reader.acquisitionParams.nShot]); 
eyeMask(:, 1, :) = []; 

eyeMask(:, :, 1:reader.acquisitionParams.nShot_off) = []; 
eyeMask = bmPointReshape(eyeMask); 


%% Run the mitosis function and compute volume elements

[y, t] = bmMitosis(y_tot, t_tot, eyeMask); 
y = bmPermuteToCol(y); 
ve  = bmVolumeElement(t, 'voronoi_full_radial3' ); 

% Save all the resulting datastructures on the disk. You are now ready
% to run your reconstruction

bmMitosius_create(mDir, y, t, ve); 
disp('Mitosius files are saved!')
disp(mDir)
























