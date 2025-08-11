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
% addpath(genpath("/media/sinf/1,0 TB Disk/Backup/Recon_fork"));
addpath(genpath('/home/debi/yiwei/forclone/Recon_scripts/'));

%% Initialize the directories and acquire the Coil
subject_num=2;

datasetDir = '/home/debi/jaime/repos/MR-EyeTrack/data/pilot';
reconDir = '/home/debi/jaime/tmp/250613_JB';

if subject_num == 1
    datasetDir = [datasetDir, '/sub-01/rawdata'];
    ETDir      = [datasetDir, '/masks_1206/Sub001'];
elseif subject_num == 2
    datasetDir = [datasetDir, '/sub-02/rawdata'];
    ETDir      = [datasetDir, '/masks_1206/Sub002'];
elseif subject_num == 3
    datasetDir = [datasetDir, '/sub-03/rawdata'];
    ETDir      = [datasetDir, '/masks_1206/Sub003'];
else
    datasetDir = [datasetDir, ' '];
    ETDir      = [datasetDir, ' '];
end

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


%% Load and Configure Data
% Read data using the library's `createRawDataReader` function
% This readers makes the usage of Siemens and ISMRMRD files equivalent for
% the library
bodyCoilreader = createRawDataReader(bodyCoilFile, false);
bodyCoilreader.acquisitionParams.nShot_off = 14;
bodyCoilreader.acquisitionParams.traj_type = 'full_radial3_phylotaxis';
%
arrayCoilReader = createRawDataReader(arrayCoilFile, false);
arrayCoilReader.acquisitionParams.nShot_off = 14;
arrayCoilReader.acquisitionParams.traj_type = 'full_radial3_phylotaxis';

% Ensure consistency in number o1f shot-off points
nShotOff = arrayCoilReader.acquisitionParams.nShot_off;

%%
quickCalC = 1; 
%1: quickly calculate C by the function organized by Mauro Leidi
%0: explore the details of C calculation in monalisa, with comments

if quickCalC
    autoFlag=false;
    nIter = 5;
    C = mlComputeCoilSensitivity(bodyCoilreader, arrayCoilReader, [48,48,48], autoFlag, nIter);
else
    % Parameters
    dK_u = [1, 1, 1] ./ arrayCoilReader.acquisitionParams.FoV;   % Cartesian grid spacing
    N_u = [48, 48, 48];             % Adjust this value as needed, low resolution is sufficient
    % Compute Trajectory and Volume Elements
    [y_body, t, ve] = bmCoilSense_nonCart_data(bodyCoilreader, N_u);
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
bmImage(C)
saveCDirList = {strcat('/Sub00',num2str(subject_num),'/T1_LIBRE_Binning/C/'),
    strcat('/Sub00',num2str(subject_num),'/T1_LIBRE_woBinning/C/')};


for idx = 1:2
    saveCDir     = [reconDir, saveCDirList{idx}];
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
end




