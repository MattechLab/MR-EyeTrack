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
addpath(genpath('/home/debi/jaime/repos/MR-EyeTrack/recon/Recon_scripts'));
addpath(genpath('/home/debi/MatTechLab/internal_monalisa'));
addpath(genpath('/Users/cag/Documents/forclone/pulseq_v15'));
addpath(genpath('/home/debi/yiwei/forclone/pulseq'));
addpath(genpath('/home/debi/yiwei/forclone/pulseqmreye'));

%% Initialize the directories and acquire the Coil

% subject_num = 3;
datatype = 2;

% subject_suffix = {'_ml', '_jb', '_yj',  '_phantom'};
mask_note_list{1}= 'ori'; mask_note_list{2}= 'ptp';
mask_note = mask_note_list{datatype};

seqFolder = '/home/debi/jaime/repos/MR-EyeTrack/data/test-pulseq/pulseq/sub_yj';
datasetDir = '/home/debi/jaime/repos/MR-EyeTrack/data/test-pulseq/pulseq/sub_yj';
reconDir = '/home/debi/jaime/repos/MR-EyeTrack/results';

seqName_list = {
    '/yj_seq2_t1w_libre_main_TR8.0ms_TE3.6ms_swap1_FA4_RF2_mreye_track_3723_traj_ptp.seq', ...
    '/yj0_seq8_t1w_libre_pre_TR6.2ms_TE3.6ms_swap1_FA4_RF2_rfmod2_trajPTP_PhNeg.seq'};

seqName = seqName_list{datatype};

bodyCoilFile    = [datasetDir, '/meas_MID00096_FID00911_track_ptp_BC.dat'];
arrayCoilFile   = [datasetDir, '/meas_MID00095_FID00910_track_ptp_HC.dat'];
measureFile     = [datasetDir, '/meas_MID00094_FID00909_track_ptp.dat'];

%% Load and Configure Data
% Read data using the library's `createRawDataReader` function
% This readers makes the usage of Siemens and ISMRMRD files equivalent for
% the library
bodyCoilReader = createRawDataReader(bodyCoilFile, true);
% bodyCoilReader.acquisitionParams.nShot_off = 14;
bodyCoilReader.acquisitionParams.traj_type = 'pulseq';
bodyCoilReader.acquisitionParams.pulseqTrajFile_name = strcat(seqFolder, seqName);
%
arrayCoilReader = createRawDataReader(arrayCoilFile, true);
% arrayCoilReader.acquisitionParams.nShot_off = 14;
arrayCoilReader.acquisitionParams.traj_type = 'pulseq';
arrayCoilReader.acquisitionParams.pulseqTrajFile_name = strcat(seqFolder, seqName);

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
saveCDirList = {strcat('/Sub',num2str(subject_num),'/T1_LIBRE_Binning/C/'),
    strcat('/Sub',num2str(subject_num),'/T1_LIBRE_woBinning/C/')};


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
