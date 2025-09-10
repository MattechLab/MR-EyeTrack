% =====================================================
% Author: Yiwei Jia
% Date: April 3
% ------------------------------------------------
% This script is used for generate mask 
% to eliminate the readouts in the non-steady state
% nShotOff * nSeg
% =====================================================
clear; clc;
addpath(genpath('/home/debi/jaime/repos/MR-EyeTrack/recon/Recon_scripts'));
addpath(genpath('/home/debi/MatTechLab/monalisa'));

%%
subject_num = 2;

% Paths
datasetDir = ['/home/debi/jaime/repos/MR-EyeTrack/data/pilot/sub-00', num2str(subject_num), '/rawdata'];
reconDir = '/home/debi/jaime/repos/MR-EyeTrack/results';
otherDir = [reconDir, '/Sub00', num2str(subject_num),'/T1_LIBRE_woBinning/other/'];
saveCDir = [reconDir, strcat('/Sub00',num2str(subject_num),'/T1_LIBRE_woBinning/C/')];

mask_note = 'woBin';
x0Dir = [reconDir, '/Sub00', num2str(subject_num), '/T1_LIBRE_woBinning/output/mask_', mask_note, '/'];

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

% Check if the directory exists
if ~isfolder(otherDir)
    % If it doesn't exist, create it
    mkdir(otherDir);
    disp(['Directory created: ', otherDir]);
else
    disp(['Directory already exists: ', otherDir]);
end

%%
autoRead = 1;
% 1:if you don't know the nSeg and nShot of the raw data, set 1, the reader
% will take care of identifying that
% 0: if you're aware of the nSeg and nShot, set 0, it will skip the
% time used for reading and parsing the raw data.

if autoRead
    reader = createRawDataReader(measureFile, autoRead);
    reader.acquisitionParams.nShot_off = 14;
    reader.acquisitionParams.traj_type = 'full_radial3_phylotaxis';
    nShotOff = 14; 
    nShot = reader.acquisitionParams.nShot;
    nSeg = reader.acquisitionParams.nSeg;
else
    %manually modify the params below
    nShotOff = 14; 
    nShot = 419;
    nSeg = 22; 
end
eMask = ones(1, nShot*nSeg);
eMask = (eMask>0);
eMask(1:nShotOff*nSeg) = 0;
% Saving data and Convert to Monalisa format
%--------------------------------------------------------------------------    
eMaskFilePath = [otherDir,'eMask_woBin.mat'];

% Save the CMask to the .mat file
save(eMaskFilePath, 'eMask');
disp('eMask has been saved here:')
disp(eMaskFilePath)
