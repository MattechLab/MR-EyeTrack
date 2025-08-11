% =====================================================
% Author: Yiwei Jia
% Date: April 3
% ------------------------------------------------
% This script is used for generate mask 
% to eliminate the readouts in the non-steady state
% nShotOff * nSeg
% =====================================================
clear; clc;
addpath(genpath('/home/debi/yiwei/forclone/Recon_scripts'));


%%


subject_num = 1;

datasetDir = '/home/debi/jaime/repos/MR-EyeTrack/data/pilot/sub-01/rawdata/';
reconDir = '/home/debi/jaime/tmp/250613_JB/';
otherDir = [reconDir, '/Sub00', num2str(subject_num),'/T1_LIBRE_woBinning/other/'];

mask_note_list={'woBin','maskNote_for_sub2', 'maskNote_for_sub3'};

mask_note = mask_note_list{subject_num};
x0Dir = [reconDir, '/Sub00', num2str(subject_num), '/T1_LIBRE_woBinning/output/mask_', mask_note, '/'];

% common for all subjects
bodyCoilFile     = [datasetDir, '/meas_MID00614_FID182868_BEAT_LIBREon_eye_BC_BC.dat'];
arrayCoilFile    = [datasetDir, '/meas_MID00615_FID182869_BEAT_LIBREon_eye_HC_BC.dat'];
measureFile = [datasetDir, '/meas_MID00605_FID182859_BEAT_LIBREon_eye_(23_09_24).dat'];

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

%%
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



%%

