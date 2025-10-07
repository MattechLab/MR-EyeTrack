% This script is used for generate the undersampled binning masks 
% according to the ET mask, where
% sampling rate of ET mask: 1ms
% sampling rate of readouts: TR=8ms
% =====================================================
clearvars; clc;
addpath(genpath('/home/debi/jaime/repos/MR-EyeTrack/recon'));
addpath(genpath('/home/debi/MatTechLab/monalisa'));
addpath(genpath('/home/debi/jaime/repos/MR-EyeTrack/recon/Recon_scripts'));

%% Config

% Variables
subject_num = 3;
mask_type = 'clean';
undersampling_factor = 0.95; % remove % of the original readouts

reconDir = '/home/debi/jaime/repos/MR-EyeTrack/results';
otherDir = [reconDir, '/Sub00', num2str(subject_num),'/T1_LIBRE_Binning/other/', mask_type, '/'];
outputDir = [reconDir, '/Sub00', num2str(subject_num),'/T1_LIBRE_Binning/other/', mask_type, '_', num2str(undersampling_factor), '/'];

% Ensure directories exist
if ~isfolder(otherDir), mkdir(otherDir); end
if ~isfolder(outputDir), mkdir(outputDir); end

%% Random undersample

region_idx_list = 0:3;  % 0:up 1:down 2:left 3:right 4:center mask

for region_idx = region_idx_list
    % Load the mask
    eMaskFilePath = fullfile(otherDir, ['eMask_th0.75_region', num2str(region_idx), '.mat']);
    mask = load(eMaskFilePath);
    eMask = mask.eMaskN;

    % Find indices of 1's
    oneIdx = find(eMask == 1);

    % Randomly choose half of them
    numToZero = floor(numel(oneIdx) * undersampling_factor);
    idxToZero = randsample(oneIdx, numToZero);

    % Set those entries to 0
    eMask(idxToZero) = 0;

    % Convert to logical type for minimal memory usage
    eMask = logical(eMask);

    % Save the modified mask
    savePath = [outputDir, 'eMask_th0.75_region', num2str(region_idx), '.mat'];
    save(savePath, 'eMask');
    disp(['Undersampled mask saved to: ', savePath]);
end
