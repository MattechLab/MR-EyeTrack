%% Clear environment
clc; clearvars; close all;
addpath(genpath('/home/debi/jaime/repos/MR-EyeTrack/analysis'));
addpath(genpath('/home/debi/jaime/repos/MR-EyeTrack/results'));
addpath(genpath('/home/debi/MatTechLab/monalisa'));  % for bmImage

%% Paths
% Define variables
subject_num = '1';  % subject number
bin = '';  % binning type (can be '' or 'wo')
mask = 'clean';  % mask type
recon = 'x0';  % reconstruction type (can be 'th8' for Steva or 'x0' for Mathilda)

% Construct base directory path
if strcmp(bin, 'wo')
    baseDir = sprintf('/home/debi/jaime/repos/MR-EyeTrack/results/Sub00%s/T1_LIBRE_%sBinning/output/%s', subject_num, bin, recon);
else
    baseDir = sprintf('/home/debi/jaime/repos/MR-EyeTrack/results/Sub00%s/T1_LIBRE_%sBinning/output/%s/%s', subject_num, bin, mask, recon);
end
ref_nii_path = '/home/debi/Desktop/tmp/webplatform/input/2022160100001.nii.gz'; % reference image

%% Find all .mat files in the directory
matFiles = dir(fullfile(baseDir, '*.mat'));

if isempty(matFiles)
    error('No .mat files found in the specified directory: %s', baseDir);
end

fprintf('Found %d .mat files to convert:\n', length(matFiles));
for i = 1:length(matFiles)
    fprintf('  %d. %s\n', i, matFiles(i).name);
end
fprintf('\n');

%% Helper function for sum of squares
function result = sumOfSquares(data, dim)
    result = sqrt(sum(abs(data).^2, dim));
end

%% Process each .mat file
for fileIdx = 1:length(matFiles)
    fprintf('Processing file %d/%d: %s\n', fileIdx, length(matFiles), matFiles(fileIdx).name);
    
    imgName = matFiles(fileIdx).name;
    imgDir = matFiles(fileIdx).folder;
    imgFilepath = fullfile(imgDir, imgName);
    
    data = load(imgFilepath);
    % Load data based on reconstruction type
    if strcmp(recon, 'th8')
        im_cs = data.x;
    elseif strcmp(recon, 'x0')
        if iscell(data.x0)
            im_cs = data.x0{1};
        else
            im_cs = data.x0;
        end
    else
        error('Unknown reconstruction type: %s', recon);
    end

    %% Compute absolute value (images are already reconstructed and in 3D, no need for SoS)
    Functional_Recon_all_lines = abs(im_cs);

    %% Adjust orientation before writing
    Functional_Recon_all_lines_flipped = permute(Functional_Recon_all_lines, [2 1 3]); % swap x/y
    Functional_Recon_all_lines_flipped = flip(Functional_Recon_all_lines_flipped, 1);   % flip L-R
    Functional_Recon_all_lines_flipped = flip(Functional_Recon_all_lines_flipped, 3);   % flip in 3rd dimension

    %% Load reference header
    ref_info = niftiinfo(ref_nii_path);

    %% Adapt header for new image
    new_info = ref_info;
    new_info.ImageSize = size(Functional_Recon_all_lines_flipped);
    
    % Generate output filename based on input .mat filename
    [~, baseFileName, ~] = fileparts(imgName);
    outputFileName = fullfile(imgDir, [baseFileName '.nii']);
    new_info.Filename = outputFileName;
    
    new_info.Filemoddate = datetime("now");
    new_info.Description = 'Reconstructed image using MR-EyeTrack pipeline';

    %% Match datatype
    new_info.Datatype = class(Functional_Recon_all_lines_flipped); % e.g. 'single' or 'uint16'
    new_info.BitsPerPixel = 8 * numel(typecast(cast(0, new_info.Datatype), 'uint8'));

    %% Option 1: Keep affine from reference (for alignment)
    % (use this if you want to overlay on the reference scan)
    % new_info.Transform.T = ref_info.Transform.T;

    %% Option 2: Reset affine (no cropping, voxel space only)
    new_info.Transform.T = eye(4);
    new_info.raw.qform_code = 0;
    new_info.raw.sform_code = 0;

    %% Write new NIfTI with copied header info
    niftiwrite(Functional_Recon_all_lines_flipped, new_info.Filename, new_info);

    %% Compress to .nii.gz
    gzip(new_info.Filename);
    delete(new_info.Filename);
    
    fprintf('✅ Saved: %s.gz\n\n', new_info.Filename);
end

fprintf('🎉 Conversion complete! Processed %d files.\n', length(matFiles));
