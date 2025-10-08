%% Clear environment
clc; clearvars; close all;

%% Paths
imgDir = '/home/debi/jaime/repos/MR-EyeTrack/results/Sub001/T1_LIBRE_woBinning/output/th8'; % directory with .mat files
imgName = 'x_steva_nIter_20_delta_1.000.mat'; % reconstructed image
ref_nii_path = '/home/debi/Desktop/tmp/webplatform/input/2022160100001.nii.gz'; % reference image

%% Select reconstructed image (.mat)
% [imgName, imgDir, ~] = uigetfile( ...
%     {'*.mat','Select image to map (*.mat)'}, ...
%     'Pick a file', 'MultiSelect', 'off', imgDir);

% if isequal(imgName,0)
%     error('No .mat file selected.');
% end

imgFilepath = fullfile(imgDir,imgName);
data = load(imgFilepath);
im_cs = data.x;   % assuming your data is stored in variable "x"

%% Compute sum of squares (if needed)
function result = sumOfSquares(data, dim)
    result = sqrt(sum(abs(data).^2, dim));
end
im_sos = sumOfSquares(im_cs,5);
Functional_Recon_all_lines = abs(sumOfSquares(im_sos,4));
Functional_Recon_all_lines_flipped = flip(Functional_Recon_all_lines,3);

%% Adjust orientation before writing
Functional_Recon_all_lines_flipped = permute(Functional_Recon_all_lines_flipped, [2 1 3]); % swap x/y
Functional_Recon_all_lines_flipped = flip(Functional_Recon_all_lines_flipped, 1);          % flip L-R

%% Load reference header
ref_info = niftiinfo(ref_nii_path);

%% Adapt header for new image
new_info = ref_info;
new_info.ImageSize = size(Functional_Recon_all_lines_flipped);
new_info.Filename = fullfile(imgDir, 'x_flipped.nii');
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
fprintf('✅ Saved new image with copied header as: %s.gz\n', new_info.Filename);
