function mat2nii_single(matFilePath, varargin)
%MAT2NII_SINGLE Convert a single .mat file to NIfTI format
%
%   mat2nii_single(matFilePath) converts the specified .mat file to NIfTI
%   
%   mat2nii_single(matFilePath, Name, Value) specifies additional options:
%       'ReconType'     - Reconstruction type: 'th8' (Steva) or 'x0' (Mathilda, default)
%       'RefNifti'      - Path to reference NIfTI file for header info
%                         (default: '/home/debi/Desktop/tmp/webplatform/input/2022160100001.nii.gz')
%       'OutputDir'     - Output directory (default: same as input file)
%       'KeepAffine'    - Keep affine from reference (true/false, default: false)
%
%   Example:
%       mat2nii_single('/path/to/your/file.mat')
%       mat2nii_single('/path/to/your/file.mat', 'ReconType', 'th8')
%       mat2nii_single('/path/to/your/file.mat', 'OutputDir', '/output/path')

%% Parse inputs
p = inputParser;
addRequired(p, 'matFilePath', @ischar);
addParameter(p, 'ReconType', 'x0', @(x) ismember(x, {'th8', 'x0'}));
addParameter(p, 'RefNifti', '/home/debi/Desktop/tmp/webplatform/input/2022160100001.nii.gz', @ischar);
addParameter(p, 'OutputDir', '', @ischar);
addParameter(p, 'KeepAffine', false, @islogical);

parse(p, matFilePath, varargin{:});

matFilePath = p.Results.matFilePath;
recon = p.Results.ReconType;
ref_nii_path = p.Results.RefNifti;
outputDir = p.Results.OutputDir;
keepAffine = p.Results.KeepAffine;

%% Add required paths
addpath(genpath('/home/debi/jaime/repos/MR-EyeTrack/analysis'));
addpath(genpath('/home/debi/jaime/repos/MR-EyeTrack/results'));
addpath(genpath('/home/debi/MatTechLab/monalisa'));  % for bmImage

%% Validate input file
if ~exist(matFilePath, 'file')
    error('Input .mat file does not exist: %s', matFilePath);
end

if ~endsWith(matFilePath, '.mat')
    error('Input file must be a .mat file: %s', matFilePath);
end

%% Validate reference NIfTI
if ~exist(ref_nii_path, 'file')
    error('Reference NIfTI file does not exist: %s', ref_nii_path);
end

%% Set output directory
[inputDir, inputName, ~] = fileparts(matFilePath);
if isempty(outputDir)
    outputDir = inputDir;
end

if ~exist(outputDir, 'dir')
    mkdir(outputDir);
    fprintf('Created output directory: %s\n', outputDir);
end

%% Load .mat file
fprintf('Loading: %s\n', matFilePath);
data = load(matFilePath);

% Determine field name based on filename
if contains(inputName, 'xrms')
    fieldName = 'xrms';
elseif contains(inputName, 'x0')
    fieldName = 'x0';
elseif contains(inputName, 'x')
    fieldName = 'x';
else
    % Fall back to reconstruction type if filename doesn't indicate field
    if strcmp(recon, 'th8')
        fieldName = 'x';
    elseif strcmp(recon, 'x0')
        fieldName = 'x0';
    else
        error('Unknown reconstruction type: %s', recon);
    end
end

% Load data based on field name
if ~isfield(data, fieldName)
    error('Field "%s" not found in .mat file', fieldName);
end

fprintf('Loading field: %s\n', fieldName);
if iscell(data.(fieldName))
    im_cs = data.(fieldName){1};
else
    im_cs = data.(fieldName);
end

fprintf('Image dimensions: %s\n', mat2str(size(im_cs)));

%% Compute absolute value
Functional_Recon_all_lines = abs(im_cs);

%% Adjust orientation - rotate axial view 90° to the left
Functional_Recon_all_lines_flipped = rot90(Functional_Recon_all_lines, -1);

%% Load reference header
fprintf('Loading reference NIfTI header from: %s\n', ref_nii_path);
ref_info = niftiinfo(ref_nii_path);

%% Adapt header for new image
new_info = ref_info;
new_info.ImageSize = size(Functional_Recon_all_lines_flipped);

% Generate output filename
outputFileName = fullfile(outputDir, [inputName '.nii']);
new_info.Filename = outputFileName;

new_info.Filemoddate = datetime("now");
new_info.Description = sprintf('Reconstructed image using MR-EyeTrack pipeline (recon: %s)', recon);

%% Match datatype
new_info.Datatype = class(Functional_Recon_all_lines_flipped);
new_info.BitsPerPixel = 8 * numel(typecast(cast(0, new_info.Datatype), 'uint8'));

%% Set affine transformation
if keepAffine
    % Keep affine from reference (for alignment)
    new_info.Transform.T = ref_info.Transform.T;
    fprintf('Using affine transformation from reference image\n');
else
    % Reset affine (no cropping, voxel space only)
    new_info.Transform.T = eye(4);
    new_info.raw.qform_code = 0;
    new_info.raw.sform_code = 0;
    fprintf('Using identity affine transformation\n');
end

%% Write new NIfTI with copied header info
fprintf('Writing NIfTI file...\n');
niftiwrite(Functional_Recon_all_lines_flipped, new_info.Filename, new_info);

%% Compress to .nii.gz
fprintf('Compressing to .nii.gz...\n');
gzip(new_info.Filename);
delete(new_info.Filename);

fprintf('✅ Successfully saved: %s.gz\n', new_info.Filename);

end
