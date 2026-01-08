function mat2nii_batch(inputFolder, varargin)
%MAT2NII_BATCH Convert all .mat files in a folder to NIfTI format
%
%   mat2nii_batch(inputFolder) converts all .mat files in the specified folder
%   
%   mat2nii_batch(inputFolder, Name, Value) specifies additional options:
%       'ReconType'     - Reconstruction type: 'th8' (Steva) or 'x0' (Mathilda, default)
%       'RefNifti'      - Path to reference NIfTI file for header info
%                         (default: '/home/debi/Desktop/tmp/webplatform/input/2022160100001.nii.gz')
%       'OutputDir'     - Output directory (default: same as input folder)
%       'KeepAffine'    - Keep affine from reference (true/false, default: false)
%       'Recursive'     - Search recursively in subfolders (true/false, default: false)
%
%   Example:
%       mat2nii_batch('/path/to/folder')
%       mat2nii_batch('/path/to/folder', 'ReconType', 'th8')
%       mat2nii_batch('/path/to/folder', 'OutputDir', '/output/path', 'Recursive', true)

%% Parse inputs
p = inputParser;
addRequired(p, 'inputFolder', @ischar);
addParameter(p, 'ReconType', 'x0', @(x) ismember(x, {'th8', 'x0'}));
addParameter(p, 'RefNifti', '/home/debi/Desktop/tmp/webplatform/input/2022160100001.nii.gz', @ischar);
addParameter(p, 'OutputDir', '', @ischar);
addParameter(p, 'KeepAffine', false, @islogical);
addParameter(p, 'Recursive', false, @islogical);

parse(p, inputFolder, varargin{:});

inputFolder = p.Results.inputFolder;
recon = p.Results.ReconType;
ref_nii_path = p.Results.RefNifti;
outputDir = p.Results.OutputDir;
keepAffine = p.Results.KeepAffine;
recursive = p.Results.Recursive;

%% Validate input folder
if ~exist(inputFolder, 'dir')
    error('Input folder does not exist: %s', inputFolder);
end

%% Find all .mat files
if recursive
    matFiles = dir(fullfile(inputFolder, '**', '*.mat'));
    fprintf('Searching recursively for .mat files in: %s\n', inputFolder);
else
    matFiles = dir(fullfile(inputFolder, '*.mat'));
    fprintf('Searching for .mat files in: %s\n', inputFolder);
end

if isempty(matFiles)
    warning('No .mat files found in: %s', inputFolder);
    return;
end

fprintf('Found %d .mat file(s) to process\n\n', length(matFiles));

%% Process each .mat file
successCount = 0;
errorCount = 0;

for i = 1:length(matFiles)
    matFilePath = fullfile(matFiles(i).folder, matFiles(i).name);
    
    fprintf('[%d/%d] Processing: %s\n', i, length(matFiles), matFiles(i).name);
    
    try
        % Prepare arguments for mat2nii_single
        args = {};
        args{end+1} = 'ReconType';
        args{end+1} = recon;
        args{end+1} = 'RefNifti';
        args{end+1} = ref_nii_path;
        args{end+1} = 'KeepAffine';
        args{end+1} = keepAffine;
        
        % Only add OutputDir if specified
        if ~isempty(outputDir)
            args{end+1} = 'OutputDir';
            args{end+1} = outputDir;
        end
        
        % Call mat2nii_single
        mat2nii_single(matFilePath, args{:});
        successCount = successCount + 1;
        
    catch ME
        errorCount = errorCount + 1;
        fprintf('  ❌ Error: %s\n', ME.message);
    end
    
    fprintf('\n');
end

%% Summary
fprintf('=====================================\n');
fprintf('Processing Complete\n');
fprintf('=====================================\n');
fprintf('Total files: %d\n', length(matFiles));
fprintf('Successful: %d\n', successCount);
fprintf('Failed: %d\n', errorCount);

if errorCount > 0
    fprintf('\n⚠️  Some files could not be processed. Check messages above.\n');
else
    fprintf('\n🎉 All files processed successfully!\n');
end

end
