function motionParams = estimate_rigid_motion_params_spm( ...
        matFile, volumeField, pixelDimensions, outputPrefix)
%ESTIMATE_RIGID_MOTION_PARAMS_SPM
%
% Estimates rigid-body motion parameters using SPM realignment. This
% function saves a 4D nifti file corresponding to the 4D timeseries input
% (Monalisa recon). And then spm is used to esimtate motion of each frame 
% w.r.t the first frame. These motion parameters can then be used to 
% motion correct the rawdata files. 
% Pipeline:
%   MAT-file (cell array of 3D volumes)
%       -> 4D NIfTI
%       -> SPM realign (estimate + reslice)
%       -> motion parameters (Tx Ty Tz Rx Ry Rz)
%
% INPUTS
%   matFile         : path to .mat file containing reconstructed volumes
%   volumeField     : variable name in MAT file (e.g. 'x_cs', 'x_test')
%   pixelDimensions : [dx dy dz dt] voxel size + TR
%   outputPrefix    : base name for output files
%
% OUTPUT
%   motionParams    : N x 6 matrix (translations [mm], rotations [rad])
%
% NOTE: SPM axis are different from Monalisa's axis. This difference results 
% in wrong rotations and translations if not correclty handled.


%% ---------------------------------------------------------
%% 0. Sanity checks
%% ---------------------------------------------------------

% Check MAT file
assert(exist(matFile,'file')==2, 'MAT file not found');

% Check SPM installation
if exist('spm','file') ~= 2
    error([ ...
        'SPM not found on the MATLAB path.\n' ...
        'You need to install SPM and add it to MATLAB:\n\n' ...
        '  addpath(genpath(''/path/to/spm''))\n\n' ...
        'SPM download: https://www.fil.ion.ucl.ac.uk/spm/' ...
    ]);
end

%% ---------------------------------------------------------
%% 0. Setup
%% ---------------------------------------------------------

[root,~,~] = fileparts(matFile);

niftiFile = fullfile(root, [outputPrefix '.nii']);
motionParamsFile = fullfile(root, ['rp_' outputPrefix '.txt']);

spm('Defaults','fmri');
spm_jobman('initcfg');

%% ---------------------------------------------------------
%% 1. Load reconstructed volumes
%% ---------------------------------------------------------
S = load(matFile);

assert(isfield(S, volumeField), ...
    'Variable "%s" not found in MAT file', volumeField);

x = S.(volumeField);           % cell array of 3D volumes
N = numel(x);

volSize = size(x{1});
vol4d = zeros([volSize N], 'single');

for k = 1:N
    vol4d(:,:,:,k) = abs(x{k});   % rigid motion works on magnitude
end

%% ---------------------------------------------------------
%% 2. Write 4D NIfTI
%% ---------------------------------------------------------
niftiwrite(vol4d, niftiFile);

niiInfo = niftiinfo(niftiFile);
niiInfo.PixelDimensions = pixelDimensions;

niftiwrite(vol4d, niftiFile, niiInfo, 'Compressed', false);

fprintf('[SPM] 4D NIfTI written: %s\n', niftiFile);

%% ---------------------------------------------------------
%% 3. SPM realignment (rigid-body)
%% ---------------------------------------------------------
matlabbatch = [];

matlabbatch{1}.spm.spatial.realign.estwrite.data = { cellstr(niftiFile) };

% --- estimation options (explicit = reproducible)
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep     = 2;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm    = 5;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm     = 0;  % to first volume
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp  = 2;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap    = [0 0 0];

% --- reslicing
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which  = [2 1];
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap   = [0 0 0];
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = 'r';

spm_jobman('run', matlabbatch);

%% ---------------------------------------------------------
%% 4. Load motion parameters
%% ---------------------------------------------------------
assert(exist(motionParamsFile,'file')==2, ...
    'Motion parameter file not found: %s', motionParamsFile);

motionParams = load(motionParamsFile);

% plot extracted params
% PLOT
figure('Color','w'); set(gca,'Color','w');   % or: gca.Color = 'w';
subplot(2,1,1);
plot(motionParams(:, 1:3), 'LineWidth',2)
xlabel('Time (Volume)'); ylabel('translation (mm/)');
legend('Tx', 'Ty', 'Tz');
title('Translation Parameters');
grid on;

subplot(2,1,2);
plot(motionParams(:, 4:6), 'LineWidth',2)
xlabel('Time (Volume)'); ylabel('rotation (rad)');
legend( 'Rx', 'Ry', 'Rz');
title('Rotation Parameters');
grid on;

end