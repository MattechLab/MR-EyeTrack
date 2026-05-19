%% check_nifti_orientation.m
% Compare affine orientation between the DICOM-derived MPRAGE and the
% pipeline reconstruction for sub-001.

baseDir = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'data', 'study');

refPath  = fullfile(baseDir, 'sub-001', 'dicom', ...
    'csTFL_mp-rage_1mm-iso_CP_acc4.6_5_MR', 'sub-001.nii.gz');
reconPath = fullfile(baseDir, 'sub-001', 'recon', 'woBin', ...
    'x_steva_nIter_20_delta_1.000.nii.gz');

assert(exist(refPath,   'file') == 2, 'Reference NIfTI not found:\n  %s', refPath);
assert(exist(reconPath, 'file') == 2, 'Recon NIfTI not found:\n  %s', reconPath);

%% Load headers
infoRef   = niftiinfo(refPath);
infoRecon = niftiinfo(reconPath);

%% Extract affine (stored as Transform.T, row-major 4x4 transposed)
% niftiinfo stores Transform.T such that world = [vox,1] * T
% For comparison we work with the standard column-major A: world = A * [vox;1]
A_ref   = infoRef.Transform.T';
A_recon = infoRecon.Transform.T';

%% Print summary
fprintf('\n===== NIfTI Orientation Check =====\n');
fprintf('  RAS convention: X = Left→Right, Y = Posterior→Anterior, Z = Inferior→Superior\n\n');

fprintf('--- Reference (DICOM MPRAGE) ---\n');
fprintf('  File       : %s\n', refPath);
fprintf('  Dimensions : %s  (voxels per axis)\n', mat2str(infoRef.ImageSize));
fprintf('  PixDim     : %s mm  (voxel size; extra values = TR for 4D)\n', mat2str(infoRef.PixelDimensions, 4));
fprintf('  Affine     : 4x4 matrix mapping voxel [i,j,k] -> mm in RAS world space\n');
fprintf('               top-left 3x3 = rotation+scaling, last column = translation (origin in mm)\n');
disp(A_ref);

fprintf('--- Reconstruction (woBin Steva) ---\n');
fprintf('  File       : %s\n', reconPath);
fprintf('  Dimensions : %s  (voxels per axis)\n', mat2str(infoRecon.ImageSize));
fprintf('  PixDim     : %s mm\n', mat2str(infoRecon.PixelDimensions, 4));
fprintf('  Affine:\n');
disp(A_recon);

%% Rotation matrices (normalise columns of the 3x3 sub-matrix)
% Dividing by voxel size strips the scaling, leaving a pure rotation matrix
% whose columns are the anatomical directions each image axis points toward
R_ref   = A_ref(1:3,1:3)   ./ infoRef.PixelDimensions(1:3);
R_recon = A_recon(1:3,1:3) ./ infoRecon.PixelDimensions(1:3);

% Relative rotation: R_rel = R_recon * R_ref^-1
% Identity = perfectly aligned; any deviation = rotation between the two
R_rel = R_recon / R_ref;

fprintf('--- Relative rotation (R_recon * R_ref^-1) ---\n');
fprintf('  Identity matrix = same orientation; off-diagonal values = misalignment\n');
disp(R_rel);

% Angle extracted via the rotation-matrix trace formula: trace(R) = 1 + 2*cos(theta)
theta_rad = acos(max(-1, min(1, (trace(R_rel) - 1) / 2)));
fprintf('  Rotation angle between frames: %.2f deg\n', rad2deg(theta_rad));
fprintf('  (0° = identical orientation)\n\n');

%% Axis direction table
% Each column of R is a unit vector [x,y,z] in RAS space showing where that
% image axis points anatomically.  A pure RAS image would give the identity:
%   axis-1 -> [1,0,0] (pure Right), axis-2 -> [0,1,0] (pure Anterior), etc.
% Deviations from identity reveal head tilt / oblique positioning.
% Tilt angle for each axis = acos(dominant component).
labels  = {'axis-1 (col)', 'axis-2 (row)', 'axis-3 (slice)'};
rasAxes = {'R→L', 'P→A', 'I→S'};

fprintf('--- Axis directions in RAS space ---\n');
fprintf('  Each image axis is expressed as a unit vector [R, A, S].\n');
fprintf('  The dominant component tells you the main anatomical direction;\n');
fprintf('  smaller components indicate oblique tilt.\n\n');
fprintf('  %-16s  %-36s  %-36s\n', 'Image axis', 'Reference (MPRAGE)', 'Reconstruction');

for ax = 1:3
    dir_ref   = R_ref(:, ax);
    dir_recon = R_recon(:, ax);

    % Dominant anatomical direction and tilt angle for the reference
    [~, dom] = max(abs(dir_ref));
    tilt_deg = rad2deg(acos(abs(dir_ref(dom))));

    fprintf('  %-16s  [R=%+.3f A=%+.3f S=%+.3f]  [R=%+.3f A=%+.3f S=%+.3f]\n', ...
        labels{ax}, ...
        dir_ref(1),   dir_ref(2),   dir_ref(3), ...
        dir_recon(1), dir_recon(2), dir_recon(3));
    fprintf('  %-16s  -> mainly %s, %.1f° oblique tilt\n', '', rasAxes{dom}, tilt_deg);
end

%% Central-slice visual comparison
fprintf('\nLoading image data for visual comparison...\n');
volRef   = niftiread(refPath);
volRecon = niftiread(reconPath);

% Normalise for display
norm_vol = @(v) double(abs(v)) / double(max(abs(v(:))) + eps);
vRef   = norm_vol(volRef);
vRecon = norm_vol(volRecon);

cRef   = round(size(vRef,   1:3) / 2);
cRecon = round(size(vRecon, 1:3) / 2);

figure('Name', 'Orientation Check — Central Slices', 'NumberTitle', 'off', ...
    'Position', [100 100 1200 700]);

views   = {'Axial (z)',   'Coronal (y)', 'Sagittal (x)'};
slRef   = {squeeze(vRef(:,:,cRef(3))),   squeeze(vRef(:,cRef(2),:)),   squeeze(vRef(cRef(1),:,:))};
slRecon = {squeeze(vRecon(:,:,cRecon(3))), squeeze(vRecon(:,cRecon(2),:)), squeeze(vRecon(cRecon(1),:,:))};

for v = 1:3
    subplot(2,3,v);
    imagesc(slRef{v}'); axis image off; colormap gray;
    title(sprintf('REF — %s', views{v}), 'FontSize', 9);

    subplot(2,3,v+3);
    imagesc(slRecon{v}'); axis image off; colormap gray;
    title(sprintf('RECON — %s', views{v}), 'FontSize', 9);
end
sgtitle('Reference (MPRAGE) vs Reconstruction (woBin Steva) — sub-001');
