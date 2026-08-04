function check_nifti_orientation(refPath, reconPath, figTitle)
%CHECK_NIFTI_ORIENTATION  Compare affine orientation between two NIfTI files.
%
% Prints:
%   Affine matrices for both images
%   Rotation matrices (direction cosines per axis in RAS space)
%   Relative rotation angle between the two images
%   Per-axis dominant anatomical direction and oblique tilt angle
%
% Displays:
%   2×3 figure of central axial / coronal / sagittal slices (ref on top,
%   recon on bottom), rotated 180° so superior is up.
%
% Usage:
%   check_nifti_orientation(refPath, reconPath)
%   check_nifti_orientation(refPath, reconPath, figTitle)
%
% Inputs:
%   refPath    Path to reference NIfTI (.nii or .nii.gz), e.g. DICOM MPRAGE
%   reconPath  Path to reconstruction NIfTI to compare against the reference
%   figTitle   (optional) Supertitle for the figure window

if nargin < 3 || isempty(figTitle)
    [~, refName]   = fileparts(refPath);
    [~, reconName] = fileparts(reconPath);
    figTitle = sprintf('%s  vs  %s', refName, reconName);
end

assert(exist(refPath,   'file') == 2, 'Reference NIfTI not found:\n  %s', refPath);
assert(exist(reconPath, 'file') == 2, 'Reconstruction NIfTI not found:\n  %s', reconPath);

% -----------------------------------------------------------------------
% Step 1 — Load headers and extract affines
%
% niftiinfo stores Transform.T in row-major form such that
%   world = [vox, 1] * T
% We transpose to get the standard column-major form:
%   world = A * [vox; 1]
% where the columns of A(1:3,1:3) are the scaled axis directions in RAS mm.
% -----------------------------------------------------------------------
infoRef   = niftiinfo(refPath);
infoRecon = niftiinfo(reconPath);

A_ref   = infoRef.Transform.T';
A_recon = infoRecon.Transform.T';

% -----------------------------------------------------------------------
% Step 2 — Print affine summaries
% -----------------------------------------------------------------------
fprintf('\n===== NIfTI Orientation Check =====\n');
fprintf('  RAS: X = Left→Right, Y = Posterior→Anterior, Z = Inferior→Superior\n\n');

fprintf('--- Reference ---\n');
fprintf('  File       : %s\n', refPath);
fprintf('  Dimensions : %s  (voxels per axis)\n', mat2str(infoRef.ImageSize));
fprintf('  PixDim     : %s mm\n', mat2str(infoRef.PixelDimensions, 4));
fprintf('  Affine (top-left 3×3 = rotation+scaling, last col = origin in mm):\n');
disp(A_ref);

fprintf('--- Reconstruction ---\n');
fprintf('  File       : %s\n', reconPath);
fprintf('  Dimensions : %s  (voxels per axis)\n', mat2str(infoRecon.ImageSize));
fprintf('  PixDim     : %s mm\n', mat2str(infoRecon.PixelDimensions, 4));
fprintf('  Affine:\n');
disp(A_recon);

% -----------------------------------------------------------------------
% Step 3 — Rotation matrices and relative rotation
%
% Dividing each column of the 3×3 block by its voxel size strips the
% scaling and leaves a pure rotation matrix whose columns are unit vectors
% pointing along the anatomical axis each image dimension corresponds to.
%
% R_rel = R_recon * R_ref^{-1} is the rotation that maps the reference
% frame onto the reconstruction frame.  Identity = same orientation.
% The rotation angle is extracted from the trace: tr(R) = 1 + 2·cos(θ).
% -----------------------------------------------------------------------
R_ref   = A_ref(1:3,1:3)   ./ infoRef.PixelDimensions(1:3);
R_recon = A_recon(1:3,1:3) ./ infoRecon.PixelDimensions(1:3);

R_rel = R_recon / R_ref;

fprintf('--- Relative rotation (R_recon · R_ref⁻¹) ---\n');
fprintf('  Identity = same orientation; off-diagonal = misalignment\n');
disp(R_rel);

theta_rad = acos(max(-1, min(1, (trace(R_rel) - 1) / 2)));
fprintf('  Rotation angle: %.2f deg\n', rad2deg(theta_rad));
fprintf('  (0° = identical orientation)\n\n');

% -----------------------------------------------------------------------
% Step 4 — Per-axis direction table
%
% Each column of R is a unit vector [Rx, Ry, Rz] in RAS space expressing
% which anatomical direction that image axis runs along.
% Dominant component → main anatomical axis; tilt = acos(dominant component).
% -----------------------------------------------------------------------
labels  = {'axis-1 (col)', 'axis-2 (row)', 'axis-3 (slice)'};
rasAxes = {'R/L', 'A/P', 'S/I'};

fprintf('--- Axis directions in RAS space ---\n');
fprintf('  %-16s  %-36s  %-36s\n', 'Image axis', 'Reference', 'Reconstruction');

for ax = 1:3
    dir_ref   = R_ref(:, ax);
    dir_recon = R_recon(:, ax);

    [~, dom] = max(abs(dir_ref));
    tilt_deg = rad2deg(acos(abs(dir_ref(dom))));

    fprintf('  %-16s  [R=%+.3f A=%+.3f S=%+.3f]  [R=%+.3f A=%+.3f S=%+.3f]\n', ...
        labels{ax}, ...
        dir_ref(1),   dir_ref(2),   dir_ref(3), ...
        dir_recon(1), dir_recon(2), dir_recon(3));
    fprintf('  %-16s  -> mainly %s, %.1f° oblique tilt\n', '', rasAxes{dom}, tilt_deg);
end

% -----------------------------------------------------------------------
% Step 5 — Central-slice visual comparison
%
% Both volumes are reoriented to canonical RAS (dim1=R, dim2=A, dim3=S)
% using their affines before slicing.  This makes the subplot assignment
% (axial / coronal / sagittal) correct regardless of how the NIfTI stores
% data internally (which varies by scanner, sequence, and dcm2niix version).
%
% Display uses flipud(slice') — equivalent to 90° CW rotation — which
% produces a neurological view (right=right, anterior=top) for any
% canonical RAS volume.  4D volumes use the first frame.
% -----------------------------------------------------------------------
fprintf('\nLoading image data...\n');
volRef   = niftiread(refPath);
volRecon = niftiread(reconPath);

if ndims(volRef)   > 3, volRef   = volRef(:,:,:,1);   end
if ndims(volRecon) > 3, volRecon = volRecon(:,:,:,1); end

norm_vol = @(v) double(abs(v)) / double(max(abs(v(:))) + eps);
vRef   = norm_vol(reorient_to_ras(volRef,   A_ref));
vRecon = norm_vol(reorient_to_ras(volRecon, A_recon));

cRef   = round(size(vRef,   1:3) / 2);
cRecon = round(size(vRecon, 1:3) / 2);

views   = {'Axial', 'Coronal', 'Sagittal'};
slRef   = {squeeze(vRef(:,:,cRef(3))),     squeeze(vRef(:,cRef(2),:)),     squeeze(vRef(cRef(1),:,:))};
slRecon = {squeeze(vRecon(:,:,cRecon(3))), squeeze(vRecon(:,cRecon(2),:)), squeeze(vRecon(cRecon(1),:,:))};

figure('Name', 'Orientation Check', 'NumberTitle', 'off', 'Position', [100 100 1200 700]);
for v = 1:3
    subplot(2,3,v);
    imagesc(flipud(slRef{v}')); axis image off; colormap gray;
    title(sprintf('REF — %s', views{v}), 'FontSize', 9);

    subplot(2,3,v+3);
    imagesc(flipud(slRecon{v}')); axis image off; colormap gray;
    title(sprintf('RECON — %s', views{v}), 'FontSize', 9);
end
sgtitle(figTitle);
end


% -----------------------------------------------------------------------
% Local helper
% -----------------------------------------------------------------------

function vol_ras = reorient_to_ras(vol, affine)
% Permute and flip vol so dim1→R, dim2→A, dim3→S.
% Works for any NIfTI storage order; only the first 3 dims are touched.
R = affine(1:3, 1:3) ./ vecnorm(affine(1:3, 1:3), 2, 1);  % unit direction cosines
[~, perm] = max(abs(R), [], 2);   % perm(d) = data dim most aligned with RAS axis d
assert(numel(unique(perm)) == 3, ...
    'Cannot determine unique axis permutation — is this an oblique acquisition?');
vol_ras = permute(vol, perm);
for d = 1:3
    if R(d, perm(d)) < 0
        vol_ras = flip(vol_ras, d);
    end
end
end
