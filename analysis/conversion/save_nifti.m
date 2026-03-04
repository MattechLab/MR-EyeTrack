function save_nifti(twixFile, dK_u, vol4d, TR, niftiFile)
%SAVE_NIFTI Save data to NIfTI with affine extracted from Twix file
%   Bovier Yannick

arguments (Input)
    twixFile
    dK_u (1,3) double                                   % In 1/mm
    vol4d {mustBeA(vol4d, ["double","single","cell"])}
    TR (1,1) {mustBePositive}                           % In sec
    niftiFile
end

fprintf("This is the updated save_nifti()\n");


% Load twix
twix = mapVBVD_JH_for_monalisa(twixFile);
if iscell(twix), twix = twix{end}; end
sa = twix.hdr.MeasYaps.sSliceArray.asSlice{1};

% Validate input type
if iscell(vol4d)
    % Expect a 1D cell array of 3D volumes
    if ~isvector(vol4d)
        error("vol4d as cell must be a 1D cell vector (1×Nt or Nt×1).");
    end
    if isempty(vol4d)
        error("vol4d cell is empty.");
    end

    % Validate each cell
    sz0 = size(vol4d{1});
    if numel(sz0) < 3
        error("Each cell element must be a 3D array (Nx×Ny×Nz).");
    end
    sz0 = sz0(1:3);

    for t = 1:numel(vol4d)
        v = vol4d{t};
        if ~isnumeric(v)
            error("vol4d{%d} is not numeric.", t);
        end
        if numel(size(v)) < 3
            error("vol4d{%d} is not 3D (Nx×Ny×Nz).", t);
        end
        if ~isequal(size(v,1), sz0(1)) || ~isequal(size(v,2), sz0(2)) || ~isequal(size(v,3), sz0(3))
            error("vol4d{%d} size mismatch. Expected [%d %d %d], got [%d %d %d].", t, sz0(1), sz0(2), sz0(3), size(v,1), size(v,2), size(v,3));
        end
    end

    % Stack into 4D (Nx×Ny×Nz×Nt)
    vol4d = cat(4, vol4d{:});
    fprintf("Converted cell array of %d volumes into 4D array [%s].\n", numel(vol4d(1,1,1,:)), num2str(size(vol4d)));

elseif isnumeric(vol4d)
    % Numeric input: accept 4D (or optionally 3D)
    nd = ndims(vol4d);
    if nd == 3
        % Optional: interpret as single timepoint
        vol4d = reshape(vol4d, size(vol4d,1), size(vol4d,2), size(vol4d,3), 1);
        fprintf("Input was 3D; promoted to 4D with Nt=1.\n");
    elseif nd ~= 4
        error("vol4d must be 3D numeric (Nx×Ny×Nz), 4D numeric (Nx×Ny×Nz×Nt) or a cell array of 3D volumes.");
    end
else
    error("vol4d must be numeric or a cell array.");
end

% NIfTI does not support complex valued data
if ~isreal(vol4d)
    fprintf("vol4d is complex, taking the magnitude\n");
    vol4d = abs(vol4d);
end

% Get matrix size
vol4d = single(vol4d);
matrix_size = size(vol4d, 1:3);  % Nx Ny Nz
fprintf("matrix_size = [%s] voxel\n", num2str(matrix_size));

% Get FOV and Twix FOV
FOV = (1 ./ dK_u);  % mm
fprintf("FOV = [%s] mm\n", num2str(FOV));

% Get voxel_size
voxel_size = FOV ./ matrix_size;  % mm
fprintf("voxel_size = [%s] mm\n", num2str(voxel_size));

%% Conventions
% Twix uses Siemens DICOM LPS
% +X = Right → Left
% +Y = Anterior → Posterior
% +Z = Foot → Head
%
% NIfTI expects RAS
% +X = Left → Right
% +Y = Posterior → Anterior
% +Z = Inferior → Superior
LPS2RAS = diag([-1 -1  1  1]);

%% Slice normal
if isfield(sa.sNormal, "dSag"), dSag = sa.sNormal.dSag; else, dSag = 0; end
if isfield(sa.sNormal, "dCor"), dCor = sa.sNormal.dCor; else, dCor = 0; end
if isfield(sa.sNormal, "dTra"), dTra = sa.sNormal.dTra; else, dTra = 0; end

n = [dSag; dCor; dTra];
n = n / norm(n);

% Pick a reference axis that is not collinear
ref = [0; 0; 1];
if abs(dot(n, ref)) > 0.9
    ref = [0; 1; 0];
end

% Construct orthonormal basis
x_dir = cross(ref, n);
x_dir = x_dir / norm(x_dir);
y_dir = cross(n, x_dir);

R_lps = [x_dir, y_dir, n];   % 3×3

%% Voxel size
pixdim = diag(voxel_size);

%% Volume center
% vol_center_vox = (matrix_size(:) + 1) / 2;
% vol_center_mm  = pixdim * (vol_center_vox - 1);

%% Slice position (LPS, mm)
if isfield(sa.sPosition, "dSag"), dSag = sa.sPosition.dSag; else, dSag = 0; end
if isfield(sa.sPosition, "dCor"), dCor = sa.sPosition.dCor; else, dCor = 0; end
if isfield(sa.sPosition, "dTra"), dTra = sa.sPosition.dTra; else, dTra = 0; end
pos_lps = [dSag;
           dCor;
           dTra];

%% Translate so that voxel center maps to slice center
% T_lps = pos_lps - R_lps * vol_center_mm;

T_lps = pos_lps ...
        - R_lps(:,1) * (FOV(1)/2) ...
        - R_lps(:,2) * (FOV(2)/2) ...
        - R_lps(:,3) * (FOV(3)/2);


%% Full affine
A_lps = eye(4);
A_lps(1:3,1:3) = R_lps * pixdim;
A_lps(1:3,4)   = T_lps;

%% LPS -> RAS
A_ras = LPS2RAS * A_lps;

%% Q-form (Rotation and translation only, not scaling)
% rot_mat = A_ras(1:3,1:3);
% quat = rotm2quat(rot_mat);
trans_vec = A_ras(1:3,4);

rot_mat = LPS2RAS(1:3,1:3) * R_lps;   % 3x3
detR   = det(rot_mat);
qfac   = sign(detR);    % +1 or -1; if detR==0, set qfac=+1
if qfac == 0, qfac = 1; end
quat = rotm2quat(rot_mat);

%% Create header
hdr = images.internal.nifti.niftiImage.niftiDefaultHeader(vol4d, true, 'NIfTI1');

% modify
hdr.pixdim = [qfac, voxel_size, TR];
hdr.descrip = "BOLD timeseries";
hdr.xyzt_units = uint8(2) + uint8(8); % mm + sec
% hdr.dim_info =

% Set sform
hdr.sform_code = 1;
hdr.srow_x = A_ras(1, :);
hdr.srow_y = A_ras(2, :);
hdr.srow_z = A_ras(3, :);

% Set qform
hdr.qform_code = 1;
hdr.quatern_b = quat(2);
hdr.quatern_c = quat(3);
hdr.quatern_d = quat(4);
hdr.qoffset_x = trans_vec(1);  % mm
hdr.qoffset_y = trans_vec(2);  % mm
hdr.qoffset_z = trans_vec(3);  % mm

% Min/max
hdr.cal_min = min(vol4d(:));
hdr.cal_max = max(vol4d(:));

% No scaling
hdr.scl_slope = 1;
hdr.scl_inter = 0;

% simplify header
NV = images.internal.nifti.niftiImage(hdr);
nifti_header = NV.simplifyStruct();
nifti_header.raw = hdr;

% Save
niftiwrite(vol4d, niftiFile, nifti_header, 'Compressed', endsWith(niftiFile, ".gz"));
fprintf("Done saving NIfTI file: %s\n", niftiFile);

end