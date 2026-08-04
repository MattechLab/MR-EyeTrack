function mat2nii_twix(matFile, twixFile, seqFile, refNifti, outputNiiGz, twixMetaFile, reorientFcn)
%MAT2NII_TWIX  Convert a reconstruction .mat volume to NIfTI.
%
% Orientation strategy — three components from three sources:
%   Direction cosines  DICOM MPRAGE reference NIfTI (empirically validated
%                      to match the raw MAT axis layout after reorientation)
%   Voxel size         Twix FOV / reconstruction matrix dimensions
%                      (acquisition-native; independent of the MPRAGE)
%   Origin             Twix sPosition (LIBRE FOV centre) shifted to the
%                      [0,0,0] corner voxel in RAS mm
%
% Usage:
%   mat2nii_twix(matFile, twixFile, seqFile, refNifti, outputNiiGz)
%   mat2nii_twix(..., twixMetaFile)
%   mat2nii_twix(..., twixMetaFile, reorientFcn)
%
% Inputs:
%   matFile       Path to .mat reconstruction file (field x, x0, or xrms)
%   twixFile      Path to raw .dat Twix acquisition file
%   seqFile       Path to pulseq .seq file (provides nShot and nSeg for the
%                 reader, needed to extract FoV from acquisitionParams)
%   refNifti      Path to DICOM-derived MPRAGE NIfTI (.nii.gz)
%   outputNiiGz   Destination path (.nii.gz)
%   twixMetaFile  (optional) Path to a .mat cache for Twix metadata
%                 (sa, FOV, TR).  Avoids re-reading the full raw .dat on
%                 repeated calls — reading raw data is slow (~minutes).
%   reorientFcn   (optional) Function handle that reorders the raw MAT volume
%                 to match the REFERENCE's axis order before the affine is
%                 applied.  Each pipeline / reconstruction has its own
%                 raw-axis convention; pass [] or omit to skip reorientation.
%
%                 Note this is the reference's storage order, which is only
%                 (+R,+A,+S) when the reference itself is stored RAS.  The
%                 output of reorientFcn is written verbatim and the direction
%                 cosines come from the reference, so the two must agree.
%                 MR-EyeTrack MPRAGEs are ('R','A','S'), but Yiwei's are
%                 ('P','I','L') — there the written array runs A-P, S-I, L-R,
%                 so the left-right axis is dim 3, not dim 1.
%
%                 Known conventions (use these in the calling test script):
%                   MR-EyeTrack LIBRE  raw=(−A,−R,+S):
%                     @(v) flip(flip(permute(v,[2 1 3 (4:ndims(v))]),1),2)
%                   Yannick AudioBOLD  raw=(−R,−A,+S):
%                     @(v) flip(flip(v,1),2)
%                   Yiwei 2.0 T1w LIBRE (MID00030) raw=(−S,+R,−A):
%                     @(v) flip(flip(permute(v,[2 3 1 (4:ndims(v))]),2),3)
%                   Yiwei 2.0 T2w LIBRE (MID00025) raw=(−R,−S,−A):
%                     @(v) flip(permute(v,[1 3 2 (4:ndims(v))]),2)
%
%                 Visual checks cannot validate the left-right sign: a mirrored
%                 brain looks plausible and a mirror changes neither axis
%                 identity nor rotation.  The T2w entry above was corrected
%                 after analysis/test_orientation/tissue_check/lr_flip_test.py
%                 found the old one mirrored.

if nargin < 6, twixMetaFile = ''; end
if nargin < 7, reorientFcn  = []; end

assert(exist(matFile,  'file') == 2, 'MAT file not found: %s',         matFile);
assert(exist(twixFile, 'file') == 2, 'Twix file not found: %s',        twixFile);
assert(exist(seqFile,  'file') == 2, 'Seq file not found: %s',         seqFile);
assert(exist(refNifti, 'file') == 2, 'Reference NIfTI not found: %s',  refNifti);

% -----------------------------------------------------------------------
% Step 1 — Load and reorient the reconstruction volume
%
% Raw MAT axes vary by reconstruction pipeline.  Pass reorientFcn in the
% calling script to map them to (+R, +A, +S) before the affine is applied.
% See the function header for known conventions.
% -----------------------------------------------------------------------
fprintf('Loading reconstruction volume:\n  %s\n', matFile);
vol = load_recon_volume(matFile);
vol = single(abs(vol));
fprintf('Raw MAT size: %s\n', mat2str(size(vol)));

if ~isempty(reorientFcn)
    vol = reorientFcn(vol);
end

fprintf('Oriented size: %s\n', mat2str(size(vol)));

% -----------------------------------------------------------------------
% Step 2 — Load Twix geometry metadata (with optional caching)
%
% Reading the full raw .dat is slow; cache sa/FOV/TR to a .mat sidecar so
% repeated conversions skip the raw-data read.
% -----------------------------------------------------------------------
[sa, FOV, TR] = load_twix_metadata(twixFile, seqFile, twixMetaFile);
if numel(FOV) == 1
    FOV = repmat(FOV, 1, 3);   % expand scalar to [Fx Fy Fz]
end

% -----------------------------------------------------------------------
% Step 3 — Build the NIfTI affine
%
% (a) Direction cosines: strip the MPRAGE voxel scaling from its affine
%     columns to recover unit-length orientation vectors.  These are
%     validated to align with the reoriented MAT array.
%
% (b) Voxel size: FOV (mm) / reconstruction grid dimensions.  Uses the
%     actual acquisition FOV, not the MPRAGE cell size.
%
% (c) Origin: Twix sPosition gives the LIBRE FOV centre in LPS (mm).
%     Convert to RAS, then shift by half the FOV to reach the [0,0,0]
%     corner voxel (using (N−1)/2 so the last voxel lands at +FOV edge).
% -----------------------------------------------------------------------
refInfo   = niftiinfo(refNifti);
refAffine = refInfo.Transform.T';          % 4×4 column-major RAS affine

% (a) Direction cosines from DICOM
R_dicom = refAffine(1:3, 1:3);
dirCos  = R_dicom ./ vecnorm(R_dicom, 2, 1);

% (b) Voxel size from Twix
Nv        = size(vol);
voxelSize = FOV(:)' ./ double(Nv(1:3));   % [vx vy vz] in mm

% (c) Origin from Twix centre position
posLps    = get_twix_vector(sa.sPosition, {'dSag', 'dCor', 'dTra'});
posRas    = [-posLps(1); -posLps(2); posLps(3)];   % LPS → RAS
originRas = posRas ...
    - dirCos(:,1) * voxelSize(1) * (Nv(1)-1)/2 ...
    - dirCos(:,2) * voxelSize(2) * (Nv(2)-1)/2 ...
    - dirCos(:,3) * voxelSize(3) * (Nv(3)-1)/2;

affine = refAffine;
affine(1:3, 1:3) = dirCos .* voxelSize;
affine(1:3, 4)   = originRas;

fprintf('DICOM voxel size: [%.4f  %.4f  %.4f] mm\n', vecnorm(R_dicom, 2, 1));
fprintf('Twix voxel size:  [%.4f  %.4f  %.4f] mm\n', voxelSize);
fprintf('Twix origin:      [%.2f  %.2f  %.2f] mm\n', originRas);

% -----------------------------------------------------------------------
% Step 4 — Write NIfTI
%
% After writing, nibabel synchronizes qform and sform to the same affine.
% MATLAB's niftiwrite can leave them inconsistent when the rotation matrix
% determinant sign is negative (left-handed frame).
% -----------------------------------------------------------------------
info = build_nifti_info_from_affine(vol, affine, TR, outputNiiGz);

if exist(outputNiiGz, 'file') == 2
    delete(outputNiiGz);
end
fprintf('Writing NIfTI:\n  %s\n', outputNiiGz);
niftiwrite(vol, outputNiiGz, info, 'Compressed', true);
synchronize_qform_with_sform(outputNiiGz);
fprintf('Done.\n');
end


% =======================================================================
% Local helpers
% =======================================================================

function vol = load_recon_volume(matFile)
% Load the reconstruction array from whichever field name was used.
% Cell arrays (e.g. BOLD timeseries) are stacked into a 4D array (Nx×Ny×Nz×Nt).
data = load(matFile);
candidateFields = {'x', 'x0', 'xrms'};
for iField = 1:numel(candidateFields)
    fieldName = candidateFields{iField};
    if isfield(data, fieldName)
        vol = data.(fieldName);
        if iscell(vol)
            vol = cat(4, vol{:});   % 1×Nt cell of 3D → 4D array
        end
        if ndims(vol) < 3 || ndims(vol) > 4
            error('Field "%s" must be a 3D or 4D array.', fieldName);
        end
        return;
    end
end
error('No supported reconstruction field found in %s.', matFile);
end


function [sa, FOV, TR] = load_twix_metadata(twixFile, seqFile, twixMetaFile)
% Return slice geometry (sa), FOV (mm), and TR (s) from the Twix header.
% If twixMetaFile is provided and contains the required variables, load from
% cache instead of re-reading the raw data.
if ~isempty(twixMetaFile) && exist(twixMetaFile, 'file') == 2
    cachedVars = who('-file', twixMetaFile);
    if all(ismember({'sa', 'FOV', 'TR'}, cachedVars))
        fprintf('Loading cached Twix metadata:\n  %s\n', twixMetaFile);
        meta = load(twixMetaFile, 'sa', 'FOV', 'TR');
        sa = meta.sa;  FOV = meta.FOV;  TR = meta.TR;
        return;
    end
    fprintf('Cache missing fields; regenerating:\n  %s\n', twixMetaFile);
end

% Full raw-data read — slow but required to populate acquisitionParams.FoV
fprintf('Reading Twix raw data (slow; will be cached):\n  %s\n', twixFile);
seqParams = extract_seq_params(seqFile);
reader = createRawDataReader(twixFile, true);
reader.acquisitionParams.traj_type          = 'pulseq';
reader.acquisitionParams.pulseqTrajFile_name = seqFile;
if isfield(seqParams, 'nshot'), reader.acquisitionParams.nShot = seqParams.nshot; end
if isfield(seqParams, 'nseg'),  reader.acquisitionParams.nSeg  = seqParams.nseg;  end

twix = mapVBVD_JH_for_monalisa(twixFile);
if iscell(twix), twix = twix{end}; end
sa  = twix.hdr.MeasYaps.sSliceArray.asSlice{1};
FOV = double(reader.acquisitionParams.FoV);
TR  = get_repetition_time_seconds(twix.hdr);

if ~isempty(twixMetaFile)
    save(twixMetaFile, 'sa', 'FOV', 'TR', '-v7.3');
    fprintf('Saved Twix metadata cache:\n  %s\n', twixMetaFile);
end
end


function v = get_twix_vector(s, names)
% Read a 3-element vector from a Twix struct, defaulting missing fields to 0.
v = zeros(3, 1);
for i = 1:3
    if isfield(s, names{i}), v(i) = double(s.(names{i})); end
end
end


function info = build_nifti_info_from_affine(vol, affine, TR, outputNiiGz)
% Build the MATLAB NIfTI info struct from a 4×4 RAS affine.
% Follows the same pattern as save_nifti.m (Yannick) which is known to work
% for both 3D and 4D: build hdr with niftiDefaultHeader, set fields as double
% (no single casts), call simplifyStruct, add raw — then let niftiwrite
% reconstruct the header internally via niftiImage(info).
voxelSize = vecnorm(affine(1:3, 1:3), 2, 1);
dirCos    = affine(1:3, 1:3) ./ voxelSize;

hdr = images.internal.nifti.niftiImage.niftiDefaultHeader(vol, true, 'NIfTI1');
% pixdim: [qfac, dx, dy, dz, dt] — 5 elements, all double (no single cast)
hdr.pixdim     = [1, voxelSize, max(TR, eps)];
hdr.xyzt_units = uint8(2) + uint8(8);   % mm + sec
hdr.descrip    = 'MR-EyeTrack Twix origin';
hdr.datatype   = int16(16);             % single float
hdr.bitpix     = int16(32);
hdr.sform_code = int16(1);
hdr.qform_code = int16(1);
hdr.srow_x     = affine(1, :);
hdr.srow_y     = affine(2, :);
hdr.srow_z     = affine(3, :);

q = rotmat_to_quat(dirCos);
hdr.quatern_b  = q(2);
hdr.quatern_c  = q(3);
hdr.quatern_d  = q(4);
hdr.qoffset_x  = affine(1, 4);
hdr.qoffset_y  = affine(2, 4);
hdr.qoffset_z  = affine(3, 4);
hdr.cal_min    = double(min(vol(:)));
hdr.cal_max    = double(max(vol(:)));
hdr.scl_slope  = 1;
hdr.scl_inter  = 0;

NV   = images.internal.nifti.niftiImage(hdr);
info = NV.simplifyStruct();
info.raw         = hdr;
info.Filename    = outputNiiGz;
info.Description = 'MR-EyeTrack Twix origin';
info.Filemoddate = char(datetime('now', 'Format', 'dd-MMM-yyyy HH:mm:ss'));
end


function TR = get_repetition_time_seconds(hdr)
% Extract TR from the Twix MeasYaps header; fall back to 1 s if absent.
TR = 1;
try
    if isfield(hdr.MeasYaps, 'alTR')
        TR = double(hdr.MeasYaps.alTR{1}) * 1e-6;
    elseif isfield(hdr.Dicom, 'alTR')
        TR = double(hdr.Dicom.alTR{1}) * 1e-6;
    end
catch
    TR = 1;
end
end


function q = rotmat_to_quat(R)
% Convert a 3×3 rotation matrix to a unit quaternion [qw qx qy qz].
% Uses Shepperd's method to avoid division by near-zero values.
traceR = trace(R);
if traceR > 0
    s  = 2 * sqrt(traceR + 1);
    qw = 0.25 * s;
    qx = (R(3,2) - R(2,3)) / s;
    qy = (R(1,3) - R(3,1)) / s;
    qz = (R(2,1) - R(1,2)) / s;
elseif (R(1,1) > R(2,2)) && (R(1,1) > R(3,3))
    s  = 2 * sqrt(1 + R(1,1) - R(2,2) - R(3,3));
    qw = (R(3,2) - R(2,3)) / s;
    qx = 0.25 * s;
    qy = (R(1,2) + R(2,1)) / s;
    qz = (R(1,3) + R(3,1)) / s;
elseif R(2,2) > R(3,3)
    s  = 2 * sqrt(1 + R(2,2) - R(1,1) - R(3,3));
    qw = (R(1,3) - R(3,1)) / s;
    qx = (R(1,2) + R(2,1)) / s;
    qy = 0.25 * s;
    qz = (R(2,3) + R(3,2)) / s;
else
    s  = 2 * sqrt(1 + R(3,3) - R(1,1) - R(2,2));
    qw = (R(2,1) - R(1,2)) / s;
    qx = (R(1,3) + R(3,1)) / s;
    qy = (R(2,3) + R(3,2)) / s;
    qz = 0.25 * s;
end
q = [qw qx qy qz];
q = q ./ norm(q);
if q(1) < 0, q = -q; end
end


function synchronize_qform_with_sform(niftiPath)
% Use nibabel to overwrite the NIfTI qform with the sform affine.
% This ensures both forms are byte-identical, which some viewers require.
pyFile = fullfile(tempdir, 'mr_eyetrack_fix_qform.py');
fid    = fopen(pyFile, 'w');
assert(fid ~= -1, 'Could not create temporary Python helper: %s', pyFile);
cleanupObj = onCleanup(@() cleanup_temp_file(fid, pyFile));

fprintf(fid, 'import sys, numpy as np, nibabel as nib\n');
fprintf(fid, 'path = sys.argv[1]\n');
fprintf(fid, 'img  = nib.load(path)\n');
fprintf(fid, 'aff  = img.get_sform()\n');
fprintf(fid, 'hdr  = img.header.copy()\n');
fprintf(fid, 'hdr.set_sform(aff, code=1)\n');
fprintf(fid, 'hdr.set_qform(aff, code=1)\n');
fprintf(fid, 'out  = nib.Nifti1Image(np.asanyarray(img.dataobj), aff, header=hdr)\n');
fprintf(fid, 'nib.save(out, path)\n');
fclose(fid);

pythonExe = '/home/debi/miniconda3/envs/mreyetrack/bin/python3';
[status, cmdout] = system(sprintf('"%s" "%s" "%s"', pythonExe, pyFile, niftiPath));
if status ~= 0
    fprintf('Skipped qform/sform sync: %s\n', strtrim(cmdout));
else
    fprintf('Synchronized qform and sform.\n');
end
end


function cleanup_temp_file(fid, pyFile)
try, fclose(fid); catch, end
if exist(pyFile, 'file') == 2, delete(pyFile); end
end
