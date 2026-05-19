%% Convert one reconstruction .mat volume to NIfTI
% Uses the DICOM reference NIfTI for direction cosines (orientation
% validated) and the Twix FOV centre for the geometrically correct origin.
% The same permute/flip validated in mat2nii_test.m is applied to the
% volume.  This fixes the translation issue of mat2nii_test.m (which
% copied the origin from the DICOM reference instead of the acquisition).
clc; clearvars; close all;

repoRoot = '/home/debi/jaime/repos/MR-EyeTrack';
addpath(genpath(fullfile(repoRoot, 'analysis')));
addpath(genpath(fullfile(repoRoot, 'recon')));
addpath(genpath('/home/debi/MatTechLab/monalisa'));
addpath(genpath('/home/debi/yiwei/forclone/pulseq'));

subjectNum = 1;
subjectStr = sprintf('sub-%03d', subjectNum);

matFile = fullfile(repoRoot, ...
    'data/study', subjectStr, 'recon/woBin/x_steva_nIter_20_delta_1.000.mat');
twixFile = fullfile(repoRoot, ...
    'data/study', subjectStr, 'rawdata', [subjectStr '_T1wLIBRE.dat']);
seqFile = fullfile(repoRoot, ...
    'data/study/pulseq/yj_seq202_t1w_libre_main_TR8.0ms_TE3.6ms_swap1_FA6_RF2_mreye_track_trajPTP_44_1872_gdsp.seq');
refNifti = fullfile(repoRoot, ...
    'data/study', subjectStr, 'dicom/csTFL_mp-rage_1mm-iso_CP_acc4.6_5_MR/sub-001.nii.gz');
outputNiiGz = fullfile(repoRoot, ...
    'data/study', subjectStr, 'recon/woBin/x_steva_nIter_20_delta_1.000_twix.nii.gz');
twixMetaFile = fullfile(repoRoot, ...
    'data/study', subjectStr, 'recon/woBin/twix_orientation_metadata.mat');

assert(exist(matFile,    'file') == 2, 'MAT file not found: %s',           matFile);
assert(exist(twixFile,   'file') == 2, 'Twix/DAT file not found: %s',      twixFile);
assert(exist(seqFile,    'file') == 2, 'Pulseq file not found: %s',         seqFile);
assert(exist(refNifti,   'file') == 2, 'Reference NIfTI not found: %s',     refNifti);

fprintf('Loading reconstruction volume:\n  %s\n', matFile);
vol = load_recon_volume(matFile);
vol = single(abs(vol));
fprintf('Raw MAT size: %s\n', mat2str(size(vol)));

% Raw MAT axes: dim1 = −A (post direction), dim2 = −R (left direction), dim3 = +S (correct).
% Target (DICOM affine cols): dim1 = +R, dim2 = +A, dim3 = +S.
% permute([2,1,3]): swap dims 1↔2 → (−R, −A, +S)
% flip(1):          negate dim1   → (+R, −A, +S)
% flip(2):          negate dim2   → (+R, +A, +S)  ✓
vol = permute(vol, [2 1 3]);
vol = flip(vol, 1);
vol = flip(vol, 2);
fprintf('Oriented size: %s\n', mat2str(size(vol)));

useCachedMetadata = false;
if exist(twixMetaFile, 'file') == 2
    cachedVars = who('-file', twixMetaFile);
    useCachedMetadata = all(ismember({'sa', 'FOV', 'TR'}, cachedVars));
end

if useCachedMetadata
    fprintf('Loading cached Twix orientation metadata:\n  %s\n', twixMetaFile);
    meta = load(twixMetaFile, 'sa', 'FOV', 'TR');
    sa = meta.sa;
    FOV = meta.FOV;
    TR = meta.TR;
else
    if exist(twixMetaFile, 'file') == 2
        fprintf('Cached metadata is missing FOV/TR; regenerating:\n  %s\n', twixMetaFile);
    end
    fprintf('Creating reader and reading raw data to match the reconstruction path...\n');
    seqParams = extract_seq_params(seqFile);
    reader = createRawDataReader(twixFile, true);
    reader.acquisitionParams.nShot_off = 14;
    reader.acquisitionParams.traj_type = 'pulseq';
    reader.acquisitionParams.pulseqTrajFile_name = seqFile;
    if isfield(seqParams, 'nshot')
        reader.acquisitionParams.nShot = seqParams.nshot;
    end
    if isfield(seqParams, 'nseg')
        reader.acquisitionParams.nSeg = seqParams.nseg;
    end
    y_tot = reader.readRawData(true, true); %#ok<NASGU>

    fprintf('Loading Twix header geometry:\n  %s\n', twixFile);
    twix = mapVBVD_JH_for_monalisa(twixFile);
    if iscell(twix)
        twix = twix{end};
    end
    sa = twix.hdr.MeasYaps.sSliceArray.asSlice{1};
    FOV = double(reader.acquisitionParams.FoV);
    TR = get_repetition_time_seconds(twix.hdr);
    save(twixMetaFile, 'sa', 'FOV', 'TR', '-v7.3');
    fprintf('Saved Twix orientation metadata:\n  %s\n', twixMetaFile);
end

if numel(FOV) == 1
    FOV = repmat(FOV, 1, 3);
end

% -------------------------------------------------------------------------
% Build affine: copy the full DICOM affine (directions + voxel sizes,
% empirically validated to match the permuted/flipped MAT array), then
% replace only the origin column with the Twix FOV centre.
%
% mat2nii_test.m copies the DICOM origin, which is the MPRAGE FOV centre,
% not the LIBRE acquisition centre.  Replacing only the origin preserves
% the validated orientation while fixing the translation.
% -------------------------------------------------------------------------

refInfo  = niftiinfo(refNifti);
refAffine = refInfo.Transform.T';   % 4×4, RAS convention

% FOV centre from Twix header (LPS → RAS).
posLps = get_twix_vector(sa.sPosition, {'dSag', 'dCor', 'dTra'});
posRas = [-posLps(1); -posLps(2); posLps(3)];

% Use the DICOM affine for directions and voxel sizes (validated orientation).
% Replace only the origin column with the Twix FOV centre.
R_dicom = refAffine(1:3, 1:3);
Nv = size(vol);
originRas = posRas ...
    - R_dicom(:,1) * (Nv(1)-1)/2 ...
    - R_dicom(:,2) * (Nv(2)-1)/2 ...
    - R_dicom(:,3) * (Nv(3)-1)/2;

affine = refAffine;
affine(1:3, 4) = originRas;

fprintf('DICOM origin:  [%.2f  %.2f  %.2f]\n', refAffine(1:3,4));
fprintf('Twix origin:   [%.2f  %.2f  %.2f]\n', originRas);
fprintf('Offset (mm):   [%.2f  %.2f  %.2f]\n', (originRas - refAffine(1:3,4))');

info = build_nifti_info_from_affine(vol, affine, TR, outputNiiGz);

tmpNii = erase(outputNiiGz, '.gz');
if exist(tmpNii, 'file') == 2
    delete(tmpNii);
end
if exist(outputNiiGz, 'file') == 2
    delete(outputNiiGz);
end

fprintf('Writing NIfTI:\n  %s\n', outputNiiGz);
niftiwrite(vol, tmpNii, info, 'Compressed', false);
gzip(tmpNii);
delete(tmpNii);
synchronize_qform_with_sform(outputNiiGz);
fprintf('Done.\n');


function vol = load_recon_volume(matFile)
data = load(matFile);
candidateFields = {'x', 'x0', 'xrms'};

for iField = 1:numel(candidateFields)
    fieldName = candidateFields{iField};
    if isfield(data, fieldName)
        vol = data.(fieldName);
        if iscell(vol)
            vol = vol{1};
        end
        if ndims(vol) ~= 3
            error('Field "%s" is not a 3D volume.', fieldName);
        end
        return;
    end
end

error('No supported reconstruction field found in %s.', matFile);
end


function v = get_twix_vector(s, names)
v = zeros(3, 1);
for i = 1:3
    if isfield(s, names{i})
        v(i) = double(s.(names{i}));
    end
end
end


function info = build_nifti_info_from_affine(vol, affine, TR, outputNiiGz)
voxelSize = vecnorm(affine(1:3, 1:3), 2, 1);
dirCos = affine(1:3, 1:3) ./ voxelSize;

hdr = images.internal.nifti.niftiImage.niftiDefaultHeader(vol, true, 'NIfTI1');
hdr.pixdim = single([1 voxelSize max(TR, eps) 0 0 0]);
hdr.xyzt_units = uint8(2) + uint8(8); % mm + sec
hdr.descrip = 'MR-EyeTrack Twix origin';
hdr.datatype = int16(16); % single
hdr.bitpix = int16(32);
hdr.sform_code = int16(1);
hdr.qform_code = int16(1);
hdr.srow_x = single(affine(1, :));
hdr.srow_y = single(affine(2, :));
hdr.srow_z = single(affine(3, :));

q = rotmat_to_quat(dirCos);
hdr.quatern_b = single(q(2));
hdr.quatern_c = single(q(3));
hdr.quatern_d = single(q(4));
hdr.qoffset_x = single(affine(1, 4));
hdr.qoffset_y = single(affine(2, 4));
hdr.qoffset_z = single(affine(3, 4));
hdr.cal_min = single(min(vol(:)));
hdr.cal_max = single(max(vol(:)));
hdr.scl_slope = 1;
hdr.scl_inter = 0;

NV = images.internal.nifti.niftiImage(hdr);
info = NV.simplifyStruct();
info.raw = hdr;
info.ImageSize = size(vol);
info.PixelDimensions = voxelSize;
info.Datatype = class(vol);
info.BitsPerPixel = 32;
info.Description = 'MR-EyeTrack Twix origin';
info.Filemoddate = char(datetime('now', 'Format', 'dd-MMM-yyyy HH:mm:ss'));
info.Filename = outputNiiGz;
info.Transform = affine3d(affine');
end


function TR = get_repetition_time_seconds(hdr)
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
traceR = trace(R);

if traceR > 0
    s = 2 * sqrt(traceR + 1);
    qw = 0.25 * s;
    qx = (R(3,2) - R(2,3)) / s;
    qy = (R(1,3) - R(3,1)) / s;
    qz = (R(2,1) - R(1,2)) / s;
elseif (R(1,1) > R(2,2)) && (R(1,1) > R(3,3))
    s = 2 * sqrt(1 + R(1,1) - R(2,2) - R(3,3));
    qw = (R(3,2) - R(2,3)) / s;
    qx = 0.25 * s;
    qy = (R(1,2) + R(2,1)) / s;
    qz = (R(1,3) + R(3,1)) / s;
elseif R(2,2) > R(3,3)
    s = 2 * sqrt(1 + R(2,2) - R(1,1) - R(3,3));
    qw = (R(1,3) - R(3,1)) / s;
    qx = (R(1,2) + R(2,1)) / s;
    qy = 0.25 * s;
    qz = (R(2,3) + R(3,2)) / s;
else
    s = 2 * sqrt(1 + R(3,3) - R(1,1) - R(2,2));
    qw = (R(2,1) - R(1,2)) / s;
    qx = (R(1,3) + R(3,1)) / s;
    qy = (R(2,3) + R(3,2)) / s;
    qz = 0.25 * s;
end

q = [qw qx qy qz];
q = q ./ norm(q);
if q(1) < 0
    q = -q;
end
end


function synchronize_qform_with_sform(niftiPath)
pyFile = fullfile(tempdir, 'mr_eyetrack_fix_qform.py');
fid = fopen(pyFile, 'w');
assert(fid ~= -1, 'Could not create temporary Python helper: %s', pyFile);

cleanupObj = onCleanup(@() cleanup_temp_file(fid, pyFile));

fprintf(fid, 'import sys\n');
fprintf(fid, 'import numpy as np\n');
fprintf(fid, 'import nibabel as nib\n');
fprintf(fid, 'path = sys.argv[1]\n');
fprintf(fid, 'img = nib.load(path)\n');
fprintf(fid, 'data = np.asanyarray(img.dataobj)\n');
fprintf(fid, 'aff = img.get_sform()\n');
fprintf(fid, 'hdr = img.header.copy()\n');
fprintf(fid, 'hdr.set_sform(aff, code=1)\n');
fprintf(fid, 'hdr.set_qform(aff, code=1)\n');
fprintf(fid, 'out = nib.Nifti1Image(data, aff, header=hdr)\n');
fprintf(fid, 'out.set_sform(aff, code=1)\n');
fprintf(fid, 'out.set_qform(aff, code=1)\n');
fprintf(fid, 'nib.save(out, path)\n');
fclose(fid);

pythonExe = '/home/debi/miniconda3/envs/mreyetrack/bin/python3';
cmd = sprintf('"%s" "%s" "%s"', pythonExe, pyFile, niftiPath);
[status, cmdout] = system(cmd);

if status ~= 0
    fprintf('Skipped qform/sform sync helper: %s\n', strtrim(cmdout));
else
    fprintf('Synchronized qform and sform using nibabel.\n');
end
end


function cleanup_temp_file(fid, pyFile)
if fid > 0
    try
        fclose(fid);
    catch
    end
end
if exist(pyFile, 'file') == 2
    delete(pyFile);
end
end
