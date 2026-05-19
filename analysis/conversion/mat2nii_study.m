%% Convert one reconstruction .mat volume to NIfTI using a reference NIfTI/JSON
clc; clearvars; close all;

repoRoot = '/home/debi/jaime/repos/MR-EyeTrack';

matFile = fullfile(repoRoot, ...
    'data/study/sub-001/recon/woBin/x_steva_nIter_20_delta_1.000.mat');
refNifti = fullfile(repoRoot, ...
    'data/study/sub-001/dicom/csTFL_mp-rage_1mm-iso_CP_acc4.6_5_MR/sub-001.nii.gz');
refJson = fullfile(repoRoot, ...
    'data/study/sub-001/dicom/csTFL_mp-rage_1mm-iso_CP_acc4.6_5_MR/sub-001.json');
outputNiiGz = fullfile(repoRoot, ...
    'data/study/sub-001/recon/woBin/x_steva_nIter_20_delta_1.000.nii.gz');

assert(exist(matFile, 'file') == 2, 'MAT file not found: %s', matFile);
assert(exist(refNifti, 'file') == 2, 'Reference NIfTI not found: %s', refNifti);
assert(exist(refJson, 'file') == 2, 'Reference JSON not found: %s', refJson);

fprintf('Loading reconstruction volume from:\n  %s\n', matFile);
vol = load_recon_volume(matFile);
fprintf('Raw MAT size: %s\n', mat2str(size(vol)));

% Reorder the MAT array so the voxel data follows the template NIfTI
% orientation instead of preserving the original MAT storage order.
vol = reorient_volume_to_template(vol);
vol = single(abs(vol));
fprintf('Oriented MAT size: %s\n', mat2str(size(vol)));

refInfo = niftiinfo(refNifti);
refJsonInfo = jsondecode(fileread(refJson));

validate_reference_orientation(refInfo, refJsonInfo);

newInfo = build_output_header(refInfo, vol, outputNiiGz);

tmpNii = erase(outputNiiGz, '.gz');
if exist(tmpNii, 'file') == 2
    delete(tmpNii);
end
if exist(outputNiiGz, 'file') == 2
    delete(outputNiiGz);
end

fprintf('Writing NIfTI:\n  %s\n', tmpNii);
niftiwrite(vol, tmpNii, newInfo, 'Compressed', false);
gzip(tmpNii);
delete(tmpNii);
synchronize_qform_with_sform(outputNiiGz);

fprintf('Done.\nSaved:\n  %s\n', outputNiiGz);


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


function vol = reorient_volume_to_template(vol)
% Raw MAT axes: dim1 = −A (post direction), dim2 = −R (left direction), dim3 = +S (correct).
% Target (DICOM affine cols): dim1 = +R, dim2 = +A, dim3 = +S.
% permute([2,1,3]): swap dims 1↔2 → (−R, −A, +S)
% flip(1):          negate dim1   → (+R, −A, +S)
% flip(2):          negate dim2   → (+R, +A, +S)  ✓
vol = permute(vol, [2 1 3]);
vol = flip(vol, 1);
vol = flip(vol, 2);
end


function infoOut = build_output_header(refInfo, vol, outputNiiGz)
refAffine = refInfo.Transform.T';
R = refAffine(1:3, 1:3);
voxelSizes = vecnorm(R, 2, 1);
dirCos = R ./ voxelSizes;

% Copy the template orientation/origin directly so the generated NIfTI
% follows the same world-orientation convention as the DICOM-derived NIfTI.
newAffine = refAffine;

hdr = images.internal.nifti.niftiImage.niftiDefaultHeader(vol, true, 'NIfTI1');
hdr.pixdim(1) = 1;
hdr.pixdim(2:4) = single(voxelSizes);
hdr.xyzt_units = refInfo.raw.xyzt_units;
hdr.descrip = 'MR-EyeTrack MAT to NIfTI';
hdr.datatype = class_to_nifti_datatype(class(vol));
hdr.bitpix = int16(8 * bytes_per_voxel(class(vol)));
hdr.sform_code = int16(1);
hdr.qform_code = int16(1);
hdr.srow_x = single(newAffine(1, :));
hdr.srow_y = single(newAffine(2, :));
hdr.srow_z = single(newAffine(3, :));

q = rotmat_to_quat(dirCos);
hdr.quatern_b = single(q(2));
hdr.quatern_c = single(q(3));
hdr.quatern_d = single(q(4));
hdr.qoffset_x = single(newAffine(1, 4));
hdr.qoffset_y = single(newAffine(2, 4));
hdr.qoffset_z = single(newAffine(3, 4));
hdr.cal_min = single(min(vol(:)));
hdr.cal_max = single(max(vol(:)));
hdr.scl_slope = 1;
hdr.scl_inter = 0;

NV = images.internal.nifti.niftiImage(hdr);
infoOut = NV.simplifyStruct();
infoOut.raw = hdr;
infoOut.ImageSize = size(vol);
infoOut.PixelDimensions = voxelSizes;
infoOut.Datatype = class(vol);
infoOut.BitsPerPixel = 8 * bytes_per_voxel(infoOut.Datatype);
infoOut.Description = 'MR-EyeTrack MAT to NIfTI';
infoOut.Filemoddate = char(datetime('now', 'Format', 'dd-MMM-yyyy HH:mm:ss'));
infoOut.Filename = outputNiiGz;
infoOut.Transform = refInfo.Transform;
infoOut.Transform.T = newAffine';
end


function validate_reference_orientation(refInfo, refJsonInfo)
if ~isfield(refJsonInfo, 'ImageOrientationPatientDICOM')
    warning('Reference JSON has no ImageOrientationPatientDICOM field. Skipping orientation check.');
    return;
end

jsonIop = double(refJsonInfo.ImageOrientationPatientDICOM(:));
rowLps = jsonIop(1:3);
colLps = jsonIop(4:6);
normalLps = cross(rowLps, colLps);

% Convert DICOM LPS directions to NIfTI RAS directions for comparison.
LpsToRas = diag([-1 -1 1]);
expectedRas = [LpsToRas * rowLps, LpsToRas * colLps, LpsToRas * normalLps];
expectedRas = normalize_columns(expectedRas);

affineRas = normalize_columns(refInfo.Transform.T(1:3, 1:3)');
alignment = abs(affineRas' * expectedRas);
fprintf('Reference orientation check, best axis alignment = [%0.4f %0.4f %0.4f]\n', ...
    max(alignment(:, 1)), max(alignment(:, 2)), max(alignment(:, 3)));
end


function A = normalize_columns(A)
for iCol = 1:size(A, 2)
    A(:, iCol) = A(:, iCol) ./ norm(A(:, iCol));
end
end


function nBytes = bytes_per_voxel(className)
switch className
    case {'uint8', 'int8', 'logical'}
        nBytes = 1;
    case {'uint16', 'int16'}
        nBytes = 2;
    case {'uint32', 'int32', 'single'}
        nBytes = 4;
    case {'uint64', 'int64', 'double'}
        nBytes = 8;
    otherwise
        error('Unsupported datatype for NIfTI export: %s', className);
end
end


function niftiDatatype = class_to_nifti_datatype(className)
switch className
    case 'uint8'
        niftiDatatype = 2;
    case 'int16'
        niftiDatatype = 4;
    case 'int32'
        niftiDatatype = 8;
    case 'single'
        niftiDatatype = 16;
    case 'double'
        niftiDatatype = 64;
    case 'int8'
        niftiDatatype = 256;
    case 'uint16'
        niftiDatatype = 512;
    case 'uint32'
        niftiDatatype = 768;
    case 'int64'
        niftiDatatype = 1024;
    case 'uint64'
        niftiDatatype = 1280;
    otherwise
        error('Unsupported MATLAB class for NIfTI datatype mapping: %s', className);
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
