%% Convert all reconstruction .mat volumes to NIfTI using a reference NIfTI/JSON
clc; clearvars; close all;

subject_num = 3;   % <-- set subject number here

repoRoot = '/home/debi/jaime/repos/MR-EyeTrack';
subID    = sprintf('sub-%03d', subject_num);
subDir   = fullfile(repoRoot, 'data/study', subID);

refNifti = fullfile(subDir, 'dicom/csTFL_mp-rage_1mm-iso_CP_acc4.6_5_MR', [subID '.nii.gz']);
refJson  = fullfile(subDir, 'dicom/csTFL_mp-rage_1mm-iso_CP_acc4.6_5_MR', [subID '.json']);

assert(exist(refNifti, 'file') == 2, 'Reference NIfTI not found: %s', refNifti);
assert(exist(refJson,  'file') == 2, 'Reference JSON not found: %s', refJson);

refInfo     = niftiinfo(refNifti);
refJsonInfo = jsondecode(fileread(refJson));
validate_reference_orientation(refInfo, refJsonInfo);

% Build list of [matFile, outputNiiGz] pairs

% woBin — all .mat files found in the folder
woBinDir  = fullfile(subDir, 'recon/woBin');
woBinMats = dir(fullfile(woBinDir, '*.mat'));
nWoBin    = numel(woBinMats);

maskTypes  = {'clean', 'clean_0.50', 'clean_0.75', 'clean_0.95'};
nJobs = nWoBin + numel(maskTypes) * 4 * 3;   % woBin files + masks × regions × (x0, x, x_joint)
jobs  = cell(nJobs, 2);
iJob  = 1;

for iFile = 1:nWoBin
    [~, baseName] = fileparts(woBinMats(iFile).name);
    jobs(iJob, :) = { ...
        fullfile(woBinDir, woBinMats(iFile).name), ...
        fullfile(woBinDir, [baseName '.nii.gz'])};
    iJob = iJob + 1;
end

% clean variants — x0 and x (steva) for each of the 4 gaze directions
for iMask = 1:numel(maskTypes)
    mask = maskTypes{iMask};
    for rIdx = 0:3
        jobs(iJob, :) = { ...
            fullfile(subDir, 'recon', mask, 'x0', sprintf('x0_regionidx%d.mat', rIdx)), ...
            fullfile(subDir, 'recon', mask, 'x0', sprintf('x0_regionidx%d.nii.gz', rIdx))};
        iJob = iJob + 1;
        jobs(iJob, :) = { ...
            fullfile(subDir, 'recon', mask, 'x',  sprintf('x_steva_regionidx_%d_nIter_20_delta_1.000.mat', rIdx)), ...
            fullfile(subDir, 'recon', mask, 'x',  sprintf('x_steva_regionidx_%d_nIter_20_delta_1.000.nii.gz', rIdx))};
        iJob = iJob + 1;
        jobs(iJob, :) = { ...
            fullfile(subDir, 'recon', mask, 'x_joint', sprintf('x_joint_regionidx_%d_nIter_20_ds_1.000_dt_0.100.mat', rIdx)), ...
            fullfile(subDir, 'recon', mask, 'x_joint', sprintf('x_joint_regionidx_%d_nIter_20_ds_1.000_dt_0.100.nii.gz', rIdx))};
        iJob = iJob + 1;
    end
end

% Process each job
nJobs = size(jobs, 1);
for iJob = 1:nJobs
    matFile     = jobs{iJob, 1};
    outputNiiGz = jobs{iJob, 2};

    if ~exist(matFile, 'file')
        fprintf('[%d/%d] skip (no MAT):  %s\n', iJob, nJobs, matFile);
        continue;
    end
    if exist(outputNiiGz, 'file') == 2
        fprintf('[%d/%d] skip (exists):  %s\n', iJob, nJobs, outputNiiGz);
        continue;
    end

    fprintf('\n[%d/%d] Converting:\n  %s\n', iJob, nJobs, matFile);
    convert_single(matFile, outputNiiGz, refInfo);
end

fprintf('\nAll done.\n');


function convert_single(matFile, outputNiiGz, refInfo)
fprintf('  Loading... ');
vol = load_recon_volume(matFile);
fprintf('raw size %s\n', mat2str(size(vol)));

vol = reorient_volume_to_template(vol);
vol = single(abs(vol));

newInfo = build_output_header(refInfo, vol, outputNiiGz);

tmpNii = erase(outputNiiGz, '.gz');
if exist(tmpNii,     'file') == 2, delete(tmpNii);     end
if exist(outputNiiGz,'file') == 2, delete(outputNiiGz); end

niftiwrite(vol, tmpNii, newInfo, 'Compressed', false);
gzip(tmpNii);
delete(tmpNii);
synchronize_qform_with_sform(outputNiiGz);
fprintf('  Saved: %s\n', outputNiiGz);
end


function vol = load_recon_volume(matFile)
data = load(matFile);
candidateFields = {'x', 'x0', 'xrms', 'x0_comp'};

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
