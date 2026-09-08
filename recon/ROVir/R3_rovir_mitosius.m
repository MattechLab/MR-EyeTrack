%% R3 - Apply the ROVir transform to the raw k-space and build the mitosius
%
% Mirrors recon/3-Mitosius/S3_mitosius_t1_woBin_coil_compression_pulseq.m, but
% recombines the coils with the ROVir matrix from R2 instead of an SVD or a
% hard coil subset.
%
% Because the coil recombination is linear and acts only on the channel index,
% it commutes with the (non-Cartesian) encoding operator. So it can be applied
% directly to the radial k-space without regridding:
%
%     y_virt = Vret.' * y          [nv x nPt]
%     C_virt = C      * Vret       [nx x ny x nz x nv]
%
% Both are PLAIN transposes, matching the A = X^H X convention in rovir_solve.
% See the long note in that file: getting this wrong is silent.
%
% Outputs -> recon/ROVir/C_rovir_<nv>.mat
%            recon/mitosius/ROVir_<nv>/...   (woBin, and optionally the 4 bins)

clc; clearvars -except subject_num nvOverride doBins variant; close all;

addpath(genpath('/home/debi/jaime/repos/MR-EyeTrack/recon'));
addpath(genpath('/home/debi/MatTechLab/monalisa'));
addpath(genpath('/home/debi/yiwei/forclone/pulseq'));

if ~exist('subject_num', 'var'); subject_num = 15;    end
% 'ROVir'  = SIR-ranked generalised eigenproblem (suppress the head)
% 'ROI-PCA'= ROI-energy-ranked (lambda -> inf); same solver, max compression.
% Outputs go to recon/<variant>/ so the two never mix.
if ~exist('variant',      'var'); variant     = 'ROVir'; end
if ~exist('nvOverride',  'var'); nvOverride  = [];    end  % else use R2's choice
if ~exist('doBins',      'var'); doBins      = false; end  % also split the 4 gaze bins

matrix_size = 240;
th_ratio    = 0.75;

baseDir    = '/home/debi/jaime/repos/MR-EyeTrack/data/study';
subjectStr = sprintf('sub-%03d', subject_num);
subjectDir = fullfile(baseDir, subjectStr);
rawDir     = fullfile(subjectDir, 'rawdata');
reconDir   = fullfile(subjectDir, 'recon');
rovirDir   = fullfile(reconDir, variant);
qcDir      = fullfile(rovirDir, 'qc');
if ~exist(qcDir, 'dir'); mkdir(qcDir); end

seqFolder   = fullfile(baseDir, 'pulseq');
seqName_list = { ...
    'yj_seq100_t1w_libre_pre_TR6.2ms_TE3.6ms_swap1_FA6_Traj1_nSeg44_nShot191_Fid0_mreye_2p0.seq', ...
    'yj_seq202_t1w_libre_main_TR8.0ms_TE3.6ms_swap1_FA6_RF2_mreye_track_trajPTP_44_1872_gdsp.seq'};

%% Load the ROVir transform

T = load(fullfile(rovirDir, 'rovir_transform.mat'), 'V', 'transform');
nv = T.transform.nv;
if ~isempty(nvOverride); nv = nvOverride; end
Vret = orth(T.V(:, 1:nv));
nChExpect = T.transform.nCh;
fprintf('ROVir transform: %d physical coils -> %d virtual coils\n', nChExpect, nv);
fprintf('  ROI energy retained %.2f%%, interference %.2f%% (gain %.1fx)\n', ...
        T.transform.roiEnergy_retain(nv), T.transform.intEnergy_retain(nv), ...
        T.transform.roiEnergy_retain(nv)/T.transform.intEnergy_retain(nv));

mask_type = sprintf('%s_%d', variant, nv);
mDirBase  = fullfile(reconDir, 'mitosius', mask_type);
if ~exist(mDirBase, 'dir'); mkdir(mDirBase); end

%% Read the raw data

bodyCoilFile  = fullfile(rawDir, dir(fullfile(rawDir, '*_BC.dat')).name);       %#ok<NASGU>
arrayCoilFile = fullfile(rawDir, dir(fullfile(rawDir, '*_HC.dat')).name);       %#ok<NASGU>
measureFile   = fullfile(rawDir, dir(fullfile(rawDir, '*_T1wLIBRE.dat')).name);
fprintf('\nMeasurement file: %s\n', measureFile);

seqFile   = fullfile(seqFolder, seqName_list{2});
seqParams = extract_seq_params(seqFile);

reader = createRawDataReader(measureFile, true);
reader.acquisitionParams.nShot_off = 14;
reader.acquisitionParams.traj_type = 'pulseq';
reader.acquisitionParams.pulseqTrajFile_name = seqFile;
if isfield(seqParams, 'nshot'); reader.acquisitionParams.nShot = seqParams.nshot; end
if isfield(seqParams, 'nseg');  reader.acquisitionParams.nSeg  = seqParams.nseg;  end

y_tot  = reader.readRawData(true, true);      % filter nShotOff and SI
t_tot  = bmTraj(reader.acquisitionParams);
ve_tot = bmVolumeElement(t_tot, 'voronoi_full_radial3');
fprintf('Raw data: %s\n', mat2str(size(y_tot)));

FoV  = reader.acquisitionParams.FoV;
N_u  = [matrix_size, matrix_size, matrix_size];
dK_u = [1, 1, 1] ./ FoV;

%% Load and rotate the coil sensitivities

load(fullfile(reconDir, 'C.mat'), 'C');
rotationMap = containers.Map([2, -1], [3, -1]);
if isKey(rotationMap, subject_num)
    C = rot90(C, rotationMap(subject_num));
    fprintf('Applied subject-specific rotation, C is now %s\n', mat2str(size(C)));
end

nCh = size(y_tot, 1);
if nCh ~= nChExpect || size(C, 4) ~= nChExpect
    error(['Channel count mismatch: transform expects %d, raw data has %d, C has %d. ' ...
           'The ROVir transform must come from the same subject.'], ...
          nChExpect, nCh, size(C, 4));
end

%% Apply the ROVir transform

szY   = size(y_tot);
y2    = reshape(y_tot, nCh, []);              % [nCh x nPt]
y_rov = reshape(Vret.' * y2, [nv, szY(2:end)]);
clear y2 y_tot;
fprintf('\nk-space recombined: %s -> %s\n', mat2str(szY), mat2str(size(y_rov)));

szC     = size(C);
C_rovir = reshape(reshape(C, [], nCh) * Vret, [szC(1:3), nv]);

CrovPath = fullfile(rovirDir, sprintf('C_rovir_%d.mat', nv));
save(CrovPath, 'C_rovir', 'Vret', 'nv', '-v7.3');
fprintf('Virtual coil maps saved: %s\n', CrovPath);

C_rovir_resized = bmImResize(C_rovir, [48, 48, 48], N_u);

%% Gridded check recon, and the normalisation constant

x0_rovir = bmMathilda(y_rov, t_tot, ve_tot, C_rovir_resized, N_u, N_u, dK_u);
x0Path   = fullfile(rovirDir, sprintf('x0_rovir_%d.mat', nv));
save(x0Path, 'x0_rovir', '-v7.3');
fprintf('Gridded ROVir recon saved: %s\n', x0Path);

% Normalise on the eye ROI rather than an interactively drawn polygon, so the
% result is reproducible and so that STEVA's delta means the same thing here as
% in a reference recon normalised the same way. delta is scale-dependent, so a
% different normalisation silently changes the regularisation strength.
M = load(fullfile(rovirDir, 'masks.mat'), 'roiMask');
roiBig = bmImResize(single(M.roiMask), [48 48 48], N_u) > 0.5;
normalize_val = mean(abs(x0_rovir(roiBig)));
fprintf('normalize_val (mean |x0| over the eye ROI) = %.6e\n', normalize_val);

y_rov_norm = y_rov / normalize_val;
clear y_rov;

%% Split into the mitosius

binsDirWo = fullfile(reconDir, 'bins', 'woBin', filesep);
eMaskPath = fullfile(binsDirWo, 'eMask_woBin.mat');
if ~exist(eMaskPath, 'file')
    eMaskPath = fullfile(binsDirWo, 'eMask_woBin');
end

writeMitosius(eMaskPath, fullfile(mDirBase, 'woBin'), y_rov_norm, t_tot, reader);

if doBins
    binsDir = fullfile(reconDir, 'bins', 'clean', filesep);
    for region_idx = 0:3
        p = fullfile(binsDir, sprintf('eMask_th%.2f_region%i.mat', th_ratio, region_idx));
        writeMitosius(p, fullfile(mDirBase, sprintf('mask_%d', region_idx)), ...
                      y_rov_norm, t_tot, reader);
    end
end

fprintf('\nR3 done. Mitosius root: %s\n', mDirBase);

%% ---------------------------------------------------------------------------

function writeMitosius(eMaskPath, mDir, y, t_tot, reader)
% Load a readout binning mask, drop the SI segment and the warm-up shots the
% reader already removed from the data, and write the mitosius.
%
% Naming caution: in this repo "eyeMask" means two different things - the
% image-space eye ROI (used by R1/R2) and this readout binning mask. This is
% the latter.

    S = load(eMaskPath);
    fn = fieldnames(S);
    binMask = S.(fn{1});

    nbins   = size(binMask, 1);
    binMask = reshape(binMask, [nbins, reader.acquisitionParams.nSeg, ...
                                reader.acquisitionParams.nShot]);
    binMask(:, 1, :) = [];                                        % SI projection
    binMask(:, :, 1:reader.acquisitionParams.nShot_off) = [];     % warm-up shots
    binMask = bmPointReshape(binMask);

    [yb, tb] = bmMitosis(y, t_tot, binMask);
    yb = bmPermuteToCol(yb);
    veb = bmVolumeElement(tb, 'voronoi_full_radial3');

    if ~exist(mDir, 'dir'); mkdir(mDir); end
    bmMitosius_create(mDir, yb, tb, veb);
    fprintf('  mitosius written: %s\n', mDir);
end
