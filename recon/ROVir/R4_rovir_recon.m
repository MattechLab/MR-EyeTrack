%% R4 - Reconstruct from the ROVir mitosius
%
% Same recon as recon/4-Recon/S4_recon_debi_4fr_pulseq.m (bmMathilda for the
% initial estimate, then bmSteva), pointed at the ROVir mitosius from R3 and
% the virtual coil maps C_rovir_<nv>.mat instead of C.mat.
%
% Outputs -> recon/ROVir/x0_ROVir_<nv>_<bin>.mat
%            recon/ROVir/x_steva_ROVir_<nv>_<bin>_nIter_<n>_delta_<d>.mat

clc; clearvars -except subject_num nv binName variant; close all;

addpath(genpath('/home/debi/jaime/repos/MR-EyeTrack/recon'));
addpath(genpath('/home/debi/MatTechLab/monalisa'));
addpath(genpath('/home/debi/yiwei/forclone/pulseq'));

if ~exist('subject_num', 'var'); subject_num = 15;      end
% 'ROVir'  = SIR-ranked generalised eigenproblem (suppress the head)
% 'ROI-PCA'= ROI-energy-ranked (lambda -> inf); same solver, max compression.
% Outputs go to recon/<variant>/ so the two never mix.
if ~exist('variant',      'var'); variant     = 'ROVir'; end
if ~exist('nv',          'var'); nv          = 20;      end
if ~exist('binName',     'var'); binName     = 'woBin'; end  % or 'mask_0' .. 'mask_3'

nIter = 20;
delta = 1;

baseDir     = '/home/debi/jaime/repos/MR-EyeTrack/data/study';
subjectStr  = sprintf('sub-%03d', subject_num);
subjectDir  = fullfile(baseDir, subjectStr);
rawDir      = fullfile(subjectDir, 'rawdata');
reconDir    = fullfile(subjectDir, 'recon');
rovirDir    = fullfile(reconDir, variant);
mDir        = fullfile(reconDir, 'mitosius', sprintf('%s_%d', variant, nv), binName);

seqFolder = fullfile(baseDir, 'pulseq');
seqName   = 'yj_seq202_t1w_libre_main_TR8.0ms_TE3.6ms_swap1_FA6_RF2_mreye_track_trajPTP_44_1872_gdsp.seq';

tag    = sprintf('%s_%d_%s', variant, nv, binName);
x0Path = fullfile(rovirDir, sprintf('x0_%s.mat', tag));
xPath  = fullfile(rovirDir, sprintf('x_steva_%s_nIter_%d_delta_%.3f.mat', tag, nIter, delta));

fprintf('Config: subject %d, nv %d, bin %s, nIter %d, delta %.3f\n', ...
        subject_num, nv, binName, nIter, delta);

if exist(xPath, 'file')
    fprintf('x already exists, skipping recon:\n  %s\n', xPath);
    return
end
if ~exist(mDir, 'dir')
    error('Mitosius not found: %s\nRun R3 first.', mDir);
end

%% Acquisition parameters (for FoV only)

arrayCoilFile = fullfile(rawDir, dir(fullfile(rawDir, '*_HC.dat')).name);
seqFile   = fullfile(seqFolder, seqName);
seqParams = extract_seq_params(seqFile);

reader = createRawDataReader(arrayCoilFile, true);
reader.acquisitionParams.nShot_off = 14;
reader.acquisitionParams.traj_type = 'pulseq';
reader.acquisitionParams.pulseqTrajFile_name = seqFile;
if isfield(seqParams, 'nshot'); reader.acquisitionParams.nShot = seqParams.nshot; end
if isfield(seqParams, 'nseg');  reader.acquisitionParams.nSeg  = seqParams.nseg;  end

FoV  = reader.acquisitionParams.FoV;
matrix_size = 240;
N_u  = [matrix_size, matrix_size, matrix_size];
n_u  = N_u;
dK_u = [1, 1, 1] ./ FoV;

%% Load mitosius and the virtual coil maps

y  = bmMitosius_load(mDir, 'y');
t  = bmMitosius_load(mDir, 't');
ve = bmMitosius_load(mDir, 've');
fprintf('Mitosius loaded from %s\n', mDir);

CrovPath = fullfile(rovirDir, sprintf('C_rovir_%d.mat', nv));
load(CrovPath, 'C_rovir');
fprintf('Virtual coil maps loaded: %s  %s\n', CrovPath, mat2str(size(C_rovir)));

% C_rovir is already in the rotated frame - R3 applied the subject rotation to
% C before recombining, so do NOT rotate again here.
C = bmImResize(C_rovir, [48, 48, 48], N_u);

if size(C, 4) ~= size(y{1}, 2)
    error('Virtual channel mismatch: C has %d, mitosius y has %d.', ...
          size(C, 4), size(y{1}, 2));
end

%% Initial gridded estimate

nFr = 1;
if exist(x0Path, 'file')
    load(x0Path, 'x0');
    fprintf('Existing x0 loaded from %s\n', x0Path);
else
    x0 = cell(nFr, 1);
    for i = 1:nFr
        x0{i} = bmMathilda(y{i}, t{i}, ve{i}, C, N_u, n_u, dK_u, [], [], [], []);
    end
    save(x0Path, 'x0', '-v7.3');
    fprintf('x0 saved: %s\n', x0Path);
end

%% STEVA

[Gu, Gut] = bmTraj2SparseMat(t, ve, N_u, dK_u);

rho    = 10 * delta;
nCGD   = 4;
ve_max = 10 * prod(dK_u(:));

x = bmSteva(x0{1}, [], [], y{1}, ve{1}, C, Gu{1}, Gut{1}, n_u, ...
            delta, rho, nCGD, ve_max, nIter, ...
            bmWitnessInfo(sprintf('steva_rovir%d', nv), []));

save(xPath, 'x', '-v7.3');
fprintf('x saved: %s\n', xPath);
close all;
