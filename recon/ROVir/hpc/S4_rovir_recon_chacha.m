%% HPC variant of R4 - reconstruct from the ROVir mitosius
%
% Same recon as recon/ROVir/R4_rovir_recon.m, with container paths and no raw
% .dat dependency, so only the mitosius and the virtual coil maps need to be
% pushed to the cluster (not the 16 GB LIBRE file).
%
% Paths are the apptainer bind points used by the other HPC scripts:
%   /usr/src/app/scripts    <- recon/
%   /usr/src/app/data/study <- data/study/
%
% Override subject_num / nv / binName by editing the block below, or generate
% one file per job with make_rovir_jobs.sh.

clc; clearvars -except subject_num nv binName FoV variant; close all;

addpath(genpath('/usr/src/app'));

if ~exist('subject_num', 'var'); subject_num = 15;      end
% 'ROVir'  = SIR-ranked generalised eigenproblem (suppress the head)
% 'ROI-PCA'= ROI-energy-ranked (lambda -> inf); same solver, max compression.
% Outputs go to recon/<variant>/ so the two never mix.
if ~exist('variant',      'var'); variant     = 'ROVir'; end
if ~exist('nv',          'var'); nv          = 20;      end
if ~exist('binName',     'var'); binName     = 'woBin'; end
% Nominal FoV in mm. Hardcoded so the raw .dat never has to be shipped; it is
% 240 for this protocol (confirmed against the reader on the workstation).
if ~exist('FoV',         'var'); FoV         = 240;     end

nIter = 20;
delta = 1;

baseDir    = '/usr/src/app/data/study';
subjectStr = sprintf('sub-%03d', subject_num);
reconDir   = fullfile(baseDir, subjectStr, 'recon');
rovirDir   = fullfile(reconDir, variant);
mDir       = fullfile(reconDir, 'mitosius', sprintf('%s_%d', variant, nv), binName);

tag    = sprintf('%s_%d_%s', variant, nv, binName);
x0Path = fullfile(rovirDir, sprintf('x0_%s.mat', tag));
xPath  = fullfile(rovirDir, sprintf('x_steva_%s_nIter_%d_delta_%.3f.mat', tag, nIter, delta));

fprintf('Config: subject %d, nv %d, bin %s, nIter %d, delta %.3f, FoV %g\n', ...
        subject_num, nv, binName, nIter, delta, FoV);

if exist(xPath, 'file')
    fprintf('x already exists, skipping recon:\n  %s\n', xPath);
    return
end
if ~exist(mDir, 'dir')
    error('Mitosius not found: %s\nPush it with push_rovir.sh first.', mDir);
end

matrix_size = 240;
N_u  = [matrix_size, matrix_size, matrix_size];
n_u  = N_u;
dK_u = [1, 1, 1] ./ FoV;

%% Load

y  = bmMitosius_load(mDir, 'y');
t  = bmMitosius_load(mDir, 't');
ve = bmMitosius_load(mDir, 've');
fprintf('Mitosius loaded from %s\n', mDir);

load(fullfile(rovirDir, sprintf('C_rovir_%d.mat', nv)), 'C_rovir');
% Already in the rotated frame - R3 rotated C before recombining. Do not rotate.
C = bmImResize(C_rovir, [48, 48, 48], N_u);
if size(C, 4) ~= size(y{1}, 2)
    error('Virtual channel mismatch: C has %d, mitosius y has %d.', size(C,4), size(y{1},2));
end

%% Initial estimate

if exist(x0Path, 'file')
    load(x0Path, 'x0');
    fprintf('Existing x0 loaded from %s\n', x0Path);
else
    x0 = cell(1, 1);
    x0{1} = bmMathilda(y{1}, t{1}, ve{1}, C, N_u, n_u, dK_u, [], [], [], []);
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
