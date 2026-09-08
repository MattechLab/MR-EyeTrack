%% R5 - Export the ROVir-compressed raw k-space, and report the compression
%
% ROVir is a *channel* compression: nCh physical coils are linearly recombined
% into nv virtual channels. The output is a normal multi-channel radial
% dataset with a fixed, smaller channel count - the signal is NOT collapsed
% into one channel, and nothing is done to the readout or the trajectory.
%
%     y_virt = Vret.' * y        [nv x nPt]   instead of [nCh x nPt]
%
% The ROVir mitosius for the 'woBin' bin already IS this dataset (woBin keeps
% every readout), so this script reads it rather than re-reading the 16 GB
% raw .dat, and writes a single self-contained file:
%
%     y   [nPt x nv]  complex single   compressed k-space
%     t   [3 x nPt]   single           trajectory (kx, ky, kz)
%     ve  [1 x nPt]   single           volume elements (density weights)
%     meta                             nv, nCh, normalisation, provenance
%
% Output -> data/study/sub-NNN/recon/ROVir/kspace_rovir_<nv>_<bin>.mat

clc; clearvars -except subject_num nv binName variant; close all;

addpath(genpath('/home/debi/jaime/repos/MR-EyeTrack/recon'));
addpath(genpath('/home/debi/MatTechLab/monalisa'));

if ~exist('subject_num', 'var'); subject_num = 15;      end
% 'ROVir'  = SIR-ranked generalised eigenproblem (suppress the head)
% 'ROI-PCA'= ROI-energy-ranked (lambda -> inf); same solver, max compression.
% Outputs go to recon/<variant>/ so the two never mix.
if ~exist('variant',      'var'); variant     = 'ROVir'; end
if ~exist('nv',          'var'); nv          = 20;      end
if ~exist('binName',     'var'); binName     = 'woBin'; end

baseDir    = '/home/debi/jaime/repos/MR-EyeTrack/data/study';
subjectStr = sprintf('sub-%03d', subject_num);
subjectDir = fullfile(baseDir, subjectStr);
reconDir   = fullfile(subjectDir, 'recon');
rovirDir   = fullfile(reconDir, variant);
mDir       = fullfile(reconDir, 'mitosius', sprintf('%s_%d', variant, nv), binName);

T = load(fullfile(rovirDir, 'rovir_transform.mat'), 'V', 'transform');
nCh = T.transform.nCh;

y  = bmMitosius_load(mDir, 'y');   y  = y{1};
t  = bmMitosius_load(mDir, 't');   t  = t{1};
ve = bmMitosius_load(mDir, 've');  ve = ve{1};

nPt = size(y, 1);
fprintf('Compressed k-space: y %s, t %s, ve %s\n', ...
        mat2str(size(y)), mat2str(size(t)), mat2str(size(ve)));

%% Compression figures

rawFile   = dir(fullfile(subjectDir, 'rawdata', '*_T1wLIBRE.dat'));
rawBytes  = rawFile(1).bytes;
bytesPerC = 8;                                   % complex single
compBytes = nPt * nv  * bytesPerC;
origBytes = nPt * nCh * bytesPerC;

fprintf('\n=== Compression (channel axis only) ===\n');
fprintf('  physical coils                : %d\n', nCh);
fprintf('  virtual channels kept         : %d\n', nv);
fprintf('  channel compression ratio     : %.2fx  (%.1f%% of channels kept)\n', ...
        nCh/nv, 100*nv/nCh);
fprintf('  readout points per channel    : %d  (unchanged)\n', nPt);
fprintf('\n  k-space payload, %2d ch        : %8.2f GiB\n', nCh, origBytes/2^30);
fprintf('  k-space payload, %2d ch        : %8.2f GiB\n', nv,  compBytes/2^30);
fprintf('  saved                         : %8.2f GiB (%.1f%%)\n', ...
        (origBytes-compBytes)/2^30, 100*(1-compBytes/origBytes));
fprintf('\n  raw Siemens .dat on disk      : %8.2f GiB (all %d coils + headers)\n', ...
        rawBytes/2^30, nCh);
fprintf('  compressed export (y+t+ve)    : %8.2f GiB\n', ...
        (compBytes + numel(t)*4 + numel(ve)*4)/2^30);

%% Write

meta = struct( ...
    'subject_num',   subject_num, ...
    'bin',           binName, ...
    'nv',            nv, ...
    'nCh',           nCh, ...
    'nPt',           nPt, ...
    'channel_ratio', nCh/nv, ...
    'roiEnergy_retain_pct', T.transform.roiEnergy_retain(nv), ...
    'intEnergy_retain_pct', T.transform.intEnergy_retain(nv), ...
    'normalised',    true, ...
    'note',          ['y is normalised by mean|x0| over the eye ROI (see R3). ' ...
                      'Trajectory and volume elements are unchanged by ROVir - ' ...
                      'only the channel axis is recombined.'], ...
    'created',       datetime('now'));

outPath = fullfile(rovirDir, sprintf('kspace_%s_%d_%s.mat', variant, nv, binName));
save(outPath, 'y', 't', 've', 'meta', '-v7.3');
fprintf('\nExported: %s\n', outPath);
d = dir(outPath);
fprintf('  on disk: %.2f GiB\n', d.bytes/2^30);
