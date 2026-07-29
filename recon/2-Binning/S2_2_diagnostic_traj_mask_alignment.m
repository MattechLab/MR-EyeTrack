% =====================================================
% Diagnostic: confirm trajectory shape and eMaskN alignment
% Run once for sub-002 / region 0 before writing the uniform-angular
% undersampling script.
%
% Expected outputs to record:
%   t_tot class + size → confirm [3, N, Nlines]
%   nShot, nSeg, nShot_off
%   numel(eMaskN) → should equal nShot*nSeg
%   (nSeg-1)*(nShot-nShot_off) → should equal Nlines (size(t_tot,3))
% =====================================================
clearvars; clc;
addpath(genpath('/home/debi/jaime/repos/MR-EyeTrack/recon'));
addpath(genpath('/home/debi/jaime/repos/MR-EyeTrack/recon/Recon_scripts'));
addpath(genpath('/home/debi/yiwei/forclone/pulseq'));
addpath(genpath('/home/debi/MatTechLab/monalisa'));

%% Config (match S3_mitosius_binning_pulseq.m exactly)
subject_num  = 2;
mask_type    = 'clean';
th_ratio     = 0.75;
region_idx   = 0;
nShot_off    = 14;

baseDir     = '/home/debi/jaime/repos/MR-EyeTrack/data/study';
subjectStr  = sprintf('sub-%03d', subject_num);
subjectDir  = fullfile(baseDir, subjectStr);
rawDir      = fullfile(subjectDir, 'rawdata');
reconDir    = fullfile(subjectDir, 'recon');
binsDir     = fullfile(reconDir, 'bins', mask_type);

seqFolder   = fullfile(baseDir, 'pulseq');
seqName_list = {
    'yj_seq100_t1w_libre_pre_TR6.2ms_TE3.6ms_swap1_FA6_Traj1_nSeg44_nShot191_Fid0_mreye_2p0.seq', ...
    'yj_seq202_t1w_libre_main_TR8.0ms_TE3.6ms_swap1_FA6_RF2_mreye_track_trajPTP_44_1872_gdsp.seq'};
seqFile     = seqFolder + "/" + seqName_list{2};
seqParams   = extract_seq_params(seqFile);

%% Build reader (trajectory only — skip readRawData to save time)
measureFile = fullfile(rawDir, dir(fullfile(rawDir, '*_T1wLIBRE.dat')).name);
autoFlag    = true;
reader      = createRawDataReader(measureFile, autoFlag);
reader.acquisitionParams.nShot_off            = nShot_off;
reader.acquisitionParams.traj_type            = 'pulseq';
reader.acquisitionParams.pulseqTrajFile_name  = strcat(seqFile);
if isfield(seqParams, 'nshot'), reader.acquisitionParams.nShot = seqParams.nshot; end
if isfield(seqParams, 'nseg'),  reader.acquisitionParams.nSeg  = seqParams.nseg;  end

nShot = reader.acquisitionParams.nShot;
nSeg  = reader.acquisitionParams.nSeg;

fprintf('\n=== Sequence parameters ===\n');
fprintf('  nShot     = %d\n', nShot);
fprintf('  nSeg      = %d\n', nSeg);
fprintf('  nShot_off = %d\n', nShot_off);

%% Compute trajectory
t_tot = bmTraj(reader.acquisitionParams);

fprintf('\n=== Trajectory (t_tot) ===\n');
fprintf('  class : %s\n', class(t_tot));
fprintf('  size  : [%s]\n', num2str(size(t_tot)));
Nlines = size(t_tot, 3);
fprintf('  Nlines (size dim-3) = %d\n', Nlines);
fprintf('  Expected (nSeg-1)*(nShot-nShot_off) = %d\n', (nSeg-1)*(nShot-nShot_off));
fprintf('  Match? %d\n', Nlines == (nSeg-1)*(nShot-nShot_off));

%% Load eMaskN for region 0
eMaskFilePath = fullfile(binsDir, sprintf('eMask_th%.2f_region%i.mat', th_ratio, region_idx));
loaded = load(eMaskFilePath, 'eMaskN');
eMaskN = loaded.eMaskN;

fprintf('\n=== eMaskN (region %d) ===\n', region_idx);
fprintf('  class        : %s\n', class(eMaskN));
fprintf('  size         : [%s]\n', num2str(size(eMaskN)));
fprintf('  numel        = %d\n', numel(eMaskN));
fprintf('  nnz (==1)    = %d\n', nnz(eMaskN == 1));
fprintf('  nShot*nSeg           = %d  (match numel? %d)\n', nShot*nSeg, numel(eMaskN)==nShot*nSeg);
fprintf('  (nShot-off)*nSeg     = %d  (match numel? %d)\n', (nShot-nShot_off)*nSeg, numel(eMaskN)==(nShot-nShot_off)*nSeg);
fprintf('  (nShot-off)*(nSeg-1) = %d  (match numel? %d  == Nlines? %d)\n', ...
    (nShot-nShot_off)*(nSeg-1), numel(eMaskN)==(nShot-nShot_off)*(nSeg-1), (nShot-nShot_off)*(nSeg-1)==Nlines);

%% Probe outermost-sample directions for first 5 spokes
fprintf('\n=== First 5 spoke directions (t_tot(:,end,l) normalised) ===\n');
for l = 1:5
    d = t_tot(:, end, l);
    d = d / norm(d);
    fprintf('  spoke %d: [%.4f  %.4f  %.4f]\n', l, d(1), d(2), d(3));
end
