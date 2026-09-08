%% mat2nii_auto_orient_test.m
% Investigates whether the NIfTI affine can be built purely from the twix
% header + pulseq trajectory, without relying on an MPRAGE reference.
%
% For each dataset (MR-EyeTrack sub-004, Yiwei T1w MID00030, Yiwei T2w
% MID00025) the script does three things:
%
%   (A) Twix-header orientation — reads sNormal / dRowDir / dColDir from
%       the MeasYaps header and converts them LPS → RAS.  This describes
%       the prescribed imaging-volume orientation.
%
%   (B) Trajectory-axis analysis — loads the pulseq .seq file, computes
%       the first full k-space line, and reports the dominant physical
%       direction of each trajectory row (kx, ky, kz) in RAS.  This tells
%       us which anatomical axis each reconstruction dim points along.
%
%   (C) Auto-affine — builds a NIfTI with the affine assembled from the
%       twix header alone (no MPRAGE), applies NO reorientFcn so dim1/2/3
%       match kx/ky/kz.  The resulting file can be compared with the
%       MPRAGE-based NIfTIs already in /home/debi/Downloads/.
%
% Outputs (all saved to /home/debi/Downloads/auto_orient/<dataset>/):
%   <dataset>_auto_orient.nii.gz   — NIfTI built from twix affine only
%   <dataset>_orientation_report.txt — printed analysis per dataset
%
% After running, open both the auto-orient NIfTI and the corresponding
% MPRAGE-based NIfTI in Mango and compare visually + via
% check_nifti_orientation.m.

clc; clearvars; close all;

repoRoot = '/home/debi/jaime/repos/MR-EyeTrack';
addpath(genpath(fullfile(repoRoot, 'analysis')));
addpath(genpath(fullfile(repoRoot, 'recon')));
addpath(genpath('/home/debi/MatTechLab/monalisa'));
addpath(genpath('/home/debi/yiwei/forclone/pulseq'));

outRoot = '/home/debi/Downloads/auto_orient';

%% =========================================================================
%  Dataset definitions
% =========================================================================

% --- MR-EyeTrack sub-004 ---
d = struct();
d.name       = 'mreyetrack';
d.matFile    = fullfile(repoRoot, 'data/study/sub-004/recon/woBin', ...
    'x_steva_nIter_20_delta_1.000.mat');
d.twixFile   = fullfile(repoRoot, 'data/study/sub-004/rawdata', ...
    'sub-004_T1wLIBRE.dat');
d.seqFile    = fullfile(repoRoot, 'data/study/pulseq', ...
    'yj_seq202_t1w_libre_main_TR8.0ms_TE3.6ms_swap1_FA6_RF2_mreye_track_trajPTP_44_1872_gdsp.seq');
d.refNifti   = fullfile(repoRoot, 'data/study/sub-004/dicom', ...
    'csTFL_mp-rage_1mm-iso_CP_acc4.6_5_MR/sub-004.nii.gz');
d.knownFcn   = '@(v) flip(flip(permute(v,[2 1 3 (4:ndims(v))]),1),2)';
d.knownRaw   = '(-A, -R, +S)';
datasets(1)  = d;

% --- Yiwei 2.0 T1w LIBRE (MID00030) ---
filerRoot = '/mnt/filer01/MatTechLab/yiwei.jia';
dataDir   = fullfile(filerRoot, 'datasets', '0005');
d = struct();
d.name       = 'yiwei_t1w';
d.matFile    = fullfile(filerRoot, 'recon_results', '0005', 'MID00030_recon', ...
    'T1_LIBRE_woBinning', 'rovir_ncoil25_ortho', 'output', ...
    'x_20260304_195803', 'x_Nx480_nIter15_delta_1.000.mat');
d.twixFile   = fullfile(dataDir, 'meas_MID00030_FID04261_t1w_libre.dat');
d.seqFile    = fullfile(filerRoot, 'datasets', '0001', ...
    'yj_seq104_t1w_libre_main_TR8.0ms_TE3.6ms_swap1_FA6_Traj1_nSeg44_nShot1055_Fid0_mreye_2p0_gdsp.seq');
d.refNifti   = fullfile(dataDir, 'dicom_0005', ...
    'csTFL_mp-rage_1mm-iso_CP_acc4.6_5_MR/mprage.nii.gz');
d.knownFcn   = '@(v) flip(flip(permute(v,[2 3 1 (4:ndims(v))]),2),3)';
d.knownRaw   = '(-S, +R, -A)';
datasets(2)  = d;

% --- Yiwei 2.0 T2w LIBRE (MID00025) ---
d = struct();
d.name       = 'yiwei_t2w';
d.matFile    = fullfile(filerRoot, 'recon_results', '0005', 'MID00025_recon', ...
    'T1_LIBRE_woBinning', 'rovir_ncoil25_ortho', 'output', ...
    'x_20260318_160957', 'x_Nx480_nIter15_delta_4.000.mat');
d.twixFile   = fullfile(dataDir, 'meas_MID00025_FID04256_t2w_libre.dat');
d.seqFile    = fullfile(dataDir, ...
    'yj_seq608_t2w_libre_main_TR8.0ms_TE3.6ms_swap1_FA6_Traj1_nSeg44_nShot1000_Fid0_t50_crusher_646p92.seq');
d.refNifti   = fullfile(dataDir, 'dicom_0005', ...
    'csTFL_mp-rage_1mm-iso_CP_acc4.6_5_MR/mprage.nii.gz');
d.knownFcn   = '@(v) flip(permute(v,[1 3 2 (4:ndims(v))]),2)';   % corrected: dim-3 flip dropped (L-R was mirrored)
d.knownRaw   = '(-R, -S, -A)';   % disputed corrected; original claim was (+R,-S,-A)
datasets(3)  = d;


%% =========================================================================
%  Process each dataset
% =========================================================================

for iD = 1:numel(datasets)
    ds      = datasets(iD);
    outDir  = fullfile(outRoot, ds.name);
    logFile = fullfile(outDir, [ds.name '_orientation_report.txt']);
    outNii  = fullfile(outDir,  [ds.name '_auto_orient.nii.gz']);

    fid = fopen(logFile, 'w');
    logPrint(fid, '==========================================================');
    logPrint(fid, 'Dataset: %s', ds.name);
    logPrint(fid, 'Known empirical reorientFcn : %s', ds.knownFcn);
    logPrint(fid, 'Known raw axis convention   : %s', ds.knownRaw);
    logPrint(fid, '==========================================================\n');

    % ------------------------------------------------------------------
    % (A) Twix header orientation
    % ------------------------------------------------------------------
    logPrint(fid, '--- (A) Twix header orientation ---');

    if ~exist(ds.twixFile, 'file')
        logPrint(fid, 'SKIP: twix file not found:\n  %s\n', ds.twixFile);
        fclose(fid); continue;
    end

    twix = mapVBVD_JH_for_monalisa(ds.twixFile);
    if iscell(twix), twix = twix{end}; end
    sa = twix.hdr.MeasYaps.sSliceArray.asSlice{1};

    % Position (LPS, mm)
    pos_lps = extract_vec3(sa, 'sPosition', {'dSag','dCor','dTra'});
    pos_ras = lps2ras(pos_lps);
    logPrint(fid, 'sPosition (LPS) [Sag Cor Tra]: [%.2f  %.2f  %.2f] mm', pos_lps);
    logPrint(fid, 'sPosition (RAS) [R   A   S  ]: [%.2f  %.2f  %.2f] mm', pos_ras);

    % Normal (LPS) — this is the slice-selection / kz direction for 2D.
    % For 3D LIBRE it describes the slab normal (which scanner gradient
    % channel is predominantly along the slab normal direction).
    norm_lps = extract_vec3(sa, 'sNormal', {'dSag','dCor','dTra'});
    norm_ras = lps2ras(norm_lps);
    logPrint(fid, '\nsNormal (LPS): [%.4f  %.4f  %.4f]', norm_lps);
    logPrint(fid, 'sNormal (RAS): [%.4f  %.4f  %.4f]  → dominant: %s', ...
        norm_ras, dominant_axis(norm_ras));

    % Row direction (readout, 1st image axis for Cartesian → Gx for LIBRE)
    if isfield(sa, 'dRowDir')
        row_lps = extract_vec3(sa, 'dRowDir', {'dSag','dCor','dTra'});
        row_ras = lps2ras(row_lps);
        logPrint(fid, '\ndRowDir (LPS): [%.4f  %.4f  %.4f]', row_lps);
        logPrint(fid, 'dRowDir (RAS): [%.4f  %.4f  %.4f]  → dominant: %s', ...
            row_ras, dominant_axis(row_ras));
    else
        row_ras = [];
        logPrint(fid, '\ndRowDir: field not present in header');
    end

    % Column direction (phase, 2nd image axis → Gy for LIBRE)
    if isfield(sa, 'dColDir')
        col_lps = extract_vec3(sa, 'dColDir', {'dSag','dCor','dTra'});
        col_ras = lps2ras(col_lps);
        logPrint(fid, '\ndColDir (LPS): [%.4f  %.4f  %.4f]', col_lps);
        logPrint(fid, 'dColDir (RAS): [%.4f  %.4f  %.4f]  → dominant: %s', ...
            col_ras, dominant_axis(col_ras));
    else
        col_ras = [];
        logPrint(fid, '\ndColDir: field not present in header');
    end

    logPrint(fid, '\n--- Predicted affine from twix header (if dim1=Row, dim2=Col, dim3=Normal) ---');
    if ~isempty(row_ras) && ~isempty(col_ras)
        logPrint(fid, '  dim1 would be: %s  (dRowDir)', dominant_axis(row_ras));
        logPrint(fid, '  dim2 would be: %s  (dColDir)', dominant_axis(col_ras));
        logPrint(fid, '  dim3 would be: %s  (sNormal)', dominant_axis(norm_ras));
    end

    % N readout points from twix header; seq definitions don't expose this directly
    N_readout = double(twix.image.NCol);

    % Seq definitions (nshot, nseg, fov, …) — safe even if file missing
    seqParams = extract_seq_params(ds.seqFile);

    % ------------------------------------------------------------------
    % (B) Trajectory-axis analysis (first projection line from seq file)
    % ------------------------------------------------------------------
    logPrint(fid, '\n--- (B) Trajectory-axis analysis (seq file) ---');

    if ~exist(ds.seqFile, 'file')
        logPrint(fid, 'SKIP: seq file not found:\n  %s\n', ds.seqFile);
    else
        try
            seq = mr.Sequence();
            seq.read(ds.seqFile);
            kpp = seq.calculateKspacePP();

            % First readout line: first N points (before any segmentation)
            % N comes from the twix header (NCol = ADC samples per readout);
            % seq.definitions rarely exposes this directly.
            N = min(N_readout, size(kpp, 2));
            firstLine = kpp(:, 1:N);   % [3 x N] in physical gradient units (1/m)

            % The midpoint of a line through the k-space origin gives the
            % projection direction.  Use the endpoint (max magnitude).
            [~, iMax] = max(vecnorm(firstLine, 2, 1));
            kdir_phys = firstLine(:, iMax);   % dominant direction in physical [Gx,Gy,Gz] = [LPS-X,LPS-Y,LPS-Z]
            kdir_phys = kdir_phys / norm(kdir_phys);

            % Convert physical gradient axes to RAS for the FIRST line
            % Physical gradient frame (Siemens): X→LPS-Sag, Y→LPS-Cor, Z→LPS-Tra
            % → RAS: R=-LPS-Sag, A=-LPS-Cor, S=+LPS-Tra
            kdir_ras = [-kdir_phys(1); -kdir_phys(2); kdir_phys(3)];

            logPrint(fid, 'First projection direction (physical Gx,Gy,Gz): [%.4f  %.4f  %.4f]', kdir_phys);
            logPrint(fid, 'First projection direction (RAS):                [%.4f  %.4f  %.4f]', kdir_ras);

            % For a 3D radial trajectory, the three trajectory ROWS
            % correspond to the three image dims.  Sample the first few
            % projections to estimate which anatomical axis each row
            % "covers" most (maximum variance in each row direction).
            logPrint(fid, '\nTrajectory row coverage (std of each row across all points):');
            nLines = min(500, size(kpp, 2));
            kSample = kpp(:, 1:nLines);
            for row = 1:3
                s = std(kSample(row,:));
                logPrint(fid, '  t row %d std = %.4f  (larger = more spatial variation along this gradient axis)', row, s);
            end

            % Estimate dominant direction per trajectory row from first line
            logPrint(fid, '\nFirst-line unit vector per trajectory row:');
            for row = 1:3
                v_lps = [firstLine(row,iMax); 0; 0];  % placeholder
            end

            % Better: look at the outer product structure of the first few lines
            nSample = min(200, floor(size(kpp, 2)/N));
            if nSample >= 1
                rowDirs_phys = zeros(3, nSample);
                for li = 1:nSample
                    seg = kpp(:, (li-1)*N+1 : li*N);
                    [~, im] = max(vecnorm(seg, 2, 1));
                    v = seg(:, im);
                    nv = norm(v);
                    if nv > 1e-9, rowDirs_phys(:,li) = v/nv; end
                end
                rowDirs_ras = rowDirs_phys .* [-1;-1;1];  % LPS→RAS sign flip

                logPrint(fid, '\nMean |projection component| per trajectory row (across %d lines):', nSample);
                logPrint(fid, '  (each column = |Gx|, |Gy|, |Gz| contribution to the projection direction)');
                % These projections are whole-line endpoint vectors stacked
                % in rowDirs.  Each row of rowDirs_ras = [R, A, S] component.
                % The row with the most variation tells us which gradient
                % axis "drives" that trajectory direction.
                covMat = rowDirs_ras * rowDirs_ras' / nSample;
                logPrint(fid, '  Trajectory-direction covariance matrix (RAS frame):');
                logPrint(fid, '  R: [%.4f  %.4f  %.4f]', covMat(1,:));
                logPrint(fid, '  A: [%.4f  %.4f  %.4f]', covMat(2,:));
                logPrint(fid, '  S: [%.4f  %.4f  %.4f]', covMat(3,:));
                logPrint(fid, '  The diagonal encodes how much each trajectory row (kx,ky,kz) ');
                logPrint(fid, '  is spread along R, A, S respectively.');
            end

        catch ME
            logPrint(fid, 'ERROR loading seq file: %s', ME.message);
        end
    end

    % ------------------------------------------------------------------
    % (C) Build auto-affine from twix and write NIfTI (no reorientFcn)
    % ------------------------------------------------------------------
    logPrint(fid, '\n--- (C) Auto-affine NIfTI (twix header only, no reorientFcn) ---');

    if ~exist(ds.matFile, 'file')
        logPrint(fid, 'SKIP: mat file not found:\n  %s\n', ds.matFile);
        fclose(fid); continue;
    end

    % Load reconstruction volume
    vol = load_recon_volume_local(ds.matFile);
    vol = single(abs(vol));
    Nv  = size(vol);
    logPrint(fid, 'Raw MAT size: %s', mat2str(Nv));

    % FOV: pulseq seq definitions store it in metres (key 'FOV' → field 'fov').
    % The seq FOV is the acquisition k-space FOV; the reconstruction doubles
    % it (2× oversampling), so multiply by 2 to get the image-space FOV.
    % Fall back to the twix Config header if not present in the seq file.
    if isfield(seqParams, 'fov') && ~isempty(seqParams.fov)
        FOV = double(seqParams.fov(:)') * 1e3 * 2;  % m → mm, ×2 for 2× oversampling
    elseif isfield(twix.hdr, 'Config') && isfield(twix.hdr.Config, 'ReadFoV')
        FOV = repmat(double(twix.hdr.Config.ReadFoV), 1, 3);
    else
        logPrint(fid, 'WARNING: FOV not found in seq params or twix header; using 240 mm');
        FOV = [240 240 240];
    end
    FOV = FOV(1:min(3,end));
    if numel(FOV) < 3, FOV = repmat(FOV(1), 1, 3); end

    TR = get_tr_seconds(twix.hdr);

    voxelSize = FOV(:)' ./ double(Nv(1:3));
    logPrint(fid, 'FOV (mm):        [%.2f  %.2f  %.2f]', FOV);
    logPrint(fid, 'Voxel size (mm): [%.4f  %.4f  %.4f]', voxelSize);

    % Build direction cosines from twix.
    % Convention: for a 3D LIBRE acquisition, the trajectory rows map as
    %   kx → image dim1 → physical Gx channel → LPS-Sag/RAS-R direction
    %   ky → image dim2 → physical Gy channel → LPS-Cor/RAS-A direction
    %   kz → image dim3 → physical Gz channel → LPS-Tra/RAS-S direction
    % This is the PHYSICAL gradient frame, BEFORE any sequence-level axis
    % swapping.  If the sequence swaps Gx↔Gy (e.g. swap1), the mapping
    % changes and the NIfTI affine will be misaligned.
    %
    % We also test the alternative: use the twix sNormal/dRowDir/dColDir
    % to build a header-based affine.

    % Approach 1: physical gradient frame (kx=Gx=-R, ky=Gy=-A, kz=Gz=+S)
    dirCos_physGrad = [-1  0  0;    % dim1 = -R  (Gx = LPS-Sag, negate for RAS)
                        0 -1  0;    % dim2 = -A  (Gy = LPS-Cor, negate for RAS)
                        0  0  1];   % dim3 = +S  (Gz = LPS-Tra)
    dirCos_physGrad = dirCos_physGrad ./ vecnorm(dirCos_physGrad, 2, 1);

    % Approach 2: header-derived (if dRowDir and dColDir are present)
    if ~isempty(row_ras) && ~isempty(col_ras)
        dirCos_hdr = [row_ras(:)/norm(row_ras), ...
                      col_ras(:)/norm(col_ras), ...
                      norm_ras(:)/norm(norm_ras)];
    else
        dirCos_hdr = [];
    end

    % Origin from sPosition
    posLps  = extract_vec3(sa, 'sPosition', {'dSag','dCor','dTra'});
    posRas  = lps2ras(posLps);

    for iApproach = 1:2
        if iApproach == 1
            dc  = dirCos_physGrad;
            tag = 'physGrad';
        else
            if isempty(dirCos_hdr), continue; end
            dc  = dirCos_hdr;
            tag = 'hdrDirCos';
        end

        originRas = posRas(:) ...
            - dc(:,1) * voxelSize(1) * (Nv(1)-1)/2 ...
            - dc(:,2) * voxelSize(2) * (Nv(2)-1)/2 ...
            - dc(:,3) * voxelSize(3) * (Nv(3)-1)/2;

        affine = eye(4);
        affine(1:3,1:3) = dc .* voxelSize;
        affine(1:3,4)   = originRas;

        outFile = fullfile(outDir, sprintf('%s_%s.nii.gz', ds.name, tag));
        info = build_nifti_info_local(vol, affine, TR, outFile);
        if exist(outFile,'file'), delete(outFile); end
        niftiwrite(vol, outFile, info, 'Compressed', true);
        logPrint(fid, '\nWrote [%s]: %s', tag, outFile);
        logPrint(fid, '  affine column 1 (dim1 direction, RAS): [%.4f  %.4f  %.4f]', dc(:,1)');
        logPrint(fid, '  affine column 2 (dim2 direction, RAS): [%.4f  %.4f  %.4f]', dc(:,2)');
        logPrint(fid, '  affine column 3 (dim3 direction, RAS): [%.4f  %.4f  %.4f]', dc(:,3)');
    end

    % ------------------------------------------------------------------
    % Compare auto-affine with MPRAGE-based affine
    % ------------------------------------------------------------------
    logPrint(fid, '\n--- (D) Comparison: auto-affine vs MPRAGE-based affine ---');

    if exist(ds.refNifti, 'file')
        refInfo   = niftiinfo(ds.refNifti);
        A_mprage  = refInfo.Transform.T';
        R_mprage  = A_mprage(1:3,1:3) ./ vecnorm(A_mprage(1:3,1:3), 2, 1);

        logPrint(fid, '\nMPRAGE affine direction cosines (each column = dim direction in RAS):');
        logPrint(fid, '  dim1: [%.4f  %.4f  %.4f]  → %s', R_mprage(:,1)', dominant_axis(R_mprage(:,1)));
        logPrint(fid, '  dim2: [%.4f  %.4f  %.4f]  → %s', R_mprage(:,2)', dominant_axis(R_mprage(:,2)));
        logPrint(fid, '  dim3: [%.4f  %.4f  %.4f]  → %s', R_mprage(:,3)', dominant_axis(R_mprage(:,3)));

        logPrint(fid, '\nPhysical-gradient-frame direction cosines:');
        logPrint(fid, '  dim1: [%.4f  %.4f  %.4f]  → %s (raw kx = Gx = LPS-Sag = -R)', ...
            dirCos_physGrad(:,1)', dominant_axis(dirCos_physGrad(:,1)));
        logPrint(fid, '  dim2: [%.4f  %.4f  %.4f]  → %s (raw ky = Gy = LPS-Cor = -A)', ...
            dirCos_physGrad(:,2)', dominant_axis(dirCos_physGrad(:,2)));
        logPrint(fid, '  dim3: [%.4f  %.4f  %.4f]  → %s (raw kz = Gz = LPS-Tra = +S)', ...
            dirCos_physGrad(:,3)', dominant_axis(dirCos_physGrad(:,3)));

        logPrint(fid, '\nKnown empirical raw axes: %s', ds.knownRaw);
        logPrint(fid, 'Physical gradient prediction: (-R, -A, +S)  ← this is the "default" before any seq-level swap');

        % Mismatch between physical-gradient prediction and empirical result
        logPrint(fid, '\nConclusion:');
        if strcmp(ds.knownRaw, '(-A, -R, +S)')
            logPrint(fid, '  MR-EyeTrack: empirical=(-A,-R,+S) vs predicted=(-R,-A,+S).');
            logPrint(fid, '  Dim1 and dim2 are SWAPPED → the sequence name "swap1" confirms');
            logPrint(fid, '  that the sequence explicitly swaps Gx↔Gy in the gradient waveforms.');
        elseif strcmp(ds.knownRaw, '(-S, +R, -A)')
            logPrint(fid, '  Yiwei T1w: empirical=(-S,+R,-A) vs predicted=(-R,-A,+S).');
            logPrint(fid, '  This is a [2,3,1] cycle permutation + sign changes.');
            logPrint(fid, '  The seq file likely cycles Gz→Gx, Gx→Gy, Gy→Gz (or similar),');
            logPrint(fid, '  resulting in a cyclic axis rotation in the gradient frame.');
        elseif strcmp(ds.knownRaw, '(-R, -S, -A)')
            logPrint(fid, '  Yiwei T2w (corrected): empirical=(-R,-S,-A) vs predicted=(-R,-A,+S).');
            logPrint(fid, '  Dim1 matches (-R); dim2↔dim3 are swapped (-S and -A).');
            logPrint(fid, '  The T2w sequence swaps ky↔kz relative to the standard LIBRE,');
            logPrint(fid, '  i.e. the phase/partition gradient channels are exchanged.');
        end
    else
        logPrint(fid, 'MPRAGE reference not found; skipping comparison.');
    end

    fclose(fid);
    type(logFile);  % print to console
end

fprintf('\nAll done. Outputs in: %s\n', outRoot);


%% =========================================================================
%  Local helpers
%% =========================================================================

function v = extract_vec3(s, fieldName, subFields)
v = zeros(3,1);
if ~isfield(s, fieldName), return; end
sub = s.(fieldName);
for i = 1:3
    f = subFields{i};
    if isfield(sub, f), v(i) = double(sub.(f)); end
end
end

function ras = lps2ras(lps)
ras = [-lps(1); -lps(2); lps(3)];
end

function label = dominant_axis(ras_vec)
ras_vec = ras_vec(:);
[~, idx] = max(abs(ras_vec));
names = {'R','A','S'};
sign_str = '';
if ras_vec(idx) < 0, sign_str = '-'; end
label = sprintf('%s%s', sign_str, names{idx});
end

function logPrint(fid, fmt, varargin)
msg = sprintf(fmt, varargin{:});
fprintf('%s\n', msg);
if fid > 0, fprintf(fid, '%s\n', msg); end
end

function vol = load_recon_volume_local(matFile)
data = load(matFile);
for f = {'x','x0','xrms'}
    if isfield(data, f{1})
        vol = data.(f{1});
        if iscell(vol), vol = cat(4, vol{:}); end
        return;
    end
end
error('No reconstruction field found in %s', matFile);
end

function TR = get_tr_seconds(hdr)
TR = 1;
try
    if isfield(hdr.MeasYaps,'alTR')
        TR = double(hdr.MeasYaps.alTR{1}) * 1e-6;
    elseif isfield(hdr.Dicom,'alTR')
        TR = double(hdr.Dicom.alTR{1}) * 1e-6;
    end
catch
end
end

function info = build_nifti_info_local(vol, affine, TR, outputPath)
voxelSize = vecnorm(affine(1:3,1:3), 2, 1);
dirCos    = affine(1:3,1:3) ./ voxelSize;
hdr = images.internal.nifti.niftiImage.niftiDefaultHeader(vol, true, 'NIfTI1');
hdr.pixdim     = [1, voxelSize, max(TR, eps)];
hdr.xyzt_units = uint8(2) + uint8(8);
hdr.descrip    = 'auto-orient from twix header';
hdr.datatype   = int16(16);
hdr.bitpix     = int16(32);
hdr.sform_code = int16(1);
hdr.qform_code = int16(1);
hdr.srow_x     = affine(1,:);
hdr.srow_y     = affine(2,:);
hdr.srow_z     = affine(3,:);
q = rotmat_to_quat_local(dirCos);
hdr.quatern_b  = q(2); hdr.quatern_c = q(3); hdr.quatern_d = q(4);
hdr.qoffset_x  = affine(1,4);
hdr.qoffset_y  = affine(2,4);
hdr.qoffset_z  = affine(3,4);
hdr.cal_min    = double(min(vol(:)));
hdr.cal_max    = double(max(vol(:)));
hdr.scl_slope  = 1; hdr.scl_inter = 0;
NV   = images.internal.nifti.niftiImage(hdr);
info = NV.simplifyStruct();
info.raw         = hdr;
info.Filename    = outputPath;
info.Description = 'auto-orient from twix header';
info.Filemoddate = char(datetime('now','Format','dd-MMM-yyyy HH:mm:ss'));
end

function q = rotmat_to_quat_local(R)
traceR = trace(R);
if traceR > 0
    s = 2*sqrt(traceR+1); qw=0.25*s; qx=(R(3,2)-R(2,3))/s; qy=(R(1,3)-R(3,1))/s; qz=(R(2,1)-R(1,2))/s;
elseif (R(1,1)>R(2,2)) && (R(1,1)>R(3,3))
    s=2*sqrt(1+R(1,1)-R(2,2)-R(3,3)); qw=(R(3,2)-R(2,3))/s; qx=0.25*s; qy=(R(1,2)+R(2,1))/s; qz=(R(1,3)+R(3,1))/s;
elseif R(2,2)>R(3,3)
    s=2*sqrt(1+R(2,2)-R(1,1)-R(3,3)); qw=(R(1,3)-R(3,1))/s; qx=(R(1,2)+R(2,1))/s; qy=0.25*s; qz=(R(2,3)+R(3,2))/s;
else
    s=2*sqrt(1+R(3,3)-R(1,1)-R(2,2)); qw=(R(2,1)-R(1,2))/s; qx=(R(1,3)+R(3,1))/s; qy=(R(2,3)+R(3,2))/s; qz=0.25*s;
end
q=[qw qx qy qz]; q=q/norm(q); if q(1)<0, q=-q; end
end
