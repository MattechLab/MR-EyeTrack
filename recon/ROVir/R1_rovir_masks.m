%% R1 - Build the ROVir region masks on the 48^3 gridding grid
%
% ROVir needs two regions:
%   roiMask - the region to KEEP  (the orbits)
%   intMask - the region to SUPPRESS ("interference": the rest of the head)
%
% Unlike the MR RawDeface abstract, which keeps the brain and suppresses only
% the face, here the ROI is small and the interference is everything else. The
% ROI therefore does NOT need to be tight: a blocky region that swallows some
% frontal lobe, nose or sinus is fine, and actually helps, because it widens
% the separation between the two regions and leaves the eigenproblem a
% better-conditioned job.
%
% The masks live on the 48^3 grid used by woBin/x0_noC_48.mat (5 mm isotropic,
% FoV 240). That is coarse for anatomy but ample for ROVir: the covariances are
% spatial integrals weighted by coil sensitivity, and coil sensitivity varies on
% centimetre scales. It also means the manually drawn 2D eyeMask needs no
% resampling at all - it is already native to this grid.
%
% Grid convention for this dataset (verified on sub-015 against xrms48):
%   dim1 = anterior -> posterior   (anterior at low index)
%   dim2 = left <-> right
%   dim3 = inferior -> superior    (superior at high index)
%
% Outputs -> data/study/sub-NNN/recon/ROVir/masks.mat  (+ QC png)

clc; close all;

if ~exist('subject_num', 'var'); subject_num = 15;    end
% 'ROVir'  = SIR-ranked generalised eigenproblem (suppress the head)
% 'ROI-PCA'= ROI-energy-ranked (lambda -> inf); same solver, max compression.
% Outputs go to recon/<variant>/ so the two never mix.
if ~exist('variant',      'var'); variant     = 'ROVir'; end
% dim3 band holding the orbits. Centred on slice 24 because that is where the
% manual 2D mask was drawn for every subject: coilSelectionEyesROI.m draws on
% xrms(:,:,round(nSlice/2)) and xrms is 48^3, so round(48/2) = 24. The eyes are
% therefore visible at slice 24 by construction, and +-6 slices (65 mm) brackets
% the orbit with margin. ROI size barely affects the result, so err generous.
if ~exist('zRange',      'var'); zRange      = 18:30; end
if ~exist('gapVox',      'var'); gapVox      = 2;     end  % voxels excluded from BOTH regions
if ~exist('headThr',     'var'); headThr     = 0.10;  end  % head mask, fraction of the 99.5th pct

addpath(genpath('/home/debi/jaime/repos/MR-EyeTrack/recon'));

baseDir    = '/home/debi/jaime/repos/MR-EyeTrack/data/study';
subjectStr = sprintf('sub-%03d', subject_num);
reconDir   = fullfile(baseDir, subjectStr, 'recon');
rovirDir   = fullfile(reconDir, variant);
if ~exist(rovirDir, 'dir'); mkdir(rovirDir); end

%% Load the RMS image that defines "where there is head"

xrmsPath = fullfile(reconDir, 'woBin', 'xrms48.mat');
X = load(xrmsPath, 'xrms');
xrms = abs(single(X.xrms));
N = size(xrms, 1);
fprintf('Loaded %s  [%s]\n', xrmsPath, num2str(size(xrms)));

%% Load the manually drawn 2D eye ROI and extrude it

eyeMaskPath = fullfile(reconDir, 'mitosius', 'woBin_comp', 'eyeMask.mat');
E = load(eyeMaskPath, 'eyeMask');
eyeMask2D = logical(E.eyeMask);
fprintf('Loaded manual 2D eye ROI  [%s]  nnz=%d\n', ...
        num2str(size(eyeMask2D)), nnz(eyeMask2D));

if ~isequal(size(eyeMask2D), [N N])
    error(['Manual eyeMask is %dx%d but the grid is %d^3. It was drawn on a ' ...
           'different matrix size - redraw it or resize before use.'], ...
          size(eyeMask2D,1), size(eyeMask2D,2), N);
end

roiMask = false(N, N, N);
roiMask(:, :, zRange) = repmat(eyeMask2D, [1 1 numel(zRange)]);

% Deliberately NOT intersected with the head mask below. The globe and the
% orbital fat are both dark in fat-suppressed LIBRE, so a signal threshold
% carves the orbit straight out of the head - exactly the voxels we want to
% keep. The extruded box is used as drawn.

%% Head mask: everything with signal

% A fixed fraction of the robust maximum, not Otsu. Otsu on the linear image
% thresholds at 0.36 and throws away ~40% of the head; Otsu on a log-compressed
% image recovers the tissue but also drags in the radial gridding ripple in the
% background. A plain percentile fraction sits between the two.
%
% Err generous. Background voxels carry almost no signal, so they contribute
% almost nothing to the covariance B and a slightly leaky mask is harmless.
% Tissue *missing* from the mask is not harmless: it is signal we then fail to
% suppress.
v = xrms(:);
lo = prctile(v, 1); hi = prctile(v, 99.5);
xn = (xrms - lo) / (hi - lo); xn(xn < 0) = 0; xn(xn > 1) = 1;
thr = headThr;
headMask = xn > thr;

% Clean up: break thin background bridges, keep the largest component, fill in
headMask = imopen(headMask, strel('sphere', 1));
headMask = imclose(headMask, strel('sphere', 2));
cc = bwconncomp(headMask, 26);
if cc.NumObjects > 1
    [~, big] = max(cellfun(@numel, cc.PixelIdxList));
    headMask = false(size(headMask));
    headMask(cc.PixelIdxList{big}) = true;
end
headMask = imfill(headMask, 'holes');

%% Interference = head, minus the ROI, minus a gap shell

roiDil  = imdilate(roiMask, strel('sphere', gapVox));
intMask = headMask & ~roiDil;

% The gap matters. If the two regions touch, the generalised eigenproblem is
% asked to suppress voxels one voxel away from the ROI, which no linear coil
% combination can do, and the solution degenerates towards noise-amplifying
% eigenvectors.
gapMask = roiDil & ~roiMask & headMask;

nVoxROI = nnz(roiMask); nVoxINT = nnz(intMask);
fprintf('\nRegion sizes (48^3 grid, 5 mm voxels):\n');
fprintf('  ROI  (keep)      : %6d voxels  (%.1f%% of head)\n', ...
        nVoxROI, 100*nVoxROI/nnz(headMask));
fprintf('  INT  (suppress)  : %6d voxels  (%.1f%% of head)\n', ...
        nVoxINT, 100*nVoxINT/nnz(headMask));
fprintf('  gap  (excluded)  : %6d voxels\n', nnz(gapMask));
fprintf('  head             : %6d voxels\n', nnz(headMask));

if nVoxROI < 200
    warning('ROI has only %d voxels - covariance estimate may be unstable.', nVoxROI);
end

%% Save

maskPath = fullfile(rovirDir, 'masks.mat');
meta = struct('subject_num', subject_num, 'zRange', zRange, 'gapVox', gapVox, ...
              'gridSize', N, 'eyeMaskSource', eyeMaskPath, 'created', datetime('now'));
save(maskPath, 'roiMask', 'intMask', 'headMask', 'gapMask', 'meta', '-v7.3');
fprintf('\nMasks saved: %s\n', maskPath);

%% QC figure: ROI red, interference blue, over the RMS image

qcDir = fullfile(rovirDir, 'qc');
if ~exist(qcDir, 'dir'); mkdir(qcDir); end

a = xrms / prctile(xrms(:), 99.5); a(a > 1) = 1;
up = 8; slices = min(zRange)-3 : max(zRange)+3;
slices = slices(slices >= 1 & slices <= N);
nC = 6; nR = ceil(numel(slices)/nC);
tile = zeros(nR*N*up, nC*N*up, 3, 'single');
for i = 1:numel(slices)
    k = slices(i); rr = floor((i-1)/nC); cc = mod(i-1,nC);
    base = kron(a(:,:,k), ones(up));
    R = base; G = base; B = base;
    eR = logical(kron(double(edgeOf(roiMask(:,:,k))), ones(up)));
    eI = logical(kron(double(edgeOf(intMask(:,:,k))), ones(up)));
    R(eI) = 0.1; G(eI) = 0.4; B(eI) = 1.0;
    R(eR) = 1.0; G(eR) = 0.1; B(eR) = 0.1;
    tile(rr*N*up+(1:N*up), cc*N*up+(1:N*up), 1) = R;
    tile(rr*N*up+(1:N*up), cc*N*up+(1:N*up), 2) = G;
    tile(rr*N*up+(1:N*up), cc*N*up+(1:N*up), 3) = B;
end
qcPath = fullfile(qcDir, sprintf('R1_masks_%s.png', subjectStr));
imwrite(tile, qcPath);
fprintf('QC figure : %s   (red = ROI, blue = interference)\n', qcPath);
fprintf('Slices shown: %d..%d, row-major, %d per row\n', slices(1), slices(end), nC);

% Zoomed view on the anterior half through the middle of the ROI band, where
% the question "did the box actually land on the globes?" is answerable.
zc = round(median(zRange));
r1 = 1:round(N*0.55); r2 = round(N*0.2):round(N*0.8);
up2 = 16; zs = zc-2 : zc+1;
h = numel(r1)*up2; w = numel(r2)*up2;
tile2 = zeros(h, numel(zs)*w, 3, 'single');
for i = 1:numel(zs)
    k = zs(i);
    base = kron(a(r1,r2,k), ones(up2));
    R = base; G = base; B = base;
    eR = logical(kron(double(edgeOf(roiMask(r1,r2,k))), ones(up2)));
    R(eR) = 1.0; G(eR) = 0.1; B(eR) = 0.1;
    tile2(:, (i-1)*w+(1:w), 1) = R;
    tile2(:, (i-1)*w+(1:w), 2) = G;
    tile2(:, (i-1)*w+(1:w), 3) = B;
end
qcPath2 = fullfile(qcDir, sprintf('R1_roi_zoom_%s.png', subjectStr));
imwrite(tile2, qcPath2);
fprintf('ROI zoom  : %s   (slices %d..%d, anterior crop)\n', qcPath2, zs(1), zs(end));

function e = edgeOf(m)
    if ~any(m(:)); e = false(size(m)); return; end
    e = m & ~imerode(m, strel('square', 3));
end
