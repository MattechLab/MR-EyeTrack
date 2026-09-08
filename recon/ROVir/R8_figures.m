%% R8 - Labelled figure set for the README
%
% Regenerates the QC images as self-explanatory panels with titles and row
% labels, and copies them into recon/ROVir/figures/ so they sit next to the
% code and render in the README (data/ is not tracked).
%
% Outputs -> recon/ROVir/figures/fig[1-5]_*.png

clc; clearvars -except subject_num nv; close all;

addpath(genpath('/home/debi/jaime/repos/MR-EyeTrack/recon'));
addpath(genpath('/home/debi/MatTechLab/monalisa'));

if ~exist('subject_num', 'var'); subject_num = 15; end
if ~exist('nv',          'var'); nv          = 20; end

repoRoot   = '/home/debi/jaime/repos/MR-EyeTrack';
subjectStr = sprintf('sub-%03d', subject_num);
reconDir   = fullfile(repoRoot, 'data/study', subjectStr, 'recon');
rovirDir   = fullfile(reconDir, 'ROVir');
figDir     = fullfile(repoRoot, 'recon/ROVir/figures');
if ~exist(figDir, 'dir'); mkdir(figDir); end

L = @(p, f) getfield(load(p, f), f);
un = @(v) subsref(v, substruct('{}', {1}));   % unwrap cell if needed
grab = @(p, f) local_grab(p, f);

M    = load(fullfile(rovirDir, 'masks.mat'));
T    = load(fullfile(rovirDir, 'rovir_transform.mat'), 'V', 'transform');
tr   = T.transform;
N48  = size(M.roiMask, 1);

%% ---------- Fig 1: the two regions ----------------------------------------
xrms = abs(single(L(fullfile(reconDir, 'woBin/xrms48.mat'), 'xrms')));
a = xrms / prctile(xrms(:), 99.5); a(a > 1) = 1;

zs = 18:3:30;
f = figure('Visible','off','Position',[100 100 1300 760]);
tl = tiledlayout(2, numel(zs), 'Padding','compact','TileSpacing','compact');
for i = 1:numel(zs)
    nexttile(i); imshow(a(:,:,zs(i)), []); title(sprintf('slice %d', zs(i)));
    if i==1, ylabel('image'); end
end
for i = 1:numel(zs)
    nexttile(numel(zs)+i);
    rgb = repmat(a(:,:,zs(i)), [1 1 3]);
    rr = M.roiMask(:,:,zs(i)); ii = M.intMask(:,:,zs(i));
    R=rgb(:,:,1); G=rgb(:,:,2); B=rgb(:,:,3);
    R(ii)=0.15*R(ii)+0.10; G(ii)=0.35*G(ii)+0.25; B(ii)=0.35*B(ii)+0.85;   % blue = suppress
    R(rr)=0.35*R(rr)+0.85; G(rr)=0.15*G(rr)+0.10; B(rr)=0.15*B(rr)+0.10;   % red  = keep
    rgb(:,:,1)=R; rgb(:,:,2)=G; rgb(:,:,3)=B;
    imshow(rgb);
end
title(tl, {'Fig 1 - the two ROVir regions on the 48^3 grid', ...
    'RED = eye ROI (kept)   BLUE = "everything else" = head mask minus dilated ROI (suppressed)   black gap between them is excluded from both'}, ...
    'FontWeight','bold');
exportgraphics(f, fullfile(figDir,'fig1_regions.png'), 'Resolution', 120); close(f);
fprintf('fig1_regions.png\n');

%% ---------- Fig 2: coil-image RSS, before vs after ------------------------
x0c = L(fullfile(reconDir, sprintf('woBin/x0_noC_%d.mat', N48)), 'x0');
nCh = numel(x0c);
X = zeros(N48^3, nCh, 'single');
for c = 1:nCh, X(:,c) = x0c{c}(:); end
clear x0c;
Vret = orth(T.V(:, 1:nv));

rssO = reshape(sqrt(sum(abs(X).^2, 2)), [N48 N48 N48]);
rssV = reshape(sqrt(sum(abs(X*Vret).^2, 2)), [N48 N48 N48]);
rssV = rssV * (mean(rssO(M.roiMask)) / mean(rssV(M.roiMask)));   % match on the ROI
hi = prctile(rssO(M.headMask), 99);

zs2 = 20:2:28;
f = figure('Visible','off','Position',[100 100 1250 620]);
tl = tiledlayout(2, numel(zs2), 'Padding','compact','TileSpacing','compact');
for i=1:numel(zs2)
    nexttile(i); imshow(min(rssO(:,:,zs2(i))/hi,1), []); title(sprintf('slice %d', zs2(i)));
    if i==1, rowlab(sprintf('all %d coils', nCh)); end
end
for i=1:numel(zs2)
    nexttile(numel(zs2)+i); imshow(min(rssV(:,:,zs2(i))/hi,1), []);
    if i==1, rowlab(sprintf('%d ROVir ch', nv)); end
end
title(tl, {'Fig 2 - coil-image energy (RSS) before vs after ROVir', ...
    'RSS = sqrt(sum_c |image_c|^2) over channels: how much signal the coils see at each voxel.', ...
    'Both rows scaled to match inside the ROI, so what you see is the head being suppressed.'}, ...
    'FontWeight','bold');
exportgraphics(f, fullfile(figDir,'fig2_coil_rss.png'), 'Resolution', 120); close(f);
fprintf('fig2_coil_rss.png\n');
clear X rssO rssV;

%% ---------- Fig 3 + 4: reconstructions -----------------------------------
tag = sprintf('ROVir_%d_woBin', nv);
xs = grab(fullfile(rovirDir, sprintf('x_steva_%s_nIter_20_delta_1.000.mat', tag)), 'x');
xr = grab(fullfile(rovirDir, sprintf('x0_%s.mat', tag)), 'x0');
xo = grab(fullfile(reconDir, 'woBin/x_steva_nIter_20_delta_1.000.mat'), 'x');
xs=abs(single(xs)); xr=abs(single(xr)); xo=abs(single(xo));
Nb = size(xs,1);
roiB = bmImResize(single(M.roiMask), [N48 N48 N48], [Nb Nb Nb]) > 0.5;
mt = @(v) v * (mean(xo(roiB))/mean(v(roiB)));
xr = mt(xr); xs = mt(xs);
hiB = prctile(xo(roiB), 99.5);

zsB = round((22:2:26) * Nb/N48);
r1 = 1:round(Nb*0.42); r2 = round(Nb*0.22):round(Nb*0.78);
rows = {xo, xr, xs};
lab  = {sprintf('52-ch STEVA (reference)'), ...
        sprintf('ROVir %d-ch, gridded x0', nv), ...
        sprintf('ROVir %d-ch, STEVA', nv)};
f = figure('Visible','off','Position',[100 100 1050 900]);
tl = tiledlayout(3, numel(zsB), 'Padding','compact','TileSpacing','compact');
for j=1:3
  for i=1:numel(zsB)
    nexttile((j-1)*numel(zsB)+i);
    imshow(min(rows{j}(r1,r2,zsB(i))/hiB,1), []);
    if i==1
        text(-0.04, 0.5, lab{j}, 'Units','normalized', 'Rotation',90, ...
             'HorizontalAlignment','center','FontWeight','bold','FontSize',10);
    end
    if j==1, title(sprintf('slice %d', zsB(i))); end
  end
end
title(tl, {'Fig 3 - the orbits: three reconstructions of the same acquisition', ...
    'Row 2 is ROVir before regularisation (noisy); row 3 is the same data after STEVA.'}, ...
    'FontWeight','bold');
exportgraphics(f, fullfile(figDir,'fig3_orbit_compare.png'), 'Resolution', 130); close(f);
fprintf('fig3_orbit_compare.png\n');

zsW = round((20:3:32) * Nb/N48);
f = figure('Visible','off','Position',[100 100 1300 560]);
tl = tiledlayout(2, numel(zsW), 'Padding','compact','TileSpacing','compact');
for i=1:numel(zsW)
    nexttile(i); imshow(min(xo(:,:,zsW(i))/hiB,1), []); title(sprintf('slice %d', zsW(i)));
    if i==1, rowlab('52-ch STEVA'); end
end
for i=1:numel(zsW)
    nexttile(numel(zsW)+i); imshow(min(xs(:,:,zsW(i))/hiB,1), []);
    if i==1, rowlab(sprintf('ROVir %d-ch STEVA', nv)); end
end
title(tl, {'Fig 4 - whole head: what the compression actually costs', ...
    'Orbits and face reconstruct cleanly; the posterior brain is discarded (it returns as noise, not signal).'}, ...
    'FontWeight','bold');
exportgraphics(f, fullfile(figDir,'fig4_wholehead.png'), 'Resolution', 120); close(f);
fprintf('fig4_wholehead.png\n');

%% ---------- Fig 6: why theirs is black and ours is noisy ------------------
%
% MR RawDeface reconstructs by root-sum-of-squares over the virtual channels
% (visualization.py: image_rsos = np.sqrt(np.sum(np.abs(image)**2, axis=3))) -
% no coil sensitivity model anywhere. In RSS the voxel value carries the coil
% sensitivity as a multiplicative factor, so where the retained virtual coils
% see nothing, the image is genuinely ~0: black.
%
% Our recon solves y = E(C_virt) m for m, so the same C_virt that suppressed
% the head is divided back out. Where it is ~0 the inversion is ill-posed and
% the region returns as amplified noise instead of zeros.
%
% Reapplying the sensitivity weighting to our own STEVA result reproduces their
% look from our data:  RSS-equivalent ~ |x| .* sqrt(sum_j |C_virt,j|^2)
load(fullfile(rovirDir, sprintf('C_rovir_%d.mat', nv)), 'C_rovir');
rssC = sqrt(sum(abs(bmImResize(C_rovir, [48 48 48], [Nb Nb Nb])).^2, 4));
xrss = xs .* rssC;

intB = bmImResize(single(M.intMask), [N48 N48 N48], [Nb Nb Nb]) > 0.5;
fprintf('\nposterior/head mean intensity, relative to the orbits:\n');
fprintf('  52-ch STEVA reference    : %.3f\n', mean(xo(intB))/mean(xo(roiB)));
fprintf('  ROVir STEVA (SENSE-like) : %.3f\n', mean(xs(intB))/mean(xs(roiB)));
fprintf('  ROVir RSS-equivalent     : %.3f\n', mean(xrss(intB))/mean(xrss(roiB)));

hiR = prctile(xrss(roiB), 99.5);
f = figure('Visible','off','Position',[100 100 1300 800]);
tl = tiledlayout(3, numel(zsW), 'Padding','compact','TileSpacing','compact');
rows6 = {xo, xs, xrss};
his   = [hiB, hiB, hiR];
lab6  = {'52-ch STEVA', sprintf('ROVir %d-ch STEVA', nv), sprintf('ROVir %d-ch, RSS-equivalent', nv)};
for j=1:3
  for i=1:numel(zsW)
    nexttile((j-1)*numel(zsW)+i);
    imshow(min(rows6{j}(:,:,zsW(i))/his(j),1), []);
    if i==1, rowlab(lab6{j}); end
    if j==1, title(sprintf('slice %d', zsW(i))); end
  end
end
title(tl, {'Fig 6 - why MR RawDeface shows black and we show noise', ...
    'Row 2 divides the coil sensitivity back out (our recon solves for magnetisation), so the discarded region returns as amplified noise.', ...
    'Row 3 re-applies that sensitivity weighting to the SAME result - the RSS view the abstract uses - and the discarded region goes black.'}, ...
    'FontWeight','bold');
exportgraphics(f, fullfile(figDir,'fig6_black_vs_noise.png'), 'Resolution', 120); close(f);
fprintf('fig6_black_vs_noise.png\n');

%% ---------- Fig 5: channel-count curves -----------------------------------
if exist(fullfile(rovirDir,'qc','R7_nv_curve.png'), 'file')
    copyfile(fullfile(rovirDir,'qc','R7_nv_curve.png'), fullfile(figDir,'fig5_nv_curve.png'));
    fprintf('fig5_nv_curve.png (copied)\n');
end

fprintf('\nFigures in %s\n', figDir);

function rowlab(txt)
    text(-0.06, 0.5, txt, 'Units','normalized', 'Rotation',90, ...
         'HorizontalAlignment','center', 'FontWeight','bold', 'FontSize',9);
end

function v = local_grab(p, f)
    S = load(p, f); v = S.(f);
    if iscell(v); v = v{1}; end
end
