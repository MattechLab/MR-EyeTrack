%% Compare left (region 2) vs right (region 3) — axial view
clc; clearvars; close all;
addpath(genpath('/home/debi/jaime/repos/MR-EyeTrack/recon'));
addpath(genpath('/home/debi/MatTechLab/monalisa'));
addpath(genpath('/home/debi/yiwei/forclone/pulseq'));

%% Config
subject_num = 1;
mask_type   = 'clean';
base_dir    = '/home/debi/jaime/repos/MR-EyeTrack';

% Paths
x_left_path  = fullfile(base_dir, sprintf('data/study/sub-%03d/recon/%s/x/x_steva_regionidx_2_nIter_20_delta_1.000.mat', subject_num, mask_type));
x_right_path = fullfile(base_dir, sprintf('data/study/sub-%03d/recon/%s/x/x_steva_regionidx_3_nIter_20_delta_1.000.mat', subject_num, mask_type));
disp(['left  path: ', x_left_path]);
disp(['right path: ', x_right_path]);

% Load and normalize together (joint scale)
x_left_raw  = load(x_left_path,  'x'); x_left  = x_left_raw.x;
x_right_raw = load(x_right_path, 'x'); x_right = x_right_raw.x;

[img_left, img_right] = norm_two_image(x_left, x_right);
img_both = cat(2, img_left, img_right);
bmImage(img_both)

%% Overlapping Axial
x_cell_ax = {norm_image(x_left, [0,0.7]), norm_image(x_right, [0,0.7])};
bmImage(x_cell_ax);

%% Overlapping Sagittal
x_cell_sag = {norm_image(rot90(permute(x_left,  [1,3,2]), 1), [0,0.7]), ...
               norm_image(rot90(permute(x_right, [1,3,2]), 1), [0,0.7])};
bmImage(x_cell_sag);

% Side-by-side axial (left | right per slice column)
ax_both = cat(2, img_left, img_right);

%% Overview — interleaved axial mosaic (left | right per column)
offset   = 8;
sl_start = 90 + offset;
inc      = 6;
sl_end   = 126 + offset;

mosaic = ax_both(:,:,sl_start);
for sl = (sl_start + inc):inc:sl_end
    mosaic = cat(2, mosaic, ax_both(:,:,sl));
end

figure('Color','white','Name','Left vs Right — Axial');
imshow(mosaic, []); title(sprintf('sub-%03d  left (left) | right (right)  axial', subject_num));

%% Difference map — axial
diff_ax = diff_volume(img_left, img_right);

diff_mosaic = diff_ax(:,:,sl_start);
for sl = (sl_start + inc):inc:sl_end
    diff_mosaic = cat(2, diff_mosaic, diff_ax(:,:,sl));
end

clim_val = max(abs(diff_ax(:))) * 0.75;
figure('Color','white','Name','Left minus Right — Axial diff');
imshow(diff_mosaic, []); colormap('redblue'); colorbar;
clim([-clim_val, clim_val]);
title(sprintf('sub-%03d  left − right  axial diff', subject_num));

%% Functions

function img_s = norm_image(img, clip_range)
if nargin < 2, clip_range = [0, 1]; end
img   = double(abs(img));
img_s = (img - min(img(:))) / (max(img(:)) - min(img(:)));
img_s = max(clip_range(1), min(clip_range(2), img_s));
end

function [img1_s, img2_s] = norm_two_image(img1, img2)
img1 = double(abs(img1));
img2 = double(abs(img2));
img1_s = (img1 - min(img1(:))) / (max(img1(:)) - min(img1(:)));
img2_s = (img2 - min(img2(:))) / (max(img2(:)) - min(img2(:)));
end

function diff_img = diff_volume(img1, img2)
[img1_s, img2_s] = norm_two_image(img1, img2);
diff_img = img1_s - img2_s;
fprintf('diff range: [%.4f, %.4f]\n', min(diff_img(:)), max(diff_img(:)));
end
