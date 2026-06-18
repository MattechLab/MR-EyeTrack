%% Compare up (region 0) vs down (region 1) — sagittal view
clc; clearvars; close all;
addpath(genpath('/home/debi/jaime/repos/MR-EyeTrack/recon'));
addpath(genpath('/home/debi/MatTechLab/monalisa'));
addpath(genpath('/home/debi/yiwei/forclone/pulseq'));

%% Config
subject_num = 1;
mask_type   = 'clean';
base_dir    = '/home/debi/jaime/repos/MR-EyeTrack';

% Paths
x_up_path   = fullfile(base_dir, sprintf('data/study/sub-%03d/recon/%s/x/x_steva_regionidx_0_nIter_20_delta_1.000.mat', subject_num, mask_type));
x_down_path = fullfile(base_dir, sprintf('data/study/sub-%03d/recon/%s/x/x_steva_regionidx_1_nIter_20_delta_1.000.mat', subject_num, mask_type));
disp(['up   path: ', x_up_path]);
disp(['down path: ', x_down_path]);

% Load and normalize together (joint scale)
x_up_raw   = load(x_up_path,   'x'); x_up   = x_up_raw.x;
x_down_raw = load(x_down_path, 'x'); x_down = x_down_raw.x;

[img_up, img_down] = norm_two_image(x_up, x_down);
img_both = cat(2, img_up, img_down);
bmImage(img_both)

%% Overlapping Axial
x_cell_ax = {norm_image(x_up, [0,0.7]), norm_image(x_down, [0,0.7])};
bmImage(x_cell_ax);

%% Overlapping Sagittal
x_cell_sag = {norm_image(rot90(permute(x_up,   [1,3,2]), 1), [0,0.7]), ...
               norm_image(rot90(permute(x_down, [1,3,2]), 1), [0,0.7])};
bmImage(x_cell_sag);

%% Permute to sagittal orientation
% permute([row, col, slice]) -> [row, slice, col], then rot90 to orient S-I axis vertically
sag_up   = rot90(permute(img_up,   [1, 3, 2]), 1);
sag_down = rot90(permute(img_down, [1, 3, 2]), 1);

sag_both = cat(2, sag_up, sag_down);

%% Overview — interleaved sagittal mosaic (up | down per column)
sl_start = 70;
inc      = 6;
sl_end   = 118;

mosaic = [sag_both(:,:,sl_start)];
for sl = (sl_start + inc):inc:sl_end
    mosaic = cat(2, mosaic, sag_both(:,:,sl));
end

figure('Color','white','Name','Up vs Down — Sagittal');
imshow(mosaic, []); title(sprintf('sub-%03d  up (left) | down (right)  sagittal', subject_num));

%% Difference map — sagittal
diff_sag = diff_volume(sag_up, sag_down);

diff_mosaic = diff_sag(:,:,sl_start);
for sl = (sl_start + inc):inc:sl_end
    diff_mosaic = cat(2, diff_mosaic, diff_sag(:,:,sl));
end

clim_val = max(abs(diff_sag(:))) * 0.75;
figure('Color','white','Name','Up minus Down — Sagittal diff');
imshow(diff_mosaic, []); colormap('redblue'); colorbar;
clim([-clim_val, clim_val]);
title(sprintf('sub-%03d  up − down  sagittal diff', subject_num));

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
