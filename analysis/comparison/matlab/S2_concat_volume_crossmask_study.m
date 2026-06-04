%% Init
clc; clearvars; close all;
addpath(genpath('/home/debi/jaime/repos/MR-EyeTrack/recon'));
addpath(genpath('/home/debi/MatTechLab/monalisa'));
addpath(genpath('/home/debi/yiwei/forclone/pulseq'));

%% Config
subject_num = 4;
region_idx  = 2;                    % 0=up 1=down 2=left 3=right
masks_to_compare = {'clean', 'filtered'};

% Paths
x_1_path = sprintf('data/study/sub-%03d/recon/%s/x/x_steva_regionidx_%d_nIter_20_delta_1.000.mat', subject_num, masks_to_compare{1}, region_idx);
x_2_path = sprintf('data/study/sub-%03d/recon/%s/x/x_steva_regionidx_%d_nIter_20_delta_1.000.mat', subject_num, masks_to_compare{2}, region_idx);
disp(['x_1_path: ', x_1_path]);
disp(['x_2_path: ', x_2_path]);

% Load images
x_1 = load(x_1_path, 'x');
x_2 = load(x_2_path, 'x');

% Prepare volumes for comparison
x1 = x_1.x;
x2 = x_2.x;

%% Normalization
[img1_trans, img2_trans] = norm_two_image(x1, x2);
img_1_2_trans = cat(2, img1_trans, img2_trans);
bmImage(img_1_2_trans)

% Sagittal permutation
img1_sag = rot90(permute(img1_trans, [1,3,2]), 1);
img2_sag = rot90(permute(img2_trans, [1,3,2]), 1);
img_1_2_sag = cat(2, img1_sag, img2_sag);

%% Overlapping
% Axial
x_cell_ax = {norm_image(x1), norm_image(x2)};
bmImage(x_cell_ax);
% Sagittal
x_cell_sag = {norm_image(rot90(permute(x1, [1,3,2]), 1)), norm_image(rot90(permute(x2, [1,3,2]), 1))};
bmImage(x_cell_sag);

%% Axial
close all;
offset = 8;
sl_start = 90+offset;
inc = 6;
sl_end = 126+offset;

show_image = img_1_2_trans(:,:,sl_start);
for slice = (sl_start+inc):inc:sl_end
    show_image = cat(2, show_image, img_1_2_trans(:,:,slice));
end
bmImage(show_image);

diff_map = diff_volume(img1_trans, img2_trans);
show_diff_image = diff_map(:,:,sl_start);
for slice_idx = sl_start+inc:inc:sl_end
    show_diff_image = cat(2, show_diff_image, diff_map(:,:,slice_idx));
end

figure('Color', 'white'); set(gca, 'Color', 'white');
imshow(show_diff_image);
colorbar; colormap('redblue'); caxis([-max(abs(diff_map(:))), max(abs(diff_map(:)))]);
title(sprintf('sub-%03d | region %d | %s vs %s | axial diff', subject_num, region_idx, masks_to_compare{1}, masks_to_compare{2}));

%% Sagittal
sl_start = 70;
inc = 6;
sl_end = 110;

show_image = img_1_2_sag(:,:,sl_start);
for slice = (sl_start+inc):inc:sl_end
    show_image = cat(2, show_image, img_1_2_sag(:,:,slice));
end
bmImage(show_image);

diff_map = diff_volume(img1_sag, img2_sag);
show_diff_image = diff_map(:,:,sl_start);
for slice_idx = sl_start+inc:inc:sl_end
    show_diff_image = cat(2, show_diff_image, diff_map(:,:,slice_idx));
end

figure('Color', 'white'); set(gca, 'Color', 'white');
imshow(show_diff_image);
colorbar; colormap('redblue'); caxis([-max(abs(diff_map(:)))/4*3, max(abs(diff_map(:)))/4*3]);
title(sprintf('sub-%03d | region %d | %s vs %s | sagittal diff', subject_num, region_idx, masks_to_compare{1}, masks_to_compare{2}));

%% Coronal
img1_cor = rot90(permute(img1_trans, [2,3,1]), 1);
img2_cor = rot90(permute(img2_trans, [2,3,1]), 1);
img_1_2_cor = cat(2, img1_cor, img2_cor);

sl_start = 94;
inc = 4;
sl_end = 110;

show_image = img_1_2_cor(:,:,sl_start);
for slice = (sl_start+inc):inc:sl_end
    show_image = cat(2, show_image, img_1_2_cor(:,:,slice));
end
bmImage(show_image);

diff_map = diff_volume(img1_cor, img2_cor);
show_diff_image = diff_map(:,:,sl_start);
for slice_idx = sl_start+inc:inc:sl_end
    show_diff_image = cat(2, show_diff_image, diff_map(:,:,slice_idx));
end

figure('Color', 'white'); set(gca, 'Color', 'white');
imshow(show_diff_image);
colorbar; colormap('redblue'); caxis([-max(abs(diff_map(:)))/4*3, max(abs(diff_map(:)))/4*3]);
title(sprintf('sub-%03d | region %d | %s vs %s | coronal diff', subject_num, region_idx, masks_to_compare{1}, masks_to_compare{2}));

%% Functions

function [img1_scaled] = norm_image(img1)
img1 = double(abs(img1));
img1_scaled = (img1 - min(img1(:))) / (max(img1(:)) - min(img1(:)));
disp('Scaling Done')
end

function [img1_scaled, img2_scaled] = norm_two_image(img1, img2)
img1 = double(abs(img1));
img2 = double(abs(img2));
img1_scaled = (img1 - min(img1(:))) / (max(img1(:)) - min(img1(:)));
img2_scaled = (img2 - min(img2(:))) / (max(img2(:)) - min(img2(:)));
disp('Scaling Done')
end

function [diff_img] = diff_volume(img1, img2)
[img1_s, img2_s] = norm_two_image(img1, img2);
diff_img = img1_s - img2_s;
disp(max(diff_img(:)))
disp(min(diff_img(:)))
end
