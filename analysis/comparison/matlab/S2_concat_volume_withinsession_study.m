%% Init
clc; clearvars; close all;
addpath(genpath('/home/debi/jaime/repos/MR-EyeTrack/results'));
addpath(genpath('/home/debi/jaime/repos/MR-EyeTrack/recon/Recon_scripts'));
addpath(genpath('/home/debi/MatTechLab/internal_monalisa'));
addpath(genpath('/home/debi/yiwei/forclone/pulseq'));

%% Config

% Paths
x_1_path = 'data/study/data/sub-004/recon/woBin/xrms.mat';
x_2_path = 'data/study/data/sub-006/recon/woBin/xrms.mat';
disp(['x_1_path: ', x_1_path]);
disp(['x_2_path: ', x_2_path]);

% Load images
x_1 = load(x_1_path, 'xrms');
x_2 = load(x_2_path, 'xrms');

% Prepare volumes for comparison
x1 = x_1.xrms;
x2 = x_2.xrms;
% x{2} = x_gdsp;
% x{3} = flip(flip(permute(x_idea, [2,1,3]),1),3);
% x2 = flip(permute(x0904, [2 1 3]), 2);

%% Normalization
% close all;

% Normalization 1 image
% x1_trans = norm_image(x1);
% x1_sag = norm_image(rot90(permute(x1_trans, [1,3,2]), 1));
% x1_coronal = norm_image(permute(x1_trans, [3,2,1]));
% bmImage(x1_trans);
% bmImage(x1_sag);
% bmImage(x1_coronal);

% x2_trans = norm_image(x2);
% x2_sag = norm_image(rot90(permute(x2_trans, [1,3,2]),1));
% x2_coronal = norm_image(permute(x2_trans, [3,2,1]));
% bmImage(x2_trans);
% bmImage(x2_sag);
% bmImage(x2_coronal);

% [img1_trans, img2_trans] = norm_two_image(x1_trans, x2_trans);
% img_1_2_trans = cat(2, img1_trans, img2_trans);
% bmImage(img_1_2_trans)

% Normalization 2 images
[img1_trans, img2_trans] = norm_two_image(x1, x2);
img_1_2_trans = cat(1, img1_trans, img2_trans);
bmImage(img_1_2_trans)

%% Axial
close all;
offset = 8;
sl_start = 90+offset;
inc=6;
sl_end=126+offset;

% show_image = img1_trans(:,:,sl_start); %90:6:126
% show_image = img2_trans(:,:,sl_start);%90+14:6:126+14
% show_image = img2_trans(:,:,sl_start);%90+14:6:126+14
show_image = img_1_2_trans(:,:,sl_start);

for slice = (sl_start+inc):inc:sl_end
    show_image = cat(2,[show_image, img_1_2_trans(:,:,slice)]);
    % show_image = cat(2,[show_image, img2_trans(:,:,slice)]);
end

bmImage(show_image);

diff_map = diff_volume(img1_trans,img2_trans);
show_diff_image = diff_map(:,:,sl_start);
for slice_idx = sl_start+inc:inc:sl_end
    show_diff_image = cat(2,[show_diff_image, diff_map(:,:,slice_idx)]);
end

figure('Color', 'white'); set(gca, 'Color', 'white'); 
imshow(show_diff_image);
colorbar; colormap('redblue'); caxis([-max(abs(diff_map(:))), max(abs(diff_map(:)))]);

% loLev = -1;upLev = 0.1;[imClip, rgb_vec] = relaxationColorMap('T1', show_diff_image, loLev, upLev);
% figure('Color', 'white'); set(gca, 'Color', 'white'); 
% imshow(imClip, 'DisplayRange', [loLev, upLev], 'InitialMagnification', 'fit'); 
% % title(strcat('difference map-', num2str(slice_idx)));
% colormap(rgb_vec); colorbar;

%% Sagittal

sl_start = 70;
inc=6;
sl_end=110;

show_image = img_1_2_sag(:,:,sl_start);

for slice = (sl_start+inc):inc:sl_end
    show_image = cat(2,[show_image, img_1_2_sag(:,:,slice)]);
end

bmImage(show_image);

diff_map = diff_volume(img1_sag,img2_sag);
show_diff_image = diff_map(:,:,sl_start);
for slice_idx = sl_start+inc:inc:sl_end
show_diff_image = cat(2,[show_diff_image, diff_map(:,:,slice_idx)]);
end

figure('Color', 'white'); set(gca, 'Color', 'white'); 
imshow(show_diff_image);
colorbar;colormap('redblue'); caxis([-max(abs(diff_map(:)))/4*3, max(abs(diff_map(:)))/4*3]);

%% Coronal
sl_start = 94;
inc=4;
sl_end=110;

show_image = img_1_2_sag(:,:,sl_start);

for slice = (sl_start+inc):inc:sl_end
    show_image = cat(2,[show_image, img_1_2_sag(:,:,slice)]);
end
bmImage(show_image);


diff_map = diff_volume(img1_sag,img2_sag);
show_diff_image = diff_map(:,:,sl_start);
for slice_idx = sl_start+inc:inc:sl_end
show_diff_image = cat(2,[show_diff_image, diff_map(:,:,slice_idx)]);
end

figure('Color', 'white'); set(gca, 'Color', 'white'); 
imshow(show_diff_image);
colorbar;colormap('redblue'); caxis([-max(abs(diff_map(:)))/4*3, max(abs(diff_map(:)))/4*3]);


%% Functions

function [img1_scaled] = norm_image(img1)
img1 = double(abs(img1));
% Scale img1 to [0, 1]
img1_scaled = (img1 - min(img1(:))) / (max(img1(:)) - min(img1(:)));
disp('Scaling Done')
end

function [img1_scaled, img2_scaled] = norm_two_image(img1, img2)
img1 = double(abs(img1));
img2 = double(abs(img2));
% Scale img1 to [0, 1]
img1_scaled = (img1 - min(img1(:))) / (max(img1(:)) - min(img1(:)));
% Scale img2 to [0, 1]
img2_scaled = (img2 - min(img2(:))) / (max(img2(:)) - min(img2(:)));
disp('Scaling Done')
end

function [diff_img]=diff_volume(img1,img2)
[img1_s,img2_s] = norm_two_image(img1,img2);
diff_img = img1_s-img2_s;
disp(max(diff_img(:)))
disp(min(diff_img(:)))
end

function V_rot = rotate_volume(V,theta)
V = double(abs(V));
% theta: degree
% positive: counterclockwise; negative: clockwise
    V_rot = zeros(size(V));
    for k = 1:size(V,3)
        V_rot(:,:,k) = imrotate(V(:,:,k), theta, 'bilinear', 'crop');  % or 'loose'
    end
end