%% Visual comparison with subplots
clc; clearvars; close all;
addpath(genpath('/home/debi/jaime/repos/MR-EyeTrack/results'));
addpath(genpath('/home/debi/MatTechLab/monalisa'));

%% Config

% Variables
subject_num = 1;
region_idx = 0; % 0:up 1:down 2:left 3:right 4:center mask
mask_type = {'clean', 'clean_0.5'};  % for the undersampled masks, the object is eMask, for the full masks, the object is eMaskN inside the struct

% Paths
resultsDir = '/home/debi/jaime/repos/MR-EyeTrack/results';
mask_1_path = [resultsDir, '/Sub00', num2str(subject_num), '/T1_LIBRE_Binning/other/', mask_type{1}, '/eMask_th0.75_region', num2str(region_idx), '.mat'];
mask_2_path = [resultsDir, '/Sub00', num2str(subject_num), '/T1_LIBRE_Binning/other/', mask_type{2}, '/eMask_th0.75_region', num2str(region_idx), '.mat'];

disp(['mask_1_path: ', mask_1_path]);
disp(['mask_2_path: ', mask_2_path]);

% Load images
m1 = load(mask_1_path, 'eMaskN');
m2 = load(mask_2_path, 'eMask');

%% Plot
figure;
plot(m1.eMaskN, '.-');
title('eMaskN 1D Plot - Clean');
xlabel('Index');
ylabel('Mask Value');

% figure;
% plot(m2.eMaskN, '.-');
% title('eMaskN 1D Plot - Filtered');
% xlabel('Index');
% ylabel('Mask Value');

%%
n = ceil(sqrt(length(m1.eMaskN)));       % closest square dimension
maskImg = zeros(1, n^2);
maskImg(1:length(m1.eMaskN)) = m1.eMaskN;   % pad to square
maskImg = reshape(maskImg, n, n);

figure;
imagesc(maskImg);
colormap(gray);
axis image off;
title('Binary Mask (reshaped view)');
