%% init
addpath(genpath('/home/debi/jaime/repos/MR-EyeTrack/recon'));
addpath(genpath('/home/debi/jaime/repos/MR-EyeTrack/recon/Recon_scripts'));

%% location
mask_center = load('/home/debi/jaime/repos/MR-EyeTrack/data/pilot/masks_1206/Sub001/subject_1_lc_mask_center.mat');
mask_left = load('/home/debi/jaime/repos/MR-EyeTrack/data/pilot/masks_1206/Sub001/subject_1_lc_mask_left.mat');
mask_right = load('/home/debi/jaime/repos/MR-EyeTrack/data/pilot/masks_1206/Sub001/subject_1_lc_mask_right.mat');
mask_up = load('/home/debi/jaime/repos/MR-EyeTrack/data/pilot/masks_1206/Sub001/subject_1_lc_mask_up.mat');
mask_down = load('/home/debi/jaime/repos/MR-EyeTrack/data/pilot/masks_1206/Sub001/subject_1_lc_mask_down.mat');

%% filtered
mask_center = load('/home/debi/jaime/repos/MR-EyeTrack/data/pilot/masks_1206/Sub001/subject_1_ft_mask_center.mat');
mask_left = load('/home/debi/jaime/repos/MR-EyeTrack/data/pilot/masks_1206/Sub001/subject_1_ft_mask_left.mat');
mask_right = load('/home/debi/jaime/repos/MR-EyeTrack/data/pilot/masks_1206/Sub001/subject_1_ft_mask_right.mat');
mask_up = load('/home/debi/jaime/repos/MR-EyeTrack/data/pilot/masks_1206/Sub001/subject_1_ft_mask_up.mat');
mask_down = load('/home/debi/jaime/repos/MR-EyeTrack/data/pilot/masks_1206/Sub001/subject_1_ft_mask_down.mat');

%% check
disp(['nnz of center mask: ' num2str(nnz(mask_center.array))]);
disp(['nnz of left mask: ' num2str(nnz(mask_left.array))]);
disp(['nnz of right mask: ' num2str(nnz(mask_right.array))]);
disp(['nnz of up mask: ' num2str(nnz(mask_up.array))]);
disp(['nnz of down mask: ' num2str(nnz(mask_down.array))]);
total_nnz = nnz(mask_center.array) + nnz(mask_left.array) + nnz(mask_right.array) + nnz(mask_up.array) + nnz(mask_down.array);
disp(['Total nnz of all masks: ' num2str(total_nnz)]);
total_elements = numel(mask_center.array);
disp(['Total elements in one mask: ' num2str(total_elements)]);
disp('-----------------------------------');
