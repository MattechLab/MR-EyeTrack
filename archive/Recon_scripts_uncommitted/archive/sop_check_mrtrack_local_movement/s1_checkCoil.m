clc; clearvars; close all;

addpath(genpath('/home/debi/jaime/repos/MR-EyeTrack/recon/Recon_scripts'));
addpath(genpath('/home/debi/MatTechLab/internal_monalisa'));
addpath(genpath('/home/debi/yiwei/forclone/pulseq'));

subject_num = 5;
C_path = sprintf('data/study/sub-%03d/recon/C.mat', subject_num);
load(C_path)

%% Rotate C
% Define rotation rules: subject_num → rotation k (rot90)
rotationMap = containers.Map( ...
    [2, -1], ...  % [subject_num, rotation] (-1 is clockwise, +1 is counter-clockwise)
    [3, -1] ...
);

% Apply rotation if subject exists in the map
if isKey(rotationMap, subject_num)
    k = rotationMap(subject_num);
    C = rot90(C, k);
    disp(size(C));
    bmImage(sqrt(sum(C.^2, 4)));
end

%%
bmImage(C)
% Calculate the rms of C
crms = cal_crms(C);
bmImage(crms);

%%
% now I need to figure out the channels which affect the eye the most
% Here I'm generating xrms with the same size as C used for creating eye mask, 
% as I'm playing with the IDEA data, I include the IDEA recon procedure in
% gen_idea_xrms
% ignore it if you can generate xrms on your own
% don't use the function if the data is from pulseq.
% baseFolder = '/Users/cag/Documents/Dataset/datasets';
%  [xrms,x0] = gen_idea_xrms(baseFolder);

%% define a mask around eye region
xrms_path = sprintf('data/study/sub-%03d/recon/woBin/xrms48.mat', subject_num);
load(xrms_path, 'xrms');
nSlice = size(xrms,3);
figure; imagesc(abs(xrms(:,:,round(nSlice/2)))); axis image
title('Draw mask around the eyes');
eyeMask = roipoly;   % binary mask
eyeMaskPath = sprintf('data/study/sub-%03d/recon/mitosius/woBin_comp/eyeMask.mat', subject_num);
mkdir(fileparts(eyeMaskPath));
save(eyeMaskPath, "eyeMask");

%%
bmImage(eyeMask);
bmImage(xrms.*eyeMask)

%%
zRange = 12:36;
numCoils = size(C, 4);
eyeMask3D = repmat(eyeMask, [1 1 length(zRange)]);

rawWeights = zeros(numCoils,1);

for coil = 1:numCoils
    coilAbs = abs(C(:,:,zRange,coil));
    rawWeights(coil) = sum(coilAbs(eyeMask3D==1), 'all');
end

% normalize across coils
weights_norm = rawWeights / sum(rawWeights);
weights_path = sprintf('data/study/sub-%03d/recon/mitosius/woBin_comp/weights_norm.mat', subject_num);
save(weights_path, 'weights_norm', '-v7.3');
disp(['weights_norm has been saved here: ', weights_path]);

%%
[sortedW, idx] = sort(weights_norm, 'descend');
idx_path = sprintf('data/study/sub-%03d/recon/mitosius/woBin_comp/idx_coilSelection.mat', subject_num);
save(idx_path, 'idx', '-v7.3');
disp(['Coil selection indices saved here: ', idx_path]);

figure;
bar( idx, sortedW, 'LineWidth', 1.2);
xlabel('Sorted coil index');
ylabel('Weight in eye region');
title('Sorted Coil Contributions to Eye Region');
grid on;

% Print coil order
disp('Coils sorted by weight (highest to lowest):');
disp(idx(:)');
%%
numTop = 10;   % change if needed

figure;
bar(weights_norm, 'FaceColor', [0.6 0.6 0.6]); hold on;
bar(idx(1:numTop), weights_norm(idx(1:numTop)), 'FaceColor', [0.9 0.3 0.3]);
xlabel('Coil index');
ylabel('Weight in eye region');
title(['Top ' num2str(numTop) ' Coil Contributions: ',  strjoin(cellstr(num2str(idx(1:numTop))), ' ')]);
legend('All coils', 'Top contributors');

grid on;

%%
for select_coil = idx(1:numTop)
    bmImage(C(:,:,:,select_coil));
end

% 7    14     8    11    13     9     3    10     1     5 
% 44         480       81906
% %============================
% Even better: ROI-optimized virtual coils (beamforming for an ROI)
% 
% There’s a line of work that does exactly what you want: synthesize virtual coils that maximize signal from an ROI 
% and suppress the rest (and can also compress channels). 
% This is often framed as maximizing a signal-to-interference ratio for a specified ROI.  
% 
% Why it’s a great fit here:
% 	•	Your ROI (eyes/orbits) is tiny compared to the whole-head FOV.
% 	•	You can get a few virtual channels that are “eye-focused,” then run your central-48-sample energy/PCA navigator on those virtual channels (instead of all 44).
% 
% If you do this, you typically don’t need to guess which physical coils matter; the method learns the best linear combinations.
% Region-optimized virtual (ROVir) coils: Localization and/or suppression of spatial regions using sensor-domain beamforming
% https://usc-mrel.github.io/Journal/2021/Kim_2021_MRM.pdf
% https://people.eecs.berkeley.edu/~mlustig/software/GCCdemo/demo_GCC.html
% demo for GCC






