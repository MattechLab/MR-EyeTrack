% =====================================================
% Batch version of S2_eyeMask_t1_binning_pulseq.m
% Runs subjects 1-15, mask_type 'location', no user prompts.
% Directions: up(0), down(1), left(2), right(3), center(4)
% =====================================================
clearvars; clc;
addpath(genpath('/home/debi/jaime/repos/MR-EyeTrack/recon'));
addpath(genpath('/home/debi/jaime/repos/MR-EyeTrack/recon/Recon_scripts'));
addpath(genpath('/home/debi/yiwei/forclone/pulseq'));
addpath(genpath('/home/debi/MatTechLab/monalisa'));

%% Config
mask_type = 'location';
baseDir   = '/home/debi/jaime/repos/MR-EyeTrack/data/study';

seqFolder    = fullfile(baseDir, 'pulseq');
seqName_list = {
    'yj_seq100_t1w_libre_pre_TR6.2ms_TE3.6ms_swap1_FA6_Traj1_nSeg44_nShot191_Fid0_mreye_2p0.seq', ...
    'yj_seq202_t1w_libre_main_TR8.0ms_TE3.6ms_swap1_FA6_RF2_mreye_track_trajPTP_44_1872_gdsp.seq'};

seqFile   = seqFolder + "/" + seqName_list{2};
seqParams = extract_seq_params(seqFile);

nShotOff        = 14;
nSeg            = seqParams.nseg;
winLen          = 10;
th_ratio        = 0.75;
costTime        = 2.5;  % Siemens timestamp scaling factor

direction_names = {'up', 'down', 'left', 'right'};
region_idx_list = 0:3;  % 0=up 1=down 2=left 3=right

%% Loop over subjects
for subject_num = 1:15
    subjectStr = sprintf('sub-%03d', subject_num);
    subjectDir = fullfile(baseDir, subjectStr);
    rawDir     = fullfile(subjectDir, 'rawdata');
    ETDir      = fullfile(subjectDir, 'eyemasks');
    reconDir   = fullfile(subjectDir, 'recon');
    binsDir    = fullfile(reconDir, 'bins', mask_type, filesep);

    fprintf('\n========== Subject %d (%s) ==========\n', subject_num, subjectStr);

    if ~isfolder(binsDir)
        mkdir(binsDir);
        fprintf('Created: %s\n', binsDir);
    end

    %% Auto-detect T1wLIBRE raw data file
    rawFiles = dir(fullfile(rawDir, '*_T1wLIBRE.dat'));
    if isempty(rawFiles)
        warning('No T1wLIBRE.dat found for %s, skipping.', subjectStr);
        continue;
    end
    rawDataPath = fullfile(rawDir, rawFiles(1).name);
    fprintf('Raw data: %s\n', rawFiles(1).name);

    %% Load twix to extract timestamps
    twix_obj_multi = mapVBVD_JH(rawDataPath);
    if iscell(twix_obj_multi)
        twix_obj = twix_obj_multi{end};
    else
        twix_obj = twix_obj_multi;
    end
    twix_obj.image.flagIgnoreSeg = true;

    PMUTimeStamp    = double(twix_obj.image.pmutime);
    TimeStamp       = double(twix_obj.image.timestamp);
    TimeStamp       = TimeStamp - min(TimeStamp);
    TimeStamp_ms    = TimeStamp * costTime;
    Timediff        = TimeStamp_ms(end) - TimeStamp_ms(1);
    NLin            = length(PMUTimeStamp);

    fprintf('Duration: %.1f ms, %d readouts\n', Timediff, NLin);

    %% Build binning matrix for all 5 directions
    binMaskMatrix = zeros(NLin, 5);
    nMeasuresOff  = nShotOff * nSeg;

    for idx_dir = 1:5
        dir_name  = direction_names{idx_dir};
        maskFiles = dir(fullfile(ETDir, sprintf('*_lc_mask_%s.mat', dir_name)));

        if isempty(maskFiles)
            warning('No lc_mask_%s found for %s, leaving zeros.', dir_name, subjectStr);
            continue;
        end
        maskPath = fullfile(ETDir, maskFiles(1).name);
        fprintf('Mask [%s]: %s\n', dir_name, maskFiles(1).name);

        raw_mask = load(maskPath);
        if isstruct(raw_mask)
            fn = fieldnames(raw_mask);
            raw_mask = raw_mask.(fn{1});
        end
        mask_method_1 = padArrayWithZeros(raw_mask, round(Timediff));

        WinWidth     = round(length(mask_method_1) / NLin) * winLen;
        HalfWinWidth = floor(WinWidth / 2);
        fprintf('  WinWidth: %d\n', WinWidth);

        for k = 1:NLin
            if k <= nMeasuresOff
                binMaskMatrix(k, idx_dir) = 0;
            else
                timeSeg   = TimeStamp_ms(k);
                win_lower = round(max(1, timeSeg - HalfWinWidth));
                win_upper = round(min(numel(mask_method_1), timeSeg + HalfWinWidth));
                window_data = mask_method_1(win_lower:win_upper);
                true_count  = sum(window_data);
                th          = win_upper - win_lower;
                if true_count >= th * th_ratio
                    binMaskMatrix(k, idx_dir) = 1;
                end
            end
            if mod(k, nSeg) == 1
                binMaskMatrix(k, idx_dir) = 0;
            end
        end

        sum_binning = sum(binMaskMatrix(:, idx_dir));
        fprintf('  [%s] preserved: %d / %d readouts\n', dir_name, sum_binning, NLin);
    end

    cMask = logical(binMaskMatrix)';  % 5 x NLin

    %% Save one .mat file per region
    binsDirStr = string(binsDir);
    for region_idx = region_idx_list
        single_eMask  = cMask(region_idx + 1, :);
        fileName      = sprintf('eMask_th%.2f_region%i.mat', th_ratio, region_idx);
        eMaskFilePath = fullfile(binsDirStr, fileName);
        eMaskN        = single_eMask;
        save(char(eMaskFilePath), 'eMaskN');
        fprintf('Saved: %s\n', char(eMaskFilePath));
    end
end

disp('Batch location binning complete.');
