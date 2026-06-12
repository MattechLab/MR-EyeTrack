% =====================================================
% Batch binning script — saccade_only mask
% Based on S2_eyeMask_t1_binning_pulseq_batch_nomo.m
% Runs subjects 1-15, mask_type 'saccade-only'.
% Generates one eMask per subject from the saccade_only ET mask.
% =====================================================
clearvars; clc;
addpath(genpath('/home/debi/jaime/repos/MR-EyeTrack/recon'));
addpath(genpath('/home/debi/jaime/repos/MR-EyeTrack/recon/Recon_scripts'));
addpath(genpath('/home/debi/yiwei/forclone/pulseq'));
addpath(genpath('/home/debi/MatTechLab/monalisa'));

%% Config
mask_type = 'saccade-only';
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

    %% Auto-detect saccade_only mask file
    maskFiles = dir(fullfile(ETDir, '*_saccade_only_mask.mat'));
    if isempty(maskFiles)
        warning('saccade_only mask not found for %s, skipping.', subjectStr);
        continue;
    end
    maskPath = fullfile(ETDir, maskFiles(1).name);
    fprintf('Mask [saccade_only]: %s\n', maskFiles(1).name);

    raw_mask = load(maskPath);
    if isstruct(raw_mask)
        fn = fieldnames(raw_mask);
        raw_mask = raw_mask.(fn{1});
    end
    mask_method_1 = padArrayWithZeros(raw_mask, round(Timediff));

    WinWidth     = round(length(mask_method_1) / NLin) * winLen;
    HalfWinWidth = floor(WinWidth / 2);
    fprintf('  WinWidth: %d\n', WinWidth);

    %% Build binning vector
    nMeasuresOff = nShotOff * nSeg;
    binMaskVec   = zeros(NLin, 1);

    for k = 1:NLin
        if k <= nMeasuresOff
            binMaskVec(k) = 0;
        else
            timeSeg   = TimeStamp_ms(k);
            win_lower = round(max(1, timeSeg - HalfWinWidth));
            win_upper = round(min(numel(mask_method_1), timeSeg + HalfWinWidth));
            window_data = mask_method_1(win_lower:win_upper);
            true_count  = sum(window_data);
            th          = win_upper - win_lower;
            if true_count >= th * th_ratio
                binMaskVec(k) = 1;
            end
        end
        if mod(k, nSeg) == 1
            binMaskVec(k) = 0;
        end
    end

    fprintf('  [saccade_only] preserved: %d / %d readouts\n', sum(binMaskVec), NLin);

    %% Save
    binsDirStr    = string(binsDir);
    eMaskFilePath = fullfile(binsDirStr, sprintf('eMask_th%.2f_winLen%i.mat', th_ratio, winLen));
    eMaskN        = logical(binMaskVec)';
    save(char(eMaskFilePath), 'eMaskN');
    fprintf('Saved: %s\n', char(eMaskFilePath));
end

disp('Batch saccade_only binning complete.');
