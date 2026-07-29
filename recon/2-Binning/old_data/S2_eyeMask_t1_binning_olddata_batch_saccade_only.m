% =====================================================
% Batch binning script (LEGACY 2018-2019 3T data) — saccade_only mask
% Based on S2_eyeMask_t1_binning_pulseq_batch_saccade_only.m
% Runs the 9 legacy subjects that have raw data, mask_type 'saccade-only'.
% Generates one eMask per subject from the saccade_only ET mask.
%
% Differences w.r.t. the study version:
%   - raw data comes from old_study/raw_data/<Subj_folder>/meas_*.dat
%     (no sub-NNN/rawdata layout, no *_T1wLIBRE.dat naming)
%   - there is no pulseq .seq file for this acquisition, so nSeg is derived
%     from the twix header instead of extract_seq_params (expected: 22)
%   - the ET masks live in the legacy filer tree, matched per subject via
%     subjectMap below (folder names encode the subject initials)
%   - outputs go to old_study/new_binning/<Subj_folder>/<mask_type>/
% The sliding-window binning itself is unchanged (winLen/th_ratio identical).
% =====================================================
clearvars; clc;
addpath(genpath('/home/debi/jaime/repos/MR-EyeTrack/recon'));
addpath(genpath('/home/debi/jaime/repos/MR-EyeTrack/recon/Recon_scripts'));

%% Config
mask_type  = 'saccade-only';
maskGlob   = '*_saccade_only_mask.mat';

repoDir = '/home/debi/jaime/repos/MR-EyeTrack';
rawRoot = fullfile(repoDir, 'old_study', 'raw_data');
outRoot = fullfile(repoDir, 'old_study', 'new_binning');
ETRoot  = ['/mnt/filer01/MatTechLab/benedetta.franceschiello/_backup_filerch/' ...
           'Benedetta_Franceschiello/1_MR_Eye_Retina/8_Results/' ...
           '1_Motion_resolved_MR_Eye_3T/1_Motion_Resolved_Anatomical_GRE/' ...
           '2_Eye_Tracker_Results_Anatomical/2. ET_3T_Data'];

% raw_data folder  <->  legacy ET folder (the two trees use different naming)
% Subjects 'mo' (20181113) and 'psr' (20190325) have ET data but no raw data.
subjectMap = {
    'Subj_20181113_1J', '20181113_exp_jr'
    'Subj_20190225_1B', '20190225_exp_bm'
    'Subj_20190225_2N', '20190225_exp_nn'
    'Subj_20190228_1E', '20190228_exp_em'
    'Subj_20190228_2N', '20190228_exp_nt'
    'subj_20190301_1G', '20190301_exp_gr'
    'Subj_20190320_1A', '20190320_exp_am'
    'Subj_20190320_2K', '20190320_exp_kg'
    'Subj_20190325_1L', '20190325_exp_la'
    };
nSubjects = size(subjectMap, 1);

nShotOff        = 14;   % warm-up shots discarded (kept from the study pipeline)
nSegFallback    = 22;   % legacy LIBRE nDeye: 22 segments x 3723 shots
winLen          = 10;
th_ratio        = 0.75;
costTime        = 2.5;  % Siemens timestamp scaling factor

summary = cell(nSubjects, 1);

%% Loop over subjects
for subject_num = 1:nSubjects
    subjectStr = subjectMap{subject_num, 1};
    etStr      = subjectMap{subject_num, 2};

    rawDir  = fullfile(rawRoot, subjectStr);
    ETDir   = fullfile(ETRoot, etStr);
    binsDir = fullfile(outRoot, subjectStr, mask_type, filesep);

    fprintf('\n========== Subject %d/%d (%s <- %s) ==========\n', ...
        subject_num, nSubjects, subjectStr, etStr);

    if ~isfolder(binsDir)
        mkdir(binsDir);
        fprintf('Created: %s\n', binsDir);
    end

    %% Auto-detect the LIBRE raw data file
    rawFiles = dir(fullfile(rawDir, 'meas_*_LIBRE_*.dat'));
    if isempty(rawFiles)
        warning('No LIBRE .dat found for %s, skipping.', subjectStr);
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

    % No .seq file for the legacy acquisition: Siemens NSeg holds the shot
    % count here, so segments-per-shot comes from NLin / nShot.
    nShot = double(twix_obj.image.NSeg);
    nSeg  = round(NLin / nShot);
    if nSeg * nShot ~= NLin
        warning('%s: NLin=%d not divisible by nShot=%d, falling back to nSeg=%d.', ...
            subjectStr, NLin, nShot, nSegFallback);
        nSeg = nSegFallback;
    end

    fprintf('Duration: %.1f ms, %d readouts (nSeg=%d, nShot=%d)\n', ...
        Timediff, NLin, nSeg, nShot);

    %% Auto-detect saccade_only mask file
    maskFiles = dir(fullfile(ETDir, maskGlob));
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
    raw_mask = double(raw_mask(:))';   % masks are saved as int64 (1 x N)
    ET_len   = numel(raw_mask);
    fprintf('  ET samples: %d (MRI duration %.0f ms)\n', ET_len, Timediff);

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

    nKept = sum(binMaskVec);
    fprintf('  [saccade_only] preserved: %d / %d readouts (%.1f%%)\n', ...
        nKept, NLin, 100 * nKept / NLin);

    %% Save
    binsDirStr    = string(binsDir);
    eMaskFilePath = fullfile(binsDirStr, sprintf('eMask_th%.2f_winLen%i.mat', th_ratio, winLen));
    eMaskN        = logical(binMaskVec)';
    save(char(eMaskFilePath), 'eMaskN');
    fprintf('Saved: %s\n', char(eMaskFilePath));

    summary{subject_num} = sprintf('%s,%s,%s,%d,%d,%d,%.1f,%d,%d,%d,%.2f', ...
        subjectStr, etStr, rawFiles(1).name, NLin, nSeg, nShot, ...
        Timediff, ET_len, WinWidth, nKept, 100 * nKept / NLin);
end

%% Cross-subject summary
summary = summary(~cellfun(@isempty, summary));
csvPath = fullfile(outRoot, sprintf('summary_%s.csv', mask_type));
fid = fopen(csvPath, 'w');
fprintf(fid, 'subject,et_folder,raw_file,NLin,nSeg,nShot,duration_ms,ET_len,WinWidth,kept,kept_pct\n');
fprintf(fid, '%s\n', summary{:});
fclose(fid);
fprintf('\nSummary: %s\n', csvPath);

disp('Batch saccade_only binning (legacy data) complete.');
