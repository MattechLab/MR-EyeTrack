% S2_batch_all_winLen - rebuild every readout bin at a chosen winLen, no prompts.
%
% The per-mask-type batch scripts (S2_eyeMask_t1_binning_pulseq_batch_*.m) each
% re-read the 16 GB raw .dat once per subject. Running all of them means ~90
% reads of ~1.4 TB. This does the same work with ONE read per subject: load the
% twix timestamps, then apply the identical windowing to every mask type.
%
% The windowing logic, the nShotOff handling and the mod(k,nSeg)==1 SI-segment
% rule are copied verbatim from the batch scripts, so only winLen (and th_ratio)
% differ from the existing bins. Outputs keep the established naming, which
% already carries the parameters:
%
%   bins/<mask_type>/eMask_th<th>_winLen<winLen>.mat        (single-class masks)
%   bins/<mask_type>/eMask_th<th>_winLen<winLen>_region<N>.mat  (4-way masks)
%
% so eMask_th0.75_winLen3.mat sits alongside the existing winLen10 files without
% overwriting them.
%
% Why winLen=3: at winLen=10 the window is 80 ms and a readout is kept only if
% 75% of it agrees, i.e. 60 ms. Saccades last 30-80 ms, so most can never be
% labelled at any position -- cohort-wide only 37.9% of expected saccade
% readouts survive, and per subject that ranges 4.5% to 58%, biased toward the
% longest events. Sustained classes keep 93-95%. winLen=3 is a 24 ms window.

clearvars; clc;
addpath(genpath('/home/debi/jaime/repos/MR-EyeTrack/recon'));
addpath(genpath('/home/debi/yiwei/forclone/pulseq'));
addpath(genpath('/home/debi/MatTechLab/monalisa'));

%% Parameters
winLen    = 3;
th_ratio  = 0.75;
applySync = true;                % apply the per-subject ET-MRI timing correction
syncFile  = '/home/debi/jaime/repos/MR-EyeTrack/data/study/sync_corrections.mat';
nShotOff  = 14;
costTime  = 2.5;                 % Siemens timestamp scaling
baseDir   = '/home/debi/jaime/repos/MR-EyeTrack/data/study';
subjects  = 1:15;

seqFolder = fullfile(baseDir, 'pulseq');
seqName   = 'yj_seq202_t1w_libre_main_TR8.0ms_TE3.6ms_swap1_FA6_RF2_mreye_track_trajPTP_44_1872_gdsp.seq';
seqParams = extract_seq_params(fullfile(seqFolder, seqName));
nSeg      = seqParams.nseg;

% mask_type -> ET mask filename pattern. Empty region list = single mask.
jobs = { ...
    'fixation-ok',   {'*_fixation_ok_mask.mat'},   {}; ...
    'saccade-only',  {'*_saccade_only_mask.mat'},  {}; ...
    'blink-event',   {'*_blink_event_mask.mat'},   {}; ...
    'tracking-loss', {'*_tracking_loss_mask.mat'}, {}; ...
    'no-mo',         {'*_nomo_mask.mat'},          {}; ...
    'clean',         {'*_clean_mask_0.3_%d.mat'},  {0,1,2,3}; ...
    'filtered',      {'*_ft_mask_%s.mat'},         {'up','down','left','right'}; ...
    'location',      {'*_lc_mask_%s.mat'},         {'up','down','left','right'} ...
};

% Per-subject timing correction (model/sync_table.py). The ET mask is indexed at
%   et_index = ratio * (offset_ms + TimeStamp_ms)
% where offset_ms is the elapsed ET time at the first readout -- about +61 ms,
% since the first readout lands ~55 ms after the trigger that released it -- and
% ratio is the EyeLink sample rate error, which spans -92 to +118 ppm across
% subjects and reaches 75 ms of accumulated shift by the end of a scan.
suffix = '';
if applySync
    S = load(syncFile);
    suffix = '_sync';
    fprintf('Applying per-subject sync correction from %s\n', syncFile);
end
fprintf('Rebuilding bins at winLen=%d, th=%.2f%s (existing files untouched)\n', ...
        winLen, th_ratio, suffix);

for subject_num = subjects
    subjectStr = sprintf('sub-%03d', subject_num);
    subjectDir = fullfile(baseDir, subjectStr);
    rawDir     = fullfile(subjectDir, 'rawdata');
    ETDir      = fullfile(subjectDir, 'eyemasks');
    reconDir   = fullfile(subjectDir, 'recon');

    fprintf('\n========== %s ==========\n', subjectStr);
    rawFiles = dir(fullfile(rawDir, '*_T1wLIBRE.dat'));
    if isempty(rawFiles)
        warning('No T1wLIBRE.dat for %s, skipping.', subjectStr); continue;
    end

    % ---- one raw read per subject, for the timestamps only ----
    twix_multi = mapVBVD_JH(fullfile(rawDir, rawFiles(1).name));
    if iscell(twix_multi); twix_obj = twix_multi{end}; else; twix_obj = twix_multi; end
    twix_obj.image.flagIgnoreSeg = true;
    TimeStamp    = double(twix_obj.image.timestamp);
    TimeStamp    = TimeStamp - min(TimeStamp);
    TimeStamp_ms = TimeStamp * costTime;
    Timediff     = TimeStamp_ms(end) - TimeStamp_ms(1);
    NLin         = numel(TimeStamp_ms);
    fprintf('  %d readouts, %.1f s\n', NLin, Timediff/1000);
    clear twix_multi twix_obj;

    offset_ms = 0; ratio = 1;
    if applySync
        ii = find(S.subjects == subject_num, 1);
        if isempty(ii)
            warning('No sync entry for %s, using uncorrected timing.', subjectStr);
        else
            offset_ms = S.offset_ms(ii); ratio = S.ratio(ii);
            fprintf('  sync: offset %+.1f ms, rate %+.0f ppm\n', ...
                    offset_ms, 1e6*(ratio-1));
        end
    end
    ET_ms = ratio * (offset_ms + TimeStamp_ms);

    for jj = 1:size(jobs,1)
        mask_type = jobs{jj,1};
        pattern   = jobs{jj,2}{1};
        regions   = jobs{jj,3};
        binsDir   = fullfile(reconDir, 'bins', mask_type);
        if ~isfolder(binsDir); mkdir(binsDir); end

        if isempty(regions); items = {[]}; else; items = regions; end
        for rr = 1:numel(items)
            if isempty(regions)
                pat = pattern; tag = mask_type; outName = ...
                    sprintf('eMask_th%.2f_winLen%i%s.mat', th_ratio, winLen, suffix);
            else
                it = items{rr};
                if ischar(it)
                    pat = sprintf(pattern, it); tag = it; ridx = rr-1;
                else
                    pat = sprintf(pattern, it); tag = num2str(it); ridx = it;
                end
                outName = sprintf('eMask_th%.2f_winLen%i%s_region%i.mat', ...
                                  th_ratio, winLen, suffix, ridx);
            end

            maskFiles = dir(fullfile(ETDir, pat));
            if isempty(maskFiles)
                fprintf('    [%s/%s] mask not found, skipped\n', mask_type, tag);
                continue;
            end
            raw_mask = load(fullfile(ETDir, maskFiles(1).name));
            if isstruct(raw_mask)
                fn = fieldnames(raw_mask); raw_mask = raw_mask.(fn{1});
            end
            m = padArrayWithZeros(raw_mask, round(Timediff));

            WinWidth     = round(numel(m) / NLin) * winLen;
            HalfWinWidth = floor(WinWidth / 2);
            nMeasuresOff = nShotOff * nSeg;
            binMaskVec   = zeros(NLin, 1);
            for k = 1:NLin
                if k <= nMeasuresOff
                    binMaskVec(k) = 0;
                else
                    timeSeg   = ET_ms(k);
                    win_lower = round(max(1, timeSeg - HalfWinWidth));
                    win_upper = round(min(numel(m), timeSeg + HalfWinWidth));
                    window_data = m(win_lower:win_upper);
                    th = win_upper - win_lower;
                    if sum(window_data) >= th * th_ratio
                        binMaskVec(k) = 1;
                    end
                end
                if mod(k, nSeg) == 1
                    binMaskVec(k) = 0;
                end
            end
            eMaskN = logical(binMaskVec)';
            save(char(fullfile(binsDir, outName)), 'eMaskN');
            fprintf('    [%s/%s] WinWidth %d ms, kept %d/%d -> %s\n', ...
                    mask_type, tag, WinWidth, sum(eMaskN), NLin, outName);
        end
    end
end

disp('Done. All mask types rebuilt at the chosen winLen.');
