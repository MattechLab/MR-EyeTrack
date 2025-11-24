%% Compute and Match Hashes between .dat and .seq files
% This script precomputes all hashes from both file types and matches them efficiently.
clearvars; clc;
dataDir = '/home/debi/jaime/repos/MR-EyeTrack/data/test-pulseq/pulseq/sub_yj';
% --- Optional toolbox paths (adjust if needed) ---
mapVBVD_path = '/home/debi/jaime/repos/mapVBVD';
pulseq_path  = '/home/debi/yiwei/forclone/pulseq';
addpath(genpath(mapVBVD_path));
addpath(genpath(pulseq_path));
%% --- Collect files ---
datFiles = dir(fullfile(dataDir, '*.dat'));
seqFiles = dir(fullfile(dataDir, '*.seq'));
fprintf('Found %d .dat and %d .seq files.\n', numel(datFiles), numel(seqFiles));
%% --- 1. Compute TWIX hashes from .dat files ---
datHashes = struct('file', {}, 'hash', {});
fprintf('\n--- Computing TWIX hashes from raw data ---\n');
for i = 1:numel(datFiles)
    f = fullfile(datFiles(i).folder, datFiles(i).name);
    try
        twixObj = mapVBVD_JB(f);
        twix_idx = 1;
        seqHash_twix = char(twixObj{twix_idx}.hdr.Dicom.tSequenceVariant);
        while numel(seqHash_twix) ~= 32
            twix_idx = twix_idx + 1;
            seqHash_twix = char(twixObj{twix_idx}.hdr.Dicom.tSequenceVariant);
        end
        datHashes(end+1).file = datFiles(i).name; %#ok<SAGROW>
        datHashes(end).hash = lower(seqHash_twix);
        fprintf('✅ %s — %s\n', datFiles(i).name, seqHash_twix);
    catch ME
        fprintf('⚠️  Could not read hash for %s: %s\n', datFiles(i).name, ME.message);
        datHashes(end+1).file = datFiles(i).name; %#ok<SAGROW>
        datHashes(end).hash = '';
    end
end
%% --- 2. Compute MD5 hashes from .seq files ---
seqHashes = struct('file', {}, 'hash', {});
fprintf('\n--- Computing MD5 hashes from sequence files ---\n');
for i = 1:numel(seqFiles)
    f = fullfile(seqFiles(i).folder, seqFiles(i).name);
    try
        rawSeq = fileread(f);
        sigPos = strfind(rawSeq, '[SIGNATURE]');
        if ~isempty(sigPos)
            rawSeq = rawSeq(1, 1:sigPos(1)-2);
        end
        md = java.security.MessageDigest.getInstance('MD5');
        md.update(uint8(rawSeq));
        hashRaw = md.digest();
        hashStr = lower(reshape(dec2hex(typecast(hashRaw,'uint8'))', 1, []));
        seqHashes(end+1).file = seqFiles(i).name; %#ok<SAGROW>
        seqHashes(end).hash = hashStr;
        fprintf('✅ %s — %s\n', seqFiles(i).name, hashStr);
    catch ME
        fprintf('⚠️  Could not compute hash for %s: %s\n', seqFiles(i).name, ME.message);
        seqHashes(end+1).file = seqFiles(i).name; %#ok<SAGROW>
        seqHashes(end).hash = '';
    end
end
%% --- 3. Match hashes ---
fprintf('\n--- Matching hashes ---\n');
matches = struct('dat', {}, 'seq', {}, 'hash', {});
for i = 1:numel(datHashes)
    dh = datHashes(i).hash;
    if isempty(dh)
        fprintf('⚠️  %s has no valid hash.\n', datHashes(i).file);
        continue;
    end
    matchIdx = find(strcmpi({seqHashes.hash}, dh), 1);
    if ~isempty(matchIdx)
        matches(end+1).dat = datHashes(i).file; %#ok<SAGROW>
        matches(end).seq = seqHashes(matchIdx).file;
        matches(end).hash = dh;
        fprintf('✅ %s ↔ %s\n', datHashes(i).file, seqHashes(matchIdx).file);
    else
        fprintf('❌ No matching .seq for %s (hash=%s)\n', datHashes(i).file, dh);
    end
end
fprintf('\n=== SUMMARY ===\n');
for i = 1:numel(matches)
    fprintf('%-35s ↔ %-50s\n', matches(i).dat, matches(i).seq);
end