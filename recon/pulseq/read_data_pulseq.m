%%
clc;
sequence_path = '/home/debi/jaime/repos/MR-EyeTrack/data/test-pulseq/pulseq/sub_yj/yj_seq2_t1w_libre_main_TR8.0ms_TE3.6ms_swap1_FA4_RF2_mreye_track_3723_traj_ptp.seq';

if ~isfile(sequence_path)
    error('File not found: %s', sequence_path);
end

% uiopen(sequence_path,1);

% Print the .seq file contents to the terminal
% txt = fileread(sequence_path);
% fprintf('%s', txt);

% If you prefer a streamed print (for very large files), use:
% fid = fopen(sequence_path, 'r');
% assert(fid ~= -1, 'Cannot open file: %s', sequence_path);
% c = onCleanup(@() fclose(fid));
% while ~feof(fid)
%     fprintf('%s', fgets(fid));
% end

% Optional: if the Pulseq MATLAB API is installed and on path:
% seq = mr.Sequence();
% seq.read(sequence_path);
% disp('Sequence loaded via Pulseq API.');

if exist('mr.Sequence','class')
    if ~exist('seq','var') || ~isa(seq,'mr.Sequence')
        seq = mr.Sequence();
        seq.read(sequence_path);
    end
    if isprop(seq,'definitions')
        fprintf('Sequence definitions:\n');
        k = keys(seq.definitions);
        v = values(seq.definitions);
        for i = 1:length(k)
            fprintf('  %s: %s\n', k{i}, mat2str(v{i}));
        end
    else
        warning('seq.definitions not available.');
    end
else
    warning('Pulseq API (mr.Sequence) not found on path.');
end