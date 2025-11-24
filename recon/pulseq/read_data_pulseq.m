%% Init
clc, clearvars;

%% Paths
addpath(genpath('/Users/cag/Documents/forclone/pulseq_v15'));
addpath(genpath('/home/debi/yiwei/forclone/pulseq'));
addpath(genpath('/home/debi/yiwei/forclone/pulseqmreye'));

sequence_path = 'data/2025-11_pulseq-p2p/sequence+prescan/yj_seq2_t1w_libre_main_TR6.2ms_TE3.6ms_swap1_FA6_RF2_mreye_track_trajPTP_44_1872.seq';

if ~isfile(sequence_path)
    error('File not found: %s', sequence_path);
end

%% Display the contents of the .seq file
uiopen(sequence_path,1);

%% Load the sequence using Pulseq MATLAB API
% if the Pulseq MATLAB API is installed and on path:
seq = mr.Sequence();
seq.read(sequence_path);
disp('Sequence loaded via Pulseq API.');

%% Display sequence definitions if Pulseq API is available
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