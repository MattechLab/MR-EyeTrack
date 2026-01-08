%% Init
clc, clearvars;

%% Paths
addpath(genpath('/home/debi/yiwei/forclone/pulseq'));

% sequence_path = 'data/study/pulseq/yj0_seq8_t1w_libre_pre_TR6.2ms_TE3.6ms_swap1_FA4_RF2_rfmod2_trajPTP_nSeg88_nShot89.seq';
% sequence_path = 'data/study/pulseq/yj_seq2_t1w_libre_main_TR6.2ms_TE3.6ms_swap1_FA6_RF2_mreye_track_trajPTP_44_1872.seq';
% sequence_path = 'data/study/pulseq/yj_seq100_t1w_libre_pre_TR6.2ms_TE3.6ms_swap1_FA6_Traj1_nSeg44_nShot191_Fid0_mreye_2p0.seq';
sequence_path = 'data/study/pulseq/yj_seq202_t1w_libre_main_TR8.0ms_TE3.6ms_swap1_FA6_RF2_mreye_track_trajPTP_44_1872_gdsp.seq';

if ~isfile(sequence_path)
    error('File not found: %s', sequence_path);
end

%% Display the contents of the .seq file
% This script opens the file in Matlab
% uiopen(sequence_path,1);

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
