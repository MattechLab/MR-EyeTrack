clc; clearvars; close all;

addpath(genpath('/home/debi/jaime/repos/MR-EyeTrack/recon'));
addpath(genpath('/home/debi/MatTechLab/monalisa'));
addpath(genpath('/home/debi/yiwei/forclone/pulseq'));
addpath(genpath('/home/debi/jaime/repos/MR-EyeTrack/analysis'));

%% Init

% Parameters
subject_num = 6;

% Base directory
baseDir = '/home/debi/jaime/repos/MR-EyeTrack/data/study';

% Construct subject folder name (zero-padded to 3 digits)
subjectStr = sprintf('sub-%03d', subject_num);

% Full path to dataset directories
subjectDir  = fullfile(baseDir, subjectStr);
rawDir      = fullfile(subjectDir, 'rawdata');
reconDir    = fullfile(subjectDir, 'recon');

% Get the list of dat files
files = dir(fullfile(rawDir, 'sub-*.dat'));

% Identify files by pattern
bodyCoilFile  = fullfile(rawDir, dir(fullfile(rawDir, '*_BC.dat')).name);
arrayCoilFile = fullfile(rawDir, dir(fullfile(rawDir, '*_HC.dat')).name);
measureFile   = fullfile(rawDir, dir(fullfile(rawDir, '*_T1wLIBRE.dat')).name);

% Display or use them
disp('Found files:');
disp(bodyCoilFile);
disp(arrayCoilFile);
disp(measureFile);

% Pulseq
seqFolder = fullfile(baseDir, 'pulseq');
seqName_list = {
    'yj_seq100_t1w_libre_pre_TR6.2ms_TE3.6ms_swap1_FA6_Traj1_nSeg44_nShot191_Fid0_mreye_2p0.seq', ...
    'yj_seq202_t1w_libre_main_TR8.0ms_TE3.6ms_swap1_FA6_RF2_mreye_track_trajPTP_44_1872_gdsp.seq'};


%% Main
path_prescan_BC  = bodyCoilFile;
path_prescan_HC  = arrayCoilFile;
path_prescan_seq = fullfile(seqFolder, seqName_list{1});

path_mainscan_HC  = measureFile;
path_mainscan_seq = fullfile(seqFolder, seqName_list{2});

measFiles = { ...
    path_prescan_BC, ...
    path_prescan_HC, ...
    path_mainscan_HC};

seqFiles = { ...
    path_prescan_seq, ...
    path_prescan_seq, ...
    path_mainscan_seq};

readers = cell(1,3);

for i = 1:3
    readers{i} = createRawDataReader(measFiles{i}, true);
    readers{i}.acquisitionParams.traj_type = 'pulseq';
    readers{i}.acquisitionParams.pulseqTrajFile_name = seqFiles{i};

    seqParams = extract_seq_params(seqFiles{i});
    if isfield(seqParams,'nshot')
        readers{i}.acquisitionParams.nShot = seqParams.nshot;
    end
    if isfield(seqParams,'nseg')
        readers{i}.acquisitionParams.nSeg = seqParams.nseg;
    end
end

% Assign meaningful names
readerBC   = readers{1};
readerHC   = readers{2};
readerMain = readers{3};

checkPrescanOrientation(readerBC, readerHC, readerMain)
