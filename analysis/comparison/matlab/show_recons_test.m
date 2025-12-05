%% Load and display xrms.mat images from different subjects
% This script loads xrms.mat from each subject's recon results and displays
% them using bmImage

clear; clc, close all;

% Define base path and subject directories
basePath = '/home/debi/jaime/repos/MR-EyeTrack/data/251204/recon_results';
subjects = {'Sub001', 'Sub002', 'Sub003', 'Sub004', 'Sub005', 'Sub006', 'Sub007'};

% Create a figure for each subject
for i = 1:length(subjects)
    outputDir = fullfile(basePath, subjects{i}, 'output');
    
    % Check if output directory exists
    if ~isfolder(outputDir)
        fprintf('Warning: output directory not found for %s\n', subjects{i});
        continue;
    end
    
    % Find subdirectories in output folder
    subdirs = dir(outputDir);
    subdirs = subdirs([subdirs.isdir] & ~startsWith({subdirs.name}, '.'));
    
    xrmsFound = false;
    
    % Search for xrms.mat in subdirectories
    for j = 1:length(subdirs)
        subDir = fullfile(outputDir, subdirs(j).name);
        xrmsFile = fullfile(subDir, 'xrms.mat');
        
        if isfile(xrmsFile)
            % Load the xrms data
            data = load(xrmsFile);
            xrms = data.xrms;
                       
            % Display using bmImage
            bmImage(xrms);
            subDirName = strrep(subdirs(j).name, '_', '-');
            title(sprintf('%s - xrms (%s)', subjects{i}, subDirName));
            
            fprintf('Loaded and displayed xrms from %s (found in %s)\n', subjects{i}, subdirs(j).name);
            xrmsFound = true;
            break;
        end
    end
    
    if ~xrmsFound
        fprintf('Warning: xrms.mat not found in any subdirectory of %s\n', outputDir);
    end
end

%%
t = 8;
t_tot = (6*t + 30*4*t + 5*t)*1000; %ms
t_tot