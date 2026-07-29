clc; clear; close all;

%% Import data: Samples
% Events indicate significant events during recordings (blink, saccades,
% etc). Samples are all the recorded (t,x,y) points.

% %Syntax to load samples file
%fid = fopen('JR_samples_20181113.txt','rt');
%datacell = textscan(fid,'%s %s        %s  %s %s  %s  %s', 'HeaderLines', 1); %Reproduce the delimiter pattern of the .txt rows
%fclose(fid);

load('JR_samples_20181113.mat');

trial_id = str2double(datacell{2}); %trial_identification
X_coord = str2double(datacell{3});
Y_coord = str2double(datacell{4});

%Provisory clean up of the NaN, which will be processed later
X_coord(isnan(X_coord)) = 0;
Y_coord(isnan(Y_coord)) = 0;
%% Import trigger and trial info from Result.txt output file of the eyetracker -in phase with the LOG

fid = fopen('JR_20181113.txt','rt');
Trigger_info = textscan(fid,'%s %s       %s     %s %s %s', 'HeaderLines', 1); %HeaderLine get rid of the first line of acquisition of a .txt file
fclose(fid);

tot_numb_trial = length(str2double(Trigger_info{2}));
trigger_time = str2double(Trigger_info{3});
trial_start = str2double(Trigger_info{4});
trial_end = str2double(Trigger_info{5});
trial_duration = str2double(Trigger_info{6}); 

%% Identification of index along the table at which the trial changes
tri_chang = trial_id(2:end) - trial_id(1:end-1);

% +1 indicates the first point of the new trial, concateneting to the
% beginning of the experiment and first trial - This is also a counter for
% the number of trials
index_change_trial = [1; find(tri_chang ~= 0) + 1];

%Insert a fourth column indicating the duration
vect_duration = zeros(length(X_coord),1);
vect_duration(index_change_trial) = trial_duration;

%% Identification of trials referring to the same presentation

mean_X = zeros(length(index_change_trial),1);
mean_Y = zeros(length(index_change_trial),1);

for i = 1:length(index_change_trial)
   
   %Take X/Y coordinate for each trial different from 0 and compute the
   %mean
   relative_indexes = find((trial_id == i) & (X_coord ~= 0));
   mean_X(i) = mean(X_coord(relative_indexes));
   mean_Y(i) = mean(Y_coord(relative_indexes));
   
   %Substitute the zero values with the mean 
   zero_indexes = find((trial_id == i) & (X_coord == 0));
   X_coord(zero_indexes) = mean_X(i);
   Y_coord(zero_indexes) = mean_Y(i);

   %assert((length(relative_indexes)+length(zero_indexes)) == trial_duration(i));
   
end

assert(length(find(X_coord == 0)) == 0);

figure; scatter(mean_X, mean_Y);
Mean_coordinates = [mean_X, mean_Y];
%% Categorize the Trials according to the center of the stimuli
% These parameters need to be taken from the protocol

%First is the center fixation, then from top left to bottom right in
%row-order
center_x = [400; 133; 309; 491; 670; 133; 309; 491; 670; 133; 309; 491; 670; 133; 309; 491; 670];
center_y = [300; 90; 90; 90; 90; 227; 227; 227; 227; 364; 364; 364; 364; 504; 504; 504; 504];

number_of_stimuli = 17;
stimuli_centers = [center_x, center_y];

indexes_clustering = kmeans(Mean_coordinates, number_of_stimuli,'Start',stimuli_centers);

err = 0;
number_of_present_per_stimulus = 6;
%We exclude the first presentation from error computation

figure;
hold on
scatter(mean_X(indexes_clustering == 1), mean_Y(indexes_clustering == 1))
text(mean_X(indexes_clustering == 1), mean_Y(indexes_clustering == 1), num2str(1),'VerticalAlignment',...
        'bottom','HorizontalAlignment','right');
    
for i=2:number_of_stimuli
    
    num = sum(indexes_clustering == i);
    if num < number_of_present_per_stimulus
        err = err + (number_of_present_per_stimulus - num);
    end
    hold on;
    scatter(mean_X(indexes_clustering == i), mean_Y(indexes_clustering == i))
    text(mean_X(indexes_clustering == i), mean_Y(indexes_clustering == i), num2str(i),'VerticalAlignment',...
        'bottom','HorizontalAlignment','right');
end

set(gca,'Ydir','reverse')
%percentage error
err = err/tot_numb_trial*100;
%% Construction of the table

%Save the mean (X,Y) coordinates
vect_Mean_X = zeros(length(X_coord),1);
vect_Mean_X(index_change_trial) = mean_X;

vect_Mean_Y = zeros(length(X_coord),1);
vect_Mean_Y(index_change_trial) = mean_Y;

%Identification of the stimulus based on k-means clustering
vect_ID_stimuli = zeros(length(X_coord),1);
vect_ID_stimuli(index_change_trial) = indexes_clustering;

%first column X coord, second column Y coord, third column Trial label
data = [X_coord, Y_coord, vect_ID_stimuli, vect_duration, vect_Mean_X, vect_Mean_Y, trial_id];
