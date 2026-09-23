addpath(genpath('/home/debi/jaime/repos/mapVBVD'));
myTwix = mapVBVD_JB();
readouts   = myTwix{end}.image.unsorted(); % size: [480, 44, 81906]

%%
center_win = 48;
r1 = max(round((size(readouts,1)/2 -center_win/2)), 1);
r2 = min((size(readouts,1)/2 +center_win/2-1), size(readouts,1));
kCenter_range = r1:r2;
readouts_kCenter = readouts(kCenter_range, :, :);
