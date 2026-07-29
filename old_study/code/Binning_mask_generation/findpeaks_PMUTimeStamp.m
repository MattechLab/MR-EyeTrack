function [RRpeaks,RRperiods] = findpeaks_PMUTimeStamp(PMUTimeStamp,TimeStamp)
% FINDPEAKS_PMUTIMESTAMP
% Corrected version of the code to find the peaks in the PMUTimeStamp:
% it takes into account the "blind-spots" we have in the non-freerunning acquisitions
%
% INPUTS
% PMUTimeSTamp: temporal period from the last ECG trigger
% TimeStamp:    temporal period from the beginning of the acquisition
%
% OUTPUTS
% RRpeaks:      temporal position of the peaks
% RRperiods:    temporal interval between two consecutive peaks


% find the "lower peaks"
[ ~, idxPeaks ] = findpeaks( -PMUTimeStamp );

% subtract from each peak the time passed from the last TRUE cardiac trigger
RRpeaks = TimeStamp( idxPeaks ) - PMUTimeStamp( idxPeaks );

% compute the interval between two consecutive corrected peaks
RRperiods = diff( RRpeaks );

end
