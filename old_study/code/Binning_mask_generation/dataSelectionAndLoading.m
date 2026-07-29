function [ twix_obj, param ] = dataSelectionAndLoading( basedir, param, new_mapVBVD_flag )
%DATASELECTIONANDLOADING
%
% INPUTS
% basedir:  starting directory for rawdata selection
%
% OUTPUTS
% twix_obj: twix object
% param:    parameters structure


if nargin < 3
    new_mapVBVD_flag = 1;
end


if isempty( param.batchParam )
    %% Select rawdata
    [rawDataName, rawDataDir, ~] = uigetfile( ...
        { '*.dat','Siemens raw data file (*.dat)'}, ...
           'Pick a file', ...
           'MultiSelect', 'off',basedir);

        if isempty(rawDataName)
            error('Not file selected');
        end
        
else
    rawDataName = param.batchParam.rawDataName;
    rawDataDir  = param.batchParam.rawDataDir; 
    
end

        filepathRawData  = fullfile(rawDataDir,rawDataName);

        param.name      = rawDataName;
        param.dir       = rawDataDir;
        param.filepath  = filepathRawData;
    
    
    
    
%% Read the twix object
if new_mapVBVD_flag
    twix_obj_multi = mapVBVD_JH(filepathRawData); % twix_obj = mapVBVDsimo(filepathRawData);
else
    twix_obj_multi = mapVBVD(filepathRawData); % twix_obj = mapVBVDsimo(filepathRawData);
end

if iscell(twix_obj_multi)
    twix_obj = twix_obj_multi{2};
else
    twix_obj = twix_obj_multi;
end

twix_obj.image.flagIgnoreSeg = true;                        % Essential to read the data correctly



%% Initialize parameters
param.Np       = double(twix_obj.image.NCol);               % Number of readout point per line
param.Nshot    = double(twix_obj.image.NSeg);               % Total number of shots
param.Nseg     = double(twix_obj.image.NLin/param.Nshot);   % Number of segments per each shot
param.Nlines   = double(twix_obj.image.NLin);               % Total number of lines
param.Ncoil    = double(twix_obj.image.NCha);               % Number of coils

param.SIprojIdx = 1:param.Nseg:(param.Nseg*param.Nshot);    % SI projections indices



end
