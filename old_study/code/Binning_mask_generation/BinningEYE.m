%--------------------------------------------------------------------------
% custom change to import mapVBVD: not sure its correct, no time to test it.
% Just need to be sure you are able to acces the mapVBVD_JH function for example.
%--------------------------------------------------------------------------
addpath('./../MATLAB/ReadRawDataSiemens/mapVBVD')
%--------------------------------------------------------------------------
% Original code from the function BinningEYE.m
%--------------------------------------------------------------------------
    param.basedir = '/Users/mauroleidi/Desktop/toDelete/EyeTracker.mat';
    basedir = param.basedir;
    param.batchParam = [];

%--------------------------------------------------------------------------
% Prompt user to select Siemens Raw Data file
%--------------------------------------------------------------------------

    [rawDataName, rawDataDir, ~] = uigetfile( ...
    { '*.dat','Siemens raw data file (*.dat)'}, ...
       'Pick a file', ...
       'MultiSelect', 'off',basedir);

    if rawDataName == 0
        warning('No file selected');
        return;
    end

    filepathRawData = fullfile(rawDataDir, rawDataName);

%--------------------------------------------------------------------------
% Initialize saving directory
%--------------------------------------------------------------------------    
  % Create saving directory
    [~,name,~] = fileparts(rawDataName);
    
  % Add number of bin to saving directory name
    name = sprintf('BinningEYE_%s',name);
    
  % Define directory name to save the results of the ICA analysis
    param.savedir = fullfile(rawDataDir,name);
    
  % Create directory
    if exist(param.savedir,'dir') ~= 7
        mkdir(param.savedir);
    end

%--------------------------------------------------------------------------
% Prompt user to enter the number of desired bins
%--------------------------------------------------------------------------           
    prompt        = {'Enter the number of bins'};
    name          = '#Bins';
    numlines      = 1;
    defaultanswer = {'4'};
    answer        = inputdlg(prompt,name,numlines,...
                           defaultanswer);
    
  % Check if the user selected cancel
    if isempty(answer)
        warning('The user selected cancel');
        return;
    end
    
  % Convert string answer to number
    nbins = str2double(answer{1});
    param.nBins = nbins;
    

%--------------------------------------------------------------------------
% Read the PMUTimeStamp triggered from external trigger
%--------------------------------------------------------------------------
    
    costTime = 2.5;
    param.batchParam.rawDataName    = rawDataName;
    param.batchParam.rawDataDir     = rawDataDir;

    [ twix_obj, param ] = dataSelectionAndLoading( basedir, param );
    
    PMUTimeStamp    = double( twix_obj.image.pmutime );
    TimeStamp       = double( twix_obj.image.timestamp );
    TimeStamp       = TimeStamp - min(TimeStamp);

    PMUTimeStamp_ms = PMUTimeStamp * costTime;
    PMUTimeStamp_s  = PMUTimeStamp_ms / 1000;
    TimeStamp_ms    = TimeStamp * costTime;
    TimeStamp_s     = TimeStamp_ms / 1000;

    param.PMUTimeStamp_ms   = PMUTimeStamp_ms;
    param.PMUTimeStamp_s    = PMUTimeStamp_s;
    param.TimeStamp_ms      = TimeStamp_ms;
    param.TimeStamp_s       = TimeStamp_s;
    
    [ triggerPeaks_ms, triggerIntervals_ms ] = findpeaks_PMUTimeStamp( param.PMUTimeStamp_ms, param.TimeStamp_ms );
    triggerPeaks_ms = [ -1 triggerPeaks_ms ];
    triggerPeaks_s      = triggerPeaks_ms       / 1000;
%     triggerIntervals_s  = triggerIntervals_ms   / 1000;
    nTrig = length( triggerPeaks_ms );
    
    
%--------------------------------------------------------------------------
% Read the sync-box data
%--------------------------------------------------------------------------
    
    % Prompt user to select respiratory binning directory
    [ syncboxDataName, syncboxDataDir, ~] = uigetfile( ...
    { '*.mat','Syncbox data file (*.mat)'}, ...
       'Pick a file', ...
       'MultiSelect', 'off',rawDataDir);

    if syncboxDataName == 0
        error('No syncbox data file were selected');
    end
    
    syncboxDataPath = fullfile(syncboxDataDir,syncboxDataName);
    
    struct = load(syncboxDataPath);
    syncboxMatrix = struct.data;
    
    
%--------------------------------------------------------------------------
% 
%--------------------------------------------------------------------------

%     x_mean = syncboxMatrix( :, 5 );
%     y_mean = syncboxMatrix( :, 6 );
%     x_mean = x_mean( x_mean~=0  );
%     y_mean = y_mean( y_mean~=0  );
    
    indexVec = syncboxMatrix( :, 3 );
    indexVec = indexVec( indexVec ~= 0 );
    
%     indexVec = indexVec(3:end);
    
    nPres = length(indexVec);
    
    indexVec_x = [];
    for k = 1:nPres
        elem = indexVec(k);
        if elem  ~= 1
            res = mod((elem-1),4);
            if res == 0
                res = 4;
            end
        else
            res = 0;
        end
        indexVec_x = [ indexVec_x res ];
    end
    
    figure, plot( indexVec_x, '*')
        ylim([0 5])
        xlim([0 nPres+1])
        xlabel('# Presentation')
        ylabel('Bin')
        title('Binning along X')
    
    indexVec_y = [];
    for k = 1:nPres
        elem = indexVec(k);
        if elem  ~= 1
            ceilVal = ceil((elem-1)/4);
        else
            ceilVal = 0;
        end
        indexVec_y = [ indexVec_y ceilVal ];
    end
    
    figure, plot( indexVec_y, 'o')
        ylim([0 5])
        title('Binning along Y')
    
    figure, plot( indexVec_x, '*')
    hold on, plot( indexVec_y, 'o')
        ylim([0 5])
        xlim([0 nPres+1])
        xlim([0 nPres+1])
        xlabel('# Presentation')
        ylabel('Bin')
        title('Binning along X and Y')
        legend('X','Y')
    
    
    meanTR = mean( diff(param.TimeStamp_ms) );
    figure, plot(diff(param.TimeStamp_ms),'.-')
        title(['Mean TR = ' num2str(meanTR) ' ms'])
    
    
%--------------------------------------------------------------------------
% Assign a bin to each line (X)
%--------------------------------------------------------------------------
    NLin    = length(param.PMUTimeStamp_ms);
    binMask = cell(nbins,1);
    binCnt  = zeros(nbins,1);
    
    binMaskMatrix = zeros([NLin,nbins]);
    
    trigCnt = 1;
    blockCounter = 0;
    
    for k = 1:NLin
        
        if trigCnt < nTrig
            timeTrig        = triggerPeaks_ms( trigCnt );
            timeTrigPlus1   = triggerPeaks_ms( trigCnt+1 );
            timeElem        = TimeStamp_ms( k );
        end
            
        if timeElem >= timeTrigPlus1
            trigCnt = trigCnt + 1;
            blockCounter = 1;
        else
            blockCounter = blockCounter + 1;
        end
        
        if trigCnt <= nPres
            idx = trigCnt;
            idx_x = indexVec_x( idx );
            if idx_x ~= 0
                if mod(k,22) ~=1
                    if blockCounter > 22
                        binMaskMatrix( k, idx_x ) = 1;
                    end
                end
            end
        end
        
    end
    
    for k = 1:nbins
        binMask{k} = binMaskMatrix(:,k);
    end
    
    param.binMask = binMask;
    
    

%--------------------------------------------------------------------------
% Saving data
%--------------------------------------------------------------------------    
  % Save data
    save(fullfile(param.savedir,'EyeBinning_X.mat'),'param');
  
    
    
    
%--------------------------------------------------------------------------
% Assign a bin to each line (Y)
%--------------------------------------------------------------------------
    NLin    = length(param.PMUTimeStamp_ms);
    binMask = cell(nbins,1);
    binCnt  = zeros(nbins,1);
    
    binMaskMatrix = zeros([NLin,nbins]);
    
    trigCnt = 1;
    blockCounter = 0;
    
    for k = 1:NLin
        
        if trigCnt < nTrig
            timeTrig        = triggerPeaks_ms( trigCnt );
            timeTrigPlus1   = triggerPeaks_ms( trigCnt+1 );
            timeElem        = TimeStamp_ms( k );
        end
            
        if timeElem >= timeTrigPlus1
            trigCnt = trigCnt + 1;
            blockCounter = 1;
        else
            blockCounter = blockCounter + 1;
        end
        
        if trigCnt <= nPres
            idx = trigCnt;
            idx_y = indexVec_y( idx );
            if idx_y ~= 0
                if mod(k,22) ~=1  
                    if blockCounter > 22
                        binMaskMatrix( k, idx_y ) = 1;
                    end
                end
            end
        end
        
    end
    
    for k = 1:nbins
        binMask{k} = binMaskMatrix(:,k);
    end
    
    param.binMask = binMask;
    
    

%--------------------------------------------------------------------------
% Saving data
%--------------------------------------------------------------------------    
  % Save data
    save(fullfile(param.savedir,'EyeBinning_Y.mat'),'param');
    
    
    
    