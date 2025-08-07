function cMask = eyeGenerateWoBinning(datasetDir, nShotOff, nSeg)
    param.basedir = datasetDir;
    basedir = param.basedir;
    param.batchParam = [];

    % Prompt user to select Siemens Raw Data file
    %--------------------------------------------------------------------------

    [rawDataName, rawDataDir, ~] = uigetfile( ...
    { '*.dat','Siemens raw data file (*.dat)'}, ...
       'Pick a raw data file', ...
       'MultiSelect', 'off', datasetDir);

    if rawDataName == 0
        warning('No file selected');
        return;
    end

% How to select the number of bins?
    prompt        = {'Enter the number of bins'};
    name          = '#Bins';
    numlines      = 1;
    defaultanswer = {'1'};
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

    costTime = 2.5;
    param.batchParam.rawDataName    = rawDataName;
    param.batchParam.rawDataDir     = rawDataDir;

    [ twix_obj, param ] = dataSelectionAndLoading( basedir, param );
    % Load the PMUTime and TimeStamp, shift the TimeStamp to the beginning
    % of 0
    PMUTimeStamp    = double( twix_obj.image.pmutime );
    TimeStamp       = double( twix_obj.image.timestamp );

    %
    TimeStamp       = TimeStamp - min(TimeStamp);
    % Do not forget to scale the times by costTime (Setting from Siemens)
    PMUTimeStamp_ms = PMUTimeStamp * costTime;
    PMUTimeStamp_s  = PMUTimeStamp_ms / 1000;
    TimeStamp_ms    = TimeStamp * costTime;
    TimeStamp_s     = TimeStamp_ms / 1000;
    % Save all the param to struct "param"
    param.PMUTimeStamp_ms   = PMUTimeStamp_ms;
    param.PMUTimeStamp_s    = PMUTimeStamp_s;
    param.TimeStamp_ms      = TimeStamp_ms;
    param.TimeStamp_s       = TimeStamp_s;

    inspect = false;
    if inspect
        figure;
        plot(param.TimeStamp_ms, param.PMUTimeStamp_ms);
        % xlim([0 1e4]);  % Set x-axis limits
        % ylim([0 180]); % Set y-axis limits
        xlabel('TimeStamp (ms)');
        ylabel('PMU TimeStamp (ms)');
        title('Time vs. PMU Stamp');
    end
    Timediff = TimeStamp_ms(end) - TimeStamp_ms(1);
    disp(['The duration of the rawdata is: ', num2str(Timediff), ...
    ' ms with data points:', num2str(length(TimeStamp_s)) ]);

    NLin    = length(param.PMUTimeStamp_ms);
    binMask = cell(nbins,1);
    binMaskMatrix = zeros([NLin,nbins]);

    for idx_bin = 1:nbins

        nMeasuresOff = nShotOff*nSeg;
        for k = 1:NLin
            %debug flag
            % Eliminate the first nShotOff
            if k<= nMeasuresOff
                 binMaskMatrix( k, idx_bin ) = 0;
            else
                binMaskMatrix( k, idx_bin) = 1;
            end

            if mod(k, 22)  == 1
                binMaskMatrix( k, idx_bin) = 0;
            end

        end
        binMask{idx_bin} = binMaskMatrix(:, idx_bin);
        sum_binning = sum(binMaskMatrix(:, idx_bin));
     
        disp(['with Binning, preserved #line: ',num2str(sum_binning), ' out of ', num2str(length(TimeStamp_s)) ])
            
    end
    param.binMask = binMask;
    cMask = logical(binMaskMatrix)';

end