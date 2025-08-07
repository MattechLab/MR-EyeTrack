function cMask = eyeGenerateBinningWin(datasetDir, nShotOff, nSeg,th_ratio, ETDir, winLen, display_binning)
    % winLen: the number of MR readouts
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

%-------------------------------------------------------------------------- 
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
        % Prompt user to select generated mask file from ET data (.mat)
    
        [maskDataName, maskDataDir, ~] = uigetfile( ...
        { '*.mat','Generated mask data file (*.mat)'}, ...
           'Pick a eye tracking mask file', ...
           'MultiSelect', 'off', ETDir);

        if maskDataName == 0
            warning('No mask file selected');
            return;
        end
        filepathMaskData = fullfile(maskDataDir, maskDataName);
        disp(filepathMaskData);

        mask_method_1 = load(filepathMaskData);
        if isstruct(mask_method_1)
            fieldNames = fieldnames(mask_method_1);
            mask_method_1 = mask_method_1.(fieldNames{1});
        end
        mask_method_1 = padArrayWithZeros(mask_method_1, round(Timediff));
        disp(sum(mask_method_1));

        % Set the window width and threshold
        WinWidth = round(length(mask_method_1)/NLin)*winLen;
        HalfWinWidth = floor(WinWidth/2);
        disp('WinWidth');
        disp(WinWidth);

        nMeasuresOff = nShotOff*nSeg;
        for k = 1:NLin
            %debug flag
            % Eliminate the first nShotOff
            if k<= nMeasuresOff
                 binMaskMatrix( k, idx_bin ) = 0;
            else
                % the timestamp of the current mri meas 
                timeSeg = TimeStamp_ms(k);
                    
                win_lower = round(max(0, timeSeg - HalfWinWidth));
                win_upper = round(min(numel(mask_method_1),timeSeg + HalfWinWidth));
                       
                window_data = mask_method_1(win_lower:win_upper);
           
                true_count = sum(window_data);
              
                th = win_upper-win_lower;
            
                if true_count >= th*th_ratio
                    binMaskMatrix( k, idx_bin) = 1;
                end
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

    if display_binning
        figure; hold on;
        colors = lines(nbins);
        for i = 1:nbins
            plot(TimeStamp_s, cMask(i, :)+i, 'Color',colors(i, :));
        end
        xlabel('Times(seconds)');
        ylabel('Mask idx');      
        title('Binning mask');
        hold off;
    end
end









