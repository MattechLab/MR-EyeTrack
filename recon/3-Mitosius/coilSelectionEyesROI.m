function weights_norm = coilSelectionEyesROI(reconDir, C)

    % Calculate the rms of C
    crms = cal_crms(C);
    bmImage(crms);

    % define a mask around eye region
    xrms_path = fullfile(reconDir, 'woBin/xrms48.mat');
    load(xrms_path, 'xrms');
    nSlice = size(xrms,3);
    figure; imagesc(abs(xrms(:,:,round(nSlice/2)))); axis image
    title('Draw mask around the eyes');
    eyeMask = roipoly;   % binary mask
    eyeMaskPath = fullfile(reconDir, 'mitosius/woBin_comp/eyeMask.mat');
    mkdir(fileparts(eyeMaskPath));
    save(eyeMaskPath, "eyeMask");

    bmImage(eyeMask);
    bmImage(xrms.*eyeMask)

    % define zRange for eye region for C image
    zRange = 12:36;
    numCoils = size(C, 4);
    eyeMask3D = repmat(eyeMask, [1 1 length(zRange)]);

    rawWeights = zeros(numCoils,1);

    for coil = 1:numCoils
        coilAbs = abs(C(:,:,zRange,coil));
        rawWeights(coil) = sum(coilAbs(eyeMask3D==1), 'all');
    end

    % normalize across coils
    weights_norm = rawWeights / sum(rawWeights);
    weights_norm_path = fullfile(reconDir, 'mitosius/woBin_comp/weights_norm.mat');
    save(weights_norm_path, "weights_norm");

    % sort coils by weights
    [sortedW, idx] = sort(weights_norm, 'descend');
    % plot
    figure;
    bar( idx, sortedW, 'LineWidth', 1.2);
    xlabel('Sorted coil index');
    ylabel('Weight in eye region');
    title('Sorted Coil Contributions to Eye Region');
    grid on;
    % Print coil order
    disp('Coils sorted by weight (highest to lowest):');
    disp(idx(:)');
    
    % top 10 coils
    numTop = 10;   % change if needed
    % plot
    figure;
    bar(weights_norm, 'FaceColor', [0.6 0.6 0.6]); hold on;
    bar(idx(1:numTop), weights_norm(idx(1:numTop)), 'FaceColor', [0.9 0.3 0.3]);
    xlabel('Coil index');
    ylabel('Weight in eye region');
    title(['Top ' num2str(numTop) ' Coil Contributions: ',  strjoin(cellstr(num2str(idx(1:numTop))), ' ')]);
    legend('All coils', 'Top contributors');
    grid on;

    for select_coil = idx(1:numTop)
        bmImage(C(:,:,:,select_coil));
    end

end
