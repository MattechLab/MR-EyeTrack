function weights_norm = coilSelectionEyesROI_intensity(reconDir, x0)
% Rank coils by integrated image intensity in the eye ROI.
% x0 : cell array (nCh x 1) of per-coil gridded reconstructions (complex),
%       OR a char/string path to a saved .mat file containing that cell array.

    if ischar(x0) || isstring(x0)
        tmp = load(x0, 'x0');
        x0 = tmp.x0;
    end

    numCoils = numel(x0);
    [nx, ny, nz] = size(x0{1});

    % RMS combination for ROI drawing reference
    sum_sq = zeros(nx, ny, nz, 'single');
    for c = 1:numCoils
        sum_sq = sum_sq + real(x0{c} .* conj(x0{c}));
    end
    xrms = sqrt(sum_sq / numCoils);

    % Draw eye ROI interactively on central slice of RMS image
    figure; imagesc(abs(xrms(:,:,round(nz/2)))); axis image; colormap gray;
    title('Draw mask around the eyes');
    eyeMask = roipoly;

    outDir = fullfile(reconDir, 'mitosius/woBin_comp');
    if ~exist(outDir, 'dir'); mkdir(outDir); end

    eyeMaskPath = fullfile(outDir, 'eyeMask_intensity.mat');
    save(eyeMaskPath, 'eyeMask');

    bmImage(eyeMask);
    bmImage(xrms .* eyeMask);

    % Score each coil: sum |x0| inside eye ROI over z-range
    zRange = 12:36;
    eyeMask3D = repmat(eyeMask, [1 1 length(zRange)]);

    rawWeights = zeros(numCoils, 1);
    for coil = 1:numCoils
        coilAbs = abs(x0{coil}(:,:,zRange));
        rawWeights(coil) = sum(coilAbs(eyeMask3D == 1), 'all');
    end

    % Normalize and sort
    weights_norm = rawWeights / sum(rawWeights);
    [sortedW, idx] = sort(weights_norm, 'descend');

    save(fullfile(outDir, 'weights_norm_intensity.mat'), 'weights_norm');
    save(fullfile(outDir, 'idx_coilSelection_intensity.mat'), 'idx', '-v7.3');
    disp(['Coil selection indices saved: ', fullfile(outDir, 'idx_coilSelection_intensity.mat')]);

    disp('Coils sorted by weight (highest to lowest):');
    disp(idx(:)');

    % Plot all coils sorted by weight
    figure;
    bar(idx, sortedW, 'LineWidth', 1.2);
    xlabel('Sorted coil index'); ylabel('Intensity in eye region');
    title('Sorted Coil Contributions to Eye Region (Intensity-based)');
    grid on;

    % Plot top 10 highlighted
    numTop = 10;
    figure;
    bar(weights_norm, 'FaceColor', [0.6 0.6 0.6]); hold on;
    bar(idx(1:numTop), weights_norm(idx(1:numTop)), 'FaceColor', [0.9 0.3 0.3]);
    xlabel('Coil index'); ylabel('Intensity in eye region');
    title(['Top ' num2str(numTop) ' Coils (Intensity): ', strjoin(cellstr(num2str(idx(1:numTop))), ' ')]);
    legend('All coils', 'Top contributors');
    grid on;

    % Display individual images for top coils
    % for k = 1:numTop
    %     bmImage(x0{idx(k)});
    % end

end
