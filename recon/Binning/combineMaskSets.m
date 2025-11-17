function [outMask] = combineMaskSets(nMaskSet, ETDir)
%COMBINEMASKSETS Summary of this function goes here
%   Detailed explanation goes here
MaskSet = cell(nMaskSet,1);
for iSet = 1:nMaskSet
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
    
        MaskSet{iSet,1} = mask_method_1;
        sum_mask = sum(mask_method_1);
        disp(['Mask-', num2str(iSet), ': sum of 1 ', num2str(sum_mask)]);
end

intersection_mask = MaskSet{1};
for i = 2:length(MaskSet)
    intersection_mask = intersection_mask & MaskSet{i};
end
outMask = intersection_mask;

disp(['The intersection of masks: sum ', num2str(sum(outMask))])


end

