clc; clear; close all;

subject_list = 1:15;

for subject_num = subject_list
    for mask_idx = 1:4
        for region_idx = 0:3
            % Rebuild runner constants because the worker script clears caller vars.
            baseDir = '/home/debi/jaime/repos/MR-EyeTrack/data/study';
            subjectStr = sprintf('sub-%03d', subject_num);
            reconDir = fullfile(baseDir, subjectStr, 'recon');
            mask_list = {'clean', 'clean_0.50', 'clean_0.75', 'clean_0.95'};
            nIter = 20;
            delta = 1.000;
            mask_type = mask_list{mask_idx};

            x0Path = fullfile(reconDir, mask_type, 'x0', sprintf('x0_regionidx%i.mat', region_idx));
            xPath = fullfile(reconDir, mask_type, 'x', ...
                sprintf('x_steva_regionidx_%i_nIter_%d_delta_%.3f.mat', region_idx, nIter, delta));

            if exist(x0Path, 'file') && exist(xPath, 'file')
                fprintf('Skipping finished task: subject=%d mask=%s region=%d\n', ...
                    subject_num, mask_type, region_idx);
                continue;
            end

            fprintf('\n=== Running task: subject=%d mask=%s region=%d ===\n', ...
                subject_num, mask_type, region_idx);
            close all force;
            run('/home/debi/jaime/repos/MR-EyeTrack/recon/4-Recon/S4_recon_debi_4fr_pulseq.m');
        end
    end
end
