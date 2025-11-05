%% Evaluate K-space fullness
clearvars; restoredefaultpath; close all;

%% PARAMETERS
addpath('/home/debi/jaime/repos/MR-EyeTrack/recon/TrajAndRCallculation/utils')
import_eMask=1;
% list of possible eMask paths
eMask_paths = {
    '/home/debi/jaime/repos/MR-EyeTrack/results/Sub001/T1_LIBRE_woBinning/other/eMask_woBin.mat'    % woBin
    '/home/debi/jaime/repos/MR-EyeTrack/results/Sub001/T1_LIBRE_Binning/other/clean/eMask_th0.75_region0.mat'    % clean
    '/home/debi/jaime/repos/MR-EyeTrack/results/Sub001/T1_LIBRE_Binning/other/clean_0.5/eMask_th0.75_region0.mat' % clean 50%
    '/home/debi/jaime/repos/MR-EyeTrack/results/Sub001/T1_LIBRE_Binning/other/clean_0.75/eMask_th0.75_region0.mat' % clean 75%
    '/home/debi/jaime/repos/MR-EyeTrack/results/Sub001/T1_LIBRE_Binning/other/clean_0.95/eMask_th0.75_region0.mat' % clean 95%
};
eMask_path = eMask_paths{5}; % Select desired eMask path
fprintf('eMask: %s\n', eMask_path);
nx  = 480;          % Number of sampled points along each spoke
nPC = 1;            % Number of phase cycles
Shots = 3723;        % Total number of spokes acquired
Segment = 22;       % Number of segments per shot
FM = 0;             % Flag to determine the type of trajectory computation (Frequency modulation)
FullRange = 1;      % Trajectory type
AngleUpdate = 1;    % Flag for using angle update in phyllotaxis computation (Multiple PC)

%% CALCULATE TRAJECTORY
if FM == 0
    if AngleUpdate == 1
        flagSelfNav = 1;
        flagPlot= 0;
        flagRosettaTraj=0;
     
        % Compute phyllotaxis trajectory with angle update
        [kx, ky, kz, azimuthal, polar, R] = computePhyllotaxis_angleUpdate(nx, Segment, Shots, flagSelfNav, flagPlot, flagRosettaTraj, FullRange, nPC);
        
        Traj3D = cat(5, kx, ky, kz);
        Traj3D = reshape(Traj3D, [nx, Segment * Shots * nPC, 3]);
    else
        % Compute standard phyllotaxis trajectory without angle update
        [kx, ky, kz, azimuthal, polar, R] = computePhyllotaxis(nx, Segment, Shots, 1, 0, 0, 0, FullRange);
        Traj3D = cat(4, kx, ky, kz);
        Traj3D = reshape(Traj3D, [nx, Segment * Shots, 3]);
        Traj3D = repmat(Traj3D, [1, nPC, 1, 1]);
    end
else
    % Compute phyllotaxis trajectory without angle update (original method)
    [kx, ky, kz, azimuthal, polar, R] = computePhyllotaxis(nx, Segment, Shots, 1, 0, 0, 0, FullRange);
    Traj3D = cat(5, kx, ky, kz);
    Traj3D = reshape(Traj3D, [nx, Segment * Shots, 3]);
end

% Plot trajectory (first line of each segment)
figure; plot3(Traj3D(:,1:Segment,1), Traj3D(:,1:Segment,2), Traj3D(:,1:Segment,3), 'LineWidth', 2);
hold on;
plot3(Traj3D(end,1:Segment,1), Traj3D(end,1:Segment,2), Traj3D(end,1:Segment,3), 'LineWidth', 2, 'Color', 'r');
grid minor;
xlabel('x-axis'); ylabel('y-axis'); zlabel('z-axis');
title('Trajectory');
xlim([-0.5 0.5]); ylim([-0.5 0.5]); zlim([-0.5 0.5]);
view(3);

%% GRIDDING (Mapping trajectory to Cartesian grid)

[nx, ntviews, ~] = size(Traj3D);
x_cord = linspace(-0.5, 0.5, nx);
dx = x_cord(2) - x_cord(1); % Grid spacing

matrix = zeros(nx, nx, nx); % Initialize empty 3D grid

% Convert trajectory coordinates to discrete grid indices
indices = round((Traj3D + 0.5) / dx) + 1;
% Example data
% import eMask
if import_eMask
    % try to load eMask or eMaskN from the .mat file
    vars = whos('-file', eMask_path);
    varnames = {vars.name};
    if ismember('eMask', varnames)
        tmp = load(eMask_path, 'eMask');
        eMask = tmp.eMask;
    elseif ismember('eMaskN', varnames)
        tmp = load(eMask_path, 'eMaskN');
        eMask = tmp.eMaskN;
    else
        warning('No variable "eMask" or "eMaskN" found in %s. Using all-ones mask.', eMask_path);
        eMask = ones(Shots * Segment, 1);
    end

    % ensure logical column vector of expected length; fall back to all-ones on mismatch
    eMask = logical(eMask(:));
    if numel(eMask) ~= Shots * Segment
        warning('eMask length (%d) does not match Shots*Segment (%d). Using all-ones mask.', numel(eMask), Shots * Segment);
        eMask = true(Shots * Segment, 1);
    end

    % Expand the mask to apply across dimensions
    eMask = reshape(eMask, [1, Shots*Segment, 1]);
    eMask = repmat(eMask, [nx, 1, 3]);
    % Apply the mask
    indices(~eMask) = 0;  % sets values to 0 where mask is 0
end

% Remove points that fall outside grid boundaries
valid = all(indices >= 1 & indices <= nx, 3);

% Populate grid with sampled k-space points
for k = 1:ntviews
    if valid(:, k) % Only consider valid points
        x_index = indices(:, k, 1);
        y_index = indices(:, k, 2);
        z_index = indices(:, k, 3);
        linear_idx = sub2ind([nx, nx, nx], x_index, y_index, z_index);
        matrix(linear_idx) = matrix(linear_idx) + 1;
    end
end
matrix(matrix > 1) = 1; % Normalize grid occupancy to binary values

%% SPHERICAL MASK (For k-space volume estimation)
[x, y, z] = ndgrid(x_cord, x_cord, x_cord); % Create 3D grid coordinates
mask = (x.^2 + y.^2 + z.^2) <= 0.5^2; % Define spherical mask

%% ACCELERATION FACTORS (Various definitions)

R1 = nx^3 / numel(find(matrix == 1)); % Cartesian (including empty corners)
R2 = (nx^3 - numel(find(mask == 0))) / numel(find(matrix == 1)); % Cartesian (excluding empty corners)
R3 = (0.2*size(matrix,1)^3)/numel(find(matrix == 1)); % 20% Nyquist with full matrix incl corners
R4 = (0.2*numel(find(mask == 1)))/numel(find(matrix == 1)); % 20% Nyquist with full matrix exl corners
R5 = (nx^2 * pi) / (size(Traj3D,2) * 8); % Radial approximation

% Visualize k-space fullness
figure
sliceViewer(matrix);
%fprintf('Cartesian (incl. corners): R=%.2f\n', R1);
fprintf('Cartesian (excl. corners): R=%.2f\n', R2);
fprintf('Matteo:  20%% Nyquist with full matrix incl corners: R=%.2f\n', R3);
fprintf('20%% Nyquist with full matrix excl corners: R=%.2f\n', R4);
fprintf('Simens:  equation for Radial: R=%.2f\n', R5);

%----------- output: wo eMask--------------------------
% Cartesian (excl. corners): R=3.87
% Matteo:  20% Nyquist with full matrix incl corners: R=1.49
% Simens:  equation for Radial: R=2.00
%----------- output: with eMask--------------------------
% Cartesian (excl. corners): R=7.27
% Matteo:  20% Nyquist with full matrix incl corners: R=2.79
% Simens:  equation for Radial: R=2.00