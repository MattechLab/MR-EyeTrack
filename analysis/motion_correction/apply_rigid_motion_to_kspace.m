function [y_corr, t_corr, meta] = apply_rigid_motion_to_kspace( ...
    y, t, lineTimeMs, motionParams, motionTimeMs, opts)
%APPLY_RIGID_MOTION_TO_KSPACE
% Applies rigid-body motion correction to raw k-space data
%
% INPUTS
%   y              : [nRead x nCoils x nLines] raw data
%   t              : [3 x nRead x nLines] k-space trajectory
%   lineTimeMs     : [nLines x 1] timestamps (ms)
%   motionParams   : [N x 6] (Tx Ty Tz Rx Ry Rz)
%   motionTimeMs   : [N x 1] motion frame times (ms)
%   opts           : options struct
%
% OUTPUTS
%   y_corr         : motion-corrected k-space
%   t_corr         : motion-corrected trajectory
%   meta           : interpolated motion parameters

arguments
    y
    t
    lineTimeMs (:,1) double
    motionParams (:,6) double
    motionTimeMs (:,1) double
    opts.interpMethod (1,:) char = 'linear'
    opts.rotationSigns (1,3) double = [1 1 -1]
    opts.translationSigns (1,3) double = [-1 1 1]
end

nLines = numel(lineTimeMs);

%% -------------------------------------------------
%% Interpolate motion to k-space lines
%% -------------------------------------------------
trans = interp1( ...
    motionTimeMs, motionParams(:,1:3), ...
    lineTimeMs, opts.interpMethod, 'extrap');

rot = interp1( ...
    motionTimeMs, motionParams(:,4:6), ...
    lineTimeMs, opts.interpMethod, 'extrap');

meta.translation_per_line = trans;
meta.rotation_per_line    = rot;

%% -------------------------------------------------
%% Apply rigid motion
%% -------------------------------------------------
y_corr = y;
t_corr = t;

for l = 1:nLines

    rx = opts.rotationSigns(1)*rot(l,1);
    ry = opts.rotationSigns(2)*rot(l,2);
    rz = opts.rotationSigns(3)*rot(l,3);

    Rx = [1 0 0; 0 cos(rx) -sin(rx); 0 sin(rx) cos(rx)];
    Ry = [cos(ry) 0 sin(ry); 0 1 0; -sin(ry) 0 cos(ry)];
    Rz = [cos(rz) -sin(rz) 0; sin(rz) cos(rz) 0; 0 0 1];

    R = Rz * Ry * Rx;

    % We shoudl apply the reverse of the detected rotation to correct for
    % them
    R = R.';
    

    t_corr(:,:,l) = R * t(:,:,l);

    tx = opts.translationSigns(1) * trans(l,1);
    ty = opts.translationSigns(2) * trans(l,2);
    tz = opts.translationSigns(3) * trans(l,3);

    phase = exp(-2*pi*1i*( ...
        t_corr(1,:,l)*tx + ...
        t_corr(2,:,l)*ty + ...
        t_corr(3,:,l)*tz ));

    y_corr(:,:,l) = y_corr(:,:,l) .* phase;
end

end
