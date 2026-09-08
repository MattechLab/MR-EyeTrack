%% ROI-PCA - maximum-compression variant of the ROVir pipeline
%
% ROI-PCA is NOT a separate method: it is the same generalised eigenproblem
% in recon/ROVir/rovir_solve.m with the Tikhonov load on B taken to infinity.
% As lambda grows, B tends to a scaled identity and
%
%       A v = lambda_gen (B + a*I) v      ->      A v = mu v
%
% i.e. the plain eigenvectors of A, ranked by ROI signal energy. Equivalently,
% the right singular vectors of X_roi (verified: subspace difference 1e-5).
%
% Why no duplicated code lives here: the two variants differ by one argument.
% Forking R1-R8 would create ~1500 lines that drift apart the first time either
% side is fixed. Instead every R script takes a `variant` string and writes to
% recon/<variant>/, so the OUTPUTS are fully separated while the code is not:
%
%   data/study/sub-NNN/recon/ROVir/      SIR-ranked   (suppress the head)
%   data/study/sub-NNN/recon/ROI-PCA/    energy-ranked (maximum compression)
%   data/study/sub-NNN/recon/mitosius/ROI-PCA_<nv>/
%
% Trade-off measured on sub-015 (ROI % / head % energy retained):
%
%        n     ROVir          ROI-PCA
%        4     5.4 / 0.09     96.6 / 44.9
%        8    10.5 / 0.21     99.4 / 66.4
%
% ROI-PCA keeps far more eye signal per channel; it does not deface. That is
% the right trade when the goal is shrinking k-space for a model rather than
% removing identifiable anatomy.

clc; clearvars; close all;

subject_num = 15;
variant     = 'ROI-PCA';
lambda      = 1e6;      % -> pure ROI-PCA (1e4 already converges)
nvManual    = 8;        % build at 8; derive 4 and 2 from it with R7
binName     = 'woBin';

run(fullfile(fileparts(mfilename('fullpath')), '..', 'ROVir', 'R1_rovir_masks.m'));
run(fullfile(fileparts(mfilename('fullpath')), '..', 'ROVir', 'R2_rovir_transform.m'));

fprintf(['\n\nNext:\n' ...
         '  variant=''ROI-PCA''; run recon/ROVir/R3_rovir_mitosius.m   (reads raw, ~40 min)\n' ...
         '  variant=''ROI-PCA''; srcNv=8; deriveTo=[4 2]; run recon/ROVir/R7_nv_table.m\n' ...
         '  bash recon/ROVir/hpc/push_rovir.sh -s 15 -n 4 -v ROI-PCA --submit\n']);
