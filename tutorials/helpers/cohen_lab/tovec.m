% --------------------------------------------------------------------------
% Copied into this repository on 2026-09-27; the code below is unchanged.
% Origin: Cohen Lab, Harvard University (2020RigControl repository,
% Analysis/Old_Analysis/Image Processing). Written by Vicente Parot.
% --------------------------------------------------------------------------
function out = tovec(mov);

% function out = tovec(mov);
% converts a movie of size [nr, nc, nt] to a matrix of size [nr*nc, nt];
%
% modified from code by Vicente Parot

[nr, nc, ~] = size(mov);
out = reshape(mov,nr*nc,[]);