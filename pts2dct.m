function coeffs = pts2dct(x, y, res, xlim, ylim)
%PTS2DCT Computes the DCT coefficients for the given point data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Patrick Clary <pclary@umail.ucsb.edu>
% 12/17/2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Sort points into bins
xbin = linspace(xlim(1), xlim(2), res);
ybin = linspace(ylim(1), ylim(2), res);
bins = hist3([x(:), y(:)], {xbin, ybin})';

% Normalize
bins = bins / (trapz(trapz(bins))/res^2);
coeffs = dct2(bins);
