function coeffs = pts2dct(xt, yt, res, xlim, ylim)
%PTS2DCT Computes the DCT coefficients for the given point data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Patrick Clary <pclary@umail.ucsb.edu>
% 12/17/2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Sort points into bins
x = linspace(xlim(1), xlim(2), res);
y = linspace(ylim(1), ylim(2), res);
bins = hist3([xt(:), yt(:)], {x, y});

% Normalize
bins = bins / trapz(trapz(bins)) * diff(xlim)*diff(ylim)/res^2;
coeffs = dct2(bins);
