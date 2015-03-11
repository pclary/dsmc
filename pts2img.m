function bins = pts2img(x, y, res, xlim, ylim)
%PTS2IMG Converts a set of points to a normalized 2D image matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Patrick Clary <pclary@umail.ucsb.edu>
% 3/10/15
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Sort points into bins
xoff = diff(xlim)/res/2;
yoff = diff(ylim)/res/2;
xbin = linspace(xlim(1) + xoff, xlim(2) - xoff, res);
ybin = linspace(ylim(1) + yoff, ylim(2) - yoff, res);
bins = hist3([x(:), y(:)], {xbin, ybin})';

% Normalize
s = sum(sum(bins));
if s == 0
    bins = bins + 1;
    s = sum(bins(:));
end
bins = bins / (s/res^2);