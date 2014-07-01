function map = mhcolormap(n)
%MHCOLORMAP Creates a mesohyperbolicity colormap

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Patrick Clary <pclary@umail.ucsb.edu>
% 5/29/2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

map = [0, 0, 1; 1, 1, 1; 0, 1, 0; 1, 1, 1; 1, 0, 0];
l1 = linspace(0, 1, 5);
l2 = linspace(0, 1, n);
map = [interp1(l1, map(:, 1),l2)', ...
    interp1(l1, map(:, 2),l2)', ...
    interp1(l1, map(:, 3),l2)'];

end

