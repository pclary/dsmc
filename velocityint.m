function v = velocityint(t, x, y, vdat, xdat, ydat, tdat)
%VELOCITYINT Interpolates velocity data in x/y and time

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Patrick Clary <pclary@umail.ucsb.edu>
% 6/14/2014
% Updated 6/30/2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Interpolate two time slices to get the velocity grid for an intermediate
% time, then use the builtin interp2 to get the velocities at (x, y) points
i = find(tdat > t, 1) - 1;
v1 = interp2(xdat, ydat, vdat(:, :, i), x, y);
v2 = interp2(xdat, ydat, vdat(:, :, i+1), x, y);
v = v1 + (v2 - v1)*(t - tdat(i))/(tdat(i+1) - tdat(i));

% Landmasses are represented by NaNs in the velocity data
v(isnan(v)) = 0;

end
