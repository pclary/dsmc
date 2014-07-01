function v = velocityint(t, x, y, vdat, xdat, ydat, tdat)
%VELOCITYINT Interpolates velocity data in x/y and time

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Patrick Clary <pclary@umail.ucsb.edu>
% 6/14/2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i = find(tdat > t, 1) - 1;
v1 = interp2(xdat, ydat, vdat(:, :, i), x, y);
v2 = interp2(xdat, ydat, vdat(:, :, i+1), x, y);
v = v1 + (v2 - v1)*(t - tdat(i))/(tdat(i+1) - tdat(i));
v(isnan(v)) = 0;

end
