function [Jdet, x0, y0, x, y] = mesohyperbolicity(vx, vy, xlim, ylim, T, ngrid)
%MESOHYPERBOLICITY Calculates the mesohyperbolicity of a 2D velocity field
%   Returns a grid of the values of det(grad(V*))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Patrick Clary <pclary@umail.ucsb.edu>
% 5/29/2014
% Updated 6/30/2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = linspace(xlim(1), xlim(2), ngrid);
y = linspace(ylim(1), ylim(2), ngrid);
[x, y] = meshgrid(x, y);
x0 = x;
y0 = y;

% Take ten steps forward in time to find the mesochronic velocity field;
% Can be adjusted to take more steps if precision is needed over long times
% in a complex velocity field
nsteps = 10;
tvals = linspace(0, T, nsteps+1);
dt = T/nsteps;

for t = tvals
    [x, y] = rk4step(t, dt, x, y, vx, vy);
end

u = (x - x0)/T;
v = (y - y0)/T;

% Compute the gradient of the velocity field using convolution
dx = (xlim(2) - xlim(1))/(ngrid-1);
dy = (xlim(2) - xlim(1))/(ngrid-1);
xker = [0, 0, 0; 1, 0, -1; 0, 0, 0];
yker = xker';
dudx = conv2(u, xker, 'same')/dx;
dudy = conv2(u, yker, 'same')/dy;
dvdx = conv2(v, xker, 'same')/dx;
dvdy = conv2(v, yker, 'same')/dy;

% Mesohyperbolicity is given by the jacobian determinant of the mesochronic
% velocity field
Jdet = dudx.*dvdy - dudy.*dvdx;
    
end

