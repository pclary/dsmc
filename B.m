function out = B(xa, ya, s, cres, xlim, ylim)
%B Computes the 'B' vector for each agent 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Patrick Clary <pclary@umail.ucsb.edu>
% 5/18/2014
% Updated 3/10/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


dx = diff(xlim)/cres;
dy = diff(ylim)/cres;

% Pad the outside of the DSMC surface to make the derivative zero at the 
% boundaries
s = [s(1, 1),   s(1, :),   s(1, end);
     s(:, 1),   s,         s(:, end);
     s(end, 1), s(end, :), s(end, end)];
cres = cres + 2;
xlim = [xlim(1) - dx, xlim(2) + dx];
ylim = [ylim(1) - dy, ylim(2) + dy];

% Numerical gradient
[dsdx, dsdy] = gradient(s, dx, dy);

% Get progerly spaced and centered axis vectors
x = linspace(xlim(1) + dx/2, xlim(2) - dx/2, cres);
y = linspace(ylim(1) + dy/2, ylim(2) - dy/2, cres);

% Interpolate numerical gradient to get B
out = [interp2(x, y, dsdx, xa, ya, 'linear', 0); 
       interp2(x, y, dsdy, xa, ya, 'linear', 0)];
