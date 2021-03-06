function [xtn, ytn] = dsmcstep(xtn0, ytn0, dt, umax, s, cres, ...
    xlim, ylim, au, spherical)
%DSMCSTEP Computes the new position of each agent using the DSMC algorithm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Patrick Clary <pclary@umail.ucsb.edu>
% 5/18/2014
% Updated 3/10/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = size(xtn0, 2);
Bdir = @(Bs) bsxfun(@rdivide, Bs, sqrt(sum((Bs).^2, 1)));

% Adjust for spherical coordinates
if ~isempty(spherical) && spherical
    usn = @(dir) -umax*bsxfun(@times, dir, sqrt(dir(2,:).^2 + ...
        bsxfun(@rdivide, dir(1,:), cosd(ytn0)).^2));
else
    usn = @(dir) -umax*dir;
end

% Take an rk4 step
Bs = B(xtn0, ytn0, s, cres, xlim, ylim);
usn0 = usn(Bdir(Bs));
xtn1 = xtn0 + dt/2 * usn0(1, :);
ytn1 = ytn0 + dt/2 * usn0(2, :);
Bs = B(xtn1, ytn1, s, cres, xlim, ylim);
usn1 = usn(Bdir(Bs));
xtn2 = xtn0 + dt/2 * usn1(1, :);
ytn2 = ytn0 + dt/2 * usn1(2, :);
Bs = B(xtn2, ytn2, s, cres, xlim, ylim);
usn2 = usn(Bdir(Bs));
xtn3 = xtn0 + dt/2 * usn2(1, :);
ytn3 = ytn0 + dt/2 * usn2(2, :);
Bs = B(xtn3, ytn3, s, cres, xlim, ylim);
usn3 = usn(Bdir(Bs));
xtn = xtn0 + dt/6 * (usn0(1, :) + 2*usn1(1, :) + 2*usn2(1, :) + usn3(1, :));
ytn = ytn0 + dt/6 * (usn0(2, :) + 2*usn1(2, :) + 2*usn2(2, :) + usn3(2, :));

% Add disturbance (au: agent uncertainty)
if (~isempty(au))
    xtn = xtn + au*(xtn-xtn0).*(rand(1, N)-0.5)*2;
    ytn = ytn + au*(ytn-ytn0).*(rand(1, N)-0.5)*2;
end

% If outside of the domain, move towards the center of the domain
if ~isempty(xlim) && ~isempty(ylim)
    outofbounds = find(xtn0 > xlim(2) | xtn0 < xlim(1) | ytn0 > ylim(2) | ytn0 < ylim(1));
    for i = outofbounds
        d = sqrt((xlim(1)+(xlim(2)-xlim(1))/2-xtn0(i)).^2 + (ylim(1)+(ylim(2)-ylim(1))/2-ytn0(i)).^2);
        xtn(i) = xtn0(i) + dt*umax*(xlim(1)+(xlim(2)-xlim(1))/2-xtn0(i))/d;
        ytn(i) = ytn0(i) + dt*umax*(ylim(1)+(ylim(2)-ylim(1))/2-ytn0(i))/d;
    end
end
