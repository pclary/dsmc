function [xt1n, xt2n] = randomwalkstep(xt1n0, xt2n0, dt, umax, xlim, ylim, spherical)
%RANDOMWALKSTEP Computes the new position of each agent using a random walk

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Patrick Clary <pclary@umail.ucsb.edu>
% 5/26/2014
% Updated 12/17/2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = size(xt1n0, 2);

% Pick a random direction
dir = 2*pi*rand(1, N);

% Adjust for spherical coordinates
if ~isempty(spherical) && spherical
    adj = sqrt(sin(dir).^2+(cos(dir)./cosd(xt2n0)).^2);
else
    adj = 1;
end

% Use an euler step to compute the new agent positions
xt1n = xt1n0 + dt.*umax.*sin(dir).*adj;
xt2n = xt2n0 + dt.*umax.*cos(dir).*adj;

% If outside of the domain, move towards the center of the domain
if ~isempty(xlim) && ~isempty(ylim)
    outofbounds = find(xt1n0 > xlim(2) | xt1n0 < xlim(1) | xt2n0 > ylim(2) | xt2n0 < ylim(1));
    for i = outofbounds
        d = sqrt((xlim(1)+(xlim(2)-xlim(1))/2-xt1n0(i)).^2 + ...
            (ylim(1)+(ylim(2)-ylim(1))/2-xt2n0(i)).^2);
        xt1n(i) = xt1n0(i) + dt*umax*(xlim(1)+(xlim(2)-xlim(1))/2-xt1n0(i))/d;
        xt2n(i) = xt2n0(i) + dt*umax*(ylim(1)+(ylim(2)-ylim(1))/2-xt2n0(i))/d;
    end
end

end