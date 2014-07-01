function [xt1n, xt2n] = randomwalkstep(xt1, xt2, h, umax, xlim, ylim, spherical)
%RANDOMWALKSTEP Computes the new position of each agent using a random walk

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Patrick Clary <pclary@umail.ucsb.edu>
% 5/26/2014
% Updated 6/14/2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = size(xt1, 2);
xt1n0 = xt1(end, :);
xt2n0 = xt2(end, :);

dir = 2*pi*rand(1, N);
if spherical
    adj = sqrt(sin(dir).^2+(cos(dir)./cosd(xt2)).^2);
else
    adj = 1;
end
xt1n = xt1n0 + h*umax*sin(dir)*adj;
xt2n = xt2n0 + h*umax*cos(dir)*adj;

outofbounds = find(xt1n0 > xlim(2) | xt1n0 < xlim(1) | xt2n0 > ylim(2) | xt2n0 < ylim(1));
for i = outofbounds
    d = sqrt((xlim(1)+(xlim(2)-xlim(1))/2-xt1n0(i)).^2 + (ylim(1)+(ylim(2)-ylim(1))/2-xt2n0(i)).^2);
    xt1n(i) = xt1n0(i) + h*umax*(xlim(1)+(xlim(2)-xlim(1))/2-xt1n0(i))/d;
    xt2n(i) = xt2n0(i) + h*umax*(ylim(1)+(ylim(2)-ylim(1))/2-xt2n0(i))/d;
end

end