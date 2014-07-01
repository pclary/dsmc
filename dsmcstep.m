function [xt1n, xt2n] = dsmcstep(xt1, xt2, muks, h, umax, hks, K, xlim, ...
    ylim, au, spherical)
%DSMCSTEP Computes the new position of each agent using the DSMC algorithm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Patrick Clary <pclary@umail.ucsb.edu>
% 5/18/2014
% Updated 6/30/2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = size(xt1, 2);
xt1n0 = xt1(end, :);
xt2n0 = xt2(end, :);
Bdir = @(Bs) bsxfun(@rdivide, Bs, sqrt(sum((Bs).^2, 1)));

% Adjust for spherical coordinates
if spherical
    usn = @(dir) -umax*bsxfun(@times, dir, sqrt(dir(2,:).^2 + ...
        bsxfun(@rdivide, dir(1,:), cosd(xt2(end,:))).^2));
else
    usn = @(dir) -umax*dir;
end

% Take an rk4 step
Bs = B(K, xt1, xt2, hks, muks, xlim, ylim);
usn0 = usn(Bdir(Bs));
xt1n1 = xt1n0 + h/2 * usn0(1, :);
xt2n1 = xt2n0 + h/2 * usn0(2, :);
Bs = B(K, [xt1; xt1n1], [xt2; xt2n1], hks, muks, xlim, ylim);
usn1 = usn(Bdir(Bs));
xt1n2 = xt1n0 + h/2 * usn1(1, :);
xt2n2 = xt2n0 + h/2 * usn1(2, :);
Bs = B(K, [xt1; xt1n2], [xt2; xt2n2], hks, muks, xlim, ylim);
usn2 = usn(Bdir(Bs));
xt1n3 = xt1n0 + h/2 * usn2(1, :);
xt2n3 = xt2n0 + h/2 * usn2(2, :);
Bs = B(K, [xt1; xt1n3], [xt2; xt2n3], hks, muks, xlim, ylim);
usn3 = usn(Bdir(Bs));
xt1n = xt1n0 + h/6 * ...
    ( usn0(1, :) + ...
    2 * usn1(1, :)  + ...
    2 * usn2(1, :)  + ...
    usn3(1, :)  );
xt2n = xt2n0 + h/6 * ...
    ( usn0(2, :) + ...
    2 * usn1(2, :)  + ...
    2 * usn2(2, :)  + ...
    usn3(2, :)  );

% Add disturbance (au: agent uncertainty)
xt1n = xt1n + au*(xt1n-xt1n0).*(rand(1, N)-0.5)*2;
xt2n = xt2n + au*(xt2n-xt2n0).*(rand(1, N)-0.5)*2;

% If outside of the domain, move towards the center of the domain
outofbounds = find(xt1n0 > xlim(2) | xt1n0 < xlim(1) | xt2n0 > ylim(2) | xt2n0 < ylim(1));
for i = outofbounds
    d = sqrt((xlim(1)+(xlim(2)-xlim(1))/2-xt1n0(i)).^2 + (ylim(1)+(ylim(2)-ylim(1))/2-xt2n0(i)).^2);
    xt1n(i) = xt1n0(i) + h*umax*(xlim(1)+(xlim(2)-xlim(1))/2-xt1n0(i))/d;
    xt2n(i) = xt2n0(i) + h*umax*(ylim(1)+(ylim(2)-ylim(1))/2-xt2n0(i))/d;
end

end