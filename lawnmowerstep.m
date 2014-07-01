function [xt1n, xt2n] = lawnmowerstep(xt1, xt2, h, umax, xlim, ylim, nlines, spherical)
%LAWNMOWERSTEP Computes the new position of each agent using a lawnmower
%algorithm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Patrick Clary <pclary@umail.ucsb.edu>
% 5/18/2014
% Updated 6/30/2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = size(xt1, 2);
xt1n = zeros(1, N);
xt2n = zeros(1, N);

for i = 1:N
    % Implements a "stateless" lawnmower trajectory; computes the new step
    % using only prior (non-advected) trajectory data

    % Find the number of times the agent has gone back and forth
    % horizontally
    lastlow = find(xt1(:, i) < xlim(1), 1, 'last');
    if numel(lastlow) == 0
        lastlow = 0;
    end
    lasthigh = find(xlim(2) - xt1(:, i) < 0, 1, 'last');
    if numel(lasthigh) == 0
        lasthigh = 0;
    end
    if lastlow < lasthigh
        dir1 = -1;
    else
        dir1 = 1;
    end
    tmp = xt1(:, i) < xlim(1);
    tmp = tmp - (xlim(2) - xt1(:, i) < 0);
    tmp = diff([0; tmp(tmp ~= 0)]);
    tripcount = sum(tmp ~= 0);
    
    % Compute a phase offset using the tripcount; prevents the agent from
    % covering the exact same path as before
    tripcount = tripcount*N + i - 1;
    n = floor(log(tripcount+1.1)/log(2));
    rem = tripcount+1-2^n;
    phase = (2*rem+1)/2^(n+1);
    
    p1 = xt1(end, i);
    p2 = xt2(end, i);
    
    % Look back in time to determine which line the agent was previously on
    ptsback = round((ylim(2)-ylim(1))/umax/2/h);
    line = max(min(round((xt1(end, i)-xlim(1))/(xlim(2)-xlim(1))*nlines-phase), nlines-1), 0);
    if size(xt1, 1) > ptsback
        lastline = max(min(round((xt1(end-ptsback, i)-xlim(1))/(xlim(2)-xlim(1))*nlines-phase), nlines-1), 0);
    else
        lastline = line;
    end
    linecenter = xlim(1) + (line+phase)*(xlim(2)-xlim(1))/nlines;
    
    % Determine which vertical direction the agent should be moving in
    lastlow = find(xt2(:, i) < ylim(1), 1, 'last');
    if numel(lastlow) == 0
        lastlow = 0;
    end
    lasthigh = find(ylim(2) - xt2(:, i) < 0, 1, 'last');
    if numel(lasthigh) == 0
        lasthigh = 0;
    end
    if lastlow < lasthigh
        dir2 = -1;
    else
        dir2 = 1;
    end
    
    % Use proportional feedback to keep the agent moving along the line
    % center even when disturbed by a velocity field
    if ((ylim(2) - p2 < 0) || (p2 < ylim(1) + 0)) && (lastline == line)
        u = [dir1, 0];
    else
        u = [0, dir2];
        kp = 2;
        u(1) = kp*(linecenter - p1)/((xlim(2)-xlim(1))/nlines/2);
    end
    
    u = u/sqrt(sum(u.^2));
    
    % Adjust for spherical coordinates
    if spherical
        adj = sqrt(u(2).^2+(u(1)./cosd(xt2(i))).^2);
    else
        adj = 1;
    end
    
    % Use an euler step to compute the new position
    xt1n(i) = p1 + h*umax*u(1)*adj;
    xt2n(i) = p2 + h*umax*u(2)*adj;
end

end