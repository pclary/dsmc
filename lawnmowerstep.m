function [xt1n, xt2n] = lawnmowerstep(xt1, xt2, h, umax, xlim, ylim, nlines, spherical)
%LAWNMOWERSTEP Computes the new position of each agent using a lawnmower
%algorithm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Patrick Clary <pclary@umail.ucsb.edu>
% 5/18/2014
% Updated 9/20/2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent currentcolumn;

N = size(xt1, 2);
xt1n = zeros(1, N);
xt2n = zeros(1, N);

% Initialize algorithm state
if isempty(currentcolumn) || size(xt1, 1) == 1
    % For each agent, set currentcolumn to the nearest column
    currentcolumn = floor((xt1(1, :) - xlim(1)) / (xlim(2) - xlim(1)) * (nlines - 1));
end

for i = 1:N
    % Compute the direction of vertical motion
    if mod(currentcolumn(i), 2) == 0
        ydir = 1;
    else
        ydir = -1;
    end
    
    % Switch directions if beyond the vertical boundaries
    if (ydir > 0 && xt2(end, i) > ylim(2)) || ...
        (ydir < 0 && xt2(end, i) < ylim(1))
        currentcolumn(i) = currentcolumn(i) + 1;
        ydir = -ydir;
    end
    
    % Determine x-position of the current column
    tripcount = floor(currentcolumn(i) / (nlines - 1));

    if mod(tripcount, 2) == 0
        modcol = mod(currentcolumn(i), nlines - 1);
    else
        modcol = nlines - mod(currentcolumn(i), nlines - 1);
    end
    
    linecenter = xlim(1) + modcol * (xlim(2) - xlim(1)) / (nlines);
    
    % Add phase offset to prevent covering the exact same path in each pass
    tripcount = tripcount*N + i - 1;
    n = floor(log(tripcount+1.1)/log(2));
    rem = tripcount+1-2^n;
    phase = (2*rem+1)/2^(n+1);
    linecenter = linecenter + phase / nlines;
    
    % Use proportional feedback to keep the agent moving along the line
    % center even when disturbed by a velocity field
    u = [0, ydir];
    kp = 2;
    u(1) = kp*(linecenter - xt1(end, i))/((xlim(2)-xlim(1))/nlines/2);
    
    u = u/sqrt(sum(u.^2));
    
    % Adjust for spherical coordinates
    if spherical
        adj = sqrt(u(2).^2+(u(1)./cosd(xt2(i))).^2);
    else
        adj = 1;
    end
    
    % Use an euler step to compute the new position
    xt1n(i) = xt1(end, i) + h*umax*u(1)*adj;
    xt2n(i) = xt2(end, i) + h*umax*u(2)*adj;
end

end
