function found = findtargets(xt, yt, tarx, tary, r, stepprob, spherical)
%FINDTARGETS Determines which targets were found during a step
%   Models target detection using a fixed detection rate for each target 
%   while it is within a distance from the searcher. 
%
%   The proper detection probability distribution for each step is the
%   convolution of a circle with a line segment, which is approximated by 
%   taking several substeps along the line and using a circular detection 
%   area.
%   
%   The number of substeps taken is automatically computer based on the
%   radius and the distance traveled in this step.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Patrick Clary <pclary@umail.ucsb.edu>
% 1/19/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

found = [];

% Get endpoint of current and previous steps
p1 = [xt(end, :)', yt(end, :)'];
if numel(xt) > 1
    p2 = [xt(end-1, :)', yt(end-1, :)'];
else
    p2 = p1;
end
validpts = ~sum(isnan(p1) | isnan(p2), 2);
if ~any(validpts)
    return
end
p1 = p1(validpts, :);
p2 = p2(validpts, :);
dp = p1 - p2;
d = sqrt(sum(dp.^2, 2));

% Compute the number of substeps necessary to adequately approximate the
% detection probability distribution
substeps = ceil(2*max(d)/r);
if substeps > 1e3
    substeps = 1e3;
    warnstring = ['Ratio of agent step size to detection radius' ...
        'is too large. Target detection is not accurate.'];
    if ~strcmp(lastwarn(), warnstring)
        warning(warnstring);
    end
end

% Proportion of the distance from p1 to p2 for each substep
x = linspace(0, 1, substeps+1);
x = x(2:end);

% Detection probability for each substep
substepprob = 1 - (1 - stepprob)^(1/substeps);

for i = 1:substeps
    ps = p2 + x(i) * dp;
    if spherical
        distx = bsxfun(@times, bsxfun(@minus, ps(:, 1)', tarx), cosd(ps(:, 2)'));
        disty = bsxfun(@minus, ps(:, 2)', tary);
    else
        distx = bsxfun(@minus, ps(:, 1)', tarx);
        disty = bsxfun(@minus, ps(:, 2)', tary);
    end
    dists = min(sqrt(distx.^2 + disty.^2), [], 2);
    inrange = find(dists <= r);
    found = [found; inrange(rand(numel(inrange), 1) < substepprob)];
end

found = unique(found);