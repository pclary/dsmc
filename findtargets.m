function found = findtargets(xt, yt, tarx, tary, r, stepprob)
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
% 5/18/2014
% Updated 12/18/2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get endpoint of current and previous steps
p1 = [xt(end, :)', yt(end, :)'];
if numel(xt) > 1
    p2 = [xt(end-1, :)', yt(end-1, :)'];
else
    p2 = p1;
end
dp = p1 - p2;
d = sqrt(sum(dp.^2, 2));

% Compute the number of substeps necessary to adequately approximate the
% detection probability distribution
substeps = ceil(2*max(d)/r);

% Proportion of the distance from p1 to p2 for each substep
x = linspace(0, 1, substeps+1);
x = x(2:end);

% Detection probability for each substep
substepprob = 1 - (1 - stepprob)^(1/substeps);

found = [];

for i = 1:substeps
    ps = p2 + x(i) * dp;
    dists = min(sqrt(bsxfun(@minus, ps(:, 1)', tarx).^2 + ...
        bsxfun(@minus, ps(:, 2)', tary).^2), [], 2);
    inrange = find(dists <= r);
    found = [found; inrange(rand(numel(inrange), 1) < substepprob)];
end

found = unique(found);