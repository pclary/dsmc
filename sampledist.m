function [xp, yp] = sampledist(dist, n, x, y, sampfun)
%SAMPLEDIST Samples a probability distribution

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Patrick Clary <pclary@umail.ucsb.edu>
% 5/18/2014
% Updated 12/19/2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Sample using the supplied sample function or the default halton sequence;
% Samples are points in the domain [0, 1] x [0, 1]
if nargin < 5
    s1 = halton(2, n);
    s2 = halton(3, n);
else
    s1 = sampfun(n);
    s2 = sampfun(n);
end

% Set output domain
if nargin < 4 || (isempty(x) && isempty(y))
    x = [0, 1];
    y = [0, 1];
end
if numel(x) == 2 && numel(y) == 2 && any(size(dist) ~= [2, 2])
    % Treat x and y as bounds, not grid vectors
    x = linspace(x(1), x(2), size(dist, 2));
    y = linspace(y(1), y(2), size(dist, 1));
end
    
% Construct an approximate 2D cumulative probability density function
% respresenting the distribution
w = cumtrapz(trapz(dist, 2));
if w(end) ~= 0
    w = w/w(end);
else
    w = linspace(0, 1, numel(w))';
end
ww = cumtrapz(dist, 2);
ww = bsxfun(@rdivide, ww, ww(:, end));
index = find(isnan(ww(:, end)));
for i = index'
    ww(i, :) = linspace(0, 1, size(ww, 2));
end

xp = zeros(n, 1);
yp = zeros(n, 1);

% Map points in [0, 1] x [0, 1] to cartesian points using the cumulative
% probability distribution function; uses interpolation
for k = 1:n
    i = find(w > s1(k), 1);
    f = (s1(k)-w(i-1))/(w(i)-w(i-1));
    yp(k) = y(i-1) + (y(i)-y(i-1))*f;
    wwi = ww(i-1,:) + (ww(i,:) - ww(i-1,:))*f;
    wwi = wwi/wwi(end);
    i = find(wwi > s2(k), 1);
    f = (s2(k)-wwi(i-1))/(wwi(i)-wwi(i-1));
    xp(k) = x(i-1) + (x(i)-x(i-1))*f;
end

end

