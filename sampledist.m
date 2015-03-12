function [xp, yp] = sampledist(dist, n, x, y, sampfun)
%SAMPLEDIST Samples a probability distribution

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Patrick Clary <pclary@umail.ucsb.edu>
% 5/18/2014
% Updated 12/19/2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Ensure that the distribution has no negative values
dist(dist < 0) = 0;

% Sample using the supplied sampling function or the default halton sequence;
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
    % Create center grid
    dx = diff(x)/size(dist, 2);
    dy = diff(y)/size(dist, 1);
    x = x(1) - dx/2:dx:x(2) + dx/2;
    y = y(1) - dy/2:dy:y(2) + dy/2;
else
    dx = diff(x(1:2));
    dy = diff(y(1:2));
end
    
% Construct an approximate 2D cumulative probability density function
% respresenting the distribution
weights = trapz(dist, 2);
w = cumtrapz(weights);
w = w + weights(1) / 2;
if w(end) ~= 0
    w = w/(w(end) + weights(end)/2);
else
    l1 = numel(q);
    w = linspace(1/l1/2, 1 - 1/l1/2, l1)';
end
weights = [0; weights; 0];
w = [-w(1); w; 2 - w(end)];
ww = cumtrapz(dist, 2);
ww = bsxfun(@plus, ww, dist(:, 1) / 2);
ww = bsxfun(@rdivide, ww, ww(:, end) + dist(:, end)/2);
index = find(isnan(ww(:, end)));
l2 = size(ww, 2);
for i = index'
    ww(i, :) = linspace(1/l2/2, 1 - 1/l2/2, l2);
end
ww = [linspace(1/l2/2, 1 - 1/l2/2, l2); ww; linspace(1/l2/2, 1 - 1/l2/2, l2)];
ww = [-ww(:, 1), ww, 2 - ww(:, end)];

xp = zeros(n, 1);
yp = zeros(n, 1);

% Map points in [0, 1] x [0, 1] to cartesian points using the cumulative
% probability distribution function; uses interpolation
for k = 1:n
    i = find(w > s1(k), 1);
    f = (s1(k)-w(i-1))/(w(i)-w(i-1));
    yp(k) = y(i-1) + (y(i)-y(i-1))*f;
    wwi = ww(i-1,:)*weights(i-1) + (ww(i,:)*weights(i) - ww(i-1,:)*weights(i-1))*f;
    wwi2 = wwi / wwi(end);
    i2 = find(wwi2 > s2(k), 1);
    f = (s2(k)-wwi2(i2-1))/(wwi2(i2)-wwi2(i2-1));
    xp(k) = x(i2-1) + (x(i2)-x(i2-1))*f;
end

end
