function [xp, yp] = sampledist(dist, x, y, n, sampfun)
%SAMPLEDIST Samples a probability distribution

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Patrick Clary <pclary@umail.ucsb.edu>
% 5/18/2014
% Updated 6/23/2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 5
    s1 = halton(2, n);
    s2 = halton(3, n);
else
    s1 = sampfun(n);
    s2 = sampfun(n);
end

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

