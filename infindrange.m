function out = infindrange(xt, yt, tarx, tary, r)

p1 = [xt(end, :)', yt(end, :)'];
if numel(xt) > 1
    p2 = [xt(end-1, :)', yt(end-1, :)'];
else
    p2 = p1;
end

dp = p1 - p2;
d = sqrt(sum(dp.^2, 2));
v = bsxfun(@rdivide, [dp(2), -dp(1)]*r, d);
rect = [p1 + v, p1 - v, p2 - v, p2 + v, ones(numel(d), 2)*NaN];
rect = reshape(rect', 2, numel(d)*5)';

dists1 = min(sqrt(bsxfun(@minus, p1(:, 1)', tarx).^2 + ...
    bsxfun(@minus, p1(:, 2)', tary).^2), [], 2);
dists2 = min(sqrt(bsxfun(@minus, p2(:, 1)', tarx).^2 + ...
    bsxfun(@minus, p2(:, 2)', tary).^2), [], 2);
if numel(dists1) == 0
    dists1 = zeros(0, 1);
    dists2 = zeros(0, 1);
end

out = dists1 <= r | (inpolygon(tarx, tary, rect(:, 1), rect(:, 2)) & dists2 > r);
