function out = B(xa, ya, Lasks, cres, xlim, ylim)
%B Computes the 'B' vector for each agent 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Patrick Clary <pclary@umail.ucsb.edu>
% 5/18/2014
% Updated 12/18/2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = numel(xa);

xa = xa(:)';
ya = ya(:)';

k1 = ((0:cres-1))*pi/diff(xlim)*cres/(cres+1);
k2 = ((0:cres-1))*pi/diff(ylim)*cres/(cres+1);
hk = ones(cres)*diff(xlim)*diff(ylim);
hk(1, :) = hk(1, :) * sqrt(2);
hk(:, 1) = hk(:, 1) * sqrt(2);

valx = bsxfun(@times, k1', xa - xlim(1) + diff(xlim)/cres/2);
valy = bsxfun(@times, k2', ya - ylim(1) + diff(ylim)/cres/2);
cx = cos(valx);
sx = sin(valx);
cy = cos(valy);
sy = sin(valy);

out = zeros(2, N);
test = zeros(1, N);
for i = 1:N
    test(i) = sum(sum(Lasks .* bsxfun(@times, cx(:, i)', cy(:, i)) ./ hk));
    qx = Lasks .* bsxfun(@times, k1, bsxfun(@times, sx(:, i)', cy(:, i))) ./ hk;
    qy = Lasks .* bsxfun(@times, k2', bsxfun(@times, cx(:, i)', sy(:, i))) ./ hk;
    out(1, i) = -sum(qx(:));
    out(2, i) = -sum(qy(:));
end
