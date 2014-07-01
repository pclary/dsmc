function out = B(K, xt1, xt2, hks, muks, xlim, ylim)
%B Computes the 'B' vector for each agent 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Patrick Clary <pclary@umail.ucsb.edu>
% 5/18/2014
% Updated 6/14/2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = size(xt1, 2);

sk = @(K1, K2) ck(K1, K2, xt1, xt2, hks(K1+1, K2+1), xlim, ylim) - muks(K1+1, K2+1);
[K2, K1] = meshgrid(0:K, 0:K);
sks = arrayfun(sk, K1, K2);
Las = 1./(1 + K1.^2 + K2.^2).^(3/2);
out = zeros(2, N);
for n1 = 1:K+1
    for n2 = 1:K+1
        k1 = K1(n1, n2)*pi/(xlim(2) - xlim(1));
        k2 = K2(n1, n2)*pi/(ylim(2) - ylim(1));
        for j = 1:N
            out(1, j) = out(1, j) + Las(n1, n2)*sks(n1, n2) * ...
                (-1/hks(n1, n2)*sin(k1*(xt1(end,j)-xlim(1))).*cos(k2*(xt2(end,j)-ylim(1))));
            out(2, j) = out(2, j) + Las(n1, n2)*sks(n1, n2) * ...
                (-1/hks(n1, n2)*cos(k1*(xt1(end,j)-xlim(1))).*sin(k2*(xt2(end,j)-ylim(1))));
        end
    end
end