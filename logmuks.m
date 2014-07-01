function muks = logmuks(x1, x2, mu1, mu2, lambda, K, fks, hks)
%LOGMUKS Computes the Fourier series coefficients correspoding to the
%logarithm of the particle distribution

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Patrick Clary <pclary@umail.ucsb.edu>
% 5/18/2014
% Updated 6/14/2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xlim = [min(x1(:)), max(x1(:))];
ylim = [min(x2(:)), max(x2(:))];

mu = zeros(size(x1));
for K1 = 0:K
    for K2 = 0:K
        mu = mu + fks{K1+1, K2+1}*muk(K1, K2, mu1, mu2, hks(K1+1, K2+1), xlim, ylim);
    end
end

logmu = max(real(log(mu/lambda*(xlim(2)-xlim(1))*(ylim(2)-ylim(1)))), 0);
logmu = logmu/trapz(x2(:, 1), trapz(x1(1, :), logmu, 2));

muks = zeros(K+1);
for i = 1:K+1
    for j = 1:K+1
        muks(i, j) = trapz(x2(:, 1), trapz(x1(1, :), fks{i, j}.*logmu, 2), 1);
    end
end

end

