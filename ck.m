function coeff = ck(K1, K2, xt1, xt2, xlim, ylim)
%CK Computes the (K1, K2)th Fourier coefficient ck for the given
%trajectory data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Patrick Clary <pclary@umail.ucsb.edu>
% 5/18/2014
% Updated 6/30/2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


N = size(xt1, 2);

% Precompute K*pi/L
k1 = K1*pi/(xlim(2) - xlim(1));
k2 = K2*pi/(ylim(2) - ylim(1));

hk = diff(xlim)*diff(ylim);
if K1 > 0
    hk = hk / 2;
end
if K2 > 0
    hk = hk / 2;
end

coeff = zeros(N, 1);

% Compute (K1, K2)th coefficient of the fourier transform of the prior
for j = 1:N
    i = find(xt1(:,j) >= xlim(1) & xt1(:,j) <= xlim(2) & ...
        xt2(:,j) >= ylim(1) & xt2(:,j) <= ylim(2));
    c = 1/hk.*cos(k1*(xt1(i,j)-xlim(1))).*cos(k2*(xt2(i,j)-ylim(1)));
    coeff(j) = coeff(j) + sum(c) - (c(1) + c(end))/2;
    coeff(j) = coeff(j)/N/numel(i);
end

coeff = sum(coeff);