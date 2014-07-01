function coeff = muk(K1, K2, mu1, mu2, hk, xlim, ylim)
%CK Computes the (K1, K2)th Fourier coefficient muk for the given
%particle locations

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Patrick Clary <pclary@umail.ucsb.edu>
% 5/18/2014
% Updated 6/30/2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Precompute K*pi/L
k1 = K1*pi/(xlim(2)-xlim(1));
k2 = K2*pi/(ylim(2)-ylim(1));

coeff = sum(1/hk.*cos(k1*(mu1-xlim(1))).*cos(k2*(mu2-ylim(1))));

coeff = coeff/size(mu1, 1);