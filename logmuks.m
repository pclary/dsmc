function muks = logmuks(muks0, lambda, res, xlim, ylim)
%LOGMUKS Computes the Fourier series coefficients correspoding to the
%logarithm of the particle distribution

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Patrick Clary <pclary@umail.ucsb.edu>
% 5/18/2014
% Updated 6/30/2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Change resolution
minres = min(size(muks0, 1), res);
muks1 = zeros(res);
muks1(minres, minres) = muks0(minres, minres);

% Reconstruct the particle distribution on a grid
mu = idct2(muks1);

% Take the log of the distribution and renormalize
logmu = max(log(max(mu, 0)/lambda*(xlim(2)-xlim(1))*(ylim(2)-ylim(1))), 0);
logmu = logmu / trapz(trapz(logmu)) * diff(xlim)*diff(ylim)/res^2;

% Compute the new fourier coefficients for the log distribution
muks = dct2(logmu);
