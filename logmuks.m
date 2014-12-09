function muks = logmuks(x1, x2, muks, lambda)
%LOGMUKS Computes the Fourier series coefficients correspoding to the
%logarithm of the particle distribution

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Patrick Clary <pclary@umail.ucsb.edu>
% 5/18/2014
% Updated 6/30/2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xlim = [min(x1(:)), max(x1(:))];
ylim = [min(x2(:)), max(x2(:))];

% Reconstruct the particle distribution on a grid
mu = idct2(muks);
mu = mu/trapz(x2(:, 1), trapz(x1(1, :), mu, 2))*diff(xlim)*diff(ylim);

% Take the log of the distribution and renormalize
logmu = max(log(max(mu, 0)/lambda*(xlim(2)-xlim(1))*(ylim(2)-ylim(1))), 0);
logmu = logmu/trapz(x2(:, 1), trapz(x1(1, :), logmu, 2));

% Compute the new fourier coefficients for the log distribution
muks = dct2(logmu);

end

