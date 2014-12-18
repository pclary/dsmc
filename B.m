function out = B(xa1, xa2, Lasks, xlim, ylim)
%B Computes the 'B' vector for each agent 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Patrick Clary <pclary@umail.ucsb.edu>
% 5/18/2014
% Updated 12/17/2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = numel(xa1);

xa1 = xa1(:)';
xa2 = xa2(:)';

% Add up fourier terms to get B vector for each agent
out = zeros(2, N);
for n1 = 1:size(sks, 1)
    for n2 = 1:size(sks, 2)
        k1 = K1(n1, n2)*pi/diff(xlim);
        k2 = K2(n1, n2)*pi/diff(ylim);
        hk = diff(xlim)*diff(ylim);
        if K1(n1, n2) > 0
            hk = hk / 2;
        end
        if K2(n1, n2) > 0
            hk = hk / 2;
        end
        
        out(1, :) = out(1, :) + Lasks(n1, n2) * ...
            (-k1/hk*sin(k1*(xa1-xlim(1))).*cos(k2*(xa2-ylim(1))));
        out(2, :) = out(2, :) + Lasks(n1, n2) * ...
            (-k2/hk*cos(k1*(xa1-xlim(1))).*sin(k2*(xa2-ylim(1))));
    end
end