function s = halton(base, n)
%HALTON Generates a Halton sequence using the specified base

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Patrick Clary <pclary@umail.ucsb.edu>
% 5/18/2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s = zeros(n, 1);

for i = 1:n
    q = i;
    x = 0;
    f = 1/base;
    
    while q > 0
        r = mod(q, base);
        q = floor(q/base);
        x = x + f*r;
        f = f/base;
    end
    s(i) = x;
end

end

