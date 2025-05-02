function h = sphankel2(m, r)
% Compute spherical Hankel function of the second kind
    h = sqrt(pi./(2*r)) .* besselh(m+0.5, 2, r);    
end