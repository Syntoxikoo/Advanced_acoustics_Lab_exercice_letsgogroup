function h = sphankel2(m, x)
% Compute spherical Hankel function of the second kind
h = sqrt(pi./(2*x)) .* besselh(m+0.5, 2, x);
end