function y = meanSqP(N,phi)
% calculate the mean square pressure


if nargin<2
    phi=rand(1,N).*2.*pi;
end

y=1/2/N*abs(sum(exp(1j.*phi)));

end