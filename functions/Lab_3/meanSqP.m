function y = meanSqP(N,phi)


phi=rand(1,N).*2.*pi;


y=1/2/N*abs(sum(exp(1j.*phi)));

end