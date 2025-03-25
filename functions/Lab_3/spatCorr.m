function [y,phi] = spatCorr(kr,N)
% function to calculate the expression for the spatial correlation p1*p2'/2


y=zeros(1,length(kr));

phi=rand(1,N).*2.*pi;
theta=asin(2.*(rand(1,N)-0.5))+pi./2;

for idx = 1:length(kr)

    wavSum = (sum( exp(1j .* (phi-kr(idx) .* cos(theta) ))));
    y(idx) = 1/(2*N) * sum(exp(1j.*phi)) * wavSum';
end

end