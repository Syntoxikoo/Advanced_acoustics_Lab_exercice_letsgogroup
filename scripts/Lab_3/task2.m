% Task 2
clear;
clc;


kr = linspace(0,8*pi,1000);

N=100;

[p1p2,phi] = spatCorr(kr,N);

 msP = meanSqP(N,phi);

 y = p1p2 ./ msP;


figure;
plot(kr,abs(y));