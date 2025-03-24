% Task 2
clear;
clc;


kr = linspace(0,8*pi,1000);
N=100;

y = spatCorr(kr,N);


figure;
plot(kr,abs(y));