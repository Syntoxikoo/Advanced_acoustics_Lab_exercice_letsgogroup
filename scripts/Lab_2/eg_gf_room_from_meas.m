% This script describe the way to compute the gren function for a room from meas


l = 0.03;
dl = 0.02; % Distance between mic
S = 0.00113; % Surface of the tube
fileN = 1; % file index
p0=20e-6;

[f,G] = measGreen(l,dl,S,fileN);

figure;
plot(f,20*log10(abs(G)/p0));
grid on
title("this is an example of the green function from measurement")

%% Dump
% [f, Hxy] = importHxy('');