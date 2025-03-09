clear;
close all;
clc;
addpath(genpath('scripts'))
addpath(genpath("functions"))
addpath(genpath("datas"))
%%
clc;

l = 0.03;
dl = 0.02;
S = 0.00113;s
fileN = 1;

[f,G] = measGreen(l,dl,S,fileN);

figure;
plot(f,20*log10(abs(G)));
xlim([100 2000])
grid on