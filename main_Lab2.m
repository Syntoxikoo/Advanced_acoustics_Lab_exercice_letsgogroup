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
S = 0.00113;
fileN = 1;
p0=20e-6;

[f, Hxy] = importHxy('');
[f,G] = measGreen(l,dl,S,fileN);

figure;
plot(f,20*log10(abs(G)/p0));
xlim([100 2000])
grid on

%% plot_room_config
plot_room_config()

%% plot cross section with green function
% This part serve just as an exemple, for understanding how the Gf simulation work
% for plotting use template !
compute_gf_cross_sec
%% plot green function for a specific pair of source receiver
% This part serve just as an exemple, for understanding how the Gf simulation work
% for plotting use template !
Plotting_n_1
%%

