clear;
close all;
clc;
addpath(genpath('scripts'))
addpath(genpath("functions"))
addpath figures/

%% Observe only certain mode
lx = 1; ly = 1;
z = linspace(0,1,10);
f = linspace(100,1000,100);

N_modes = [40,40];
xM = 1; yM = 1;
rS = [0.,0.,0];

G = Gf_duct([xM,yM],rS, z, [lx,ly], f, N_modes,[1;0 ],1);

figure;
plot(f,20*log10(abs(G)/2e-5))




% -------------------------------------------
%% Observe the cross section pressure
lx = 1; ly = 1;
N_modes = [5,5];
fmn = compute_modes(lx,ly,N_modes);


% f= 500; the result is clean

z=1;

% xM = 1; yM = 1;
xM = linspace(0,lx,100);
yM = linspace(0,ly,100);
rS = [0.,0.,0];
f = diag(fmn);
f(1) = 100;
%% Compute
G = Gf_duct([xM;yM]',rS, z, [lx,ly], f, N_modes);

%% Plot

idx= [1,2,3,4];
plot_cross_sec

%%

