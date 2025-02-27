clear;
close all;
clc;
addpath(genpath('scripts'))
addpath(genpath("functions"))
addpath figures/
addpath datas/

%% Observe only certain mode
lx = 1; ly = 1;
% z = linspace(0,1,10);
z=0.1;
f = linspace(1,1000,1000-1) ;

N_modes = [40,40];
xM = 1; yM = 1; % receiver
rS = [0.,0.,0]; %source

%G = Gf_duct([xM,yM],rS, z, [lx,ly], f, N_modes,[1;0],1); % use either N_modes array and empty 
G = Gf_duct([xM,yM],rS, z, [lx,ly], f, N_modes,[],2); % use either N_modes array and empty 


plot(f,20*log10(abs(G)/2e-5))


% -------------------------------------------
%% Observe the cross section pressure
lx = 0.7; ly = 1;
N_modes = [5,5];
fmn = compute_modes(lx,ly,N_modes);

z=0;

xM = linspace(0,lx,100);
yM = linspace(0,ly,100);
rS = [0.35,0.5,0];
f = diag(fmn);
f(1) = 10;
% Compute
G = Gf_duct([xM;yM]',rS, z, [lx,ly], f, N_modes);

% Plot
idx= [1,2,3,4];
plot_cross_sec

%% plot green function convergence
compute_convergence
plot_convergence

