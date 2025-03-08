clear;
close all;
clc;
addpath(genpath('scripts'))
addpath(genpath("functions"))
addpath figures/
addpath datas/



% -------------------------------------------
%% Observe the cross section pressure
lx = 0.7; ly = 1;
N_modes = [5,5];
fmn = compute_modes(lx,ly,N_modes);

z=0;

xM = linspace(0,lx,100);
yM = linspace(0,ly,100);
rS = [0.,0.,0];
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


%% plot velocities
f = 1:10:1000;
[Cp,Cg,nx] = calc_modes_vel([0.7,1],f,[9,9]);

n = linspace(0,length(nx),length(nx));
figure;
surf(f,n,20*log10(abs(Cp)))
xlabel("frequency (Hz)")
ylabel("mode number")
title("Phase velocity in function of frequency and number of modes")

figure;
surf(f,n,20*log10(abs(Cg)))
xlabel("frequency (Hz)")
ylabel("mode number")
title("Group velocity in function of frequency and number of modes")
