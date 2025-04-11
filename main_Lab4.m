clear;
close all;
clc;
addpath(genpath('scripts'))
addpath(genpath("functions"))
addpath(genpath("datas"))

%% ---- Import Data

calcSP;
plotSP(fc,P_dipole,P_dipole_wall,'Dipole');
plotSP(fc,P_highImp,P_highImp_wall,'High Impedance Source');
plotSP(fc,P_lowImp,P_lowImp_wall,'Low Impedance Source');
plotRT(fc,RT,'Reverberation time T60');