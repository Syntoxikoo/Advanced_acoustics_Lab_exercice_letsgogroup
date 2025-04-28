clear;
close all;
clc;
addpath(genpath('scripts'))
addpath(genpath("functions"))
addpath(genpath("datas"))

%% ---- Import Data
load('Lab_4_data.mat'); 
calcSP;
% plotSP(fc,P_dipole,P_dipole_wall,'Dipole');
% plotSP(fc,P_highImp,P_highImp_wall,'High Impedance Source');
% plotSP(fc,P_lowImp,P_lowImp_wall,'Low Impedance Source');
% RT = 20 * log10(RT./2e-5);
plotRT(fc,RT,'Reverberation_time_T60');