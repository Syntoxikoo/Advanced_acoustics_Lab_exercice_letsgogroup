clear;
close all;
clc;
addpath(genpath('scripts'))
addpath(genpath("functions"))
addpath(genpath("datas"))


%% --- SETTINGS ---
load('Lab_4_data.mat'); 
load('Lab_5_data.mat'); 
datas = [Dipole_IS; Dipole_fan_IS; Fan_IS; Omni_IS; Omni_1_IS]; % gather data for lab 5
calc_SP_diff;

surf = 1;
calc_SP_I

% Account for windscreen
windscreen = zeros(1,length(Dipole_IS));
windscreen(end - 10:end) = [1.4 1.4 1.4 1.7 1.8 1.8 2.2 2.4 3 3.3 3.6];

%% -- Plot -- 


plotSP2(fc, SP_dipole,"SoundPowerDipole",SP_dipole_wall,["Diff. Field","Diff. Field wall", "Intensity scan"],'P3',SP_datas(1,:));
plotSP2(fc, SP_highImp,"SoundPowerOmni",SP_highImp_wall,["Diff. Field","Diff. Field wall", "Intensity scan"],'P3',SP_datas(4,:));
plotSP2(fc, SP_lowImp+windscreen,"SoundPowerFan",SP_lowImp_wall+windscreen,["Diff. Field","Diff. Field wall", "Intensity scan"],'P3',SP_datas(3,:)+windscreen);
plotSP2(fc, SP_datas(4,:),"SoundPower_for_2_operator",SP_datas(5,:),["Operator A","Operator B"]);
plotSP3(fc, SP_datas(1,:),"SoundPower_pI_noisy_meas",SP_datas(2,:),["Dipole","Dipole with noise source"],"PI1",Dipole_PI,"PI2",Dipole_fan_PI);

