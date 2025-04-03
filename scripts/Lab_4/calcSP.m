% calc sound power
load('Lab_4_data.mat'); % load data

avgDipole=mean(Dipole,1);
avghighImp=mean(highImp,1);
avglowImp=mean(lowImp,1);

P_dipole=diffFieldSP(avgDipole,RT,fc);
P_highImp=diffFieldSP(avghighImp,RT,fc);
P_lowImp=diffFieldSP(avglowImp,RT,fc);

P_dipole_wall=diffFieldSP(Dipole_wall,RT,fc);
P_highImp_wall=diffFieldSP(highImp_wall,RT,fc);
P_lowImp_wall=diffFieldSP(lowImp_wall,RT,fc);


