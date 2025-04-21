% calc sound power

% Avg over the 3 microphone position
avgDipole=mean(Dipole,1);
avghighImp=mean(highImp,1);
avglowImp=mean(lowImp,1);

% Compute Sound power from pressure datas
P_dipole=diffFieldSP(avgDipole,RT,fc);
P_highImp=diffFieldSP(avghighImp,RT,fc);
P_lowImp=diffFieldSP(avglowImp,RT,fc);

% Compute Sound power from pressure datas near the wall
P_dipole_wall=diffFieldSP(Dipole_wall,RT,fc);
P_highImp_wall=diffFieldSP(highImp_wall,RT,fc);
P_lowImp_wall=diffFieldSP(lowImp_wall,RT,fc);


