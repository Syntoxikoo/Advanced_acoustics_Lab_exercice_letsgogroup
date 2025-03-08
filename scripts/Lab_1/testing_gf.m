addpath(genpath('scripts'))
addpath(genpath("functions"))
addpath figures/


lx = 1; ly = 1;
z = linspace(0,1,100);
f = linspace(0,8000,100);
N_modes = [20,30];
xM = 1; yM = 1;
rS = [0.,0.,0];

% Results with method 2
G = Gf_duct([xM,yM],rS, z, [lx,ly], f, N_modes,[],2);
figure(1);
surf(f,z,20*log10(abs(G)/2e-5),'LineStyle','none')
set(gca, 'XScale', 'log');
colorbar;


% Results with method 1
G = Gf_duct([xM,yM],rS, z, [lx,ly], f, N_modes,[],1);
figure(2);
surf(f,z,20*log10(abs(G)/2e-5))
set(gca, 'XScale', 'log');
colorbar;
