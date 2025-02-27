% Silas Plot
clear;
close all;
clc;
addpath(genpath('scripts'))
addpath(genpath("functions"))
addpath figures/
addpath datas/
% constants
c=344;
rho=1.2;
p0=20e-6;

% set duct dimensions
a = 1;
b = 0.7;
duct=[a, b];

% Eliot constraint
constrain=false;

% max m and n
M = 60;
N = 60;
N_modes=[N, M];


% source vector
r01=[0, 0, 0];
r02=[a/2, 0, 0];

% receiver vector
rec1=[0, 0];
rec2=[a/2, 0];
z1= 0.1;
z2= 10;

f=linspace(1,600,1000);
w=2*pi*f;
% calculate Green's function
Gf11 = Gf_duct(rec1,r01,z1,duct,f,N_modes,[],1,constrain);
Gf12 = Gf_duct(rec1,r01,z2,duct,f,N_modes,[],1,constrain);

LGf11=20*log10(abs(1j.*w.*rho.*Gf11)/p0);
LGf12=20*log10(abs(1j.*w.*rho.*Gf12)/p0);

Gf21 = Gf_duct(rec2,r02,z1,duct,f,N_modes,[],1,constrain);
Gf22 = Gf_duct(rec2,r02,z2,duct,f,N_modes,[],1,constrain);

LGf21=20*log10(abs(1j.*w.*rho.*Gf21)/p0);
LGf22=20*log10(abs(1j.*w.*rho.*Gf22)/p0);

% Define the number of rows and columns for tiled layout
nrows = 1;
ncols = 2;
heightScale = 0.75; % Adjust height scaling if needed

% ------------------------------------------------------------------------------------------------------------
[columnwidth, ~] = get_widths();
height = get_height() * heightScale; 
fig = figure("Position", [0, 0, columnwidth*1.5, height], "Units", "points");
tiled = tiledlayout(nrows, ncols, "TileSpacing", "tight", "Padding", "loose");
corder = colororder;

limy=[100 200];

% ------------------------------------- First Tile ------------------------------------------------------------
nexttile

% Set X Y...
Leg(1) = plot(f, LGf11,"LineStyle",'-',"LineWidth", 1.0, "Color", corder(1,:)); 

grid on; hold on;

Leg(2) = plot(f, LGf12,"LineStyle",'--',"LineWidth", 1.0, "Color", corder(1,:)); 

% Set axis limits 
ylim(limy) 
% xlim([10, 1000]) 

% Add title 
title('Source and Receiver placed in corner');

% Uncomment these if individual tile labels are preferred
xlabel('frequency in Hz');
% ylabel('Sound pressure level in dBSPL');

% ------------------------------------- Second Tile ---------------------------------------------------------
nexttile 

Leg(1) = plot(f, LGf21,"LineStyle",'-',"LineWidth", 1.0, "Color", corder(1,:)); 

grid on; hold on;

Leg(2) = plot(f, LGf22,"LineStyle",'--',"LineWidth", 1.0, "Color", corder(1,:)); 

% Set axis limits 
ylim(limy) 
% xlim([10, 1000]) 

% Add title 
title('Source and receiver placed at nodal line');

% Uncomment these if individual tile labels are preferred
xlabel('frequency in Hz');
% ylabel('Sound pressure level in dBSPL');

% ------------------------------------- Misc for Figure -----------------------------------------------------------

leg = legend(Leg, ['z = ', num2str(z1), 'm'],['z = ', num2str(z2), 'm'],'NumColumns', 2); 
leg.Layout.Tile = 'south'; 

% Add common X and Y axis labels for all tiles
%xlabel(tiled,'frequency in Hz');
ylabel(tiled,'Level in dBSPL');

% Save the figure in EPS format (modify file name)
saveas(gcf, 'figures/PoverF_S+R_v2.eps', "epsc");


% OI
