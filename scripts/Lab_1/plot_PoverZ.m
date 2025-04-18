% Silas Plot


% set duct dimensions
a = 0.7;
b = 1;
duct=[a, b];

% Eliot constraint
constrain=false;

% max m and n
M = 8;
N = 8;
N_modes=[N, M];

% set frequency
f1=160;
w1=2*pi*f1;

f2=590;
w2=2*pi*f2;

% constants
c=344;
rho=1.2;
p0=20e-6;

% source vector
x0 = 0;
y0 = 0;
r0=[x0, y0, 0];

% receiver vector
x= 0;
y= 0;
rec=[x, y];
z= linspace(0,10,1000);


% Calculate Green's function
Gz1 = Gf_duct(rec,r0,z,duct,f1,N_modes,[],1,constrain);
Gz2 = Gf_duct(rec,r0,z,duct,f2,N_modes,[],1,constrain);
% log of Green's function
LGz1=20*log10(abs(Gz1)/p0);
LGz2=20*log10(abs(Gz2)/p0);

% Define the number of rows and columns for tiled layout
nrows = 1;
ncols = 1;
heightScale =0.75; % Adjust height scaling if needed

% ------------------------------------------------------------------------------------------------------------
[columnwidth, ~] = get_widths();
height = get_height() * heightScale; 
fig = figure("Position", [0, 0, columnwidth*1.5, height], "Units", "points");
tiled = tiledlayout(nrows, ncols, "TileSpacing", "tight", "Padding", "loose");
corder = colororder;

% ------------------------------------- First Tile ------------------------------------------------------------
nexttile

% Set X Y...
Leg(1) = plot(z, LGz1,"LineStyle",'-',"LineWidth", 1.0, "Color", corder(1,:)); 

grid on; hold on;

Leg(2) = plot(z, LGz2,"LineStyle",'--',"LineWidth", 1.0, "Color", corder(1,:)); 

% Set axis limits 
% ylim([0, 1]) 
% xlim([10, 1000]) 

% Add title 
title('Sound pressure over distance');

% Uncomment these if individual tile labels are preferred
xlabel('distance in m');
ylabel('Level in dBSPL');

% ------------------------------------- Second Tile ---------------------------------------------------------
% nexttile 
% 
% % Set X Y...
% Leg(1) = plot(X, Y,"LineStyle",'-',"LineWidth", 1.0, "Color", corder(1,:)); 
% 
% grid on; hold on;
% 
% Leg(2) = plot(X, Y,"LineStyle",'-',"LineWidth", 1.0, "Color", corder(2,:)); 
% 
% % Set axis limits 
% ylim([0, 1]) 
% xlim([10, 1000]) 
% 
% % Add title (optional)
% title('Second Tile Title');

% Uncomment these if individual tile labels are preferred
% xlabel('X-axis Label');
% ylabel('Y-axis Label');

% ------------------------------------- Misc for Figure -----------------------------------------------------------

leg = legend(Leg, ['f = ', num2str(f1), 'Hz'],['f = ', num2str(f2), 'Hz'],'NumColumns', 2); 
leg.Layout.Tile = 'south'; 

% Add common X and Y axis labels for all tiles
% xlabel(tiled, 'Common X-axis Label');
% ylabel(tiled, 'Common Y-axis Label');

% Save the figure in EPS format (modify file name)
saveas(gcf, 'figures/PoverZ_S+R-Corner.eps', "epsc");


% OI
