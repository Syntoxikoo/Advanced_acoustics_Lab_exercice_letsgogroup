% Silas Plot


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

% constants
c=344;
rho=1.2;
p0=20e-6;

% source vector
x0 = a/2;
y0 = 0;
r0=[x0, y0, 0];

% receiver vector
x= a/2;
y= 0;
rec=[x, y];
z1= 0.1;
z2= 10;

f=linspace(1,600,1000);

% calculate Green's function
Gf1 = Gf_duct(rec,r0,z1,duct,f,N_modes,[],1,constrain);
Gf2 = Gf_duct(rec,r0,z2,duct,f,N_modes,[],1,constrain);

LGf1=20*log10(abs(Gf1)/p0);
LGf2=20*log10(abs(Gf2)/p0);

% Define the number of rows and columns for tiled layout
nrows = 1;
ncols = 1;
heightScale =0.5; % Adjust height scaling if needed

% ------------------------------------------------------------------------------------------------------------
[columnwidth, ~] = get_widths();
height = get_height() * 1.0; 
fig = figure("Position", [0, 0, columnwidth, height], "Units", "points");
tiled = tiledlayout(nrows, ncols, "TileSpacing", "tight", "Padding", "loose");
corder = colororder;

% ------------------------------------- First Tile ------------------------------------------------------------
nexttile

% Set X Y...
Leg(1) = plot(f, LGf1,"LineStyle",'-',"LineWidth", 1.0, "Color", corder(1,:)); 

grid on; hold on;

Leg(2) = plot(f, LGf2,"LineStyle",'--',"LineWidth", 1.0, "Color", corder(1,:)); 

% Set axis limits 
% ylim([0, 1]) 
% xlim([10, 1000]) 

% Add title 
title('Sound pressure over frequency');

% Uncomment these if individual tile labels are preferred
xlabel('frequency in Hz');
ylabel('Sound pressure level in dBSPL');

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

leg = legend(Leg, ['z = ', num2str(z1), 'm'],['f = ', num2str(z2), 'm'],'NumColumns', 2); 
leg.Layout.Tile = 'south'; 

% Add common X and Y axis labels for all tiles
% xlabel(tiled, 'Common X-axis Label');
% ylabel(tiled, 'Common Y-axis Label');

% Save the figure in EPS format (modify file name)
saveas(gcf, 'figures/PoverF_S+R-halfA.eps', "epsc");


% OI
