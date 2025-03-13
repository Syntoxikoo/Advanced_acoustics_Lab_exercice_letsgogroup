% Gf SIMULATION visualization
% Room dimensions
room =[3.14, 4.38, 3.27]; 

% Source position:

%rS = [0, 0, 0];
rS = [0.16, 0.155, 0.155];

% Observation point #n.1 coordinates:
x = 0;
y = room(2); 
z = 0;
rM = [x,y,z];
[G1,~] = green_func_room(rM,rS,room, 'absorption', true); 


% % Observation point #n.2 coordinates:
% x = room(1)/2; 
% y = room(2)/2; 
% z = room(3)/2;
% rM = [x,y,z];
% [G2,f] = green_func_room(rM,rS,room, 'absorption', true); 



% Gf MEASUREMENT visualization
l = 0.03;
dl = 0.02;
S = 0.00113;
fileN = 1;

[f_m,G_m] = measGreen(l,dl,S,fileN);
f_m = f_m';
G_m = G_m';

%%
% Define the number of rows and columns for tiled layout
nrows = 1;
ncols = 1;
heightScale = 0.5; % Adjust height scaling if needed

% ------------------------------------------------------------------------------------------------------------
[columnwidth, ~] = get_widths();
height = get_height() * heightScale; 
fig = figure("Position", [0, 0, columnwidth, height], "Units", "points");
tiled = tiledlayout(nrows, ncols, "TileSpacing", "tight", "Padding", "loose");
corder = colororder;

% ------------------------------------- First Tile ------------------------------------------------------------
nexttile

% Set X Y...
Leg(1) = plot(f,20*log10(abs(G1)/2e-5),"LineStyle",'-',"LineWidth", 1.0, "Color", corder(1,:)); 

grid on; hold on;

%Leg(2) = plot(f,20*log10(abs(G2)/2e-5),"LineStyle",'-',"LineWidth", 1.0, "Color", corder(2,:)); 
Leg(2) = plot(f,20*log10(abs(G_m)/2e-5),"LineStyle",'-',"LineWidth", 1.0, "Color", corder(2,:)); 


% Set axis limits 
%ylim([0, 1]) 
xlim([min(f) max(f)])


% Add title 
title('First Tile Title');

% Uncomment these if individual tile labels are preferred
% xlabel('X-axis Label');
% ylabel('Y-axis Label');

% ------------------------------------- Second Tile ---------------------------------------------------------
% nexttile 
% 
% % Set X Y...
% Leg(1) = plot(f_m,20*log10(abs(G_m)/2e-5),"LineStyle",'-',"LineWidth", 1.0, "Color", corder(1,:)); 
% 
% grid on; hold on;
% 
% %Leg(2) = plot(X, Y,"LineStyle",'-',"LineWidth", 1.0, "Color", corder(2,:)); 
% 
% % Set axis limits 
% %ylim([0, 1]) 
% xlim([min(f) max(f)])
% 
% % Add title (optional)
% title('Second Tile Title');
% 
% Uncomment these if individual tile labels are preferred
% xlabel('X-axis Label');
% ylabel('Y-axis Label');

% ------------------------------------- Misc for Figure -----------------------------------------------------------

leg = legend(Leg, 'Simulation', 'Measurement', 'NumColumns', 2); 
leg.Layout.Tile = 'north'; 

% Add common X and Y axis labels for all tiles
xlabel(tiled, 'Frequency - Linear / Hz');
ylabel(tiled, 'G - Magnitude');

% Save the figure in EPS format (modify file name)
%saveas(gcf, 'figures/namefig.eps', "epsc");


% Alpha

















