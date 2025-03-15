% Gf SIMULATION visualization
% Room dimensions
room =[3.14, 4.38, 3.27]; 

% Source position:
rS_1 = [0.16, 0.155, 0.155];
rS_2 = [0.16, 0.155, 0.155];
rS_3 = [0.16, 0.155, 0.155];
rS_4 = [0.16, 0.155, 0.155];
rS_5 = [0.035, 2.19, 0.155];
rS_6 = [0.13, 0.13, 1.64];

% Receiver position:
rM_1 = [0.03, room(2) - 0.03, 0.05];
[G1,~] = green_func_room(rM_1,rS_1,room, 'absorption', true); 

rM_2 = [0.025, room(2)/2, 0.05];
[G2,~] = green_func_room(rM_2,rS_2,room, 'absorption', true); 

rM_3 = [0.015, room(2)/2, 1.64];
[G3,~] = green_func_room(rM_3,rS_3,room, 'absorption', true); 

rM_4 = [1.58, room(2)/2, 1.64];
[G4,~] = green_func_room(rM_4,rS_4,room, 'absorption', true); 

rM_5 = [0.16, 0.155, 0.155];
[G5,~] = green_func_room(rM_5,rS_5,room, 'absorption', true); 

rM_6 = [0.015, 1.46, 1.64];
[G6,f] = green_func_room(rM_6,rS_6,room, 'absorption', true); 



% Gf MEASUREMENT visualization
l = 0.03;
dl = 0.02;
S = 0.00113;

[f_1,G_1] = measGreen(l,dl,S,1);
f_1 = f_1';
G_1 = G_1';

[f_2,G_2] = measGreen(l,dl,S,2);
f_2 = f_2';
G_2 = G_2';

[f_3,G_3] = measGreen(l,dl,S,3);
f_3 = f_3';
G_3 = G_3';

[f_4,G_4] = measGreen(l,dl,S,4);
f_4 = f_4';
G_4 = G_4';

[f_5,G_5] = measGreen(l,dl,S,5);
f_5 = f_5';
G_5 = G_5';

[f_6,G_6] = measGreen(l,dl,S,6);
f_6 = f_6';
G_6 = G_6';
%%
% Define the number of rows and columns for tiled layout
nrows = 1;
ncols = 1;
heightScale = 0.6; % Adjust height scaling if needed

% ------------------------------------------------------------------------------------------------------------
[columnwidth, ~] = get_widths();
height = get_height() * heightScale; 
fig = figure("Position", [0, 0, columnwidth, height], "Units", "points");
tiled = tiledlayout(nrows, ncols, "TileSpacing", "tight", "Padding", "loose");
corder = colororder;

% ------------------------------------- First Tile ------------------------------------------------------------
nexttile

% Set X Y...
Leg(1) = plot(f,20*log10(abs(G_3)/2e-5),"LineStyle",'-',"LineWidth", 1.0, "Color", corder(1,:)); 
grid on; hold on;
Leg(2) = plot(f,20*log10(abs(G_5)/2e-5),"LineStyle",'-',"LineWidth", 1.0, "Color", corder(2,:)); 


% Set axis limits 
%ylim([0, 1]) 
xlim([min(f) max(f)])


% Add title 
title('Reciprocity between source and receiver positions');
%title('Source and Receiver placed in corner');

% Uncomment these if individual tile labels are preferred
% xlabel('X-axis Label');
% ylabel('Y-axis Label');

% ------------------------------------- Second Tile ---------------------------------------------------------
% nexttile 
% % Set X Y...
% Leg(1) = plot(f,20*log10(abs(G_3)/2e-5),"LineStyle",'-',"LineWidth", 1.0, "Color", corder(1,:)); 
% 
% grid on; hold on;
% 
% %Leg(2) = plot(f,20*log10(abs(G2)/2e-5),"LineStyle",'-',"LineWidth", 1.0, "Color", corder(2,:)); 
% Leg(2) = plot(f,20*log10(abs(G_5)/2e-5),"LineStyle",'-',"LineWidth", 1.0, "Color", corder(2,:)); 
% 
% 
% % Set axis limits 
% %ylim([0, 1]) 
% xlim([min(f) max(f)])
% 
% 
% % Add title 
% 
% title('Source placed in corner and Receiver at (x ≈ 0, y/2, z/2)');
% % 
% % Uncomment these if individual tile labels are preferred
% % xlabel('X-axis Label');
% % ylabel('Y-axis Label');

% ------------------------------------- Misc for Figure -----------------------------------------------------------

leg = legend(Leg, 'Measurement 3', 'Measurement 5', 'NumColumns', 2); 
leg.Layout.Tile = 'north'; 

% Add common X and Y axis labels for all tiles
xlabel(tiled, 'Frequency in Hz');
ylabel(tiled, 'Level in dB SPL');

% Save the figure in EPS format (modify file name)
saveas(gcf, '/Users/matijamacus/Desktop/SR_M3vs5.eps', 'epsc');



% Alpha

















