
% ---------------------- Gf simulation -------------------
% Room dimensions
room =[3.14, 4.38, 3.27]; 

% Source position:
rS_3 = [0.16, 0.155, 0.155];
rS_5 = [0.035, 2.19, 0.155];
% Receiver position:
rM_3 = [0.015, room(2)/2, 1.64];
rM_5 = [0.16, 0.155, 0.155];


[G3_sim,~] = green_func_room(rM_3,rS_3,room, 'absorption', true); 

[G5_sim,~] = green_func_room(rM_5,rS_5,room, 'absorption', true); 

% ---------------------- Gf Measurement -------------------
l = 0.03;
dl = 0.02;
S = 0.00113;

[f_3,G_3] = measGreen(l,dl,S,3);
f_3 = f_3';
G_3 = G_3';


[f_5,G_5] = measGreen(l,dl,S,5);
f_5 = f_5';
G_5 = G_5';


%% Plotting
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
Leg(1) = plot(f,20*log10(abs(G_3)/2e-5),"LineStyle",'-',"LineWidth", 1.0, "Color", corder(2,:)); 
grid on; hold on;
plot(f,20*log10(abs(G3_sim)/2e-5),"LineStyle",'--',"LineWidth", 1.0, "Color", corder(1,:)); 

Leg(2) = plot(f,20*log10(abs(G_5)/2e-5),"LineStyle",'-',"LineWidth", 1.0, "Color", corder(3,:)); 
plot(f,20*log10(abs(G_5)/2e-5),"LineStyle",'--',"LineWidth", 1.0, "Color", corder(4,:)); 

xlim([min(f) max(f)])


% ------------------------------------- Misc for Figure -----------------------------------------------------------

leg = legend(Leg, 'Measurement 3', 'Measurement 5', 'NumColumns', 2); 
leg.Layout.Tile = 'north'; 

% Add common X and Y axis labels for all tiles
xlabel(tiled, 'Frequency in Hz');
ylabel(tiled, 'Level in dB SPL');

% Save the figure in EPS format (modify file name)

saveas(gcf, 'figures/reciprocity_src_rec.eps', 'epsc');




















