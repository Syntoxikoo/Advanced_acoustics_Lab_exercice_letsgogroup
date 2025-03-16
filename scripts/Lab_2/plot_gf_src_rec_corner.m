
% ---------------------- Gf simulation -------------------
% Room dimensions
room =[3.14, 4.38, 3.27]; 

% Source position:
rS_1 = [0.16, 0.155, 0.155];


% Receiver position:
rM_1 = [0.03, room(2) - 0.03, 0.05];
[G1_sim,~] = green_func_room(rM_1,rS_1,room, 'absorption', true); 

% ---------------------- Gf Measurement -------------------
l = 0.03;
dl = 0.02;
S = 0.00113;

[f_1,G_1] = measGreen(l,dl,S,1);
f_1 = f_1';
G_1 = G_1';

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
Leg(1) = plot(f,20*log10(abs(G1_sim)/2e-5),"LineStyle",'-',"LineWidth", 1.0, "Color", corder(1,:)); 
grid on; hold on;
Leg(2) = plot(f,20*log10(abs(G_1)/2e-5),"LineStyle",'-',"LineWidth", 1.0, "Color", corder(2,:)); 

xlim([min(f) max(f)])
ylim([30 110])


% ------------------------------------- Misc for Figure -----------------------------------------------------------

leg = legend(Leg, 'Simulation', 'Measurement 1', 'NumColumns', 2); 
leg.Layout.Tile = 'north'; 

% Add common X and Y axis labels for all tiles
xlabel(tiled, 'Frequency in Hz');
ylabel(tiled, 'Level in dB SPL');

% Save the figure in EPS format (modify file name)

saveas(gcf, 'figures/gf_room_src_rec_corner.eps', 'epsc');




















