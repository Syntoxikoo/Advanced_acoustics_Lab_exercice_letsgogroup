% Define the number of rows and columns for tiled layout
nrows = 2;
ncols = 2;
heightScale = 1.0; % Adjust height scaling if needed

% ------------------------------------------------------------------------------------------------------------
[columnwidth, ~] = get_widths();
height = get_height() * 1.0; 
fig = figure("Position", [0, 0, columnwidth, height], "Units", "points");
tiled = tiledlayout(nrows, ncols, "TileSpacing", "tight", "Padding", "loose");
corder = colororder;

% ------------------------------------- First Tile ------------------------------------------------------------
nexttile

contourf(xM,yM,(abs(G(:,:,idx(1)))))
% Set axis limits 
ylim([0 1]) 
xlim([0 1]) 

% Add title 
title("Modes: n_x=0,n_y=0");


% ------------------------------------- Second Tile ---------------------------------------------------------
nexttile 

contourf(xM,yM,(abs(G(:,:,idx(2)))))
% Set axis limits 
ylim([0 1]) 
xlim([0 1]) 

% Add title 
title("Modes: n_x=1,n_y=1");



% ------------------------------------- third Tile ---------------------------------------------------------
nexttile 

contourf(xM,yM,(abs(G(:,:,idx(4)))))
% Set axis limits 
ylim([0 1]) 
xlim([0 1]) 

% Add title 
title("Modes: n_x=2,n_y=2");



% ------------------------------------- fourth Tile ---------------------------------------------------------
nexttile 

contourf(xM,yM,(abs(G(:,:,idx(4)))))
% Set axis limits 
ylim([0 1]) 
xlim([0 1]) 

% Add title 
title("Modes: n_x=3,n_y=3");



% ------------------------------------- Misc for Figure -----------------------------------------------------------

% leg = legend(Leg, 'Dataset 1', 'Dataset 2', 'NumColumns', 2); 
% leg.Layout.Tile = 'north'; 

% Add common X and Y axis labels for all tiles
xlabel(tiled, 'x (m)');
ylabel(tiled, 'y (m)');

% Save the figure in EPS format (modify file name)
saveas(gcf, 'figures/cross_sec_4modes.eps', "epsc");
