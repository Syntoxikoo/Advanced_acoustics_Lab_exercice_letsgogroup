% Define the number of rows and columns for tiled layout
nrows = 1;
ncols = 2;
heightScale = 0.5; % Adjust height scaling if needed

% ------------------------------------------------------------------------------------------------------------
[columnwidth, ~] = get_widths();
height = get_height() * heightScale; 
fig = figure("Position", [0, 0, columnwidth, height], "Units", "points");
tiled = tiledlayout(nrows, ncols, "TileSpacing", "tight", "Padding", "loose");
corder = colororder;

% ------------------------------------- First Tile ------------------------------------------------------------
nexttile
G1 = abs(G(:,:,idx(2))/max(G(:,:,idx(2)),[],'all'));
contourf(xM,yM,abs(G(:,:,idx(2))))

% Set axis limits 
% ylim([0 ly]) 
% xlim([0 lx]) 


% Add title 
title("Modes: n_x=1,n_y=1");


% ------------------------------------- Second Tile ---------------------------------------------------------
ax2 = nexttile ;
% G2 = abs(G(:,:,idx(3))/max(G(:,:,idx(3)),[],'all'));
contourf(xM,yM,abs(G(:,:,idx(3))))

% Set axis limits 
ylim([0 ly]) 
xlim([0 lx]) 


% Add title 
title("Modes: n_x=2,n_y=2");





% ------------------------------------- Misc for Figure -----------------------------------------------------------

% leg = legend(Leg, 'Dataset 1', 'Dataset 2', 'NumColumns', 2); 
% leg.Layout.Tile = 'north'; 

% Add common X and Y axis labels for all tiles
xlabel(tiled, 'x (m)');
ylabel(tiled, 'y (m)');
ftsize = get(tiled.XLabel, 'FontSize');
c = colorbar(ax2);
c.Label.String = '|P| norm';
c.Label.FontSize = ftsize;
% Save the figure in EPS format (modify file name)
saveas(gcf, 'figures/cross_sec_2modes_srcMid.eps', "epsc");
