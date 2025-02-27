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
Leg(1) = plot(X, Y,"LineStyle",'-',"LineWidth", 1.0, "Color", corder(1,:)); 

grid on; hold on;

Leg(2) = plot(X, Y,"LineStyle",'-',"LineWidth", 1.0, "Color", corder(2,:)); 

% Set axis limits 
ylim([0, 1]) 
xlim([10, 1000]) 

% Add title 
title('First Tile Title');

% Uncomment these if individual tile labels are preferred
% xlabel('X-axis Label');
% ylabel('Y-axis Label');

% ------------------------------------- Second Tile ---------------------------------------------------------
nexttile 

% Set X Y...
Leg(1) = plot(X, Y,"LineStyle",'-',"LineWidth", 1.0, "Color", corder(1,:)); 

grid on; hold on;

Leg(2) = plot(X, Y,"LineStyle",'-',"LineWidth", 1.0, "Color", corder(2,:)); 

% Set axis limits 
ylim([0, 1]) 
xlim([10, 1000]) 

% Add title (optional)
title('Second Tile Title');

% Uncomment these if individual tile labels are preferred
% xlabel('X-axis Label');
% ylabel('Y-axis Label');

% ------------------------------------- Misc for Figure -----------------------------------------------------------

leg = legend(Leg, 'Dataset 1', 'Dataset 2', 'NumColumns', 2); 
leg.Layout.Tile = 'north'; 

% Add common X and Y axis labels for all tiles
xlabel(tiled, 'Common X-axis Label');
ylabel(tiled, 'Common Y-axis Label');

% Save the figure in EPS format (modify file name)
saveas(gcf, 'figures/namefig.eps', "epsc");


% OI
