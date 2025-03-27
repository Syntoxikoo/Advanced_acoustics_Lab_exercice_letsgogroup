% Define the number of rows and columns for tiled layout
nrows = 1;
ncols = 2;
heightScale = 0.6; % Adjust height scaling if needed

% ------------------------------------------------------------------------------------------------------------
[columnwidth, ~] = get_widths();
height = get_height() * heightScale; 
fig = figure("Position", [0, 0, columnwidth, height], "Units", "points");
tiled = tiledlayout(nrows, ncols, "TileSpacing", "tight", "Padding", "loose");
corder = colororder;

% ------------------------------------- First Tile ------------------------------------------------------------
nexttile


bins = 30;
Leg(1) = histogram(X_11, bins,'Normalization', 'pdf','FaceAlpha',0.4); hold on;
plot(x_X11, pd_X11, 'LineWidth', 2,"Color",corder(1,:))

Leg(2) = histogram(X_1, bins,'Normalization', 'pdf','FaceAlpha',0.4); 
plot(x_X1, pd_X1, 'LineWidth', 2,"Color",corder(2,:))

% Add title 
title('Distrib pressure rms');

% Uncomment these if individual tile labels are preferred
% xlabel('X-axis Label');
% ylabel('Y-axis Label');
xlabel('p_{rms}^2');
% ------------------------------------- Second Tile ---------------------------------------------------------
nexttile 

bins = 30;
Leg(1) = histogram(X_21, bins,'Normalization', 'pdf','FaceAlpha',0.4); hold on;
plot(x_X21, pd_X21, 'LineWidth', 2,"Color",corder(1,:))

Leg(2) = histogram(X_2, bins,'Normalization', 'pdf','FaceAlpha',0.4); 
plot(x_X2, pd_X2, 'LineWidth', 2,"Color",corder(2,:))
xlim([50 100])
% Add title 
title('Distrib pressure level');
xlabel('Lp (dB re 20\mu Pa)');
% ------------------------------------- Misc for Figure -----------------------------------------------------------

leg = legend(Leg, "Green's f estimate", 'Pwave exp estimate', 'NumColumns', 2); 
leg.Layout.Tile = 'north'; 

% Add common X and Y axis labels for all tiles

ylabel(tiled, 'Probability density');

% Save the figure in EPS format (modify file name)
saveas(gcf, 'figures/namefig.eps', "epsc");


