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


bins = 30;
Leg(1) = histogram(X_11, bins,'Normalization', 'pdf','FaceAlpha',0.4); hold on;
plot(x_X11, pd_X11, 'LineWidth', 2,"Color",corder(1,:))
Leg(2) = histogram(X_12, bins,'Normalization', 'pdf','FaceAlpha',0.4); 


% Add title 
title('Mean square pressure distribution');

xlabel('p_{rms}^2');
% ------------------------------------- Second Tile ---------------------------------------------------------
nexttile 

bins = 30;
Leg(1) = histogram(X_21, bins,'Normalization', 'pdf','FaceAlpha',0.4); hold on;
plot(x_X21, pd_X21, 'LineWidth', 2,"Color",corder(1,:))
Leg(2) = histogram(X_22, bins,'Normalization', 'pdf','FaceAlpha',0.4); 

xlim([50 100])
% Add title 
title('Pressure level distribution');
xlabel('Lp (dB re 20\mu Pa)');
% ------------------------------------- Misc for Figure -----------------------------------------------------------

leg = legend(Leg, "Independant", 'correlated', 'NumColumns', 2); 
leg.Layout.Tile = 'north'; 



ylabel(tiled, 'Probability density');


saveas(gcf, 'figures/pressure_and_LP_distribution_part4Sf.eps', "epsc");


