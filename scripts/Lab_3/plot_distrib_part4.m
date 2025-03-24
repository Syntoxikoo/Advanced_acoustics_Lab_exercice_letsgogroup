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
x = linspace(0,max(X_11),100);
pd = fitdist(X_11, 'Beta');
y = pdf(pd, x);
plot(x, y, 'LineWidth', 2)

Leg(2) = histogram(X_1, bins,'Normalization', 'pdf','FaceAlpha',0.4); hold on;
% x = linspace(0,max(X_1),100);
% pd = fitdist(X_1, 'Beta');
% y = pdf(pd, x);
% plot(x, y, 'LineWidth', 2)


% Add title 
title('Distrib pressure rms');

% Uncomment these if individual tile labels are preferred
% xlabel('X-axis Label');
% ylabel('Y-axis Label');

% ------------------------------------- Second Tile ---------------------------------------------------------
nexttile 

bins = 30;
Leg(1) = histogram(X_21, bins,'Normalization', 'pdf','FaceAlpha',0.4); hold on;
x = linspace(0,max(X_21),100);
pd = fitdist(X_21, 'Normal');
y = pdf(pd, x);
plot(x, y, 'LineWidth', 2)

Leg(2) = histogram(X_2, bins,'Normalization', 'pdf','FaceAlpha',0.4); hold on;
% x = linspace(0,max(X_2),100);
% pd = fitdist(X_2, 'Beta');
% y = pdf(pd, x);
% plot(x, y, 'LineWidth', 2)


% Add title 
title('Distrib pressure level');

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


