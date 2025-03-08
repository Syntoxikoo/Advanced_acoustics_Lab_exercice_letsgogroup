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
Leg(1) = plot(nx*8,convf,"LineStyle",'-',"LineWidth", 1.0, "Color", corder(1,:)); 
% set(gca, 'YScale', 'log')
% xlim([])

grid on; hold on;

Leg(2) = plot(nx*8, convz,"LineStyle",'-',"LineWidth", 1.0, "Color", corder(2,:)); 
hold off
% Set axis limits 
% ylim([0, 1]) 
xlim([0, nx(end)*8]) 

% Add title 
% title('First Tile Title');

% Uncomment these if individual tile labels are preferred
xlabel('Number of modes (n \cdot m)');
ylabel('C (Pa^2)');

% ------------------------------------- Misc for Figure -----------------------------------------------------------

leg = legend('z fixed', 'f fixed', 'NumColumns', 2); 

saveas(gcf, 'figures/convergence.eps', "epsc");
