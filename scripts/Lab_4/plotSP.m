%% plotting sound power

% Define the number of rows and columns for tiled layout
nrows = 1;
ncols = 3;
heightScale = 0.5; % Adjust height scaling if needed

% ------------------------------------------------------------------------------------------------------------
[columnwidth, ~] = get_widths();
height = get_height() * heightScale; 
fig = figure("Position", [0, 0, columnwidth, height], "Units", "points");
tiled = tiledlayout(nrows, ncols, "TileSpacing", "tight", "Padding", "loose");
corder = colororder;

% ------------------------------------- First Tile ------------------------------------------------------------
nexttile

% Create an index-based x-axis for equal spacing
x = 1:length(fc); % Generates equally spaced x-values

% Set X Y...
Leg(1) = plot(x, P_highImp,"LineStyle",'-',"LineWidth", 1.0, "Color", corder(1,:),"DisplayName","average"); 

grid on; hold on;

Leg(2) = plot(x, P_highImp_wall,"LineStyle",'-',"LineWidth", 1.0, "Color", corder(2,:),"DisplayName","at wall"); 

% Set axis limits 
% ylim([0, 1]) 
% xlim([10, 1000]) 

% Adjust x-axis labels to show original frequency values
xticks(x);           % Set ticks to match x indices
xticklabels(fc);   % Label them with the actual frequency values

% Add title 
title('High impedance source');

% Uncomment these if individual tile labels are preferred
% xlabel('X-axis Label');
% ylabel('Y-axis Label');

% ------------------------------------- Second Tile ---------------------------------------------------------
nexttile 

% Set X Y...
Leg(1) = plot(x, P_lowImp,"LineStyle",'-',"LineWidth", 1.0, "Color", corder(1,:),"DisplayName","average"); 

grid on; hold on;

Leg(2) = plot(x, P_lowImp_wall,"LineStyle",'-',"LineWidth", 1.0, "Color", corder(2,:),"DisplayName","at wall"); 

% Set axis limits 
% ylim([0, 1]) 
% xlim([10, 1000]) 

% Adjust x-axis labels to show original frequency values
xticks(x);           % Set ticks to match x indices
xticklabels(fc);   % Label them with the actual frequency values

% Add title (optional)
title('Low impedance source');

% Uncomment these if individual tile labels are preferred
% xlabel('X-axis Label');
% ylabel('Y-axis Label');

% ------------------------------------- third Tile ---------------------------------------------------------
nexttile 

% Set X Y...
Leg(1) = plot(x, P_dipole,"LineStyle",'-',"LineWidth", 1.0, "Color", corder(1,:),"DisplayName","average"); 

grid on; hold on;

Leg(2) = plot(x, P_dipole_wall,"LineStyle",'-',"LineWidth", 1.0, "Color", corder(2,:),"DisplayName","at wall"); 

% Set axis limits 
% ylim([0, 1]) 
% xlim([10, 1000]) 

% Adjust x-axis labels to show original frequency values
xticks(x);           % Set ticks to match x indices
xticklabels(fc);   % Label them with the actual frequency values

% Add title (optional)
title('Dipole');

legend("location","southeast");

% Uncomment these if individual tile labels are preferred
% xlabel('X-axis Label');
% ylabel('Y-axis Label');

% ------------------------------------- Misc for Figure -----------------------------------------------------------

% leg = legend(Leg, 'Dataset 1', 'Dataset 2', 'NumColumns', 2); 
% leg.Layout.Tile = 'north'; 

% Add common X and Y axis labels for all tiles
xlabel(tiled, 'Common X-axis Label');
ylabel(tiled, 'Common Y-axis Label');

% Save the figure in EPS format (modify file name)
saveas(gcf, 'figures/namefig.eps', "epsc");


