function []=plotRT(fc,RT,name)
% function to plot the sound power of two measurements. one regular and one
% at a wall
% fc = 1/3 octave band center frequencies
% P1 = first Power spectrum
% P2 = second power spectrum
% name = title of graph

% plotting sound power

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

% change measurement to only display 100 to 10kHz
idx = find (fc == 100);


% Create an index-based x-axis for equal spacing
x = 1:length(fc(idx:end)); % Generates equally spaced x-values

RT_x = RT(idx:end);

% Set X Y...
Leg(1) = plot(x, RT_x,"LineStyle",'-',"LineWidth", 1.0, "Color", corder(1,:),"DisplayName","average"); 

grid on; hold on;

%Leg(2) = plot(x, P2(idx:end),"LineStyle",'-',"LineWidth", 1.0, "Color", corder(2,:),"DisplayName","at wall"); 

% Adjust x-axis labels to show original frequency values
xticks(x);           % Set ticks to match x indices
xticklabels(fc(idx:end));   % Label them with the actual frequency values
xlim([x(1)-0.5, x(end)+0.5]);

% Set axis limits 
% ylim([0, 60]) 
%xlim([100, 1000]) 

% Add title 
title(name);

% Uncomment these if individual tile labels are preferred
xlabel('Frequency in Hz');
ylabel('T60 in s');

% ------------------------------------- Misc for Figure -----------------------------------------------------------

% leg = legend(Leg, 'Dataset 1', 'Dataset 2', 'NumColumns', 2); 
% leg.Layout.Tile = 'north'; 

% Add common X and Y axis labels for all tiles
% xlabel(tiled, 'Common X-axis Label');
% ylabel(tiled, 'Common Y-axis Label');

% Save the figure in EPS format (modify file name)
saveas(gcf, sprintf('figures/%s_SP.eps',name), "epsc");
end


