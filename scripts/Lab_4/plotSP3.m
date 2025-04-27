function []=plotSP3(fc,P1,name,P2,legend_list,varargin)
% function to plot the sound power of two measurements. one regular and one
% at a wall
% fc = 1/3 octave band center frequencies
% P1 = first Power spectrum
% P2 = second power spectrum
% name = title of graph

p = inputParser;
addRequired(p, 'fc');
addRequired(p, 'P1');
addOptional(p, 'P2', NaN(1,24));
addRequired(p, 'name');
addOptional(p, 'legend_list', ["P1","P2", "P3"]);
addOptional(p, 'ylimit', [50 87]);
addParameter(p, 'PI1', NaN(1,24));
addParameter(p, 'PI2', NaN(1,24));


parse(p, fc, P1, name,P2, varargin{:});

PI1  = p.Results.PI1;
PI2  = p.Results.PI2;
ylimit = p.Results.ylimit;

% Define the number of rows and columns for tiled layout
nrows = 2;
ncols = 1;
heightScale = 1.0; % Adjust height scaling if needed
Nleg = 1;

% ------------------------------------------------------------------------------------------------------------
[columnwidth, ~] = get_widths();
height = get_height() * heightScale; 
fig = figure("Position", [0, 0, columnwidth, height], "Units", "points");
tiled = tiledlayout(nrows, ncols, "TileSpacing", "tight", "Padding", "loose");
corder = colororder;

% ------------------------------------- First Tile ------------------------------------------------------------
ax1 = nexttile;

idx = find (fc == 100);

x = 1:length(fc(idx:end));  %space equally

leg(1) = plot(x, P1(idx:end),"LineStyle",'-',"LineWidth", 1.0,"marker",".","MarkerSize",10, "Color", corder(1,:)); 


grid on; hold on;
if ~any(isnan(P2))
    Nleg = Nleg+1;
    leg(2) = plot(x, P2(idx:end),"LineStyle",'-',"LineWidth", 1.0,"marker",".","MarkerSize",10, "Color", corder(2,:)); 
end


xticks(x);     
xticklabels(" ");
xlim([x(1) x(end)])
ylim(ylimit)



% Uncomment these if individual tile labels are preferred
% xlabel('Frequency (Hz)');
ylabel('Sound Power (dB re 1pW)');
title("Sound power from Intensity")


nexttile
ax2 = gca;
idx = find (fc == 100);

x = 1:length(fc(idx:end));  %space equally


grid on; hold on;
if ~any(isnan(PI1))
    plot(x, PI1(idx:end),"LineStyle",'-',"LineWidth", 1.0,"marker",".","MarkerSize",10, "Color", corder(1,:)); 
end

if ~any(isnan(PI2))
    plot(x, PI2(idx:end),"LineStyle",'-',"LineWidth", 1.0,"marker",".","MarkerSize",10, "Color", corder(2,:)); 
end

xticks(x);     
xticklabels(fc(idx:end));
xlim([x(1) x(end)])
% ylim(ylimit)



% Uncomment these if individual tile labels are preferred
xlabel('Frequency (Hz)');
ylabel('Pressure-Intensity index');
title("Pressure-Intensity Index")

% ------------------------------------- Misc for Figure -----------------------------------------------------------
linkaxes([ax1, ax2], 'x');
leg = legend(leg, legend_list(1:Nleg), 'NumColumns', Nleg); 
leg.Layout.Tile = 'north'; 

% Save the figure in EPS format (modify file name)
saveas(gcf, sprintf('figures/%s_SP.eps',name), "epsc");
end


