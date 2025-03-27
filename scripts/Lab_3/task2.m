% Task 2
clear;
clc;


kr = linspace(0,8*pi,1000);

N=100; % Number of wawes
M=1000;  % Number of Monte carlo simulations for averaging

p1p2=zeros(M,length(kr));
phi=zeros(M,N);
msP=zeros(1,M);


for i=1:M
    [p1p2(i,:),phi(i,:)] = spatCorr(kr,N);
    msP(i) = meanSqP(N,phi);
end
    
avgmsP=mean(msP);
avgp1p2=mean(p1p2,1);

 y = avgp1p2 ./ avgmsP;


%% plotting

% Define the number of rows and columns for tiled layout
nrows = 1;
ncols = 1;
heightScale = 0.75; % Adjust height scaling if needed

% ------------------------------------------------------------------------------------------------------------
[columnwidth, ~] = get_widths();
height = get_height() * heightScale; 
fig = figure("Position", [0, 0, columnwidth, height], "Units", "points");
tiled = tiledlayout(nrows, ncols, "TileSpacing", "tight", "Padding", "loose");
corder = colororder;

% ------------------------------------- First Tile ------------------------------------------------------------
nexttile

% Set X Y...
Leg(1) = plot(kr, abs(y),"LineStyle",'-',"LineWidth", 1.0, "Color", corder(1,:)); 

grid on; hold on;

% Set x-axis ticks at multiples of π
xticks(0:pi:8*pi); 
xticklabels(compose('%g\\pi', 0:8)); % Auto-generate labels as 0π, 1π, ..., 8π


%Leg(2) = plot(X, Y,"LineStyle",'-',"LineWidth", 1.0, "Color", corder(2,:)); 

% Set axis limits 
% ylim([0, 1]) 
% xlim([10, 1000]) 

% Add title 
%title('First Tile Title');

% Uncomment these if individual tile labels are preferred
% xlabel('X-axis Label');
% ylabel('Y-axis Label');

% ------------------------------------- Misc for Figure -----------------------------------------------------------

% leg = legend(Leg, 'Dataset 1', 'Dataset 2', 'NumColumns', 2); 
% leg.Layout.Tile = 'north'; 

% Add common X and Y axis labels for all tiles
xlabel(tiled, 'normalized distance between points, kr');
ylabel(tiled, 'normalized spatial covariance');

% Save the figure in EPS format (modify file name)
saveas(gcf, 'figures/normSpatCov.eps', "epsc");