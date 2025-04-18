%% Task 2
clear;
clc

krMAX= 8*pi; % maximum of kr for simulations

N = [1 10 100]; % Number of Waves
M = [10 100 1000]; % Number of Monte carlo simulations for averaging

kr = linspace(0,krMAX,1000); % kr Vector

% initializing matrices
yM=zeros(length(M),length(kr)); 
yN=zeros(length(N),length(kr));

% Calculate for fixed Number of waves and different number of MC-Simulations
for m=1:length(M) 

    % initializing matrices
    p1p2=zeros(M(m),length(kr)); 
    phi=zeros(M(m),N(length(N)));
    msP=zeros(1,M(m));

    % calculate array of MC-Simulations
    for i=1:M(m)
        [p1p2(i,:),phi(i,:)] = spatCorr(kr,N(length(N)));
        msP(i) = meanSqP(N(length(N)),phi(i,:));
    end
    
    % average over number of MC-Simulations
    avgmsP=mean(msP);
    avgp1p2=mean(p1p2,1);
    
    % Calculate normalized spatial correlation
    yM(m,:) = avgp1p2 ./ avgmsP;
end

% Calculate for fixed Number of MC-Simulations  and different number of waves
for n=1:length(N) 

    % initializing matrices
    p1p2=zeros(M(length(M)),length(kr));
    phi=zeros(M(length(M)),N(n));
    msP=zeros(1,M(length(M)));

    % calculate array of MC-Simulations
    for i=1:M(length(M))
        [p1p2(i,:),phi(i,:)] = spatCorr(kr,N(n));
        msP(i) = meanSqP(N(n),phi(i,:));
    end
    
    % average over number of MC-Simulations
    avgmsP=mean(msP);
    avgp1p2=mean(p1p2,1);
   
    % Calculate normalized spatial correlation
    yN(n,:) = avgp1p2 ./ avgmsP;
end

%% calculate the sinc function
sinc=sin(kr)./kr;

%% plotting

% Define the number of rows and columns for tiled layout
nrows = 1;
ncols = 2;
heightScale = 0.75; % Adjust height scaling if needed

% ------------------------------------------------------------------------------------------------------------
[columnwidth, ~] = get_widths();
columnwidth=800;
height = get_height() * heightScale; 
fig = figure("Position", [0, 0, columnwidth, height], "Units", "points");
tiled = tiledlayout(nrows, ncols, "TileSpacing", "tight", "Padding", "loose");
corder = colororder;

% ------------------------------------- First Tile ------------------------------------------------------------
nexttile

% Set X Y...
Leg(1) = plot(kr, sinc,"LineStyle",'-',"LineWidth", 1.0, "Color", corder(1,:),"DisplayName",'sin(kr)/kr'); 

hold on;
    for l=1:(length(M))
        Leg(l) = plot(kr, real(yM(l,:)),"LineStyle",'--',"LineWidth", 1.0, "Color", corder(l+1,:),"DisplayName",sprintf('M = %d', M(l)));
    end
hold off;

grid on; hold on;

% Set x-axis ticks at multiples of π
xticks(0:pi:8*pi); 
xticklabels(compose('%g\\pi', 0:8)); % Auto-generate labels as 0π, 1π, ..., 8π

legend;

% Set axis limits 
% ylim([0, 1]) 
% xlim([10, 1000]) 

% Add title 
title(sprintf('N = %d',N(length(N))));

% Uncomment these if individual tile labels are preferred
% xlabel('X-axis Label');
% ylabel('Y-axis Label');

% ------------------------------------- Second Tile ------------------------------------------------------------
nexttile

% Set X Y...
Leg(1) = plot(kr, sinc,"LineStyle",'-',"LineWidth", 1.0, "Color", corder(1,:),"DisplayName","sin(kr)/kr"); 


hold on;
    for l=1:(length(N))
        Leg(l) = plot(kr, real(yN(l,:)),"LineStyle",'--',"LineWidth", 1.0, "Color", corder(l+1,:),"DisplayName",sprintf('N = %d', N(l)));
    end
hold off;

grid on; hold on;

% Set x-axis ticks at multiples of π
xticks(0:pi:8*pi); 
xticklabels(compose('%g\\pi', 0:8)); % Auto-generate labels as 0π, 1π, ..., 8π

legend;

% Set axis limits 
% ylim([0, 1]) 
% xlim([10, 1000]) 

% Add title 
title(sprintf('M = %d',M(length(M))));

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
