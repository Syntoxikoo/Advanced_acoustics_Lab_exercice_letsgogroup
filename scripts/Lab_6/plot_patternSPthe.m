
nrows = 1;
ncols = 1;
heightScale = 0.7; % Adjust height scaling if needed
Nleg = 1;

[columnwidth, ~] = get_widths();
height = get_height() * heightScale; 
fig = figure("Position", [0, 0, columnwidth*0.7, height], "Units", "points");
tiled = tiledlayout(nrows, ncols, "TileSpacing", "tight", "Padding", "loose");
corder = colororder;


tile1 = nexttile;


Leg(1)=polarplot(theta, SPL(p0), 'LineWidth', 1.5,Color=corder(1,:)); hold on;
Leg(2)=polarplot(theta, SPL(p1), 'LineWidth', 1.5,Color=corder(2,:)); 
Leg(3)=polarplot(theta, SPL(p2), 'LineWidth', 1.5,Color=corder(3,:)); 


grid on;
pax = gca;
pax.ThetaZeroLocation = "top";
% rlim([-60, 0]);
% rticks(-60:6:0);
thetaticks(-180:15:180);
thetalim([-180 180]);

hold off;
title('Normalized Pressure Level (dB SPL) in function of \theta (in °)'); % might want to change title
leg1 = legend(Leg, 'm = 0, Monopole', 'm = 1, Dipole', 'm = 2, Quadrupole',NumColumns=3);
leg1.Layout.Tile = 'north';
hold off;

saveas(gcf,'figures/6plot_patternSPthe.svg')