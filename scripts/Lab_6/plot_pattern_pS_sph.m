
nrows = 1;
ncols = 1;
heightScale = 0.7; % Adjust height scaling if needed
Nleg = 1;

[columnwidth, ~] = get_widths();
height = get_height() * heightScale; 
fig = figure("Position", [0, 0, columnwidth*0.7, height], "Units", "points");
tiled = tiledlayout(nrows, ncols, "TileSpacing", "tight", "Padding", "loose");
corder = ["#181748","#810100","#575579","#7F7F7F"];


tile1 = nexttile;

for ii = 1: length(p0)
    Leg(ii)=polarplot(theta, SPL(p(:,ii)/p0(ii)), 'LineWidth', 1.5,Color=corder(ii)); hold on;
end

grid on;
pax = gca;
pax.ThetaZeroLocation = "top";
% rlim([-60, 0]);
% rticks(-60:6:0);
thetaticks(-180:15:180);
thetalim([-180 180]);

hold off;
title({'Normalized Pressure Level of a point source (dB SPL)',' in far field function of \theta (in °) for various ka'}); % might want to change title
leg1 = legend(Leg, 'ka=0.1', 'ka=1', 'ka=5', 'ka=10',NumColumns=2);
leg1.Layout.Tile = 'north';
hold off;

saveas(gcf,'figures/6plot_pattern_pS_sph.svg')

