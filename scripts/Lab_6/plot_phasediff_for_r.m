nrows = 1;
ncols = 1;
heightScale = 0.55; % Adjust height scaling if needed
Nleg = 1;

[columnwidth, ~] = get_widths();
height = get_height() * heightScale; 
fig = figure("Position", [0, 0, columnwidth, height], "Units", "points");
tiled = tiledlayout(nrows, ncols, "TileSpacing", "tight", "Padding", "loose");
corder = ["#181748","#810100","#575579","#7F7F7F"];

tile1 = nexttile;
for ii = 1 : length(m)
    semilogx(r/a, X(:,ii), 'LineWidth', 1.5,"Color", corder(ii));hold on
end
grid on ; hold off


xlabel("r / a (m)")
ylabel("Phase angle (rad)")
ylim([0 pi/2])
yticks([0 pi/4, pi/2])
yticklabels(["0","\pi/4","\pi/2"])

leg1 = legend('m = 0, Monopole', 'm = 1, Dipole', 'm = 2, Quadrupole',NumColumns=3);
leg1.Layout.Tile = 'north';
title("Phase angle between sound pressure and radial component of the particle velocity")
hold off
saveas(gcf,'figures/phasediff_for_r_m123.svg')