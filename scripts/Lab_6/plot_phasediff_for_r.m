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
ylim([-pi pi])
yticks([-pi, -pi/2, 0, pi/2, pi])
yticklabels(["-\pi","-\pi/2","0","\pi/2", "\pi"])

leg1 = legend('m = 0, Monopole', 'm = 1, Dipole', 'm = 2, Quadrupole',NumColumns=3);
leg1.Layout.Tile = 'north';
title("Pressure on axis in function of distance")
hold off
saveas(gcf,'figures/phasediff_for_r_m123.svg')