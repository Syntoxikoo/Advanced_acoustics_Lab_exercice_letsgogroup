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
    semilogx(r/a, 20*log10(abs(p(:,ii))/2e-5), 'LineWidth', 1.5,"Color", corder(ii));hold on
end
grid on ; hold off


xlabel("r / a (m)")
ylabel("Normalized Pressure (dB re 20 \mu Pa)")
% ylim([40 80])
leg1 = legend('m = 0, Monopole', 'm = 1, Dipole', 'm = 2, Quadrupole',NumColumns=3);
leg1.Layout.Tile = 'north';
title("Pressure on axis in function of distance")
hold off
saveas(gcf,'figures/p_for_r_m123.svg')