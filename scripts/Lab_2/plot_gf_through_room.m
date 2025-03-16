room =[3.14, 4.38, 3.27];
rS = [0,0,0];
x = linspace(0,room(1),100);
y = linspace(0,room(2),100);
z = linspace(0,room(3),100);
rM = [x',y',z'];
[idx,fm] = find_f_modes(1,1,1);
[G,f] = green_func_room(rM,rS,room, 'absorption', true,'no_const',true, 'f', fm);

[~,fm2] = find_f_modes(2,1,1);
fm = mean([fm,fm2]);
[G2,f2] = green_func_room(rM,rS,room, 'absorption', true,'no_const',true, 'f', fm);
% G_1 = abs(G(:,idx));


% Define the number of rows and columns for tiled layout
nrows = 1;
ncols = 1;
heightScale = 0.6; % Adjust height scaling if needed

% ------------------------------------------------------------------------------------------------------------
[columnwidth, ~] = get_widths();
height = get_height() * heightScale; 
fig = figure("Position", [0, 0, columnwidth, height], "Units", "points");
tiled = tiledlayout(nrows, ncols, "TileSpacing", "tight", "Padding", "loose");
corder = colororder;

% ------------------------------------- First Tile ------------------------------------------------------------

nexttile



Leg(1) = plot(x,20*log10(abs(G)/2e-5),"LineStyle",'-',"LineWidth", 1.0, "Color", corder(1,:)); 
grid("on"); hold on;
Leg(2) = plot(x,20*log10(abs(G2)/2e-5),"LineStyle",'--',"LineWidth", 1.0, "Color", corder(1,:)); 

hold off
xlabel( tiled,'Position');
xlim([0 max(x)])
xticks([0 max(x)])
xticklabels( ["(0,0,0)", "(l_x,l_y,l_z)"])

ylabel(tiled,'Level in dB SPL');

% ------------------------------------- Misc for Figure -----------------------------------------------------------

leg = legend(Leg, 'f \approx f_m [1,1,1]', 'f \approx f_m [2,1,1]', 'NumColumns', 2); 
leg.Layout.Tile = 'north'; 


%---- Misc for Figure -----------------------------------------------------------


% Save the figure in EPS format (modify file name)

saveas(gcf, 'figures/gf_through_room.eps', "epsc");


