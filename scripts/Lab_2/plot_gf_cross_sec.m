room =[3.14, 4.38, 3.27];
rS = [0,0,0];
x = linspace(0,room(1),100);
y = linspace(0,room(2),100);
z = ones(size(x)) * room(3);
rM = [x',y',z'];
[G,f] = green_func_room(rM,rS,room, 'absorption', true);
idx = find_f_modes(1,1,1);
G_1 = abs(G(:,:,idx));



y2 = linspace(0,room(2),100);
z2 = linspace(0,room(3),100);
x2 = ones(size(y)) * room(1);
rM = [x2',y2',z2'];
[G,f] = green_func_room(rM,rS,room, 'absorption', true);
idx = find_f_modes(1,1,1);
G_2 = abs(G(:,:,idx));

% Define the number of rows and columns for tiled layout
nrows = 1;
ncols = 2;
heightScale = 0.5; % Adjust height scaling if needed

% ------------------------------------------------------------------------------------------------------------
[columnwidth, ~] = get_widths();
height = get_height() * heightScale; 
fig = figure("Position", [0, 0, columnwidth, height], "Units", "points");
tiled = tiledlayout(nrows, ncols, "TileSpacing", "tight", "Padding", "loose");
corder = colororder;

% ------------------------------------- First Tile ------------------------------------------------------------

nexttile
% G1 = abs(G(:,:,idx(2))/max(G(:,:,idx(2)),[],'all'));
contourf(x,y,G_1)

title("z=0, Modes: n_x=1,n_y=1,n_z=1");

xlabel( 'X (m)');
ylabel( 'Y (m)');
% ------------------------------------- Second Tile ---------------------------------------------------------
ax2 = nexttile ;
% G2 = abs(G(:,:,idx(3))/max(G(:,:,idx(3)),[],'all'));
contourf(z2,y2,G_2)

title("x=0, Modes: n_x=1,n_y=1,n_z=1");
xlabel( 'Z (m)');
ylabel( 'Y (m)');


%---- Misc for Figure -----------------------------------------------------------


% Save the figure in EPS format (modify file name)
saveas(gcf, 'figures/gf_cross_sec_111.eps', "epsc");


