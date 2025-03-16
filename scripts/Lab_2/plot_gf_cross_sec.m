room =[3.14, 4.38, 3.27];
rS = [0,0,0];

nx=3;
ny=2;
nz=1;

x = linspace(0,room(1),100);
y = linspace(0,room(2),100);
z = ones(size(x)) * room(3);
rM = [x',y',z'];

[G,f] = green_func_room(rM,rS,room, 'absorption', true);%, 'mode', [nx,ny,nz]);
idx = find_f_modes(nx,ny,nz);
G_1 = abs(G(:,:,idx));


y2 = linspace(0,room(2),100); 
x2 = ones(size(y2))*room(1); 
z2 = linspace(0,room(3),100); 
rM = [x2',y2',z2']; 

[G,f] = green_func_room(rM,rS,room, 'absorption', true;%, 'mode', [nx,ny,nz]);
idx = find_f_modes(nx,ny,nz);
G_2 = abs(G(:,:,idx));

x3 = linspace(0,room(1),100);
z3 = linspace(0,room(3),100);
y3 = ones(size(x3))*room(2);
rM = [x3',y3',z3'];

[G,f] = green_func_room(rM,rS,room, 'absorption', true);%, 'mode', [nx,ny,nz]);
idx = find_f_modes(nx,ny,nz);
G_3 = abs(G(:,:,idx));


% Define the number of rows and columns for tiled layout
nrows = 1;
ncols = 3;
heightScale = 0.6; % Adjust height scaling if needed

% ------------------------------------------------------------------------------------------------------------
[columnwidth, ~] = get_widths();
height = get_height() * heightScale; 
fig = figure("Position", [0, 0, columnwidth*1.5, height], "Units", "points");
tiled = tiledlayout(nrows, ncols, "TileSpacing", "tight", "Padding", "loose");
corder = colororder;

% ------------------------------------- First Tile ------------------------------------------------------------

nexttile
% G1 = abs(G(:,:,idx(2))/max(G(:,:,idx(2)),[],'all'));
contourf(x,y,G_1)
title("z=0, Modes: n_x=3,n_y=2,n_z=1");


xlabel( 'X (m)');
ylabel( 'Y (m)');

% ------------------------------------- Second Tile ---------------------------------------------------------
ax2 = nexttile ;

contourf(y2,y2,G_2)

title("x=0, Modes: n_x=3,n_y=2,n_z=1");
xlabel( 'Y (m)');
ylabel( 'Z (m)');


% ------------------------------------- Third Tile ---------------------------------------------------------
nexttile ;

contourf(x3,z3,G_3)

title("y=0, Modes: n_x=3,n_y=2,n_z=1");
xlabel( 'X (m)');
ylabel( 'Z (m)');
colorbar()
%---- Misc for Figure -----------------------------------------------------------


% Save the figure in EPS format (modify file name)
saveas(gcf, 'figures/gf_cross_sec_321_whole_sum.eps', "epsc");


