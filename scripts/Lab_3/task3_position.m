room = [6.93, 4.78, 5.55]; % Arbitrary room dimensions
rS = [0, 0, 0];

num_val_r = 50;
r_values = linspace(0, room(3), num_val_r); % Generate r values from 0 to room(3)

[G_temp, ~] = green_func_room_lab3([0, 0, r_values(1)], rS, room, 'absorption', true);
Nfreq = length(G_temp); % Number of frequency bins
G_values = zeros(num_val_r, Nfreq); % Preallocate matrix for storing G values

for i = 1:num_steps
    r = r_values(i);
    rM = [0, 0, r];
    
    [G, ~] = green_func_room_lab3(rM, rS, room, 'absorption', true);
    G_values(i, :) = G;
end

% Plot results
figure;
imagesc(r_values, f, 20*log10(abs(G_values)/2e-5)); % Color plot of G vs frequency and r
set(gca, 'YDir', 'normal'); % Ensure frequency axis is in correct order
xlabel("r (height from floor)");
ylabel("Frequency (Hz)");
title("Green's Function Variation with r");
colorbar;
colormap jet;


%%

room = [6.93, 4.78, 5.55]; % Arbitrary room dimensions
rS = [0,0,0];
x = zeros(100);
y = zeros(100);
z = linspace(0,room(3),100);
rM = [x',y',z'];
%[idx,fm] = find_f_modes_lab3(1,1,1);
%[G,f] = green_func_room_lab3(rM,rS,room, 'absorption', true,'no_const',true, 'f', fm);
[G, ~] = green_func_room_lab3(rM, rS, room, 'absorption', true, 'no_const',true);


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
% xticks([0 max(x)])
% xticklabels( ["(0,0,0)", "(l_x,l_y,l_z)"])
% 
% ylabel(tiled,'Level in dB SPL');
% 
% % ------------------------------------- Misc for Figure -----------------------------------------------------------
% 
% leg = legend(Leg, 'f \approx f_m [1,1,1]', 'f \approx f_m [2,1,1]', 'NumColumns', 2); 
% leg.Layout.Tile = 'north'; 
% 
% 
% %---- Misc for Figure -----------------------------------------------------------
% 
% 
% % Save the figure in EPS format (modify file name)
% 
% saveas(gcf, 'figures/gf_through_room.eps', "epsc");


%% average pressure for 1 position

N = 10^2;
f = 700; % Hz
c = 343 ; % m/s
omega = 2*pi*f;
k = omega/c;
l_z = 7; % m

sum=0;
for i = 1:N
    sum=sum+exp(-1j*(k*rcos(theta_i+phi_i)));
end

p_avg_sqr = 1/N * abs(sum)^2;
