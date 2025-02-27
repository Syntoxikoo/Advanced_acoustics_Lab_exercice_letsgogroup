clear;
close all;
clc;
addpath(genpath('scripts'))
addpath(genpath("functions"))
addpath figures/

%% PLOT VARYING FREQ AND VARYING Z POSITION
%VARY F
c = 343;

lx = 0.7; ly = 1;
if lx >= ly
    lowest_fc = c / (2 * lx);
else
    lowest_fc = c / (2 * ly);
end
N_modes = [4,4];
xM = lx; yM = ly; % receiver
rS = [0,0,0]; % source
fmn = compute_modes(lx, ly, N_modes);

z = linspace(5,0.01,2);
f = linspace(1,1000,1000-1);

G = Gf_duct([xM, yM], rS, z, [lx, ly], f, N_modes, [], 1);

figure(1);
clf;
hold on;

for n = 1:size(G, 1)
    plot(f, 20*log10(abs(G(n, :)) / 2e-5), 'LineWidth', 2);
end

%add vertical dotted lines at frequencies closest to fmn
for k = 1:numel(fmn)
    [~, idx] = min(abs(f - fmn(k))); % Find closest frequency in f
    xline(f(idx), 'k:', 'LineWidth', 1.5); % Black dotted line
end

grid on;
xlabel('Frequency (Hz)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Pressure Level (dB SPL)', 'FontSize', 14, 'FontWeight', 'bold');

title(sprintf('Pressure = f(frequency), rR(x,y,z)=[%.2f, %.2f, zR], rS(x,y,z)=[%.2f,%.2f,%.2f], N_{modes}(x,y)=[%d,%d], cross-section=[lx=%.2f x ly=%.2f]', ...
    xM, yM, rS(1), rS(2), rS(3), N_modes(1), N_modes(2), lx, ly), ...
    'FontSize', 16, 'FontWeight', 'bold');

legendStrings = arrayfun(@(n, zVal) sprintf('z = %.2f', zVal), 1:size(G, 1), z, 'UniformOutput', false);
legend(legendStrings, 'FontSize', 14, 'Location', 'best');

set(gca, 'FontSize', 14);
hold off;

%
% VARY Z
z = linspace(0,10,1000-1) ;
f = linspace(lowest_fc-20,10*lowest_fc,2) ;

G = Gf_duct([xM,yM],rS, z, [lx,ly], f, N_modes,[],1);
G=G';

figure(2);
clf;
hold on;

for n = 1:size(G, 1)
    plot(z, 20*log10(abs(G(n, :)) / 2e-5), 'LineWidth', 2);
end

grid on;
xlabel('z (m)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Pressure Level (dB SPL)', 'FontSize', 14, 'FontWeight', 'bold');

title(sprintf('Pressure = f(z), rR(x,y,z)=[%.2f, %.2f, zR], rS(x,y,z)=[%.2f,%.2f,%.2f], N_{modes}(x,y)=[%d,%d], cross-section=[lx=%.2f x ly=%.2f]', ...
    xM, yM, rS(1), rS(2), rS(3), N_modes(1), N_modes(2), lx, ly), ...
    'FontSize', 16, 'FontWeight', 'bold');

legendStrings = arrayfun(@(n, zVal) sprintf('f = %.2f', zVal), 1:size(G, 1), f, 'UniformOutput', false);
legend(legendStrings, 'FontSize', 14, 'Location', 'best');

set(gca, 'FontSize', 14);
hold off;

%% add the 3D distance/freq plot + the 3D surface plot

z = linspace(0,5,2500-1) ;
f = linspace(1,1000,2500-1) ;
N_modes = [10,10];

G = Gf_duct([xM,yM],rS, z, [lx,ly], f, N_modes,[],1);

figure(3);
clf;
hold on;

% Plot using imagesc
imagesc(z, f, 20*log10(abs(G') / 2e-5)); 
set(gca, 'YDir', 'normal'); % Ensure correct orientation

% Labels and title
xlabel('z (m)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Frequency (Hz)', 'FontSize', 14, 'FontWeight', 'bold');
xlim([min(z) max(z)]);
ylim([min(f) max(f)]);

% Colorbar for reference
c = colorbar;
colormap('default'); % Use a colormap for better visualization
c.Label.String = 'Level (dB SPL)';  % Add label to colorbar
c.Label.FontSize = 14;
c.Label.FontWeight = 'bold';

title(sprintf('Pressure = f(frequency,z), rR(x,y,z)=[%.2f, %.2f, zR], rS(x,y,z)=[%.2f,%.2f,%.2f], N_{modes}(x,y)=[%d,%d], cross-section=[lx=%.2f x ly=%.2f]', ...
    xM, yM, rS(1), rS(2), rS(3), N_modes(1), N_modes(2), lx, ly), ...
    'FontSize', 16, 'FontWeight', 'bold');

set(gca, 'FontSize', 14);
hold off;
