clear;
close all;
clc;
addpath(genpath('scripts'))
addpath(genpath("functions"))
addpath figures/

%% PLOT VARYING FREQ AND VARYING Z POSITION
lx = 1; ly = 1;
z = linspace(0,1,10);
f = linspace(100,1000,100);

N_modes = [40,40];
xM = 1; yM = 1; % receiver
rS = [0.,0.,0]; %source

G = Gf_duct([xM,yM],rS, z, [lx,ly], f, N_modes,[1;0 ],1);

figure(1);
hold on;

for n = 1:size(G, 1)
% for n = 1:2
    plot(f, 20*log10(abs(G(n, :)) / 2e-5)+3*n, 'LineWidth', 2);
end

grid on;
xlabel('Frequency (Hz)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Pressure Level (dB SPL)', 'FontSize', 14, 'FontWeight', 'bold');

title(sprintf('Pressure = f(frequency), xM=%.2f, yM=%.2f, rS=[%.2f,%.2f,%.2f], N_{modes}=[%d,%d], lx=%.2f, ly=%.2f', ...
    xM, yM, rS(1), rS(2), rS(3), N_modes(1), N_modes(2), lx, ly), ...
    'FontSize', 16, 'FontWeight', 'bold');

legendStrings = arrayfun(@(n) sprintf('Position %d', n), 1:size(G, 1), 'UniformOutput', false);
legend(legendStrings, 'FontSize', 14, 'Location', 'best');

set(gca, 'FontSize', 14);
hold off;
