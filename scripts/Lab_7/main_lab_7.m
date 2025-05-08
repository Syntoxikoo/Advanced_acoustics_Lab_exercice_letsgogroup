%% Wavenumber transform, Lab 7

% This script solves tasks A to E from the lab exercise
% It is divided into sections and contains result discussion in comments
% Authors: s243674, s243592, s243590, s242987

clear; close all; clc;

%% Constants

c = 343;            % Speed of sound (m/s)
rho = 1.21;           % Air density (kg/m^3)
Q = 1e-4;           % Monopole volume velocity (m^3/s)
L = 1;               % Mic array size in x and y axis (m)
Lx=L;
Ly=L;
M = 15;            % Number of microphones per x/y dimension
N = 512;              % Number of lines per dimension of the FFTs (for zero padding)
f = 1000;             % Frequency (Hz)
omega = 2 * pi * f;   
k = omega / c;        

%% Define microphone grid

x = linspace(-Lx/2, Lx/2, M);   % x-Mics creations
y = linspace(-Ly/2, Ly/2, M);   % y-Mics creations
[X, Y] = meshgrid(x, y);
dx = x(2) - x(1);             % space between each mics
dy = y(2) - y(1);             % dx and dy are both computed but dx=dy as Lx=Ly

figure(1);
plot(X, Y, 'ko');
title('Microphone Array Mesh');
xlim([-0.7*L 0.7*L]); ylim([-0.7*L 0.7*L]);   % allows to see mics nicely
xlabel('x (m)'); ylabel('y (m)');
grid on; axis equal;

%% Define wavenumber grid

kx = (-N/2:N/2-1) * (2 * pi / (N * dx));    % Using N to make size match P_k size
ky = (-N/2:N/2-1) * (2 * pi / (N * dy));
kx = -kx;                  % Correct FFT sign convention???
ky = -ky;                  % Correct FFT sign convention???
[KX, KY] = meshgrid(kx, ky);









%% TASK A - Source at (0,0,-0.1)
z_source = -0.1;                                                 % Normal (z) distance of the monopole to the array plane (x,y) 
R = sqrt(X.^2 + Y.^2 + (0 - z_source).^2);                      % Distance from the monopole to each microphone position (15x15)
p = 1j * omega * rho * Q ./ (4 * pi * R) .* exp(-1j * k * R);   % Pressure at each measurment microphone of the array (15x15) (eq.9.2 from the book)

% SPL calculation
p_rms = abs(p)/ sqrt(2);
SPL = 20 * log10(p_rms / 2e-5);

x0 =0; y0 = 0;     

figure(2);
hold on;
plot3(X(:), Y(:), zeros(size(X(:))), 'bo', 'MarkerFaceColor', 'b');
plot3(x0, y0, -0.1, 'rp', 'MarkerSize', 15, 'MarkerFaceColor', 'r');
title('Task C: Microphone Array and Source Position (3D view)');
xlabel('x [m]'); ylabel('y [m]'); 
zlabel('z [m]');
grid on;
legend('Microphones','Source');
view(3); axis equal;

figure(3);
pcolor(X, Y, SPL); shading flat; colorbar;
title('Task A: SPL (dB SPL) for monopole at (0,0,-0.1)');
xlabel('x (m)'); ylabel('y (m)');

figure(4);
pcolor(X, Y, angle(p)); shading flat; colorbar;
title('Task A: Phase (rad) of pressure for monopole at (0,0,-0.1)');
xlabel('x (m)'); ylabel('y (m)');

% ANALYSIS : The sound pressure level figure 1 shows an expected maxima of
% pressure for the microphone which distance to the monopole is the
% smallest (the microphone at coordonates (0,0,0)). The sound pressure
% level then decreases with R, so in each x-y direction to the center 
% microphone similarly, creating a radially symetrical figure. Those points
% being further and further away from the center microphone also receive a
% wave that travelled during more time, therefore the phase also follows a
% sinusoidal behavior in all radial directions to the center microphone. 





%% TASK B - Wavenumber spectrum at 1kHz and 500Hz

% Same methodology but for two different frequencies
f = 1000;
omega = 2 * pi * f;
k = omega / c;
R = sqrt(X.^2 + Y.^2 + (0 - z_source).^2);
p = 1j * omega * rho * Q ./ (4 * pi * R) .* exp(-1j * k * R);
P_k = fftshift(fft2(p, N, N));
kz = sqrt(k^2 - KX.^2 - KY.^2);


figure(5);
pcolor(KX, KY, abs(P_k)); 
shading flat; 
colorbar;
title('Task B: Wavenumber spectrum (rad/m) at 1kHz');
xlabel('k_x (rad/m)'); ylabel('k_y (rad/m)');

figure(6);
pcolor(KX, KY, angle(P_k)); 
shading flat; 
colorbar;
title('Task B: Wavenumber spectrum phase (rad) at 1kHz');
xlabel('k_x (rad/m)'); ylabel('k_y (rad/m)');

figure(7);
pcolor(KX, KY, real(kz)); shading flat; colorbar;
title('Task B: Real part of k_z at 1kHz');
xlabel('k_x [rad/m]'); ylabel('k_y [rad/m]');

figure(8);
pcolor(KX, KY, imag(kz)); shading flat; colorbar;
title('Task B: Imag part of k_z at 1kHz');
xlabel('k_x [rad/m]'); ylabel('k_y [rad/m]');

% Changing f=100 to f500
f = 500;
omega = 2 * pi * f;
k = omega / c;
R = sqrt(X.^2 + Y.^2 + (0 - z_source).^2);
p = 1j * omega * rho * Q ./ (4 * pi * R) .* exp(-1j * k * R);
P_k = fftshift(fft2(p, N, N));
kz = sqrt(k^2 - KX.^2 - KY.^2);

figure(9);
pcolor(KX, KY, abs(P_k)); 
shading flat; 
colorbar;
title('Task B: Wavenumber spectrum (rad/m) at 500Hz');
xlabel('k_x (rad/m)'); ylabel('k_y (rad/m)');

figure(10);
pcolor(KX, KY, angle(P_k)); 
shading flat; 
colorbar;
title('Task B: Wavenumber spectrum phase (rad) at 500Hz');
xlabel('k_x (rad/m)'); ylabel('k_y (rad/m)');

figure(11);
pcolor(KX, KY, real(kz)); shading flat; colorbar;
title('Task B: Real part of k_z at 500 Hz');
xlabel('k_x [rad/m]'); ylabel('k_y [rad/m]');

figure(12);
pcolor(KX, KY, imag(kz)); shading flat; colorbar;
title('Task B: Imag part of k_z at 500 Hz');
xlabel('k_x [rad/m]'); ylabel('k_y [rad/m]');

% HAS TO BE ANALYSED





%% TASK C - Source at corner (L/2, L/2, -0.1)

f = 1000;
omega = 2 * pi * f;
k = omega / c;
x0 = L/2; y0 = L/2;                                       % New source coordinates
R = sqrt((X - x0).^2 + (Y - y0).^2 + (0 - z_source).^2);
p = 1j * omega * rho * Q ./ (4 * pi * R) .* exp(-1j * k * R);
P_k = fftshift(fft2(p, N, N));
kz = sqrt(k^2 - KX.^2 - KY.^2);

figure(13);
hold on;
plot3(X(:), Y(:), zeros(size(X(:))), 'bo', 'MarkerFaceColor', 'b');
plot3(x0, y0, -0.1, 'rp', 'MarkerSize', 15, 'MarkerFaceColor', 'r');
title('Task C: Microphone Array and Source Position (3D view)');
xlabel('x [m]'); ylabel('y [m]'); 
zlabel('z [m]');
grid on;
legend('Microphones','Source');
view(3); axis equal;

figure(14);
pcolor(KX, KY, abs(P_k)); 
shading flat; 
colorbar;
title('Task C: Wavenumber magnitude spectrum at corner (L/2, L/2)');
xlabel('k_x [rad/m]'); ylabel('k_y [rad/m]');

figure(15);
pcolor(KX, KY, angle(P_k)); 
shading flat; 
colorbar;
title('Task C: Wavenumber phase spectrum at corner (L/2, L/2)');
xlabel('k_x [rad/m]'); ylabel('k_y [rad/m]');

figure(16);
pcolor(KX, KY, real(kz)); shading flat; colorbar;
title('Task B: Real part of k_zz for cornered source');
xlabel('k_x [rad/m]'); ylabel('k_y [rad/m]');

figure(17);
pcolor(KX, KY, imag(kz)); shading flat; colorbar;
title('Task B: Imag part of k_z for cornered source');
xlabel('k_x [rad/m]'); ylabel('k_y [rad/m]');

% HAS TO BE ANALYSED

