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
view(3);
zlim([-0.2, 0.2])
%  axis equal;

figure(3);
pcolor(X, Y, SPL); shading flat; colorbar;
title('Task A: SPL (dB SPL) for monopole at (0,0,-0.1)');
xlabel('x (m)'); ylabel('y (m)');

figure(4);
pcolor(X, Y, angle(p)); shading flat; colorbar;
title('Task A: Phase (rad) of pressure for monopole at (0,0,-0.1)');
xlabel('x (m)'); ylabel('y (m)');

% ANALYSIS : The sound pressure level figure 1 shows an expected maxima of
% pressure for the microphone which has the smallest distance to the monopole
% (the microphone at coordonates (0,0,0)). The sound pressure
% level then decreases with R, so in each x-y direction to the center 
% microphone similarly, creating a radially symetrical figure. Those points
% being further and further away from the center microphone also receive a
% wave that travelled during more time, therefore the phase also follows a
% sinusoidal behavior in all radial directions to the center microphone.
% completing a full wavelength at the corner of the array.

%% TASK B - Wavenumber spectrum at 1kHz and 500Hz

% Same methodology but for two different frequencies
f = 1000;
omega = 2 * pi * f;
k = omega / c;
R = sqrt(X.^2 + Y.^2 + (0 - z_source).^2);
p = 1j * omega * rho * Q ./ (4 * pi * R) .* exp(-1j * k * R);
P_k = fftshift(fft2(p, N, N));
kz = sqrt(k^2 - KX.^2 - KY.^2);


fig5 = figure(5);
set(fig5,"position",[0,0,1000,400]);
subplot(1,2,1);
pcolor(KX, KY, abs(P_k)); 
shading flat; 
colorbar;
title('Task B: Wavenumber spectrum (rad/m) at 1kHz');
xlabel('k_x (rad/m)'); ylabel('k_y (rad/m)');

fig6 = figure(6);
set(fig6,"position",[0,0,1000,400]);
subplot(1,2,1)
pcolor(KX, KY, angle(P_k)); 
shading flat; 
colorbar;
title('Task B: Wavenumber spectrum phase (rad) at 1kHz');
xlabel('k_x (rad/m)'); ylabel('k_y (rad/m)');

fig7 = figure(7);
set(fig7,"position",[0,0,1000,400]);
subplot(1,2,1)
pcolor(KX, KY, real(kz)); shading flat; colorbar;
title('Task B: Real part of k_z at 1kHz');
xlabel('k_x [rad/m]'); ylabel('k_y [rad/m]');


fig8 = figure(8);
set(fig8,"position",[0,0,1000,400]);
subplot(1,2,1)
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

figure(5);
subplot(1,2,2);
pcolor(KX, KY, abs(P_k)); 
shading flat; 
colorbar;
title('Task B: Wavenumber spectrum of pressure at 500Hz');
xlabel('k_x (rad/m)'); ylabel('k_y (rad/m)');

figure(6);
subplot(1,2,2)
pcolor(KX, KY, angle(P_k)); 
shading flat; 
colorbar;
title('Task B: Wavenumber phase spectrum of pressure at 500Hz');
xlabel('k_x (rad/m)'); ylabel('k_y (rad/m)');

figure(7);
subplot(1,2,2)
pcolor(KX, KY, real(kz)); shading flat; colorbar;
title('Task B: Real part of k_z at 500 Hz');
xlabel('k_x [rad/m]'); ylabel('k_y [rad/m]');

figure(8);
subplot(1,2,2)
pcolor(KX, KY, imag(kz)); shading flat; colorbar;
title('Task B: Imag part of k_z at 500 Hz');
xlabel('k_x [rad/m]'); ylabel('k_y [rad/m]');

% B Analysis:
% Since the microphone array is quite close to the source we are able to
% observe the evanescent waves as well as the traveling waves. Since the
% source is a point source positioned in the middle of the microphone array (x,y) % plane, we can clearly see a circle-like 
% shape of the wavenumber spectrum.
% The dot with high pressure directly in the middle represents the plane wave traveling normal to the surface we are observing. 
% The circle around this represents the spherical behavior from the source propagating waves as all the waves travel 
% from the center towards the ends of our array. It can be seen that the
% circle has a larger radius for the 1 kHz signal compared to the 500 Hz
% signal because the larger wavelength of 500 Hz results in a smaller value
% for the wavenumbers since k^2=kx^2+ky^2+kz^2. This equation illustrates the
% link between the radius of the seen circle to the frequency of the
% signal. This can even better be observed in the plot of the real and
% imaginary part of the wave spectrum. The imaginary part of kz being zero inside
% the circle and growing outside the circle. This shows that outside of the
% radius kx+ky>k indicating evanescent waves and inside kx+ky<k indicating
% propagating waves. The periodic artefacts seen towards larger values of kx and ky are results of the truncation of the sum of modes. 
% Increasing the number of microphone withing the grid and therefore reducing the spacing between each one of them would increase the 
% wave number resolution, an allow to capture the full characteristic of the sound field at higher frequency and with more precision
% on the specific direction of propagation of each wave.




%% TASK C - Source at corner (L/2, L/2, -0.1)

f = 1000;
omega = 2 * pi * f;
k = omega / c;
x0 = L/2; y0 = L/2;                                       % New source coordinates
R = sqrt((X - x0).^2 + (Y - y0).^2 + (0 - z_source).^2);
p = 1j * omega * rho * Q ./ (4 * pi * R) .* exp(-1j * k * R);
P_k = fftshift(fft2(p, N, N));
kz = sqrt(k^2 - KX.^2 - KY.^2);

figure(9);
hold on;
plot3(X(:), Y(:), zeros(size(X(:))), 'bo', 'MarkerFaceColor', 'b');
plot3(x0, y0, -0.1, 'rp', 'MarkerSize', 15, 'MarkerFaceColor', 'r');
title('Task C: Microphone Array and Source Position (3D view)');
xlabel('x [m]'); ylabel('y [m]'); 
zlabel('z [m]');
grid on;
legend('Microphones','Source');
view(3);

fig10 = figure(10);
set(fig10,"position",[0,0,1000,400]);
subplot(1,2,1)
pcolor(KX, KY, abs(P_k)); 
shading flat; 
colorbar;
title('Task C: Wavenumber magnitude spectrum of pressure at corner (L/2, L/2)');
xlabel('k_x [rad/m]'); ylabel('k_y [rad/m]');

figure(10);
subplot(1,2,2)
pcolor(KX, KY, angle(P_k)); 
shading flat; 
colorbar;
title('Task C: Wavenumber phase spectrum of pressure at corner (L/2, L/2)');
xlabel('k_x [rad/m]'); ylabel('k_y [rad/m]');


% C Analysis:
% Now in the case of the source being positioned in a corner of our array,
% all the plane waves are traveling with an initial angle to the array, and 
% therefore travel towards the negative kx and ky. The wavenumber 
% specturm can here be analysed like a vector pointing towards the direction of 
% propagation of the wavefront comming from the source.
% Compared to the first wavenumber spectrum in figure 5, it is just one
% corner of the circle without the center dot. Showing a similar radius
% from 0 but only propagating in the negative x/y direction.
% Again, the artifacts which look like a repetition of the shape are due to
% the truncation error. 



%% TASK D - Source at (0,0,-3)
x0 = 0; y0 = 0;  
z_source = -3;
f = 1000;
omega = 2 * pi * f;
k = omega / c;
R = sqrt(X.^2 + Y.^2 + (0 - z_source).^2);
p = 1j * omega * rho * Q ./ (4 * pi * R) .* exp(-1j * k * R);

figure(12);
hold on;
plot3(X(:), Y(:), zeros(size(X(:))), 'bo', 'MarkerFaceColor', 'b');
plot3(x0, y0, z_source, 'rp', 'MarkerSize', 15, 'MarkerFaceColor', 'r');
title('Task C: Microphone Array and Source Position (3D view)');
xlabel('x (m)'); ylabel('y (m)'); 
zlabel('z (m)');
grid on;
legend('Microphones','Source');
view(3);


% SPL calculation
p_rms = abs(p) / sqrt(2);
SPL = 20 * log10(p_rms / 2e-5);

fig13 = figure(13);
set(fig13,"position",[0,0,1000,400]);
subplot(1,2,1)
pcolor(X, Y, SPL); shading flat; colorbar;
title('Task D: SPL (dB SPL) for source at (0,0,-3)');
xlabel('x (m)'); ylabel('y (m)');

figure(13);
subplot(1,2,2)
pcolor(X, Y, angle(p)); shading flat; colorbar;
title('Task D: Phase (rad) of pressure for source at (0,0,-3)');
xlabel('x (m)'); ylabel('y (m)');

P_k = fftshift(fft2(p, N, N));
kz = sqrt(k^2 - KX.^2 - KY.^2);

fig14 = figure(14);
set(fig14,"position",[0,0,1000,400]);
subplot(1,2,1)
pcolor(KX, KY, abs(P_k)); shading flat; colorbar;
title('Task D: Wavenumber magnitude spectrum of pressure at (0,0,-3)');
xlabel('k_x (rad/m)'); ylabel('k_y (rad/m)');

figure(14);
subplot(1,2,2)
pcolor(KX, KY, angle(P_k)); shading flat; colorbar;
title('Task D: Wavenumber phase spectrum of pressure at (0,0,-3)');
xlabel('k_x (rad/m)'); ylabel('k_y (rad/m)');


% D Analysis% In this case, the source is positioned in the center again but 3 m away 
% from the microphone array indicating that the array is more or less
% in the farfield of the source. Far from the source, all the evanescent waves magnitudes have
% already decayed a lot. This indicates that, for the microphone array, it looks in
% the wavenumber spectrum like a plane wave propagating in the normal direction, this is shown by a dot in the
% middle (kx,ky ->0). By comparing the SPL values from figure 3 and figure 13 we can
% also see the level decay caused by the distance and the 1/r decrease.



%% TASK E - Wavenumber spectrum of normal velocity (Source at (0,0,-0.1))
z_source = -0.1;
f = 1000;
omega = 2 * pi * f;
k = omega / c;
R = sqrt(X.^2 + Y.^2 + (0 - z_source).^2);
p = 1j * omega * rho * Q ./ (4 * pi * R) .* exp(-1j * k * R);

uz = p ./ (1j * omega * rho) .* (0 - z_source) ./ R;
Uz_k = fftshift(fft2(uz, N, N));
kz = sqrt(k^2 - KX.^2 - KY.^2);

fig15 = figure(15);
set(fig15,"position",[0,0,1000,400]);
subplot(1,2,1)
pcolor(KX, KY, abs(Uz_k)); shading flat; colorbar;
title('Task E: Wavenumber magnitude spectrum of normal velocity');
xlabel('k_x (rad/m)'); ylabel('k_y (rad/m)');

figure(15);
subplot(1,2,2)
pcolor(KX, KY, angle(Uz_k)); shading flat; colorbar;
title('Task E: Wavenumber phase spectrum of normal velocity');
xlabel('k_x (rad/m)'); ylabel('k_y (rad/m)');

fig16 = figure(16);
set(fig16,"position",[0,0,1000,400]);
subplot(1,2,1)
pcolor(KX, KY, real(kz)); shading flat; colorbar;
title('Task E: Real part of k_z for normal velocity');
xlabel('k_x (rad/m)'); ylabel('k_y (rad/m)');

figure(16);
subplot(1,2,2)
pcolor(KX, KY, imag(kz)); shading flat; colorbar;
title('Task E: Imag part of k_z for normal velocity');
xlabel('k_x (rad/m)'); ylabel('k_y (rad/m)');

% HAS TO BE ANALYSED
% here I'm a bit lost but i think the point is that they are somewhat out
% of phase I guess ??

% E Analysis :
% Where, in figure 5 the pressure shower a maxima on a circle, in this case the particle velocity 
% exibit high value also at the center.
% As the particle velocity decay faster than ??
% It tends also to have less artifact than the wavenumber transform of the pressure.