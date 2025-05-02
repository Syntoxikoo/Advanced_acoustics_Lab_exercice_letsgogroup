% Plot pressure magnitude =f( theta) for m = 0, 1, 2

theta = linspace(0, 2*pi, 500);
cth = cos(theta); % cos(theta) for Legendre

%polynoms
P0 = legendre(0, cth); P0 = P0(1,:);
P1 = legendre(1, cth); P1 = P1(1,:);
P2 = legendre(2, cth); P2 = P2(1,:);

% normalization by axial pressure P_i(1)
P0_norm = abs(P0) / abs(P0(1));
P1_norm = abs(P1) / abs(P1(1));
P2_norm = abs(P2) / abs(P2(1));

P0_dB = 20 * log10(P0_norm);
P1_dB = 20 * log10(P1_norm);
P2_dB = 20 * log10(P2_norm);

% eps = 1e-6; % avoid error
% P0_dB = 20 * log10(P0_norm + eps);
% P1_dB = 20 * log10(P1_norm + eps);
% P2_dB = 20 * log10(P2_norm + eps);

figure;
polarplot(theta, P0_dB, 'LineWidth', 2); 
hold on;
polarplot(theta, P1_dB, 'LineWidth', 2);
polarplot(theta, P2_dB, 'LineWidth', 2);
%rlim([-40 0]);

legend('m = 0 (Monopole)', 'm = 1 (Dipole)', 'm = 2 (Quadrupole)', 'Location', 'southoutside');
title('Normalized Pressure Magnitude (in dB) in function of \theta (in °)');
grid on;
