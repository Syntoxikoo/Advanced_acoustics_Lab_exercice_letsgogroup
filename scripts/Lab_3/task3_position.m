clc;
clear;

Nwaves = 5000; % amount of plane waves for each point
Npoints = 100; % space points
Lz = 7;        % z length room (m)
c = 343;       % (m/s)
f = 700;       % (Hz)
omega = 2 * pi * f ;
k = omega / c;

z_values = linspace(0, Lz, Npoints); % space points
p2_values = zeros(size(z_values));

% calculate p_2 for each position as a sum of every random wave
for idx = 1:length(z_values)
    r = z_values(idx);
    phi_i = rand(1, Nwaves) * 2 * pi; % random phases
    theta_i = acos(2 * rand(1, Nwaves) - 1); % random angles
    sum_exp = sum(exp(-1j * (k * r * cos(theta_i) + phi_i))); 
    p2_values(idx) = (1/Nwaves) * abs(sum_exp)^2;
end

% Print array size info
fprintf('Size of p2_values array: %d x %d\n', size(p2_values,1), size(p2_values,2));

p_ref = 2 * 10^(-5);


% Plot relative sound pressure level (SPL)
figure;
plot(z_values, 10*log10(p2_values / (p_ref)^2), 'LineWidth', 1.5);
grid on;
title('Relative Sound Pressure Level vs Distance at 700 Hz');
xlabel('z (m)');
ylabel('L_p (dB SPL)');

% stats
mean_p2 = mean(p2_values);
std_p2 = std(p2_values);
rel_std_p2 = std_p2 / mean_p2 * 100; % relative std (%)

mean_p2_dB = mean(10*log10(p2_values / (p_ref)^2));
std_p2_dB = std(10*log10(p2_values / (p_ref)^2));
rel_std_p2_dB = std_p2_dB / mean_p2_dB * 100; % relative std (%)

% Histogram (p2 values)
N_bins = 30;
figure;
hist(p2_values, N_bins);
grid on;
title(['Histogram of p_2 values across spatial points ' ...
    '(Mean = ' num2str(mean_p2, '%.2f') ' Pa, Std = ' num2str(std_p2, '%.2f') ...
    ' Pa, Rel. Std = ' num2str(rel_std_p2, '%.2f') ' %)']);
xlabel('p_2 values');
ylabel('Counts');

% Optional: plot histogram in dB scale
figure;
hist(10*log10(p2_values / (p_ref)^2), N_bins);
grid on;
title(['Histogram of SPL across spatial points ' ...
    '(Mean = ' num2str(10*log10(mean_p2_dB / (p_ref)^2), '%.2f') ' dB, ' ...
    'Std = ' num2str(10*log10(std_p2_dB / (p_ref)^2), '%.2f') ' dB, ' ...
    'Rel. Std = ' num2str(rel_std_p2_dB, '%.2f') '%)']);
xlabel('L_p (dB SPL)');
ylabel('Counts');
