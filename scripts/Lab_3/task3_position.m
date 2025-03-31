%% function

Nwaves = 10000; % amount of plane waves for each point
Npoints = 10000; % space points
Lz = 7;        % z length room (m)
c = 343;       % (m/s)
f = 700;       % (Hz)
omega = 2 * pi * f ;
k = omega / c;

z_values = linspace(0, Lz, Npoints); % space points
p2_values = zeros(size(z_values));

phi_random = [];
theta_random = [];

rand_idx = randi([1, Npoints]);

for idx = 1:length(z_values)
    r = z_values(idx);
    phi_i = rand(1, Nwaves) * 2 * pi; 
    theta_i = acos(2 * rand(1, Nwaves) - 1); 
    
    if idx == rand_idx
        phi_random = phi_i;
        theta_random = theta_i;
    end
    
    sum_exp = sum(exp(-1j * (k * r * cos(theta_i) + phi_i))); 
    p2_values(idx) = (1/Nwaves) * abs(sum_exp)^2;
end

fprintf('Size of p2_values array: %d x %d\n', size(p2_values,1), size(p2_values,2));

p_ref = 2 * 10^(-5);
pref_sq = p_ref^2;

% Plot SPL in dB
figure;
plot(z_values, 10*log10(p2_values / pref_sq), 'LineWidth', 1.5);
grid on;
title('Relative Sound Pressure Level vs Distance at 700 Hz');
xlabel('z (m)');
ylabel('L_p (dB SPL)');

%% Compute pressure-related statistics
p2_rms = p2_values; % p^2_{rms} is simply p2_values

mean_p2 = mean(p2_rms); % Expectation of mean square pressure
mean_p = sqrt(mean_p2); % mean of p_rms
std_p = std(sqrt(p2_rms)); % std of p_rms
rel_std_p = std_p / mean_p * 100; 

Lp_values = 20 * log10(sqrt(p2_rms) / p_ref); % Convert p_rms to SPL in dB
mean_Lp = mean(Lp_values);
std_Lp = std(Lp_values);
rel_std_Lp = std_Lp / mean_Lp * 100;  % Relative standard deviation for Lp
L0 = 10 * log10(mean_p2 / (p_ref^2)); % Reference pressure level

% Compute standard deviation and relative standard deviation for p^2_rms
std_p2_rms = std(p2_rms);
rel_std_p2_rms = std_p2_rms / mean_p2 * 100;  % Relative standard deviation for p^2_rms

N_bins = 30;

% Compute histograms for PDFs
[p2_rms_hist, p2_rms_bins] = hist(p2_rms, N_bins); % Using p^2_rms
p2_rms_hist_normalized = p2_rms_hist / sum(p2_rms_hist);

[Lp_hist, Lp_bins] = hist(Lp_values, N_bins);
Lp_hist_normalized = Lp_hist / sum(Lp_hist);

% Theoretical distributions
p2_rms_pdf_theoretical = (1 / mean_p2) * exp(-p2_rms_bins / mean_p2);
p2_rms_pdf_theoretical_scaled = (p2_rms_pdf_theoretical / max(p2_rms_pdf_theoretical)) * max(p2_rms_hist_normalized);

Lp_pdf_theoretical = (log(10)/10) * exp((log(10)/10) * (Lp_bins - L0) - exp((log(10)/10) * (Lp_bins - L0)));
Lp_pdf_theoretical_scaled = (Lp_pdf_theoretical / max(Lp_pdf_theoretical)) * max(Lp_hist_normalized);

% Create a figure for probability functions
figure;

% Plot p^2_rms probability function
subplot(1, 2, 1);
bar(p2_rms_bins, p2_rms_hist_normalized, 'hist');
hold on;
plot(p2_rms_bins, p2_rms_pdf_theoretical_scaled, 'r', 'LineWidth', 2);
grid on;
title(['PDF of p^2_{rms} (Mean = ' num2str(mean_p2, '%.2f') ' Pa^2, Std = ' num2str(std_p2_rms, '%.2f') ' Pa^2, Rel. Std = ' num2str(rel_std_p2_rms, '%.2f') '%)']);
xlabel('p^2_{rms} (Pa^2)', 'FontSize', 14);
ylabel('Probability Density', 'FontSize', 14);
legend('Observed Distribution', 'Scaled Theoretical Distribution', 'FontSize', 12);
set(gca, 'FontSize', 14);

% Plot SPL probability function
subplot(1, 2, 2);
bar(Lp_bins, Lp_hist_normalized, 'hist');
hold on;
plot(Lp_bins, Lp_pdf_theoretical_scaled, 'r', 'LineWidth', 2);
grid on;
title(['PDF of SPL (Mean = ' num2str(mean_Lp, '%.2f') ' dB, Std = ' num2str(std_Lp, '%.2f') ' dB, Rel. Std = ' num2str(rel_std_Lp, '%.2f') '%)']);
xlabel('L_p (dB SPL)', 'FontSize', 14);
ylabel('Probability Density', 'FontSize', 14);
legend('Observed Distribution', 'Scaled Theoretical Distribution', 'FontSize', 12);
set(gca, 'FontSize', 14);

% Compute statistics on phi and theta
theta_bins = linspace(0, pi, N_bins);
theta_hist = hist(theta_random, theta_bins);
theta_hist_normalized = theta_hist / sum(theta_hist);
max_theta_hist = max(theta_hist_normalized);
theta_pdf_theoretical = sin(theta_bins) / 2;
theta_pdf_theoretical_scaled = (theta_pdf_theoretical / max(theta_pdf_theoretical)) * max_theta_hist;

phi_bins = linspace(0, 2*pi, N_bins);
phi_hist = hist(phi_random, phi_bins);
phi_hist_normalized = phi_hist / sum(phi_hist);
max_phi_hist = max(phi_hist_normalized);
phi_pdf_theoretical = ones(size(phi_bins)) / (2*pi);
phi_pdf_theoretical_scaled = (phi_pdf_theoretical / max(phi_pdf_theoretical)) * max_phi_hist;

% Create a figure for theta and phi probability functions
figure;

% Plot theta histogram
subplot(1, 2, 1);
bar(theta_bins, theta_hist_normalized, 'hist');
hold on;
plot(theta_bins, theta_pdf_theoretical_scaled, 'r', 'LineWidth', 2);
grid on;
title('Histogram of theta at Random z point with theoretical curve');
xlabel('Theta Angle (radians)', 'FontSize', 14);
xlim([0 pi]);
ylabel('Probability Density', 'FontSize', 14);
legend('Observed Distribution', 'Scaled Theoretical Distribution', 'FontSize', 12);
set(gca, 'FontSize', 14);

% Plot phi histogram
subplot(1, 2, 2);
bar(phi_bins, phi_hist_normalized, 'hist');
hold on;
plot(phi_bins, phi_pdf_theoretical_scaled, 'r', 'LineWidth', 2);
grid on;
title('Histogram of Phi at random z Point with theoretical curve');
xlabel('Phase Angle (radians)', 'FontSize', 14);
xlim([0 2*pi]);
ylabel('Probability Density', 'FontSize', 14);
legend('Observed Distribution', 'Scaled Theoretical Distribution', 'FontSize', 12);
set(gca, 'FontSize', 14);

%% OLD VERSION OF STATS PART

% p_ref = 2 * 10^(-5);
% 
% % Relative sound pressure level (SPL)
% figure;
% plot(z_values, 10*log10(p2_values / (p_ref)^2), 'LineWidth', 1.5);
% grid on;
% title('Relative Sound Pressure Level vs Distance at 700 Hz');
% xlabel('z (m)');
% ylabel('L_p (dB SPL)');
% 
% % stats
% mean_p2 = mean(p2_values);
% std_p2 = std(p2_values);
% rel_std_p2 = std_p2 / mean_p2 * 100; % relative std (%)
% 
% mean_p2_dB = mean(10*log10(p2_values / (p_ref)^2));
% std_p2_dB = std(10*log10(p2_values / (p_ref)^2));
% rel_std_p2_dB = std_p2_dB / mean_p2_dB * 100; % relative std (%)
% 
% % Histograms (p2 and dB) for one sound field
% N_bins = 30;
% figure;
% hist(p2_values, N_bins);
% grid on;
% title(['Histogram of p_2 values across spatial points ' ...
%     '(Mean = ' num2str(mean_p2, '%.2f') ' Pa, Std = ' num2str(std_p2, '%.2f') ...
%     ' Pa, Rel. Std = ' num2str(rel_std_p2, '%.2f') ' %)']);
% xlabel('p_2 values');
% ylabel('Counts');
% 
% figure;
% hist(10*log10(p2_values / (p_ref)^2), N_bins);
% grid on;
% title(['Histogram of SPL across spatial points ' ...
%     '(Mean = ' num2str(10*log10(mean_p2_dB / (p_ref)^2), '%.2f') ' dB, ' ...
%     'Std = ' num2str(10*log10(std_p2_dB / (p_ref)^2), '%.2f') ' dB, ' ...
%     'Rel. Std = ' num2str(rel_std_p2_dB, '%.2f') '%)']);
% xlabel('L_p (dB SPL)');
% ylabel('Counts');


%% plot 2 Npoints config

% Define parameters
Lz = 7;        % z length room (m)
c = 343;       % speed of sound (m/s)
f = 700;       % frequency (Hz)
omega = 2 * pi * f;
k = omega / c;

% Reference pressure
p_ref = 2 * 10^(-5);
pref_sq = p_ref^2;

% Create a figure
figure;
hold on;

% Define the two combinations of Npoints and Nwaves
Npoints_values = [10000, 1000 , 100];  % Npoints values for plotting
Nwaves = 10000;  % Keep Nwaves constant for simplicity

% Loop through both Npoints values and plot
for i = 1:length(Npoints_values)
    Npoints = Npoints_values(i);  % Current Npoints value
    
    % Define the z-values and p2_values array for the current combination
    z_vals = linspace(0, Lz, Npoints);
    p2_vals = zeros(size(z_vals));
    
    % Loop over each point to calculate pressure values
    for idx = 1:length(z_vals)
        r = z_vals(idx);
        phi_i = rand(1, Nwaves) * 2 * pi; 
        theta_i = acos(2 * rand(1, Nwaves) - 1); 

        sum_exp = sum(exp(-1j * (k * r * cos(theta_i) + phi_i))); 
        p2_vals(idx) = (1/Nwaves) * abs(sum_exp)^2;
    end
    
    % Plot the results in dB SPL
    plot(z_vals, 10*log10(p2_vals / pref_sq), 'LineWidth', 2, ...
         'DisplayName', sprintf('Npoints = %d', Npoints));
end

% Add labels, title, and legend
grid on;
%title('Relative Sound Pressure Level vs Distance at 700 Hz');
set(gca, 'GridLineStyle', '-', 'LineWidth', 1); % Thicker grid lines and axis lines
ax = gca;
ax.XAxis.FontSize = 14;
ax.YAxis.FontSize = 14;
xlabel('z (m)', 'FontSize', 14);
ylabel('L_p (dB SPL)', 'FontSize', 14);
legend('show', 'FontSize', 14);
hold off;



















%% function V2

Nwaves = 100; % Amount of plane waves for each point
Npoints = 10000; % Space points
Lz = 7;        % z length room (m)
c = 343;       % (m/s)
f = 700;       % (Hz)
omega = 2 * pi * f;
k = omega / c;
M = 100;       % Number of Monte Carlo simulations per point

z_values = linspace(0, Lz, Npoints); % Space points
p2_values = zeros(size(z_values));

phi_random = [];
theta_random = [];

rand_idx = randi([1, Npoints]);

for idx = 1:length(z_values)
    r = z_values(idx);
    p2_temp = zeros(1, M);
    
    for m = 1:M
        phi_i = rand(1, Nwaves) * 2 * pi; 
        theta_i = acos(2 * rand(1, Nwaves) - 1);
        
        if idx == rand_idx && m == 1  % Store a single realization for reference
            phi_random = phi_i;
            theta_random = theta_i;
        end
        
        sum_exp = sum(exp(-1j * (k * r * cos(theta_i) + phi_i))); 
        p2_temp(m) = (1/Nwaves) * abs(sum_exp)^2;
    end
    
    p2_values(idx) = mean(p2_temp); % Average over Monte Carlo simulations
end

fprintf('Size of p2_values array: %d x %d\n', size(p2_values,1), size(p2_values,2));

p_ref = 2 * 10^(-5);
pref_sq = p_ref^2;

figure;
plot(z_values, 10*log10(p2_values / pref_sq), 'LineWidth', 1.5);
grid on;
title('Relative Sound Pressure Level vs Distance at 700 Hz');
xlabel('z (m)');
ylabel('L_p (dB SPL)');



%% plot 2 Npoints config with Monte Carlo averaging

% Define parameters
Lz = 7;        % z length room (m)
c = 343;       % speed of sound (m/s)
f = 700;       % frequency (Hz)
omega = 2 * pi * f;
k = omega / c;

% Reference pressure
p_ref = 2 * 10^(-5);
pref_sq = p_ref^2;

% Number of Monte Carlo simulations per point
M = 100;

% Create a figure
figure;
hold on;

% Define the combinations of Npoints
Npoints_values = [10000, 1000, 100];  % Npoints values for plotting
Nwaves = 100;  % Keep Nwaves constant for simplicity

% Loop through different Npoints values and plot
for i = 1:length(Npoints_values)
    Npoints = Npoints_values(i);  % Current Npoints value
    
    % Define the z-values and p2_values array for the current combination
    z_vals = linspace(0, Lz, Npoints);
    p2_vals = zeros(size(z_vals));
    
    % Loop over each point to calculate pressure values with Monte Carlo averaging
    for idx = 1:length(z_vals)
        r = z_vals(idx);
        p2_temp = zeros(1, M);
        
        for mc = 1:M
            phi_i = rand(1, Nwaves) * 2 * pi;
            theta_i = acos(2 * rand(1, Nwaves) - 1);
            
            sum_exp = sum(exp(-1j * (k * r * cos(theta_i) + phi_i)));
            p2_temp(mc) = (1/Nwaves) * abs(sum_exp)^2;
        end
        
        p2_vals(idx) = mean(p2_temp); % Average over Monte Carlo simulations
    end
    
    % Plot the results in dB SPL
    plot(z_vals, 10*log10(p2_vals / pref_sq), 'LineWidth', 2, ...
         'DisplayName', sprintf('Npoints = %d', Npoints));
end

% Add labels, title, and legend
grid on;
set(gca, 'GridLineStyle', '-', 'LineWidth', 1); % Thicker grid lines and axis lines
ax = gca;
ax.XAxis.FontSize = 14;
ax.YAxis.FontSize = 14;
xlabel('z (m)', 'FontSize', 14);
ylabel('L_p (dB SPL)', 'FontSize', 14);
legend('show', 'FontSize', 14);
hold off;