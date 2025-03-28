%% function

Nwaves = 5000; % amount of plane waves for each point
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

figure;
plot(z_values, 10*log10(p2_values / pref_sq), 'LineWidth', 1.5);
grid on;
title('Relative Sound Pressure Level vs Distance at 700 Hz');
xlabel('z (m)');
ylabel('L_p (dB SPL)');

%% stats

p_rms = sqrt(p2_values/2);

mean_p = mean(p_rms);
std_p = std(p_rms);
rel_std_p = std_p / mean_p * 100; 

Lp_values = 20 * log10(sqrt(p2_values)  / p_ref);
mean_Lp = mean(Lp_values);
std_Lp = std(Lp_values);

N_bins = 30;
figure;
hist(p_rms, N_bins);
grid on;
title(['Histogram of p_rms values across spatial points ' ...
    '(Mean = ' num2str(mean_p, '%.2f') ' Pa, Std = ' num2str(std_p, '%.2f') ...
    ' Pa, Rel. Std = ' num2str(rel_std_p, '%.2f') ' %)']);
xlabel('p_2 values');
ylabel('Counts');

figure;
hist(Lp_values, N_bins);
grid on;
title(['Histogram of SPL across spatial points ' ...
    '(Mean = ' num2str(mean_Lp, '%.2f') ' dB, ' ...
    'Std = ' num2str(std_Lp, '%.2f') ' dB)']);
xlabel('L_p (dB SPL)');
ylabel('Counts');


% stats on phi and theta
theta_bins = linspace(0, pi, N_bins);
theta_hist = hist(theta_random, theta_bins);
theta_pdf_theoretical = sin(theta_bins) / 2;

theta_hist_normalized = theta_hist / sum(theta_hist); 
max_theta_hist = max(theta_hist_normalized); 
theta_pdf_theoretical_scaled = (theta_pdf_theoretical / max(theta_pdf_theoretical)) * max_theta_hist;

% Phi distribution
phi_bins = linspace(0, 2*pi, N_bins);
phi_hist = hist(phi_random, phi_bins);
phi_pdf_theoretical = ones(size(phi_bins)) / (2*pi); 

% Normalization
phi_hist_normalized = phi_hist / sum(phi_hist); 
max_phi_hist = max(phi_hist_normalized); 
phi_pdf_theoretical_scaled = (phi_pdf_theoretical / max(phi_pdf_theoretical)) * max_phi_hist;

% Plot theta histogram with theoretical curve
figure;
bar(theta_bins, theta_hist_normalized, 'hist'); 
hold on;
plot(theta_bins, theta_pdf_theoretical_scaled, 'r', 'LineWidth', 2);
grid on;
title('Histogram of theta at Random z point with theoretical curve');
xlabel('Theta Angle (radians)');
ylabel('Probability density function');
legend('Observed Distribution', 'Scaled Theoretical Distribution');

% Plot phi histogram with theoretical curve
figure;
bar(phi_bins, phi_hist_normalized, 'hist'); 
hold on;
plot(phi_bins, phi_pdf_theoretical_scaled, 'r', 'LineWidth', 2);
grid on;
title('Histogram of Phi at random z Point with theoretical curve');
xlabel('Phase Angle (radians)');
ylabel('Probability density function');
legend('Observed Distribution', 'Scaled Theoretical Distribution');

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
