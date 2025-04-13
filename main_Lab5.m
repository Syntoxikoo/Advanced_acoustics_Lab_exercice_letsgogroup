%%

clear;
close all;
clc;
addpath(genpath('scripts'))
addpath(genpath("functions"))
addpath(genpath("datas"))

%% ---- Import Data

load('Lab_5_data.mat'); % load data

%% --- SETTINGS ---
lineWidth = 2;
fontSize = 14;
figureSize = [100, 100, 1000, 600];

% Frequency tick labels for 1/3 octave bands
freq_labels = f_IS;

%% --- CONVERT INTENSITY SPECTRUM DATA TO dB re 1 pW ---
Dipole_IS      = 10 * log10(Dipole_IS      / 1e-12);
Dipole_fan_IS  = 10 * log10(Dipole_fan_IS  / 1e-12);
Fan_IS         = 10 * log10(Fan_IS         / 1e-12);
Omni_IS        = 10 * log10(Omni_IS        / 1e-12);
Omni_1_IS      = 10 * log10(Omni_1_IS      / 1e-12);
% and PI data is already in dB

%% --- MICROPHONE FREE-FIELD CORRECTION (for IS spectra only) ---

% Given correction values in dB
corr_freqs = [1250, 1600, 2000, 2500, 3150, 4000, 5000, 6300, 8000, 10000]; % Hz
% corr_values = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0];          % dB
corr_values = [0.1, 0.2, 0.4, 0.6, 1.0, 1.7, 2.5, 3.8, 5.8, 8.3];          % dB
% Interpolate correction for your frequency vector
mic_corr = interp1(corr_freqs, corr_values, f_IS, 'linear', 0); % extrapolated as 0 for <1250 Hz

% Apply correction to all IS spectra (in dB)
Dipole_IS      = Dipole_IS     + mic_corr;
Dipole_fan_IS  = Dipole_fan_IS + mic_corr;
Fan_IS         = Fan_IS        + mic_corr;
Omni_IS        = Omni_IS       + mic_corr;
Omni_1_IS      = Omni_1_IS     + mic_corr;


%% --- PLOT: INTENSITY SPECTRUM (IS) ---
figure('Position', figureSize);
semilogx(f_IS, Dipole_IS,     'DisplayName','Dipole (Operator A)', 'LineWidth', lineWidth);
hold on;
semilogx(f_IS, Dipole_fan_IS, 'DisplayName','Dipole + Fan',        'LineWidth', lineWidth);
semilogx(f_IS, Fan_IS,        'DisplayName','Fan only',            'LineWidth', lineWidth);
semilogx(f_IS, Omni_IS,       'DisplayName','Omni Mic',            'LineWidth', lineWidth);
semilogx(f_IS, Omni_1_IS,     'DisplayName','Omni Mic (Op B)',     'LineWidth', lineWidth);
hold off;
grid on;
title('Intensity Spectrum', 'FontSize', fontSize + 2);
xlabel('Frequency [Hz]', 'FontSize', fontSize);
ylabel('Sound Power Level [dB re 1 pW]', 'FontSize', fontSize);
legend('FontSize', fontSize);
set(gca, 'FontSize', fontSize);
xticks(freq_labels);
xticklabels(string(freq_labels));
xlim([50, 10000]);
ylim auto;

%% --- PLOT: PI INDEX SPECTRUM ---
figure('Position', figureSize);
semilogx(f_PI, Dipole_PI,     'DisplayName','Dipole (Op A)',       'LineWidth', lineWidth);
hold on;
semilogx(f_PI, Dipole_fan_PI, 'DisplayName','Dipole + Fan',        'LineWidth', lineWidth);
semilogx(f_PI, Fan_PI,        'DisplayName','Fan only',            'LineWidth', lineWidth);
semilogx(f_PI, Omni_PI,       'DisplayName','Omni Mic',            'LineWidth', lineWidth);
semilogx(f_PI, Omni_1_PI,     'DisplayName','Omni Mic (Op B)',     'LineWidth', lineWidth);
hold off;
grid on;
title('PI Index Spectrum (Signal 1, Signal 2)', 'FontSize', fontSize + 2);
xlabel('Frequency [Hz]', 'FontSize', fontSize);
ylabel('PI Index [dB]', 'FontSize', fontSize);
legend('Location','northeast', 'FontSize', fontSize);
set(gca, 'FontSize', fontSize);
xticks(freq_labels);
xticklabels(string(freq_labels));
xlim([50, 10000]);
ylim auto;
