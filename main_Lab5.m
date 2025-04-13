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
