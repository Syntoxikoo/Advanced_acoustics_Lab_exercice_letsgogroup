clc;
clear;

addpath(genpath("functions"))

N = 1000;  % Number of plane waves
M = 1000; % Number of Monte Carlo realizations

[mean_p_hat_sq, rel_std_p_hat_sq, rel_std_p_hat, std_Lp] = mean_square_p(N, M);
