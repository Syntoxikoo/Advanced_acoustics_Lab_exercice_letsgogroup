N = 10000;        % Number of plane waves
M = 10000;       % Number of Monte Carlo iterations

p_hat_sq = zeros(M,1);  % Mean square pressure
p_hat = zeros(M,1);     % Sound pressure amplitude
Lp = zeros(M,1);        % Sound pressure level in dB

for m = 1:M
    phi = 2 * pi * rand(N,1);% phases (uniform in [0, 2*pi])
    p_hat_m = (1/sqrt(2*N)) * sum(exp(1j * phi)); % Compute pressure sum and normalize
    p_hat_sq(m) = abs(p_hat_m)^2; % Compute mean square pressure, amplitude, and SPL
    p_hat(m) = abs(p_hat_m);
    Lp(m) = 10 * log10(p_hat_sq(m));  % SPL in dB
end


mean_p_hat_sq = mean(p_hat_sq);
std_p_hat_sq = std(p_hat_sq);
rel_std_p_hat_sq = std_p_hat_sq / mean_p_hat_sq;  % Should be 1

mean_p_hat = mean(p_hat);
std_p_hat = std(p_hat);
rel_std_p_hat = std_p_hat / mean_p_hat;  % Should be 0.52

std_Lp = std(Lp);  % Should be 5.6 dB

% Display results
fprintf('Mean of |p̂|^2: %.3f\n', mean_p_hat_sq);
fprintf('Relative Std Dev of |p̂|^2: %.3f (Expected: 1)\n', rel_std_p_hat_sq);
fprintf('Relative Std Dev of |p̂|: %.3f (Expected: 0.52)\n', rel_std_p_hat);
fprintf('Std Dev of Lp: %.3f dB (Expected: 5.6 dB)\n', std_Lp);

%%
N = 1000;  % Number of plane waves
M_values = round(linspace(1,10e3,100)); % Different Monte Carlo iterations

rel_std_p_hat_sq = zeros(length(M_values),1);
rel_std_p_hat = zeros(length(M_values),1);
std_Lp = zeros(length(M_values),1);


for idx = 1:length(M_values)
    M = M_values(idx);
   
    p_hat_sq = zeros(M,1);
    p_hat = zeros(M,1);
    Lp = zeros(M,1);
    
    for m = 1:M
        phi = 2 * pi * rand(N,1); % phases (uniform in [0, 2*pi])
        
        p_hat_m = (1/sqrt(2*N)) * sum(exp(1j * phi)); % Compute pressure sum and normalize
       
        % Compute mean square pressure, amplitude, and SPL
        p_hat_sq(m) = abs(p_hat_m)^2;
        p_hat(m) = abs(p_hat_m);
        Lp(m) = 10 * log10(p_hat_sq(m));  % SPL in dB
    end
    
    rel_std_p_hat_sq(idx) = std(p_hat_sq) / mean(p_hat_sq);
    rel_std_p_hat(idx) = std(p_hat) / mean(p_hat);
    std_Lp(idx) = std(Lp);
end


figure;
subplot(3,1,1);
semilogx(M_values, rel_std_p_hat_sq, 'o-', 'LineWidth', 2);
xlabel('Number of Monte Carlo Iterations (M)');
ylabel('Relative Std Dev of |p̂|^2');
title('Convergence of Relative Std Dev of |p̂|^2');
grid on;

subplot(3,1,2);
semilogx(M_values, rel_std_p_hat, 's-', 'LineWidth', 2);
xlabel('Number of Monte Carlo Iterations (M)');
ylabel('Relative Std Dev of |p̂|');
title('Convergence of Relative Std Dev of |p̂|');
grid on;

subplot(3,1,3);
semilogx(M_values, std_Lp, 'd-', 'LineWidth', 2);
xlabel('Number of Monte Carlo Iterations (M)');
ylabel('Std Dev of L_p (dB)');
title('Convergence of Std Dev of L_p');
grid on;

