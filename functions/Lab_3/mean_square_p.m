function [mean_p_hat_sq, rel_std_p_hat_sq, rel_std_p_hat, std_Lp] = mean_square_p(N, M)

    p_hat_sq = zeros(M,1);
    p_hat = zeros(M,1);
    Lp = zeros(M,1);

    % Monte Carlo simulation
    for m = 1:M
        phi = 2 * pi * rand(N,1);  % Random phases in [0, 2*pi]
        p_hat_m = (1/sqrt(2*N)) * sum(exp(1j * phi)); % Compute normalized pressure sum
        p_hat_sq(m) = abs(p_hat_m)^2; % Mean square pressure
        p_hat(m) = abs(p_hat_m);     % Pressure amplitude
        Lp(m) = 10 * log10(p_hat_sq(m)); % SPL in dB
    end

    mean_p_hat_sq = mean(p_hat_sq);
    std_p_hat_sq = std(p_hat_sq);
    rel_std_p_hat_sq = std_p_hat_sq / mean_p_hat_sq;  % Expected: 1

    mean_p_hat = mean(p_hat);
    std_p_hat = std(p_hat);
    rel_std_p_hat = std_p_hat / mean_p_hat;  % Expected: 0.52

    std_Lp = std(Lp);  % Expected: 5.6 dB

    % Display results
    fprintf('Mean of |p̂|^2: %.3f\n', mean_p_hat_sq);
    fprintf('Relative Std Dev of |p̂|^2: %.3f (Expected: 1)\n', rel_std_p_hat_sq);
    fprintf('Relative Std Dev of |p̂|: %.3f (Expected: 0.52)\n', rel_std_p_hat);
    fprintf('Std Dev of Lp: %.3f dB (Expected: 5.6 dB)\n', std_Lp);
end

    