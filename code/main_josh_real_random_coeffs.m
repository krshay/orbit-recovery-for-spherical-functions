clear; clc;

rng(10);

% --- Parameters ---
d = 3;  % Bandlimit
n_trials = 3;
shell_values = 1:5;
avg_errors = zeros(size(shell_values));
std_errors = zeros(size(shell_values));

for idx = 1:length(shell_values)
    num_shells = shell_values(idx);
    rel_errors = zeros(1, n_trials);
    fprintf('\n=== Running for %d shells ===\n', num_shells);

    for trial = 1:n_trials
        n_coeffs_per_shell = (d+1)^2;
        true_coeffs = randn(num_shells, n_coeffs_per_shell);

        fixed_bands = [0,1];
        fixed_indices = [];
        for l = fixed_bands
            for m = -l:l
                idx_sh = sh_index(l, m);
                fixed_indices = [fixed_indices, idx_sh];
            end
        end

        est_coeffs = zeros(size(true_coeffs));
        est_coeffs(:, fixed_indices) = true_coeffs(:, fixed_indices);

        options = optimoptions('fminunc', ...
             'Algorithm', 'quasi-newton', ...
            'Display', 'off', ...
            'MaxIterations', 500, ...
            'MaxFunctionEvaluations', 2e5, ...
            'OptimalityTolerance', 1e-10, ...
            'StepTolerance', 1e-10, ...
            'FunctionTolerance', 1e-10, ...
            'SpecifyObjectiveGradient', true);

        for ell = 2:d
            current_indices = [];
            for m = -ell:ell
                current_indices = [current_indices, sh_index(ell, m)];
            end

            for shell = 1:num_shells
                x_init = randn(size(current_indices));
                obj = @(x_shell) specialized_loss_and_fd_grad(x_shell, est_coeffs, ...
                    shell, current_indices, ell, d, true_coeffs);

                [x_sol, ~] = fminunc(obj, x_init, options);
                est_coeffs(shell, current_indices) = x_sol;
            end
        end

        final_err = norm(est_coeffs - true_coeffs, 'fro');
        rel_errors(trial) = final_err / norm(true_coeffs, 'fro');
        fprintf('Trial %d: Relative error = %.4f%%\n', trial, 100 * rel_errors(trial));
    end

    avg_errors(idx) = mean(rel_errors);
    std_errors(idx) = std(rel_errors);
    fprintf('\nAverage error for %d shells: %.4f%% (std = %.4f%%)\n', ...
        num_shells, 100 * avg_errors(idx), 100 * std_errors(idx));
end

figure;
errorbar(shell_values, 100 * avg_errors, 100 * std_errors, 'o-', 'LineWidth', 2);
xlabel('Number of Shells');
ylabel('Average Relative Error (%)');
title('Frequency Marching Recovery Error vs. Number of Shells');
grid on;
