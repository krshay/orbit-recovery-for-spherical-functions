clear; clc;

% --- Parameters ---
d = 3;  % Bandlimit
n_trials = 10;
shell_values = 5;
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
            'Display', 'iter-detailed', ...
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
                x_init = true_coeffs(shell, current_indices) + 1* randn(size(current_indices));
                obj = @(x_shell) specialized_loss_and_grad(x_shell, est_coeffs, ...
                    shell, current_indices, ell, d, num_shells, true_coeffs);
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
        num_shells, 100*avg_errors(idx), 100*std_errors(idx));
end

% Plot results
figure;
errorbar(shell_values, 100*avg_errors, 100*std_errors, 'o-', 'LineWidth', 2);
xlabel('Number of Shells');
ylabel('Average Relative Error (%)');
title('Frequency Marching Recovery Error vs. Number of Shells');
grid on;

% --- Helper functions ---
function idx = sh_index(l, m)
    idx = l^2 + m + l + 1;
end

function [loss, grad] = specialized_loss_and_grad(x_shell, est_coeffs, shell_idx, ...
        current_indices, ell, d, num_shells, true_coeffs)

    est_coeffs(shell_idx, current_indices) = x_shell;
    B_target = compute_SO3_bispectrum(true_coeffs, d);

    loss = scalar_loss_wrapper(x_shell, est_coeffs, shell_idx, current_indices, ell, d, num_shells, B_target);
    grad = analytic_gradient(x_shell, est_coeffs, shell_idx, current_indices, ell, d, num_shells, B_target);
end

function loss = scalar_loss_wrapper(x_shell, est_coeffs, shell_idx, current_indices, ell, d, num_shells, B_target)
    est_coeffs(shell_idx, current_indices) = x_shell;
    B_current = compute_SO3_bispectrum(est_coeffs, d);
    B_diff = B_current - B_target;

    loss = sum(abs(B_diff(:)).^2);
end

function grad = analytic_gradient(x_shell, est_coeffs, shell_idx, current_indices, ell, d, num_shells, B_target)
    grad = zeros(size(x_shell));
    est_coeffs(shell_idx, current_indices) = x_shell;
    B_current = compute_SO3_bispectrum(est_coeffs, d);
    B_diff = B_current - B_target;

    for i = 1:length(current_indices)
        l = ell;
        m = i - (l^2 + l);
        idx_i = current_indices(i);
        grad_val = 0;

        for l1 = 0:d
            for l2 = 0:d
                for l3 = abs(l1 - l2):min(l1 + l2, d)
                    lvals = [l1, l2, l3];
                    idx_ell = find(lvals == ell);
                    if length(idx_ell) ~= 1 || any(lvals(setdiff(1:3, idx_ell)) >= ell)
                        continue;
                    end
                    for s1 = 1:num_shells
                        for s2 = 1:num_shells
                            for s3 = 1:num_shells
                                shells = [s1, s2, s3];
                                if shells(idx_ell) ~= shell_idx
                                    continue;
                                end
                                for m1 = -l1:l1
                                    for m2 = -l2:l2
                                        m3 = -m1 - m2;
                                        if abs(m3) > l3, continue; end
                                        idx1 = sh_index(l1, m1);
                                        idx2 = sh_index(l2, m2);
                                        idx3 = sh_index(l3, m3);
                                        cg = clebsch_gordan(l1, m1, l2, m2, l3, m3);
                                        phase = (-1)^m1;
                                        coeffs = [est_coeffs(s1,idx1), est_coeffs(s2,idx2), est_coeffs(s3,idx3)];
                                        if idx_ell == 1 && idx1 == idx_i
                                            term = coeffs(2) * coeffs(3);
                                        elseif idx_ell == 2 && idx2 == idx_i
                                            term = coeffs(1) * coeffs(3);
                                        elseif idx_ell == 3 && idx3 == idx_i
                                            term = coeffs(1) * coeffs(2);
                                        else
                                            continue;
                                        end
                                        grad_val = grad_val + 2 * real(B_diff(s1,s2,s3,l1+1,l2+1,l3+1)) * term * cg * phase;
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
        grad(i) = grad_val;
    end
end

function B = compute_SO3_bispectrum(a_lm_matrix, d)
    num_shells = size(a_lm_matrix, 1);
    B = zeros(num_shells, num_shells, num_shells, d+1, d+1, d+1);
    for l1 = 0:d
        for l2 = 0:d
            for l3 = abs(l1 - l2):min(l1 + l2, d)
                for s1 = 1:num_shells
                    for s2 = 1:num_shells
                        for s3 = 1:num_shells
                            sum_m = 0;
                            for m1 = -l1:l1
                                for m2 = -l2:l2
                                    m3 = -m1 - m2;
                                    if abs(m3) > l3, continue; end
                                    idx1 = sh_index(l1, m1);
                                    idx2 = sh_index(l2, m2);
                                    idx3 = sh_index(l3, m3);
                                    cg = clebsch_gordan(l1, m1, l2, m2, l3, m3);
                                    a1 = a_lm_matrix(s1, idx1);
                                    a2 = a_lm_matrix(s2, idx2);
                                    a3 = a_lm_matrix(s3, idx3);
                                    sum_m = sum_m + (-1)^m1 * a1 * a2 * a3 * cg;
                                end
                            end
                            B(s1,s2,s3,l1+1,l2+1,l3+1) = sum_m;
                        end
                    end
                end
            end
        end
    end
end

function cg = clebsch_gordan(j1, m1, j2, m2, j3, m3)
    if m1 + m2 ~= m3 || j3 < abs(j1 - j2) || j3 > j1 + j2
        cg = 0; return;
    end
    t_min = max([0, j2 - m1 - j3, j1 + m2 - j3]);
    t_max = min([j1 + j2 - j3, j1 - m1, j2 + m2]);
    sum_term = 0;
    for t = t_min:t_max
        denom = factorial(t) * factorial(j1 + j2 - j3 - t) * ...
                factorial(j1 - m1 - t) * factorial(j2 + m2 - t) * ...
                factorial(j3 - j2 + m1 + t) * factorial(j3 - j1 - m2 + t);
        sum_term = sum_term + (-1)^t / denom;
    end
    prefactor = sqrt(...
        (2*j3 + 1) * factorial(j1 + j2 - j3) * factorial(j1 - j2 + j3) * ...
        factorial(-j1 + j2 + j3) / factorial(j1 + j2 + j3 + 1) * ...
        factorial(j3 + m3) * factorial(j3 - m3) * ...
        factorial(j1 + m1) * factorial(j1 - m1) * ...
        factorial(j2 + m2) * factorial(j2 - m2));
    cg = prefactor * sum_term;
end
