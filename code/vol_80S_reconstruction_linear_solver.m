close all;
clear variables;

warning('off');
addpath(genpath('~/Documents/MATLAB/aspire'))
addpath('../../easyspin-5.2.33/easyspin')
addpath('../../SphericalHarmonics')

options = optimoptions('fminunc', ...
             'Algorithm', 'quasi-newton', ...
            'Display', 'off', ...
            'MaxIterations', 500, ...
            'MaxFunctionEvaluations', 20000, ...
            'OptimalityTolerance', 1e-12, ...
            'StepTolerance', 1e-12, ...
            'FunctionTolerance', 1e-10, ...
            'SpecifyObjectiveGradient', true);

L = 31;


%% Volume generation
vol = load('../data/emd_2660.mat');
vol = vol.emd_2660;
vol = vol(50:300, 50:300, 50:300);

vol_true_downsampled = cryo_downsample(vol, L);
vol_true_downsampled = 100 * vol_true_downsampled / norm(vol_true_downsampled, "fro");
%% Volume reconstruction
rng(1);
n_trials = 1;
vol_inits = zeros(L, L, L, n_trials);
for trial = 1:n_trials
    tmp = randn(L,L,L);
    vol_inits( :, :, :, trial) = 100 * tmp / norm(tmp, "fro");
end

ell_max = 10;
shell_values = 8;
avg_errors = zeros(size(shell_values));
std_errors = zeros(size(shell_values));
rel_errors = zeros(length(shell_values), n_trials);

for idx = 1:length(shell_values)
    num_shells = shell_values(idx);

    [x_true_NUM_SHELLS, vol_trunc_NUM_SHELLS, Psilms_NUM_SHELLS, jball_NUM_SHELLS, rad_size] = expand_vol_spherical_harmonics_shells(vol_true_downsampled, L, ell_max, num_shells);
    x_true_NUM_SHELLS_real = complex2real_spherical_harmonics(x_true_NUM_SHELLS, ell_max);

    true_coeffs = [];
    for ell=0:ell_max
        true_coeffs = [true_coeffs x_true_NUM_SHELLS_real{ell + 1}];
    end
    fprintf('\n=== Running for %d shells ===\n', num_shells);

    for trial = 1:n_trials
        vol_init = vol_inits( :, :, :, trial);
        [x_init_NUM_SHELLS, vol_init_trunc_NUM_SHELLS] = expand_vol_spherical_harmonics_shells(vol_init, L, ell_max, num_shells);
        x_init_NUM_SHELLS_real = complex2real_spherical_harmonics(x_init_NUM_SHELLS, ell_max);

        init_coeffs = [];
        for ell=0:ell_max
            init_coeffs = [init_coeffs x_init_NUM_SHELLS_real{ell + 1}];
        end

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

        
        for ell = 2:ell_max
            current_indices = [];
            for m = -ell:ell
                current_indices = [current_indices, sh_index(ell, m)];
            end

            for shell = 1:num_shells
                B_target = compute_SO3_bispectrum_restricted(true_coeffs, ell_max, ell, shell);
                [A, b] = construct_linear_system(est_coeffs, shell, current_indices, ell, ell_max, B_target);
                est_coeffs(shell, current_indices) = A \ b;
            end
        end

        rel_errors(idx, trial) = norm(est_coeffs - true_coeffs, "fro") / norm(true_coeffs, "fro");
        fprintf('Trial %d: Relative error = %.4f \n', trial, rel_errors(idx, trial));
    end

    avg_errors(idx) = mean(rel_errors(idx, :));
    std_errors(idx) = std(rel_errors(idx, :));
    fprintf('\nAverage error for %d shells: %.4f (std = %.4f) \n', num_shells, avg_errors(idx), std_errors(idx));
end

save("rel_errors_80S_shells_ell5.mat", "rel_errors")
figure;
errorbar(shell_values, avg_errors, std_errors, 'o-', 'LineWidth', 2);
xlabel('Number of Shells');
ylabel('Average Relative Error');
title('Frequency Marching Recovery Error vs. Number of Shells');
grid on;

x_est_real = coeffs_vec2cell(est_coeffs, ell_max);
x_est = real2complex_spherical_harmonics(x_est_real, ell_max, num_shells);
vol_rec = real(expand_vol_psilms(x_est, rad_size, jball_NUM_SHELLS, Psilms_NUM_SHELLS, L));
figure('Name', 'True Volume (Ground Truth)');
volumeViewer(real(vol_trunc_NUM_SHELLS));
vol_80S_true = real(vol_trunc_NUM_SHELLS);
save("vol_80S_est.mat", "vol_rec");
save("vol_80S_true.mat", "vol_80S_true");
figure('Name', 'Reconstructed Volume');
volumeViewer(vol_rec);

WriteMRC(vol_rec, 1, "vol_80S_est.mrc");
WriteMRC(vol_rec, 1, "vol_80S_true.mrc");

% === Helper functions ===
function idx = sh_index(l, m)
    idx = l^2 + m + l + 1;
end

function [A, b] = construct_linear_system(est_coeffs, shell_idx, current_indices, ell, d, B_target)
    num_shells = size(est_coeffs, 1);
    num_coeffs = length(current_indices);
    max_rows = 10000;
    A = zeros(max_rows, num_coeffs);
    b = zeros(max_rows, 1);
    row = 1;

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
                            coeff_row = zeros(1, num_coeffs);
                            sum_known = 0;
                            for m1 = -l1:l1
                                for m2 = -l2:l2
                                    m3 = -m1 - m2;
                                    if abs(m3) > l3, continue; end
                                    idx1 = sh_index(l1, m1);
                                    idx2 = sh_index(l2, m2);
                                    idx3 = sh_index(l3, m3);
                                    cg = clebsch_gordan(l2, m2, l3, m3, l1, -m1);
                                    idxs = [idx1, idx2, idx3];
                                    s_idx = shells;
                                    if ismember(idxs(idx_ell), current_indices)
                                        i = find(current_indices == idxs(idx_ell));
                                        val1 = est_coeffs(s_idx(mod(idx_ell,3)+1), idxs(mod(idx_ell,3)+1));
                                        val2 = est_coeffs(s_idx(mod(idx_ell+1,3)+1), idxs(mod(idx_ell+1,3)+1));
                                        coeff_row(i) = coeff_row(i) + (-1)^m1 * val1 * val2 * cg;
                                    else
                                        a1 = est_coeffs(s_idx(1), idx1);
                                        a2 = est_coeffs(s_idx(2), idx2);
                                        a3 = est_coeffs(s_idx(3), idx3);
                                        sum_known = sum_known + (-1)^m1 * a1 * a2 * a3 * cg;
                                    end
                                end
                            end
                            A(row, :) = coeff_row;
                            b(row) = B_target(s1, s2, s3, l1+1, l2+1, l3+1) - sum_known;
                            row = row + 1;
                        end
                    end
                end
            end
        end
    end

    A = A(1:row-1, :);
    b = b(1:row-1);
end

function B = compute_SO3_bispectrum_restricted(a_lm_matrix, d, ell, shell_idx)
    num_shells = size(a_lm_matrix, 1);
    B = zeros(num_shells, num_shells, num_shells, d+1, d+1, d+1);
    for l1 = 0:d
        for l2 = 0:d
            for l3 = abs(l1 - l2):min(l1 + l2, d)
                lvals = [l1, l2, l3];
                idx_ell = find(lvals == ell);
                if length(idx_ell) ~= 1
                    continue;
                end
                for s1 = 1:num_shells
                    for s2 = 1:num_shells
                        for s3 = 1:num_shells
                            shells = [s1, s2, s3];
                            if shells(idx_ell) ~= shell_idx
                                continue;
                            end
                            sum_m = 0;
                            for m1 = -l1:l1
                                for m2 = -l2:l2
                                    m3 = -m1 - m2;
                                    if abs(m3) > l3, continue; end
                                    idx1 = sh_index(l1, m1);
                                    idx2 = sh_index(l2, m2);
                                    idx3 = sh_index(l3, m3);
                                    cg = clebsch_gordan(l2, m2, l3, m3, l1, -m1);
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
