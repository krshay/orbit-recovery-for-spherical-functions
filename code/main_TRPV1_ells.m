close all;
clear variables;

warning('off');
addpath(genpath('../../ASPIRE'))
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

% Length of each dimension of the downsampled volume.
L = 21;

%% Volume generation
vol = load("../data/TRPV1.mat");
vol = vol.TRPV1_vol;

LL = size(vol, 1);

vol = LL * vol / sqrt(sum(vol.^2, 'all'));

[X, Y, Z] = meshgrid(-LL/2:LL/2-1, -LL/2:LL/2-1, -LL/2:LL/2-1);
tmp = sqrt((X.^2 + Y.^2 + Z.^2)) <= LL/2;
vol_cropped = zeros(size(vol));
for i=1:LL
    for j=1:LL
        for k=1:LL
            if tmp(i, j, k)
                vol_cropped(i, j, k) = vol(i, j, k);
            end
        end
    end
end
vol = vol_cropped;
vol_true_downsampled = cryo_downsample(vol, L);
vol_true_downsampled = 100 * vol_true_downsampled / norm(vol_true_downsampled, "fro");


%% Volume reconstruction
rng(1);
n_trials = 5;
vol_inits = zeros(L, L, L, n_trials);
for trial = 1:n_trials
    tmp = normrnd(0, 1, [L, L, L]);
    vol_inits( :, :, :, trial) = 100 * tmp / norm(tmp, "fro");
end
% --- Parameters ---
ell_maxs = 2:7;  % Bandlimit
num_shells = 5;
avg_errors = zeros(size(ell_maxs));
std_errors = zeros(size(ell_maxs));

    rel_errors = zeros(length(ell_maxs), n_trials);
for idx = 1:length(ell_maxs)
    ell_max = ell_maxs(idx);
    rel_errors = zeros(1, n_trials);
    % Volume expansion in 3-D Fourier-Bessel basis

    [x_true_NUM_SHELLS, vol_trunc_NUM_SHELLS, Psilms_NUM_SHELLS, jball_NUM_SHELLS, rad_size] = expand_vol_spherical_harmonics_shells(vol_true_downsampled, L, ell_max, num_shells);

    % Volume expansion in real spherical harmonics
    x_true_NUM_SHELLS_real = complex2real_spherical_harmonics(x_true_NUM_SHELLS, ell_max);

    true_coeffs = [];
    for ell=0:ell_max
        true_coeffs = [true_coeffs x_true_NUM_SHELLS_real{ell + 1}];
    end
    fprintf('\n=== Running for maximum frequency %d ===\n', ell_max);
    for trial = 1:n_trials
        vol_init = vol_inits( :, :, :, trial);
        % Volume expansion in 3-D Fourier-Bessel basis
        [x_init_NUM_SHELLS, vol_init_trunc_NUM_SHELLS] = expand_vol_spherical_harmonics_shells(vol_init, L, ell_max, num_shells);

        % Volume expansion in real spherical harmonics
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

        % Frequency marching
        for ell = 2:ell_max
            current_indices = [];
            for m = -ell:ell
                current_indices = [current_indices, sh_index(ell, m)];
            end

            for shell = 1:num_shells
                x_init = init_coeffs(shell, current_indices);
                obj = @(x_shell) specialized_loss_and_fd_grad(x_shell, est_coeffs, ...
                    shell, current_indices, ell, ell_max, true_coeffs);

                [x_sol, ~] = fminunc(obj, x_init, options);
                est_coeffs(shell, current_indices) = x_sol;
            end
        end
        rel_errors(idx, trial) = norm(est_coeffs - true_coeffs, "fro") / norm(true_coeffs, "fro");
        fprintf('Trial %d: Relative error = %.4f%%\n', trial, 100 * rel_errors(trial));
    end

    avg_errors(idx) = mean(rel_errors(idx, :));
    std_errors(idx) = std(rel_errors(idx, :));
    fprintf('\nAverage error for maximum frequency %d: %.4f%% (std = %.4f%%)\n', ...
        ell_max, 100 * avg_errors(idx), 100 * std_errors(idx));
end
save("rel_errors_TRPV1_ells_shells5.mat", "rel_errors")
figure;
errorbar(ell_maxs, avg_errors, std_errors, 'o-', 'LineWidth', 2);
xlabel('Maximum Frequency');
ylabel('Average Relative Error');
title('Frequency Marching Recovery Error vs. Maximum Frequency');
grid on;


x_est_real = coeffs_vec2cell(est_coeffs, ell_max);
x_est = real2complex_spherical_harmonics(x_est_real, ell_max, num_shells);
vol_rec = real(expand_vol_psilms(x_est, rad_size, jball_NUM_SHELLS, Psilms_NUM_SHELLS, L));
volumeViewer(real(vol_trunc_NUM_SHELLS))

volumeViewer(vol_rec)
