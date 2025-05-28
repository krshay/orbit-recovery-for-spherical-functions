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
n_trials = 1;
vol_inits = zeros(L, L, L, n_trials);
for trial = 1:n_trials
    tmp = normrnd(0, 1, [L, L, L]);
    vol_inits( :, :, :, trial) = 100 * tmp / norm(tmp, "fro");
end
% --- Parameters ---
ell_max = 5;  % Bandlimit
shell_values = 1:5;
avg_errors = zeros(size(shell_values));
std_errors = zeros(size(shell_values));

rel_errors = zeros(length(shell_values), n_trials);
rotated_errors = zeros(length(shell_values), n_trials, size(angles, 2));

for idx = 1:length(shell_values)
    num_shells = shell_values(idx);
    % Volume expansion in 3-D Fourier-Bessel basis

    [x_true_NUM_SHELLS, vol_trunc_NUM_SHELLS, Psilms_NUM_SHELLS, jball_NUM_SHELLS, rad_size] = expand_vol_spherical_harmonics_shells(vol_true_downsampled, L, ell_max, num_shells);

    % Volume expansion in real spherical harmonics
    x_true_NUM_SHELLS_real = complex2real_spherical_harmonics(x_true_NUM_SHELLS, ell_max);

    true_coeffs = [];
    for ell=0:ell_max
        true_coeffs = [true_coeffs x_true_NUM_SHELLS_real{ell + 1}];
    end
    fprintf('\n=== Running for %d shells ===\n', num_shells);
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
            fprintf('Trial %d: Relative error = %.4f \n', trial, rel_errors(idx, trial));
    end

    avg_errors(idx) = mean(rel_errors(idx, :));
    std_errors(idx) = std(rel_errors(idx, :));
    fprintf('\nAverage error for %d shells: %.4f (std = %.4f) \n', ...
        num_shells, avg_errors(idx), std_errors(idx));
end
save("rel_errors_TRPV1_shells_ell5.mat", "rel_errors")
figure;
errorbar(shell_values, avg_errors, std_errors, 'o-', 'LineWidth', 2);
xlabel('Number of Shells');
ylabel('Average Relative Error');
title('Frequency Marching Recovery Error vs. Number of Shells');
grid on;


x_est_real = coeffs_vec2cell(x_rotated, ell_max);
x_est = real2complex_spherical_harmonics(x_est_real, ell_max, num_shells);
vol_rec = real(expand_vol_psilms(x_est, rad_size, jball_NUM_SHELLS, Psilms_NUM_SHELLS, L));
volumeViewer(real(vol_trunc_NUM_SHELLS))
vol_TRPV1_true = real(vol_trunc_NUM_SHELLS);
save("vol_TRPV1_est.mat", "vol_rec");
save("vol_TRPV1_true.mat", "vol_TRPV1_true");
save("vol_TRPV1_init.mat", "vol_init_trunc_NUM_SHELLS");
volumeViewer(vol_rec)
