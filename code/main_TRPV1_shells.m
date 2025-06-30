close all;
clear variables;

warning('off');
addpath(genpath('~/Documents/MATLAB/aspire'))
addpath('../../easyspin-5.2.33/easyspin')
addpath('../../SphericalHarmonics')

L = 31;


%% Volume generation
vol = load('../data/TRPV1.mat');
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
vol_true_downsampled = vol_true_downsampled / norm(vol_true_downsampled, "fro");

%% Volume reconstruction
rng(1);
noise_var = 0.5;   % variance  σ²  of the complex Gaussian noise
num_avg   = 500;    % number of noisy realisations to average

ell_max = 10;
shell_values = 3:8;
errors = zeros(size(shell_values));
ell_start   = 2;                            % first band you solve
n_bands     = ell_max - ell_start + 1;
n_shellsets = numel(shell_values);

cond_band   = NaN(n_shellsets, n_bands);    % rows = shell counts, cols = ℓ
% -------------------------------------------------
% 1) PRE-COMPUTE Clebsch–Gordan coefficients
% -------------------------------------------------
cg_tab = precompute_cg_table(ell_max);
%   cg_tab(l1+1 , m1+shift , l2+1 , m2+shift , l3+1) = C^{l1,m1}_{l2,m2,l3,-m1-m2}
shift  = ell_max;   % max |m|  (used for index shift)

% -------------------------------------------------
% 2) Plot preparation
% -------------------------------------------------
colors = lines(length(noise_var));
figure; hold on;
noisy_volumes = zeros(L, L, L, num_avg);
for k = 1:num_avg
    % add Gaussian noise ( √σ²·(N(0,1) )
    tmp = randn(size(vol_true_downsampled));
    tmp = tmp / norm(tmp, 'fro');
    noise_term = sqrt(noise_var) * tmp;
    noisy_volumes(:, :, :, k) = vol_true_downsampled + noise_term;
end

noisy_coeffs = cell(num_avg, 1);

for idx = 1:length(shell_values)
    num_shells = shell_values(idx);

    [x_true_NUM_SHELLS, vol_trunc_NUM_SHELLS, Psilms_NUM_SHELLS, jball_NUM_SHELLS, rad_size] = expand_vol_spherical_harmonics_shells(vol_true_downsampled, L, ell_max, num_shells);
    x_true_NUM_SHELLS_real = complex2real_spherical_harmonics(x_true_NUM_SHELLS, ell_max);
    for k = 1:num_avg
        noisy_coeffs_complex = expand_vol_spherical_harmonics_shells(noisy_volumes( :, :, :, k), L, ell_max, num_shells);
        noisy_coeffs_real = complex2real_spherical_harmonics(noisy_coeffs_complex, ell_max);
        noisy_coeffs_k = [];
        for ell=0:ell_max
            noisy_coeffs_k = [noisy_coeffs_k noisy_coeffs_real{ell + 1}];
        end
        noisy_coeffs{k} = noisy_coeffs_k;
    end
    true_coeffs = [];
    for ell=0:ell_max
        true_coeffs = [true_coeffs x_true_NUM_SHELLS_real{ell + 1}];
    end
    fprintf('\n=== Running for %d shells ===\n', num_shells);


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
            % ------------------------------------------------------------
            % NEW: build B_target by averaging over 'num_avg' noisy samples
            % ------------------------------------------------------------
            B_target = zeros(num_shells,num_shells,num_shells,ell_max+1,ell_max+1,ell_max+1);

            for k = 1:num_avg
                % bispectrum of this noisy realisation
                B_k = compute_SO3_bispectrum_restricted(noisy_coeffs{k}, ell_max, ell, shell, cg_tab,shift);
                % accumulate
                B_target = B_target + B_k;
            end

            % average
            B_target = B_target / num_avg;

            [A, b] = construct_linear_system(est_coeffs, shell, current_indices, ell, ell_max, B_target, cg_tab,shift);
            cond_number = cond(A);                      % condition of current system

            % -- record only once per band ℓ (take the first shell) --
            if shell == 1
                row = idx;                              % 1→3 shells, 2→4 shells, ...
                col = ell - ell_start + 1;              % ℓ=2 → col 1, ℓ=3 → col 2 …
                cond_band(row, col) = cond_number;
            end


            est_coeffs(shell, current_indices) = A \ b;
        end
    end
    rel_errors(idx) = norm(est_coeffs - true_coeffs, "fro") / norm(true_coeffs, "fro");
    fprintf('Relative error = %.4f \n', rel_errors(idx));

end

% ========  Condition-number table  =====================================
row_headers = arrayfun(@(s) sprintf('%d shells', s), shell_values, 'uni', false);
col_headers = arrayfun(@(ell) sprintf('ell = %d', ell), ell_start:ell_max, 'uni', false);

fprintf('\nCondition-number table (first shell of each band):\n');
fprintf('%12s', '');    fprintf('%12s', col_headers{:});  fprintf('\n');
for r = 1:n_shellsets
    fprintf('%12s', row_headers{r});
    for c = 1:n_bands
        fprintf('%12.3e', cond_band(r,c));
    end
    fprintf('\n');
end

% ---------- optional: save CSV for LaTeX / Excel ------------------------
fid = fopen('cond_table.csv','w');
fprintf(fid, ',%s\n', strjoin(col_headers, ','));
for r = 1:n_shellsets
    fprintf(fid, '%s,%s\n', row_headers{r}, ...
        strjoin(cellfun(@(x) sprintf('%.3e', x), ...
        num2cell(cond_band(r,:)), 'uni', false), ','));
end
fclose(fid);
disp('Saved  cond_table.csv');



save("rel_errors_80S_shells_ell10.mat", "rel_errors")
figure;
plot(shell_values, rel_errors, 'o-', 'LineWidth', 2);
xlabel('Number of Shells');
ylabel('Relative Error');
xticks([3 4 5 6 7 8])        % Set the tick positions
xticklabels({'3','4','5','6','7','8'}) 

grid on;

x_est_real = coeffs_vec2cell(est_coeffs, ell_max);
x_est = real2complex_spherical_harmonics(x_est_real, ell_max, num_shells);
vol_rec = real(expand_vol_psilms(x_est, rad_size, jball_NUM_SHELLS, Psilms_NUM_SHELLS, L));
vol_80S_true = real(vol_trunc_NUM_SHELLS);