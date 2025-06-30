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
noise_var = 1;   % variance  σ²  of the complex Gaussian noise
num_avg   = 1;    % number of noisy realisations to average
vol_inits = zeros(L, L, L, n_trials);
for trial = 1:n_trials
    tmp = randn(L,L,L);
    vol_inits( :, :, :, trial) = 100 * tmp / norm(tmp, "fro");
end

ell_max = 10;
shell_values = 8;
errors = zeros(size(shell_values));
std_errors = zeros(size(shell_values));
rel_errors = zeros(length(shell_values), n_trials);
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
tic;  % Start timing
     
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
    % add complex Gaussian noise ( √σ²·(N(0,1)+iN(0,1)) )
    noise_term = sqrt(noise_var) * (randn(size(true_coeffs)) + 1i*randn(size(true_coeffs)));
    noisy_coeffs = true_coeffs + 0*norm(true_coeffs)*noise_term/norm(noise_term);

    % bispectrum of this noisy realisation
    B_k = compute_SO3_bispectrum_restricted(noisy_coeffs, ell_max, ell, shell, cg_tab,shift);
    % accumulate
    B_target = B_target + B_k;
end

% average
B_target = B_target / num_avg;

                [A, b] = construct_linear_system(est_coeffs, shell, current_indices, ell, ell_max, B_target, cg_tab,shift);
               cond_number = cond(A);                      % condition of current system

% -- record only once per band ℓ (take the first shell) --
if shell == 1 && trial == 1
    row = idx;                              % 1→3 shells, 2→4 shells, ...
    col = ell - ell_start + 1;              % ℓ=2 → col 1, ℓ=3 → col 2 …
    cond_band(row, col) = cond_number;
end


                est_coeffs(shell, current_indices) = A \ b;
            end
        end
elapsedTime = toc;  % End timing

        rel_errors(idx, trial) = norm(est_coeffs - true_coeffs, "fro") / norm(true_coeffs, "fro");
        fprintf('Trial %d: Relative error = %.4f \n', trial, rel_errors(idx, trial));
    end

    avg_errors(idx) = mean(rel_errors(idx, :));
    std_errors(idx) = std(rel_errors(idx, :));
    fprintf('\nAverage error for %d shells: %.4f (std = %.4f) \n', num_shells, avg_errors(idx), std_errors(idx));
    
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
WriteMRC(vol_rec, 1, "vol_TRPV1_true.mrc");
volumeViewer(real(vol_trunc_NUM_SHELLS));
vol_80S_true = real(vol_trunc_NUM_SHELLS);
save("vol_80S_est.mat", "vol_rec");
save("vol_80S_true.mat", "vol_80S_true");
save("vol_80S_init.mat", "vol_init_trunc_NUM_SHELLS");
WriteMRC(vol_rec, 1, "vol_TRPV1_est.mrc");
volumeViewer(vol_rec);

fprintf('Total runtime: %.2f seconds\n', elapsedTime);


% === Helper functions ===
function idx = sh_index(l, m)
    idx = l^2 + m + l + 1;
end

function cg_tab = precompute_cg_table(d)
    shift = d;                              % max |m|
    cg_tab = zeros(d+1,2*shift+1,d+1,2*shift+1,d+1);  % l1,m1,l2,m2,l3 (m3 implicit)
    for l1 = 0:d
        for m1 = -l1:l1
            for l2 = 0:d
                for m2 = -l2:l2
                    for l3 = abs(l1-l2):min(l1+l2,d)
                        m3 = -m1-m2;
                        if abs(m3)>l3, continue; end
                        cg = clebsch_gordan(l2,m2,l3,m3,l1,-m1);
                        cg_tab(l1+1,m1+shift+1,l2+1,m2+shift+1,l3+1) = cg;
                    end
                end
            end
        end
    end
end


function [A,b] = construct_linear_system(est_coeffs,shell_idx,curr_idx,ell,d,B_target,cg_tab,shift)
    num_shells = size(est_coeffs,1); num_coeffs = numel(curr_idx);
    A = []; b = [];
    for l1 = 0:d
      for l2 = 0:d
        for l3 = abs(l1-l2):min(l1+l2,d)
          lvals = [l1,l2,l3]; idx_ell = find(lvals==ell,1); if isempty(idx_ell), continue; end
          if any(lvals(setdiff(1:3,idx_ell)) >= ell), continue; end

          for s1 = 1:num_shells
            for s2 = 1:num_shells
              for s3 = 1:num_shells
                shells = [s1,s2,s3]; if shells(idx_ell)~=shell_idx, continue; end
                coeff_row = zeros(1,num_coeffs); sum_known = 0;

                for m1=-l1:l1
                  for m2=-l2:l2
                    m3=-m1-m2; if abs(m3)>l3, continue; end
                    idx1=sh_index(l1,m1); idx2=sh_index(l2,m2); idx3=sh_index(l3,m3);
                    cg = cg_tab(l1+1,m1+shift+1,l2+1,m2+shift+1,l3+1);

                    idxs=[idx1,idx2,idx3]; sidx=shells;
                    if ismember(idxs(idx_ell),curr_idx)
                        k = find(curr_idx==idxs(idx_ell));
                        val1 = est_coeffs(sidx(mod(idx_ell,3)+1),idxs(mod(idx_ell,3)+1));
                        val2 = est_coeffs(sidx(mod(idx_ell+1,3)+1),idxs(mod(idx_ell+1,3)+1));
                        coeff_row(k) = coeff_row(k) + (-1)^m1*val1*val2*cg;
                    else
                        a1=est_coeffs(sidx(1),idx1); a2=est_coeffs(sidx(2),idx2); a3=est_coeffs(sidx(3),idx3);
                        sum_known = sum_known + (-1)^m1*a1*a2*a3*cg;
                    end
                  end
                end
                A(end+1,:) = coeff_row; %#ok<AGROW>
                b(end+1,1) = B_target(s1,s2,s3,l1+1,l2+1,l3+1) - sum_known;
              end
            end
          end
        end
      end
    end
end

function B = compute_SO3_bispectrum_restricted(a_lm,d,ell,shell_idx,cg_tab,shift)
    num_shells=size(a_lm,1);
    B=zeros(num_shells,num_shells,num_shells,d+1,d+1,d+1);
    for l1=0:d
      for l2=0:d
        for l3=abs(l1-l2):min(l1+l2,d)
          lvals=[l1,l2,l3]; idx_ell=find(lvals==ell,1); if isempty(idx_ell),continue; end
          for s1=1:num_shells,for s2=1:num_shells,for s3=1:num_shells
            shells=[s1,s2,s3]; if shells(idx_ell)~=shell_idx,continue; end
            acc=0;
            for m1=-l1:l1,for m2=-l2:l2
              m3=-m1-m2; if abs(m3)>l3, continue; end
              cg=cg_tab(l1+1,m1+shift+1,l2+1,m2+shift+1,l3+1);
              a1=a_lm(s1,sh_index(l1,m1));
              a2=a_lm(s2,sh_index(l2,m2));
              a3=a_lm(s3,sh_index(l3,m3));
              acc = acc + (-1)^m1 * a1*a2*a3 * cg;
            end,end
            B(s1,s2,s3,l1+1,l2+1,l3+1)=acc;
          end,end,end
    end,end,end
end

function cg = clebsch_gordan(j1,m1,j2,m2,j3,m3)
    if m1+m2~=m3 || j3<abs(j1-j2) || j3>j1+j2, cg=0; return; end
    tmin=max([0 j2-m1-j3 j1+m2-j3]);
    tmax=min([j1+j2-j3 j1-m1 j2+m2]);
    sumTerm=0;
    for t=tmin:tmax
        sumTerm=sumTerm+ (-1)^t / ...
          (factorial(t)*factorial(j1+j2-j3-t)*factorial(j1-m1-t)* ...
           factorial(j2+m2-t)*factorial(j3-j2+m1+t)*factorial(j3-j1-m2+t));
    end
    pref=sqrt((2*j3+1)*factorial(j1+j2-j3)*factorial(j1-j2+j3)* ...
              factorial(-j1+j2+j3)/factorial(j1+j2+j3+1)* ...
              factorial(j3+m3)*factorial(j3-m3)*factorial(j1+m1)* ...
              factorial(j1-m1)*factorial(j2+m2)*factorial(j2-m2));
    cg=pref*sumTerm;
end
