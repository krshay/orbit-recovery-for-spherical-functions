close all;
clear variables;

warning('off');
addpath(genpath('../../ASPIRE'))
addpath('../../easyspin-5.2.33/easyspin')
addpath('../../SphericalHarmonics')

rng(10);


% Length of each dimension of the volume.
LL = 240;

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
%
vol_true_downsampled = randn(L, L, L);
vol_true_downsampled = 100 * vol_true_downsampled / norm(vol_true_downsampled, "fro");

% % Volume expansion in 3-D Fourier-Bessel basis
ell_max = 2;
R_max = L / 2;
NUM_SHELLS = 5;
shells = round(linspace(2, R_max, NUM_SHELLS));
r_cut = 1 / 2;
rad_size = floor(L / 2);
% Generate 3-D Fourier-Bessel basis
s_lens = gen_s_list(ell_max, r_cut, rad_size);
lms_list = calc_lms_list(ell_max, s_lens);

[Psilms_NUM_SHELLS, ~, jball_NUM_SHELLS, ~] = precompute_spherical_basis_psilms_NUM_SHELLS(rad_size, ...
    r_cut, ell_max, L, NUM_SHELLS);

[~, vol_trunc_NUM_SHELLS] = expand_vol_spherical_basis(vol_true_downsampled, rad_size, ...
    ell_max, L, Psilms_NUM_SHELLS, jball_NUM_SHELLS);
x_true_NUM_SHELLS = expand_vol_spherical_basis(vol_trunc_NUM_SHELLS, rad_size, ell_max, ...
    L, Psilms_NUM_SHELLS, jball_NUM_SHELLS);

for ell=0:ell_max
    x_true_NUM_SHELLS{ell + 1} = 100 * randn(size(x_true_NUM_SHELLS{ell+1}));
end

vol_init = normrnd(0, 1, [L, L, L]); vol_init = 100 * vol_init / norm(vol_init, "fro");
% vol_init = vol_true_downsampled +  randn(L, L, L);
[~, vol_init_trunc] = expand_vol_spherical_basis(vol_init, rad_size, ...
    ell_max, L, Psilms_NUM_SHELLS, jball_NUM_SHELLS);
vol_init_trunc = real(vol_init_trunc);
vol_init_trunc =  vol_init_trunc / norm(vol_init_trunc, "fro");

x_init = expand_vol_spherical_basis(vol_init_trunc, rad_size, ell_max, ...
    L, Psilms_NUM_SHELLS, jball_NUM_SHELLS);
for ell=0:ell_max
    x_init{ell + 1} = 100 * randn(size(x_true_NUM_SHELLS{ell + 1}));% + x_true_NUM_SHELLS{ell + 1} + ;
end
x_true_vec = vec_cell(x_true_NUM_SHELLS);
x_true_vec = [real(x_true_vec); imag(x_true_vec)];

x_init_vec = vec_cell(x_init);
x_init_vec = [real(x_init_vec); imag(x_init_vec)];
% 
% tic
% B_old = calc_bispectrum(ell_max, NUM_SHELLS, x_true_NUM_SHELLS);
% toc

tic
B = calc_bispectrum_and_grad_real(ell_max, NUM_SHELLS, x_true_vec);
toc


% tic
% B_init = calc_bispectrum(ell_max, NUM_SHELLS, x_init);
% toc
%x_true_vec + 1e-4*randn(size(x_true_vec))
% x_init_vec = x_true_vec;
Aeq = zeros(length(x_true_vec), length(x_true_vec));
for ii=1:(NUM_SHELLS + 3 * NUM_SHELLS)
    Aeq(ii, ii) = 1;
    Aeq(ii + length(x_true_vec) / 2, ii + length(x_true_vec) / 2) = 1;
end
beq = Aeq * x_true_vec;
[z_con, cst_con] = optimize_constrained(x_init_vec, B, ell_max, NUM_SHELLS, Aeq, beq);
N = length(x_true_vec) / 2;
norm(x_true_vec - z_con) / norm(x_true_vec)
norm(x_true_vec - x_init_vec) / norm(x_true_vec)

% 
% [z_unc, cst_unc] = optimize(x_init_vec, B, ell_max, NUM_SHELLS);
