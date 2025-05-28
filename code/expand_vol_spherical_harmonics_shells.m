function [x_true_NUM_SHELLS, vol_trunc_NUM_SHELLS, Psilms_NUM_SHELLS, jball_NUM_SHELLS, rad_size] = expand_vol_spherical_harmonics_shells(vol, L, ell_max, num_shells)
% Volume expansion in 3-D Fourier-Bessel basis
r_cut = 1 / 2;
rad_size = floor(L / 2);
% Generate 3-D Fourier-Bessel basis
s_lens = gen_s_list(ell_max, r_cut, rad_size);
lms_list = calc_lms_list(ell_max, s_lens);

[Psilms_NUM_SHELLS, ~, jball_NUM_SHELLS, ~] = precompute_spherical_basis_psilms_NUM_SHELLS(rad_size, ...
    r_cut, ell_max, L, num_shells);

[~, vol_trunc_NUM_SHELLS] = expand_vol_spherical_basis(vol, rad_size, ...
    ell_max, L, Psilms_NUM_SHELLS, jball_NUM_SHELLS);
x_true_NUM_SHELLS = expand_vol_spherical_basis(vol_trunc_NUM_SHELLS, rad_size, ell_max, ...
    L, Psilms_NUM_SHELLS, jball_NUM_SHELLS);
vol_trunc_NUM_SHELLS = expand_vol_psilms(x_true_NUM_SHELLS, rad_size, jball_NUM_SHELLS, Psilms_NUM_SHELLS, L);
end