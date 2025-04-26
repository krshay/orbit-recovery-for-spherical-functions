close all;
clear variables;

warning('off');
addpath(genpath('../../ASPIRE'))
addpath('../../easyspin-5.2.33/easyspin')
addpath('../../SphericalHarmonics')

rng(1);


% Length of each dimension of the volume.
LL = 240;

% Length of each dimension of the downsampled volume.
L = 120;

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
vol_true_downsampled = vol_true_downsampled / norm(vol_true_downsampled, "fro");

%% Initializations


%% ell_max = 6
% Volume expansion in 3-D Fourier-Bessel basis
ell_max = 10;
R_max = L / 2;
NUM_SHELLS = 59;
shells = round(linspace(2, R_max, NUM_SHELLS));
r_cut = 1 / 2;
rad_size = floor(L / 2);
% Generate 3-D Fourier-Bessel basis
% [Psilms, Psilms_2D, jball, jball_2D] = precompute_spherical_basis_old(rad_size, ...
%     r_cut, ell_max, L);

[Ylm, jball, jball_2D] = precompute_spherical_basis(rad_size, ...
    r_cut, ell_max, LL);


[X, Y, Z] = meshgrid(-L/2:L/2-1, -L/2:L/2-1, -L/2:L/2-1);
xs = cell(NUM_SHELLS, 1);
vol_truncs = cell(NUM_SHELLS, 1);

vol_spherical_harmonics = zeros(L, L, L);
for i=1:NUM_SHELLS
    R_curr = shells(i);
    tmp = (sqrt((X.^2 + Y.^2 + Z.^2)) <= R_curr + 0.01) & (sqrt((X.^2 + Y.^2 + Z.^2)) >= R_curr - 0.01);
    shell_curr = zeros(L, L, L);
    for ii=1:L
        for jj=1:L
            for kk=1:L
                if tmp(ii, jj, kk)
                    shell_curr(ii, jj, kk) = vol_true_downsampled(ii, jj, kk);
                end
            end
        end
    end
    shell_curr = shell_curr(round(L/2) - R_curr + 1:round(L/2) + R_curr, round(L/2) - R_curr + 1:round(L/2) + R_curr, round(L/2) - R_curr + 1:round(L/2) + R_curr);
    rad_size = R_curr;
    [Ylm, jball, jball_2D] = precompute_spherical_basis(rad_size, ...
        r_cut, ell_max, L);
    [~, vol_trunc] = expand_vol_spherical_basis(shell_curr, rad_size, ...
        ell_max, mod(rad_size, 2) + 1, Ylm, jball);
    max(imag(vol_trunc(:)))
    vol_trunc = real(vol_trunc);
    vol_truncs{i} = vol_trunc;
    xs{i} = expand_vol_spherical_basis(vol_trunc, rad_size, ell_max, ...
        L, Ylm, jball);

end


x = cell(ell_max + 1, 1);
for ell=1:ell_max+1
    x{ell} = xs{2}{ell} * 2;
end
for i=3:NUM_SHELLS
    for ell=1:ell_max+1
        x{ell} = x{ell} + xs{i}{ell} * i;
    end
end
rad_size = floor(L / 2);
vol_sh = expand_vol_psilms(x, L/2, jball, Ylm, mod(rad_size, 2) + 1);
vol_sh = real(vol_sh);