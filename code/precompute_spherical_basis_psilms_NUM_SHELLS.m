function [Psilms, Psilms_2D, jball, jball_2D, jjorigin] = ...
    precompute_spherical_basis_psilms_NUM_SHELLS(N, r_cut, ell_max, L, num_shells)

%  Precompute spherical harmonics.
[jball, jjorigin, Y, dems_Y, dels_Y, xyplane, jball_2D] = ...
    prepare_psilms(N, r_cut, ell_max, L);

r_select_ratio = 1; % auxiliary factor in sampling criterion

Psilms = cell(ell_max + 1, 1);
for ii=0:ell_max
    Psilms{ii+1} = generate_psilms_NUM_SHELLS(jball, jjorigin, Y, dems_Y, dels_Y,...
        r_cut, N, r_select_ratio, ii, ell_max, L, num_shells);
end

Psilms_2D = cellfun(@(x) x(xyplane, :), Psilms, 'UniformOutput', 0);
end
