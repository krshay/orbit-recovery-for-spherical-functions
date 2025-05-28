function [Ylm, jball, jball_2D, jjorigin] = ...
    precompute_spherical_basis(N, r_cut, ell_max, L)

%  Precompute spherical harmonics.
[jball, jjorigin, Y, ~, dels_Y, ~, jball_2D] = ...
    prepare_psilms(N, r_cut, ell_max, L);

Ylm = cell(ell_max + 1, 1);
for ii=0:ell_max
    Ylm{ii+1, 1} = generate_rYlm(jball, Y, dels_Y,...
        N, ii, L);
end

end
