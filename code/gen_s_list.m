function s_list = gen_s_list(maxL, r_cut, N)
% Generate list with number of radial frequencies for volume expansion for
% each L-subspace of spherical harmonics
% 
% Inputs:
%   * maxL: cutoff for spherical harmonics expansion
%   * r_cut: assumed bandlimit for the volume (and the projections)
%   * N: radius of the support of the volume (N = floor(L/2))
% 
% Outputs:
%   * s_list: list of number of radial frequencies for each order of
%   spherical harmonics
% 
% Eitan Levin, August 2018

% Load table of Bessel zeros
fname='../data/Besselj_L200_S500.mat';
load(fname, 'B') 

r_select_ratio = 1;
 
B = B( B(:, 3)< 2*pi*N*r_cut * r_select_ratio & B(:,1) <= maxL, :); %Nyquist criterion
ell_grid = B(:, 1);

s_list = zeros(max(ell_grid)+1,1);
for l = 0:max(ell_grid)
    s_list(l+1) = length(find(ell_grid == l));
end
end
