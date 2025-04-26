function [f, g] = cost_and_grad_bispectrum(x_vec, B, ell_max, num_shells)
%cost_Bi Calculates the cost function for the optimization
%using the Bispectrum; eq. (4.2)

%   z is the currentf guess for the fft of the signal
%   y is the measurement
%   A is the sensing matrix
%   k1k2k3_map is a N^2*3 matrix containing the triplets of frequencies

% x_cell = cell(ell_max + 1, 1);
% s_lens = num_shells * ones(ell_max + 1, 1);
N = length(x_vec);
% x_cell = cell_vec(x_vec(1:N / 2) + 1j * x_vec(N / 2 + 1:end), x_cell, ell_max, s_lens);
% B_x = calc_bispectrum(ell_max, num_shells, x_cell);
% [B_x_vec, gB_x_vec] = calc_bispectrum_and_grad_vec(ell_max, num_shells, x_vec(1:N / 2) + 1j * x_vec(N / 2 + 1:end));
[B_x, gB_x] = calc_bispectrum_and_grad_real(ell_max, num_shells, x_vec);

f = norm(B - B_x, "fro") ^ 2;
% g = - 2 * reshape(gB_x_vec, N / 2, []) * (B(:) - B_x_vec(:));
g = - 2 * reshape(real(gB_x), N, []) * real(B(:) - B_x(:)) - 2 * reshape(imag(gB_x), N, []) * imag(B(:) - B_x(:));

% g = [real(g); imag(g)];
% display("Gradient norm");
% norm(g, "fro")
end