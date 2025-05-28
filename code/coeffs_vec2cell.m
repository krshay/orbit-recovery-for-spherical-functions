function x_cell = coeffs_vec2cell(x_vec, ell_max)
x_cell = cell(ell_max + 1, 1);
current_idx = 0;
for ell=0:ell_max
    num_m = 2 * ell + 1;
    x_cell{ell + 1} = x_vec(:, current_idx + 1:current_idx + num_m);
    current_idx = current_idx + num_m;
end
end