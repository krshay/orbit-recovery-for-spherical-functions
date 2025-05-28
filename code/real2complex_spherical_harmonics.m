function x_complex = real2complex_spherical_harmonics(x_real, ell_max, num_shells)
x_complex = cell(size(x_real));
for ell=0:ell_max
    x_complex{ell+1} = zeros(num_shells, 2 * ell + 1);
    for m=-ell:ell
        if m < 0
            x_complex{ell+1}(:, m + ell + 1) = (1 / sqrt(2)) * (x_real{ell + 1}(:, abs(m) + ell + 1) - 1j * x_real{ell + 1}(:, -abs(m) + ell + 1));
        elseif m > 0
            x_complex{ell+1}(:, m + ell + 1) = (((-1) ^ m) / sqrt(2)) * (x_real{ell + 1}(:, abs(m) + ell + 1) + x_real{ell + 1}(:, -abs(m) + ell + 1));
        elseif m == 0
            x_complex{ell+1}(:, m + ell + 1) = x_real{ell + 1}(:, m + ell + 1);
        end
    end
end
end