function x_real = complex2real_spherical_harmonics(x_complex, ell_max)
% Volume expansion in real spherical harmonics
x_real = cell(size(x_complex));
num_shells = size(x_complex{1}, 1);
for ell=0:ell_max
    x_real{ell+1} = zeros(num_shells, 2 * ell + 1);
    for m=-ell:ell
        if m < 0
            x_real{ell+1}(:, m + ell + 1) = real((1j / sqrt(2)) * (x_complex{ell + 1}(:, -abs(m) + ell + 1) - ((-1) ^ m) * x_complex{ell + 1}(:, abs(m) + ell + 1)));
        elseif m > 0
            x_real{ell+1}(:, m + ell + 1) = real((1 / sqrt(2)) * (x_complex{ell + 1}(:, -abs(m) + ell + 1) + ((-1) ^ m) * x_complex{ell + 1}(:, abs(m) + ell + 1)));
        elseif m == 0
            x_real{ell+1}(:, m + ell + 1) = real(x_complex{ell + 1}(:, m + ell + 1));
        end
    end
end
end