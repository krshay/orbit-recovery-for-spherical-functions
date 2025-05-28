function rYlm = generate_rYlm(jball, Y, dels_Y, N, ell, L)


siz_lms = 2*ell+1;
siz_grid = numel(jball);
rYlm = zeros(siz_grid, siz_lms);

 
% Y_{l, m} computation.
Y_lm = Y(dels_Y == ell, :).';

% Fourier-Bessel basis functions computation.
indcol = 0;
for ii=1:(2 * ell + 1)
     rYlm(:, indcol+1) = Y_lm(:, ii);
    indcol = indcol + 1;
end
end
