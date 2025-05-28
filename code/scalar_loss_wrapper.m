function loss = scalar_loss_wrapper(x_shell, est_coeffs, shell_idx, current_indices, ell, d, B_target)
    est_coeffs(shell_idx, current_indices) = x_shell;
    B_current = compute_SO3_bispectrum_restricted(est_coeffs, d, ell, shell_idx);
    B_diff = B_current - B_target;
    loss = sum(abs(B_diff(:)).^2);
end