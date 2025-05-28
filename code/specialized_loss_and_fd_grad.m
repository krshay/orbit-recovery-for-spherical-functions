function [loss, grad] = specialized_loss_and_fd_grad(x_shell, est_coeffs, shell_idx, ...
        current_indices, ell, d, true_coeffs)

    est_coeffs(shell_idx, current_indices) = x_shell;
    B_target = compute_SO3_bispectrum_restricted(true_coeffs, d, ell, shell_idx);

    loss_fun = @(x) scalar_loss_wrapper(x, est_coeffs, shell_idx, current_indices, ell, d, B_target);
    loss = loss_fun(x_shell);
    grad = finite_difference_gradient(loss_fun, x_shell);
end