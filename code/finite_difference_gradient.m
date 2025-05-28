function grad_fd = finite_difference_gradient(fun, x)
    epsilon = 1e-6;
    grad_fd = zeros(size(x));
    for i = 1:length(x)
        x_fwd = x; x_fwd(i) = x_fwd(i) + epsilon;
        x_bwd = x; x_bwd(i) = x_bwd(i) - epsilon;
        grad_fd(i) = (fun(x_fwd) - fun(x_bwd)) / (2 * epsilon);
    end
end