function [z, cst] = optimize(x_init, B, ell_max, num_shells)
%optimize Optimizes the cost function (based on the provided map - for the 
%case of bispectrum or trispectrum) using fminunc

%   z_init is the initial guess for the Fourier-Bessel coefficients
%   B is the true bispectrum
%   A is the sensing matrix

fun = @(x) cost_and_grad_bispectrum(x, B, ell_max, num_shells);
options = optimoptions(@fminunc,'Display', ...
    'iter-detailed','SpecifyObjectiveGradient',true, ...
        'MaxFunctionEvaluations', 50, ...
    'CheckGradients',false);
    % 'OptimalityTolerance', 1e-15, ...
    % 'StepTolerance', 1e-15, ...

% options = optimset('PlotFcns',@optimplotfval,'Display','iter-detailed','MaxFunEvals', 10);

[z, cst, exitflag,output,grad,hessian] =  fminunc(fun, x_init, options);
% [z, cst] =  fminsearch(fun, x_init, options);

end