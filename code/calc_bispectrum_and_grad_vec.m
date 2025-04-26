function [B, grad_B] = calc_bispectrum_and_grad_vec(ell_max, num_shells, x)
B = zeros(ell_max + 1, ell_max + 1, ell_max + 1, num_shells, num_shells, num_shells);
grad_B = zeros(length(x), ell_max + 1, ell_max + 1, ell_max + 1, num_shells, num_shells, num_shells);
for ell1=0:ell_max
    for ell2=0:ell_max
        for ell3=0:ell_max
            for k1=-ell1:ell1
                for k2=-ell2:ell2
                    for k3=-ell3:ell3
                        for s1=1:num_shells
                            for s2=1:num_shells
                                for s3=1:num_shells
                                    if k1 + k2 + k3 == 0
                                        cg_coeff = (-1)^k1 * ClebschGordan( ...
                                            ell2, ell3, ell1, k2, k3, -k1);
                                        B(ell1 + 1, ell2 + 1, ell3 + 1, ...
                                            s1, s2, s3) = cg_coeff * x(sum((2 * (0:ell1-1) + 1) * num_shells) + (k1 + ell1) * num_shells + s1) * x(sum((2 * (0:ell2-1) + 1) * num_shells) + (k2 + ell2) * num_shells + s2) * x(sum((2 * (0:ell3-1) + 1) * num_shells) + (k3 + ell3) * num_shells + s3);
                                        grad_B(sum((2 * (0:ell1-1) + 1) * num_shells) + (k1 + ell1) * num_shells + s1, ell1 + 1, ell2 + 1, ell3 + 1, ...
                                            s1, s2, s3) = grad_B(sum((2 * (0:ell1-1) + 1) * num_shells) + (k1 + ell1) * num_shells + s1, ell1 + 1, ell2 + 1, ell3 + 1, ...
                                            s1, s2, s3) + cg_coeff * x(sum((2 * (0:ell2-1) + 1) * num_shells) + (k2 + ell2) * num_shells + s2) * x(sum((2 * (0:ell3-1) + 1) * num_shells) + (k3 + ell3) * num_shells + s3);
                                        grad_B(sum((2 * (0:ell2-1) + 1) * num_shells) + (k2 + ell2) * num_shells + s2, ell1 + 1, ell2 + 1, ell3 + 1, ...
                                            s1, s2, s3) = grad_B(sum((2 * (0:ell2-1) + 1) * num_shells) + (k2 + ell2) * num_shells + s2, ell1 + 1, ell2 + 1, ell3 + 1, ...
                                            s1, s2, s3) + cg_coeff * x(sum((2 * (0:ell1-1) + 1) * num_shells) + (k1 + ell1) * num_shells + s1) * x(sum((2 * (0:ell3-1) + 1) * num_shells) + (k3 + ell3) * num_shells + s3);
                                        grad_B(sum((2 * (0:ell3-1) + 1) * num_shells) + (k3 + ell3) * num_shells + s3, ell1 + 1, ell2 + 1, ell3 + 1, ...
                                            s1, s2, s3) = grad_B(sum((2 * (0:ell3-1) + 1) * num_shells) + (k3 + ell3) * num_shells + s3, ell1 + 1, ell2 + 1, ell3 + 1, ...
                                            s1, s2, s3) + cg_coeff * x(sum((2 * (0:ell1-1) + 1) * num_shells) + (k1 + ell1) * num_shells + s1) * x(sum((2 * (0:ell2-1) + 1) * num_shells) + (k2 + ell2) * num_shells + s2);
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
end
