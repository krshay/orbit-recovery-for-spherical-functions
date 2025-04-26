function [B, grad_B] = calc_bispectrum_and_grad_real(ell_max, num_shells, x)
N = length(x);
B = zeros(ell_max + 1, ell_max + 1, ell_max + 1, num_shells, num_shells, num_shells);
grad_B = zeros(N, ell_max + 1, ell_max + 1, ell_max + 1, num_shells, num_shells, num_shells);
for ell2=0:ell_max
    for ell3=0:ell_max
        for ell1=abs(ell2-ell3):min(ell2+ell3, ell_max)
            for k1=-ell1:ell1
                for k2=-ell2:ell2
                    for k3=-ell3:ell3
                        for s1=1:num_shells
                            for s2=1:num_shells
                                for s3=1:num_shells
                                    if k1 - k2 - k3 == 0
                                        cg_coeff = ClebschGordan( ...
                                            ell2, ell3, ell1, k2, k3, k1);
                                        idx1 = sum((2 * (0:ell1-1) + 1) * num_shells) + (k1 + ell1) * num_shells + s1;
                                        idx2 = sum((2 * (0:ell2-1) + 1) * num_shells) + (k2 + ell2) * num_shells + s2;
                                        idx3 = sum((2 * (0:ell3-1) + 1) * num_shells) + (k3 + ell3) * num_shells + s3;
                                        x1 = x(idx1);
                                        x2 = x(idx2);
                                        x3 = x(idx3);
                                        y1 = -x(idx1 + N / 2);
                                        y2 = x(idx2 + N / 2);
                                        y3 = x(idx3 + N / 2);

                                        B_re = x1 * x2 * x3 - x1 * y2 * y3 - y1 * y2 * x3 - y1 * x2 * y3;
                                        B_im = x1 * x2 * y3 + x1 * y2 * x3 + y1 * x2 * x3 - y1 * y2 * y3;
                                        B(ell1 + 1, ell2 + 1, ell3 + 1, ...
                                            s1, s2, s3) = cg_coeff * (B_re + 1i * B_im);
                                        grad_B(idx1, ell1 + 1, ell2 + 1, ell3 + 1, s1, s2, s3) = grad_B(idx1, ell1 + 1, ell2 + 1, ell3 + 1, ...
                                            s1, s2, s3) + cg_coeff * ((x2 * x3 - y2 * y3) + 1i * (x2 * y3 + y2 * x3));
                                        grad_B(idx2, ell1 + 1, ell2 + 1, ell3 + 1, s1, s2, s3) = grad_B(idx2, ell1 + 1, ell2 + 1, ell3 + 1, ...
                                            s1, s2, s3) + cg_coeff * ((x1 * x3 - y1 * y3) + 1i * (x1 * y3 + y1 * x3));
                                        grad_B(idx3, ell1 + 1, ell2 + 1, ell3 + 1, s1, s2, s3) = grad_B(idx3, ell1 + 1, ell2 + 1, ell3 + 1, ...
                                            s1, s2, s3) + cg_coeff * ((x1 * x2 - y1 * y2) + 1i * (y1 * x2 + x1 * y2));
                                        grad_B(idx1 + N / 2, ell1 + 1, ell2 + 1, ell3 + 1, s1, s2, s3) = grad_B(idx1 + N / 2, ell1 + 1, ell2 + 1, ell3 + 1, ...
                                            s1, s2, s3) + cg_coeff * ((-y2 * x3 - x2 * y3) + 1i * (x2 * x3 - y2 * y3));
                                        grad_B(idx2 + N / 2, ell1 + 1, ell2 + 1, ell3 + 1, s1, s2, s3) = grad_B(idx2 + N / 2, ell1 + 1, ell2 + 1, ell3 + 1, ...
                                            s1, s2, s3) + cg_coeff * ((-x1 * y3 - y1 * x3) + 1i * (x1 * x3 - y1 * y3));
                                        grad_B(idx3 + N / 2, ell1 + 1, ell2 + 1, ell3 + 1, s1, s2, s3) = grad_B(idx3 + N / 2, ell1 + 1, ell2 + 1, ell3 + 1, ...
                                            s1, s2, s3) + cg_coeff * ((-x1 * y2 - y1 * x2) + 1i * (x1 * x2 - y1 * y2));
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
