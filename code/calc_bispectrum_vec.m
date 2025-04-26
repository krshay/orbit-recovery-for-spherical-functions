function B = calc_bispectrum_vec(ell_max, num_shells, x)
B = zeros(ell_max + 1, ell_max + 1, ell_max + 1, num_shells, num_shells, num_shells);
idx_ell1 = 0;
for ell1=0:ell_max
    idx_ell2 = 0;
    for ell2=0:ell_max
        idx_ell3 = 0;
        for ell3=0:ell_max
            idx_k1 = 0;
            for k1=-ell1:ell1
                idx_k2 = 0;
                for k2=-ell2:ell2
                    idx_k3 = 0;
                    for k3=-ell3:ell3

                        for s1=1:num_shells
                            for s2=1:num_shells
                                for s3=1:num_shells
                                    if k1 + k2 + k3 == 0
                                        B(ell1 + 1, ell2 + 1, ell3 + 1, ...
                                            s1, s2, s3) = (-1)^k1 * ClebschGordan( ...
                                            ell2, ell3, ell1, k2, k3, -k1) * x(sum((2 * (0:ell1-1) + 1) * num_shells) + (k1 + ell1) * num_shells + s1)* x(sum((2 * (0:ell2-1) + 1) * num_shells) + (k2 + ell2) * num_shells + s2) * x(sum((2 * (0:ell3-1) + 1) * num_shells) + (k3 + ell3) * num_shells + s3);
                                    end
                                end
                                
                            end
                                                

                        end
                                        

idx_k3 = idx_k3 + num_shells;
                    end
idx_k2 = idx_k2 + num_shells;
                end

                idx_k1 = idx_k1 + num_shells;
            end
            idx_ell3 = idx_ell3 + (2 * ell3 + 1) * num_shells;
        end
        idx_ell2 = idx_ell2 + (2 * ell2 + 1) * num_shells;
    end
    idx_ell1 = idx_ell1 + (2 * ell1 + 1) * num_shells;
end
end
