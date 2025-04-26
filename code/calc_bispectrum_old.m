function B = calc_bispectrum_old(ell_max, num_shells, x)
B = zeros(ell_max+1, ell_max+1, ell_max+1, num_shells, num_shells, num_shells);
for ell1=0:ell_max
    for ell2=0:ell_max
        for ell3=0:ell_max
            for s1=1:num_shells
                for s2=1:num_shells
                    for s3=1:num_shells
                        for k1=-ell1:ell1
                            for k2=-ell2:ell2
                                for k3=-ell3:ell3
                                    if k1 + k2 + k3 == 0
                                        B(ell1 + 1, ell2 + 1, ell3 + 1, ...
                                            s1, s2, s3) = (-1)^k1 * ClebschGordan( ...
                                            ell2, ell3, ell1, k2, k3, -k1 ...
                                            ) * x{ell1 + 1}( ...
                                            s1, k1 + ell1 + 1 ...
                                            ) * x{ell2 + 1}( ...
                                            s2, k2 + ell2 + 1 ...
                                            ) * x{ell3 + 1}( ...
                                            s3, k3 + ell3 + 1 ...
                                            );
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
