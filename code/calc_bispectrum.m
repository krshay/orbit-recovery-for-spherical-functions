function B = calc_bispectrum(ell_max, num_shells, x)
B = zeros(ell_max+1, ell_max+1, ell_max+1, num_shells, num_shells, num_shells);
for ell2=0:ell_max
    for ell3=0:ell_max
        for ell1=abs(ell2-ell3):min(ell2+ell3, ell_max)
            for k1=-ell1:ell1
                for k2=-ell2:ell2
                    for k3=-ell3:ell3
                        if k1 - k2 - k3 == 0
                            for s1=1:num_shells
                                for s2=1:num_shells
                                    for s3=1:num_shells
                                        B(ell1 + 1, ell2 + 1, ell3 + 1, ...
                                            s1, s2, s3) = ClebschGordan( ...
                                            ell2, ell3, ell1, k2, k3, k1 ...
                                            ) * conj(x{ell1 + 1}( ...
                                            s1, k1 + ell1 + 1 ...
                                            )) * x{ell2 + 1}( ...
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