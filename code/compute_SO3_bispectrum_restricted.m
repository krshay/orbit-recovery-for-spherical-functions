function B = compute_SO3_bispectrum_restricted(a_lm_matrix, d, ell, shell_idx)
    num_shells = size(a_lm_matrix, 1);
    B = zeros(num_shells, num_shells, num_shells, d+1, d+1, d+1);
    for l1 = 0:d
        for l2 = 0:d
            for l3 = abs(l1 - l2):min(l1 + l2, d)
                lvals = [l1, l2, l3];
                idx_ell = find(lvals == ell);
                if length(idx_ell) ~= 1
                    continue;
                end
                for s1 = 1:num_shells
                    for s2 = 1:num_shells
                        for s3 = 1:num_shells
                            shells = [s1, s2, s3];
                            if shells(idx_ell) ~= shell_idx
                                continue;
                            end
                            sum_m = 0;
                            for m1 = -l1:l1
                                for m2 = -l2:l2
                                    m3 = -m1 - m2;
                                    if abs(m3) > l3, continue; end
                                    idx1 = sh_index(l1, m1);
                                    idx2 = sh_index(l2, m2);
                                    idx3 = sh_index(l3, m3);
                                    cg = clebsch_gordan(l2, m2, l3, m3, l1, -m1);
                                    a1 = a_lm_matrix(s1, idx1);
                                    a2 = a_lm_matrix(s2, idx2);
                                    a3 = a_lm_matrix(s3, idx3);
                                    sum_m = sum_m + (-1)^m1 * a1 * a2 * a3 * cg;
                                end
                            end
                            B(s1,s2,s3,l1+1,l2+1,l3+1) = sum_m;
                        end
                    end
                end
            end
        end
    end
end
