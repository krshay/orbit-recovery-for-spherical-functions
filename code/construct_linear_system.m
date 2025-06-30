function [A,b] = construct_linear_system(est_coeffs,shell_idx,curr_idx,ell,d,B_target,cg_tab,shift)
num_shells = size(est_coeffs,1); num_coeffs = numel(curr_idx);
A = []; b = [];
for l1 = 0:d
    for l2 = 0:d
        for l3 = abs(l1-l2):min(l1+l2,d)
            lvals = [l1,l2,l3]; idx_ell = find(lvals==ell,1); if isempty(idx_ell), continue; end
            if any(lvals(setdiff(1:3,idx_ell)) >= ell), continue; end

            for s1 = 1:num_shells
                for s2 = 1:num_shells
                    for s3 = 1:num_shells
                        shells = [s1,s2,s3]; if shells(idx_ell)~=shell_idx, continue; end
                        coeff_row = zeros(1,num_coeffs); sum_known = 0;

                        for m1=-l1:l1
                            for m2=-l2:l2
                                m3=-m1-m2; if abs(m3)>l3, continue; end
                                idx1=sh_index(l1,m1); idx2=sh_index(l2,m2); idx3=sh_index(l3,m3);
                                cg = cg_tab(l1+1,m1+shift+1,l2+1,m2+shift+1,l3+1);

                                idxs=[idx1,idx2,idx3]; sidx=shells;
                                if ismember(idxs(idx_ell),curr_idx)
                                    k = find(curr_idx==idxs(idx_ell));
                                    val1 = est_coeffs(sidx(mod(idx_ell,3)+1),idxs(mod(idx_ell,3)+1));
                                    val2 = est_coeffs(sidx(mod(idx_ell+1,3)+1),idxs(mod(idx_ell+1,3)+1));
                                    coeff_row(k) = coeff_row(k) + (-1)^m1*val1*val2*cg;
                                else
                                    a1=est_coeffs(sidx(1),idx1); a2=est_coeffs(sidx(2),idx2); a3=est_coeffs(sidx(3),idx3);
                                    sum_known = sum_known + (-1)^m1*a1*a2*a3*cg;
                                end
                            end
                        end
                        A(end+1,:) = coeff_row; 
                        b(end+1,1) = B_target(s1,s2,s3,l1+1,l2+1,l3+1) - sum_known;
                    end
                end
            end
        end
    end
end
end