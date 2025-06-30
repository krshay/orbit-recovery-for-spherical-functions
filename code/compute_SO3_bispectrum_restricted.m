function B = compute_SO3_bispectrum_restricted(a_lm,d,ell,shell_idx,cg_tab,shift)
num_shells=size(a_lm,1);
B=zeros(num_shells,num_shells,num_shells,d+1,d+1,d+1);
for l1=0:d
    for l2=0:d
        for l3=abs(l1-l2):min(l1+l2,d)
            lvals=[l1,l2,l3]; idx_ell=find(lvals==ell,1); if isempty(idx_ell),continue; end
            for s1=1:num_shells
                for s2=1:num_shells
                    for s3=1:num_shells
                        shells=[s1,s2,s3]; if shells(idx_ell)~=shell_idx,continue; end
                        acc=0;
                        for m1=-l1:l1
                            for m2=-l2:l2
                                m3=-m1-m2; if abs(m3)>l3, continue; end
                                cg=cg_tab(l1+1,m1+shift+1,l2+1,m2+shift+1,l3+1);
                                a1=a_lm(s1,sh_index(l1,m1));
                                a2=a_lm(s2,sh_index(l2,m2));
                                a3=a_lm(s3,sh_index(l3,m3));
                                acc = acc + (-1)^m1 * a1*a2*a3 * cg;
                            end
                        end
                        B(s1,s2,s3,l1+1,l2+1,l3+1)=acc;
                    end
                end
            end
        end
    end
end
end