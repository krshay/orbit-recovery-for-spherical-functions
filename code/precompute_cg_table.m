function cg_tab = precompute_cg_table(d)
shift = d;                              % max |m|
cg_tab = zeros(d+1,2*shift+1,d+1,2*shift+1,d+1);  % l1,m1,l2,m2,l3 (m3 implicit)
for l1 = 0:d
    for m1 = -l1:l1
        for l2 = 0:d
            for m2 = -l2:l2
                for l3 = abs(l1-l2):min(l1+l2,d)
                    m3 = -m1-m2;
                    if abs(m3)>l3, continue; end
                    cg = clebsch_gordan(l2,m2,l3,m3,l1,-m1);
                    cg_tab(l1+1,m1+shift+1,l2+1,m2+shift+1,l3+1) = cg;
                end
            end
        end
    end
end
end