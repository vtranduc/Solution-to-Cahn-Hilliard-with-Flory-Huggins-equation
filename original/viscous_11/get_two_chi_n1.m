function Two_chi_n1=get_two_chi_n1(ne,coef_T,entropy,n1)
Two_chi_n1=zeros(ne,3,3);
for e=1:1:ne
    for ix=1:1:3
        for iy=1:1:3
            Two_chi_n1(e,ix,iy)=...
                get_Two_chi_n1_with_order(coef_T(e,ix,iy,1),entropy,n1,0);
        end
    end
end
end