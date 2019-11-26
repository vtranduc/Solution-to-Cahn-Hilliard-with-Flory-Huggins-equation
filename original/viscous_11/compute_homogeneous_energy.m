function energy_homogeneous=compute_homogeneous_energy(...
    ci_ave,diff,coef_T,n1,n2,Two_chi_n1,xlen,ylen,ne,dxdy,grad_T)

% return
if grad_T==0
    c2=1-ci_ave;
    energy_homogeneous=ci_ave*log(ci_ave)/n1+c2*log(c2)/n2...
        +ci_ave*c2*Two_chi_n1/2;
    energy_homogeneous=2*coef_T*energy_homogeneous;
    energy_homogeneous=energy_homogeneous*xlen*ylen;
    
elseif grad_T==1
    
    chi_n1=Two_chi_n1/2;
    c2=1-ci_ave;
    
    z1=2*diff;
    z2=ci_ave*log(ci_ave)/n1...
        +c2*log(c2)/n2;
    z3=ci_ave*c2;
    
    w=[5/18 4/9 5/18];
    w_set=zeros(3,3);
    for m=1:1:3
        for n=1:1:3
            w_set(m,n)=w(m)*w(n);
        end
    end
    
    energy_homogeneous=0;
    
    for e=1:1:ne
        for iny=1:1:3
            for inx=1:1:3
                
                energy_homogeneous=energy_homogeneous+w_set(inx,iny)*(...
                    z1*coef_T(e,inx,iny,1)...
                    *(z2+z3*chi_n1(e,inx,iny)));
                
            end
        end
    end
    energy_homogeneous=energy_homogeneous*dxdy;
    
end


end