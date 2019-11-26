function total_energy=compute_total_energy_2d(...
    coef_T,n1,n2,ne,ney,Two_chi_n1,dxdy,...
    ny,c,weights,...
    grad_T,diff)
w=[5/18 4/9 5/18];
w_set=zeros(3,3);
for m=1:1:3
    for n=1:1:3
        w_set(m,n)=w(m)*w(n);
    end
end
[conc_abs,abs_gradient_squared]=get_conc_abs(ne,ney,ny,c,weights);
total_energy=0;
if grad_T==0
    two_diffT=2*coef_T;
    chi_n1=Two_chi_n1/2;
    for e=1:1:ne
        for iny=1:1:3
            for inx=1:1:3
                c2=1-conc_abs(inx,iny,e);
                total_energy=total_energy+w_set(inx,iny)*(...
                    two_diffT...
                    *(conc_abs(inx,iny,e)*log(conc_abs(inx,iny,e))/n1...
                    +c2*log(c2)/n2...
                    +conc_abs(inx,iny,e)*c2*chi_n1)...
                    +abs_gradient_squared(inx,iny,e));
            end
        end
    end
elseif grad_T==1 %============================================
    two_diffT=2*coef_T*diff;
    chi_n1=Two_chi_n1/2;
    inxoo=0;
    for e=1:1:ne
        if mod(e,ney)==1
            inxoo=inxoo+1;
        end
        for iny=1:1:3
            for inx=1:1:3
                c2=1-conc_abs(inx,iny,e);
                total_energy=total_energy+w_set(inx,iny)*(...
                    two_diffT(inxoo,inx)...
                    *(conc_abs(inx,iny,e)*log(conc_abs(inx,iny,e))/n1...
                    +c2*log(c2)/n2...
                    +conc_abs(inx,iny,e)*c2*chi_n1(inxoo,inx))...
                    +abs_gradient_squared(inx,iny,e));
            end
        end
    end
end
total_energy=dxdy*total_energy;
end

function [conc_abs,abs_gradient_squared]=get_conc_abs(ne,ney,ny,c,weights)
conc_abs=zeros(3,3,ne);
abs_gradient_squared=zeros(3,3,ne);
for element=1:1:ne
    [inx,iny]=inxiny_elemental(element,ney);
    gbfs=elemental_gbf(inx,iny,ny);
    cElemental=get_cElemental(gbfs,c);
    conc_abs(:,:,element)=conc(cElemental,weights,1);
    abs_gradient_squared(:,:,element)=conc(cElemental,weights,2).^2+conc(cElemental,weights,3).^2;
end
end