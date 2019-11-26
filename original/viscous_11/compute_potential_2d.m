function P=compute_potential_2d(c_nodal,nx,ny,n1,n2,coef_T,diff,grad_T,Two_chi_n1,entropy,T)
P=zeros(ny,nx);
if grad_T==0
    two_diffT=coef_T*2;
    chi_n1=Two_chi_n1/2;
    for inx=1:1:nx
        for iny=1:1:ny
            P(iny,inx)=two_diffT*((log(c_nodal(iny,inx))+1)/n1...
                -(log((1-c_nodal(iny,inx)))+1)/n2...
                +chi_n1*(1-2*c_nodal(iny,inx)));
        end
    end
elseif grad_T==1
    two_diffT=linspace(T(1),T(2),nx);
    chi_n1=get_Two_chi_n1_with_order(two_diffT,entropy,n1,0)/2;
    two_diffT=two_diffT*diff*2;
    for inx=1:1:nx
        for iny=1:1:ny
            P(iny,inx)=two_diffT(inx)*((log(c_nodal(iny,inx))+1)/n1...
                -(log((1-c_nodal(iny,inx)))+1)/n2...
                +chi_n1(inx)*(1-2*c_nodal(iny,inx)));
        end
    end
end
end