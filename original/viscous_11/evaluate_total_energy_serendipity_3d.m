function total_energy=evaluate_total_energy_serendipity_3d(...
    c,diff,weights,nex,ney,nx,ny,ne,n1,n2,T,entropy,dxyz)

total_energy=0;

w=[5/18 4/9 5/18];
% w_set=zeros(3,3);
% for m=1:1:3
%     for n=1:1:3
%         w_set(m,n)=w(m)*w(n);
%     end
% end

conc_=compute_elemental_elemental_conc_serendipity_2_orders(...
    c,weights,nex,nex*ney,nx*ny,nx,ne);

if length(T)==1
    grad_T=0;
else
    error('This function has not been designed for varying temperature condition')
end

if grad_T==0
    
    coef_T=2*diff*T;
    chi_n1=dimensionless_chi(T,entropy)/n1;
    
    for e=1:1:ne
        for inz=1:1:3
            for iny=1:1:3
                for inx=1:1:3
                    
                    c2=1-conc_(e,1,inx,iny,inz);
                    
                    total_energy=total_energy+w(inx)*w(iny)*w(inz)*(...
                        ...
                        ...
                        coef_T...
                        *(conc_(e,1,inx,iny,inz)*log(conc_(e,1,inx,iny,inz))/n1...
                        +c2*log(c2)/n2...
                        +conc_(e,1,inx,iny,inz)*c2*chi_n1)...
                        ...
                        ...
                        +conc_(e,2,inx,iny,inz)^2+...
                        conc_(e,3,inx,iny,inz)^2+...
                        conc_(e,4,inx,iny,inz)^2);
                    
                end
            end
        end
    end
    
    total_energy=total_energy*dxyz;
    
elseif grad_T==1
    
    
end


end

% function conc_=compute_conc_serendipity_2_orders_3d(ne,c,weights,nex,nexney,nxny,nx)
% 
% conc_=zeros(ne,7,3,3,3);
% 
% for e=1:1:ne %APPARENTLY, MAKING THIS PARFOR WILL MAKE IT SLOWER
%     conc_(e,:,:,:,:)=compute_elemental_elemental_conc_serendipity(e,c,weights,nex,nexney,nxny,nx);
% end
% 
% end

function elemental_conc_=compute_elemental_elemental_conc_serendipity_2_orders(...
    c,weights,nex,nexney,nxny,nx,ne)

elemental_conc_=zeros(ne,4,3,3,3);

for e=1:1:ne
    gbfs=get_gbfs_of_element_serendipity(e,nex,nexney,nxny,nx);

    for order=1:1:4
        for ix=1:1:3
            for iy=1:1:3
                for iz=1:1:3
                    index=0;
                    for orientation=1:1:8
                        for type=1:1:4
                            index=index+1;
                            elemental_conc_(e,order,ix,iy,iz)=...
                                elemental_conc_(e,order,ix,iy,iz)...
                                +c(gbfs(index))...
                                *weights(orientation,type,order,ix,iy,iz);
                        end
                    end
                end
            end
        end
    end
end

end