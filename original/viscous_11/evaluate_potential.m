function potential=evaluate_potential(c,nex,ney,ny,n1,n2,Two_chi_n1,...
    weights,grad_T,coef_T)

potential=zeros(ney,nex);

% test1=potential;
% test2=potential;

%--------------------------

% solution=conc(c_,weights,order_type)

% ie=0;

% c_=ones(1,16);
if grad_T==0
    e=0;
    for exth=1:1:nex
        for eyth=1:1:ney
            e=e+1;

            [inx,iny]=inxiny_elemental(e,ney);
            gbfs=elemental_gbf(inx,iny,ny);
            cElemental=get_cElemental(gbfs,c);

            u=0;

            laplacian=0;
            for ilocal=1:1:16
                laplacian=laplacian+cElemental(ilocal)*weights(2,2,ilocal,4);
                laplacian=laplacian+cElemental(ilocal)*weights(2,2,ilocal,5);

                u=u+cElemental(ilocal)*weights(2,2,ilocal,1);

            end
            
            
            potential(eyth,exth)=coef_T*(...
                log(u)/n1-log(1-u)/n2+1/n1-1/n2+Two_chi_n1*(1-2*u)/2)...
                -laplacian;
            
                
        end
    end
    
elseif grad_T==1
    
    e=0;
    for exth=1:1:nex
        for eyth=1:1:ney
            e=e+1;
            
            [inx,iny]=inxiny_elemental(e,ney);
            gbfs=elemental_gbf(inx,iny,ny);
            cElemental=get_cElemental(gbfs,c);

            u=0;

            laplacian=0;
            
            
            
            
            for ilocal=1:1:16
                laplacian=laplacian+cElemental(ilocal)*weights(2,2,ilocal,4);
                laplacian=laplacian+cElemental(ilocal)*weights(2,2,ilocal,5);

                u=u+cElemental(ilocal)*weights(2,2,ilocal,1);

            end
            
            potential(eyth,exth)=coef_T(e,2,2,1)*(...
                log(u)/n1-log(1-u)/n2+1/n1-1/n2+Two_chi_n1(e,2,2)*(1-2*u)/2)...
                -laplacian;
                
        end
    end
end

return

% error('rinchan')
% 
% %------------------------
% 
% dx=zeros(ny,nx);
% dy=zeros(ny,nx);
% 
% u=dx;
% 
% index=-2;
% for xth=1:1:nx
%     for yth=1:1:ny
%         index=index+4;
%         dx(yth,xth)=c(index);
%         dy(yth,xth)=c(index+1);
%         
%         u(yth,xth)=c(index-1);
%         
%     end
% end
% 
% % [~,dxx]=gradient(dx);
% % [dyy,~]=gradient(dy);
% 
% 
% [dxx,dyy]=gradient(u);
% 
% % 
% % c'
% % 
% % error('fdaafd')

end

% function solution=get_conc_1_3(ne,ney,ny,c,weights)
% 
% solution=zeros(3,3,5,ne);
% 
% parfor element=1:1:ne
%     [inx,iny]=inxiny_elemental(element,ney);
%     gbfs=elemental_gbf(inx,iny,ny);
%     cElemental=get_cElemental(gbfs,c);
%     
%     
%     conc_=zeros(3,3,5);
%     for order_type=1:1:5
%         conc_(:,:,order_type)=conc(cElemental,weights,order_type);
%     end
%     solution(:,:,:,element)=conc_;
% end
% 
% end