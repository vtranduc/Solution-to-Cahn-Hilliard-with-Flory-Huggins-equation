function sf=compute_sf_3d(neight,conc_,conco_,nx,ny,nz,nex,ney,weights,dt,n1,n2,dxyz,diffT,chi,wTerms)

sf=zeros(1,neight);

%----
% diff=1800;
% T=0.6;
% chi=(0.5-1*(1-1/T))*ones(1,3);
% diffT=diff*T*ones(1,3);
% % size(weights)
% % size(conc)
%-------
if isscalar(diffT)
    parfor gbf=1:1:neight %PARFOR MUST BE FEASIBLE!!!
        sf(gbf)=compute_sfth_sf_3d(gbf,conc_,conco_,nx,ny,nz,nex,ney,weights,dt,chi,diffT,n1,n2,dxyz,wTerms);
    end
else
    parfor gbf=1:1:neight %PARFOR MUST BE FEASIBLE!!!
        sf(gbf)=compute_sfth_sf_gradient_3d(gbf,conc_,conco_,nx,ny,nz,nex,ney,weights,dt,chi,diffT,n1,n2,dxyz,wTerms);
    end
end

end

function solution=compute_sfth_sf_gradient_3d(sfth,conc_,conco_,nx,ny,nz,nex,ney,weights,dt,chi,diffT,n1,n2,dxyz,wTerms)

solution=0;

w=[5/18 4/9 5/18];

[node,type]=analyze_gbs_3d(sfth);
[xth,yth,zth]=get_xyzth_3d(node,nx,ny);
elements=get_elements_3d(xth,yth,zth,nex,ney,nx,ny,nz);

for i_positional_element=1:1:8
    if elements(i_positional_element)~=0
        
        if mod(i_positional_element,2)==1
            exth=xth-1;
        else
            exth=xth;
        end
        
        con=dim_reducer2_3d(conc_(elements(i_positional_element),1,:,:,:));
        cono=dim_reducer1_3d(conco_(elements(i_positional_element),:,:,:));
        orientation=9-i_positional_element;
        cont=(con-cono)/dt;
        
        for ix=1:1:3
            for iy=1:1:3
                for iz=1:1:3
                    solution=solution+w(ix)*w(iy)*w(iz)*(...
                        weights(orientation,type,1,ix,iy,iz)*cont(ix,iy,iz)-diffT(exth,ix)*weights(orientation,type,1,ix,iy,iz)*...
                        ((-1/(con(ix,iy,iz)^2*n1)+1/((1-con(ix,iy,iz))^2*n2))*...
                        (conc_(elements(i_positional_element),2,ix,iy,iz)^2+conc_(elements(i_positional_element),3,ix,iy,iz)^2+conc_(elements(i_positional_element),4,ix,iy,iz)^2)+...
                        (1/(con(ix,iy,iz)*n1)+1/((1-con(ix,iy,iz))*n2)-2*chi(exth,ix))*...
                        (conc_(elements(i_positional_element),5,ix,iy,iz)+conc_(elements(i_positional_element),6,ix,iy,iz)+conc_(elements(i_positional_element),7,ix,iy,iz)))+...
                        (conc_(elements(i_positional_element),5,ix,iy,iz)+conc_(elements(i_positional_element),6,ix,iy,iz)+conc_(elements(i_positional_element),7,ix,iy,iz))*...
                        wTerms(orientation,type,ix,iy,iz)...
                        );
                end
            end
        end
    end
end
solution=dxyz*solution;
end

function solution=compute_sfth_sf_3d(sfth,conc_,conco_,nx,ny,nz,nex,ney,weights,dt,chi,diffT,n1,n2,dxyz,wTerms)

solution=0;

w=[5/18 4/9 5/18];

[node,type]=analyze_gbs_3d(sfth);
[xth,yth,zth]=get_xyzth_3d(node,nx,ny);
elements=get_elements_3d(xth,yth,zth,nex,ney,nx,ny,nz);

for i_positional_element=1:1:8
    if elements(i_positional_element)~=0
        
        con=dim_reducer2_3d(conc_(elements(i_positional_element),1,:,:,:));
        cono=dim_reducer1_3d(conco_(elements(i_positional_element),:,:,:));
        orientation=9-i_positional_element;
        cont=(con-cono)/dt;
        
        for ix=1:1:3
            for iy=1:1:3
                for iz=1:1:3
                    solution=solution+w(ix)*w(iy)*w(iz)*(...
                        weights(orientation,type,1,ix,iy,iz)*cont(ix,iy,iz)-diffT*weights(orientation,type,1,ix,iy,iz)*...
                        ((-1/(con(ix,iy,iz)^2*n1)+1/((1-con(ix,iy,iz))^2*n2))*...
                        (conc_(elements(i_positional_element),2,ix,iy,iz)^2+conc_(elements(i_positional_element),3,ix,iy,iz)^2+conc_(elements(i_positional_element),4,ix,iy,iz)^2)+...
                        (1/(con(ix,iy,iz)*n1)+1/((1-con(ix,iy,iz))*n2)-2*chi)*...
                        (conc_(elements(i_positional_element),5,ix,iy,iz)+conc_(elements(i_positional_element),6,ix,iy,iz)+conc_(elements(i_positional_element),7,ix,iy,iz)))+...
                        (conc_(elements(i_positional_element),5,ix,iy,iz)+conc_(elements(i_positional_element),6,ix,iy,iz)+conc_(elements(i_positional_element),7,ix,iy,iz))*...
                        wTerms(orientation,type,ix,iy,iz)...
                        );
                end
            end
        end
    end
end
solution=dxyz*solution;
end