function sf=test62(...
    neightTotal,conc_,conco_,nx,ny,nz,nex,ney,weights,dt,determinants,...
    diffT,wTerms,terms,n,ne,n_eSurf,n_nSurf,nTotal,nx_,...
    determinants_surface,terms_surface,normals,weights_surface)

sf=zeros(1,neightTotal);

w=[5/18 4/9 5/18];

if isscalar(diffT)
    parfor gbf=1:1:neightTotal %PARFOR MUST BE FEASIBLE!!!
        sf(gbf)=compute_sfth_sf_with_terms_sph(...
            gbf,conc_,conco_,nx,ny,nz,nex,ney,weights,dt,diffT,determinants,wTerms,terms,...
            n,ne,n_eSurf,n_nSurf,nTotal,nx_,w,...
            determinants_surface,terms_surface,normals,weights_surface);
    end
% else
%     parfor gbf=1:1:neightTotal %PARFOR MUST BE FEASIBLE!!!
%         sf(gbf)=compute_sfth_sf_gradient_with_terms_sph(gbf,conc_,conco_,nx,ny,nz,nex,ney,weights,dt,diffT,determinants,wTerms,terms);
%     end
end

end

% function solution=compute_sfth_sf_gradient_with_terms_sph(sfth,conc_,conco_,nx,ny,nz,nex,ney,weights,dt,diffT,determinants,wTerms,terms)
% 
% solution=0;
% 
% w=[5/18 4/9 5/18];
% 
% [node,type]=analyze_gbs_3d(sfth);
% [xth,yth,zth]=get_xyzth_3d(node,nx,ny);
% elements=get_elements_3d(xth,yth,zth,nex,ney,nx,ny,nz);
% 
% for i_positional_element=1:1:8
%     if elements(i_positional_element)~=0
%         
%         con=dim_reducer2_3d(conc_(elements(i_positional_element),1,:,:,:));
%         cono=dim_reducer1_3d(conco_(elements(i_positional_element),:,:,:));
%         orientation=9-i_positional_element;
%         cont=(con-cono)/dt;
%         
%         for ix=1:1:3
%             for iy=1:1:3
%                 for iz=1:1:3
%                     solution=solution+w(ix)*w(iy)*w(iz)...
%                         *determinants(elements(i_positional_element),ix,iy,iz)*(...
%                         weights(elements(i_positional_element),orientation,type,1,ix,iy,iz)*cont(ix,iy,iz)-diffT(elements(i_positional_element),ix,iy,iz)...
%                         *weights(elements(i_positional_element),orientation,type,1,ix,iy,iz)*...
%                         (terms(elements(i_positional_element),1,ix,iy,iz)+...
%                         terms(elements(i_positional_element),2,ix,iy,iz))+...
%                         terms(elements(i_positional_element),3,ix,iy,iz)*...
%                         wTerms(elements(i_positional_element),orientation,type,ix,iy,iz)...
%                         );
%                 end
%             end
%         end
%     end
% end
% % solution=dxyz*solution;
% end

function solution=compute_sfth_sf_with_terms_sph(...
    sfth,conc_,conco_,nx,ny,nz,nex,ney,weights,dt,diffT,determinants,wTerms,terms,...
    n,ne,n_eSurf,n_nSurf,nTotal,nx_,w,...
    determinants_surface,terms_surface,normals,weights_surface)

solution=0;

[node,type]=analyze_gbs_3d(sfth);

[elements,eExtruding_x,eExtruding_y,eExtruding_z,...
    xOrientations,yOrientations,zOrientations]=...
    get_adjacent_elements_sphere(...
    node,n,nx,ny,nz,ne,nex,ney,n_eSurf,n_nSurf,nTotal,nx_);


%--------------------------------------------------------
if node>n
    if ~isnan(eExtruding_x)
%         if x_dir==-1
%             orientations=[8 6 4 2];
%         elseif x_dir==1
%             orientations=[7 5 3 1];
%         end
        for i_positional_element=1:1:4
            if eExtruding_x(i_positional_element)~=0
                e_surface=eExtruding_x(i_positional_element)-ne;
                for dim1=1:1:3
                    for dim2=1:1:3
                        solution=solution-w(dim1)*w(dim2)...
                            *determinants_surface(e_surface,dim1,dim2)...
                            *terms_surface(e_surface,dim1,dim2)...
                            *(normals(e_surface,dim1,dim2,1)*weights_surface(e_surface,xOrientations(i_positional_element),type,2,dim1,dim2)...
                            +normals(e_surface,dim1,dim2,2)*weights_surface(e_surface,xOrientations(i_positional_element),type,3,dim1,dim2)...
                            +normals(e_surface,dim1,dim2,3)*weights_surface(e_surface,xOrientations(i_positional_element),type,4,dim1,dim2));
                    end
                end
            end
        end
    end
    if ~isnan(eExtruding_y)
%         if y_dir==-1
%             orientations=[8 7 4 3];
%         elseif y_dir==1
%             orientations=[6 5 2 1];
%         end
        for i_positional_element=1:1:4
            if eExtruding_y(i_positional_element)~=0
                e_surface=eExtruding_y(i_positional_element)-ne;
                for dim1=1:1:3
                    for dim2=1:1:3
                        solution=solution-w(dim1)*w(dim2)...
                            *determinants_surface(e_surface,dim1,dim2)...
                            *terms_surface(e_surface,dim1,dim2)...
                            *(normals(e_surface,dim1,dim2,1)*weights_surface(e_surface,yOrientations(i_positional_element),type,2,dim1,dim2)...
                            +normals(e_surface,dim1,dim2,2)*weights_surface(e_surface,yOrientations(i_positional_element),type,3,dim1,dim2)...
                            +normals(e_surface,dim1,dim2,3)*weights_surface(e_surface,yOrientations(i_positional_element),type,4,dim1,dim2));
                    end
                end
            end
        end
    end
    if ~isnan(eExtruding_z)
%         if z_dir==-1
%             orientations=[8 7 6 5];
%         elseif z_dir==1
%             orientations=[4 3 2 1];
%         end
        for i_positional_element=1:1:4
            if eExtruding_z(i_positional_element)~=0
                e_surface=eExtruding_z(i_positional_element)-ne;
                for dim1=1:1:3
                    for dim2=1:1:3
                        solution=solution-w(dim1)*w(dim2)...
                            *determinants_surface(e_surface,dim1,dim2)...
                            *terms_surface(e_surface,dim1,dim2)...
                            *(normals(e_surface,dim1,dim2,1)*weights_surface(e_surface,zOrientations(i_positional_element),type,2,dim1,dim2)...
                            +normals(e_surface,dim1,dim2,2)*weights_surface(e_surface,zOrientations(i_positional_element),type,3,dim1,dim2)...
                            +normals(e_surface,dim1,dim2,3)*weights_surface(e_surface,zOrientations(i_positional_element),type,4,dim1,dim2));
                    end
                end
            end
        end
    end
end

% --------------------------------------------------------------


if ~isnan(elements)
    for i_positional_element=1:1:8
        if elements(i_positional_element)~=0
            orientation=9-i_positional_element;
            cont=get_cont_curved_3d(elements(i_positional_element),dt,conc_,conco_);
            for ix=1:1:3
                for iy=1:1:3
                    for iz=1:1:3
                        solution=solution+w(ix)*w(iy)*w(iz)...
                            *determinants(elements(i_positional_element),ix,iy,iz)*(...
                            weights(elements(i_positional_element),orientation,type,1,ix,iy,iz)*cont(ix,iy,iz)...
                            -diffT*weights(elements(i_positional_element),orientation,type,1,ix,iy,iz)*...
                            (terms(elements(i_positional_element),1,ix,iy,iz)+...
                            terms(elements(i_positional_element),2,ix,iy,iz))+...
                            terms(elements(i_positional_element),3,ix,iy,iz)*...
                            wTerms(elements(i_positional_element),orientation,type,ix,iy,iz)...
                            );
                    end
                end
            end
        end
    end
end
if ~isnan(eExtruding_x)
    for i_positional_element=1:1:4
        if eExtruding_x(i_positional_element)~=0
            orientation=xOrientations(i_positional_element);
            cont=get_cont_curved_3d(eExtruding_x(i_positional_element),dt,conc_,conco_);
            for ix=1:1:3
                for iy=1:1:3
                    for iz=1:1:3
                        solution=solution+w(ix)*w(iy)*w(iz)...
                            *determinants(eExtruding_x(i_positional_element),ix,iy,iz)*(...
                            weights(eExtruding_x(i_positional_element),orientation,type,1,ix,iy,iz)*cont(ix,iy,iz)...
                            -diffT*weights(eExtruding_x(i_positional_element),orientation,type,1,ix,iy,iz)*...
                            (terms(eExtruding_x(i_positional_element),1,ix,iy,iz)+...
                            terms(eExtruding_x(i_positional_element),2,ix,iy,iz))+...
                            terms(eExtruding_x(i_positional_element),3,ix,iy,iz)*...
                            wTerms(eExtruding_x(i_positional_element),orientation,type,ix,iy,iz)...
                            );
                    end
                end
            end
        end
    end
end
if ~isnan(eExtruding_y)
    for i_positional_element=1:1:4
        if eExtruding_y(i_positional_element)~=0
            orientation=yOrientations(i_positional_element);
            cont=get_cont_curved_3d(eExtruding_y(i_positional_element),dt,conc_,conco_);
            for ix=1:1:3
                for iy=1:1:3
                    for iz=1:1:3
                        solution=solution+w(ix)*w(iy)*w(iz)...
                            *determinants(eExtruding_y(i_positional_element),ix,iy,iz)*(...
                            weights(eExtruding_y(i_positional_element),orientation,type,1,ix,iy,iz)*cont(ix,iy,iz)...
                            -diffT*weights(eExtruding_y(i_positional_element),orientation,type,1,ix,iy,iz)*...
                            (terms(eExtruding_y(i_positional_element),1,ix,iy,iz)+...
                            terms(eExtruding_y(i_positional_element),2,ix,iy,iz))+...
                            terms(eExtruding_y(i_positional_element),3,ix,iy,iz)*...
                            wTerms(eExtruding_y(i_positional_element),orientation,type,ix,iy,iz)...
                            );
                    end
                end
            end
        end
    end
end
if ~isnan(eExtruding_z)
    for i_positional_element=1:1:4
        if eExtruding_z(i_positional_element)~=0
            orientation=zOrientations(i_positional_element); 
            cont=get_cont_curved_3d(eExtruding_z(i_positional_element),dt,conc_,conco_);
            for ix=1:1:3
                for iy=1:1:3
                    for iz=1:1:3
                        solution=solution+w(ix)*w(iy)*w(iz)...
                            *determinants(eExtruding_z(i_positional_element),ix,iy,iz)*(...
                            weights(eExtruding_z(i_positional_element),orientation,type,1,ix,iy,iz)*cont(ix,iy,iz)...
                            -diffT*weights(eExtruding_z(i_positional_element),orientation,type,1,ix,iy,iz)*...
                            (terms(eExtruding_z(i_positional_element),1,ix,iy,iz)+...
                            terms(eExtruding_z(i_positional_element),2,ix,iy,iz))+...
                            terms(eExtruding_z(i_positional_element),3,ix,iy,iz)*...
                            wTerms(eExtruding_z(i_positional_element),orientation,type,ix,iy,iz)...
                            );
                    end
                end
            end
        end
    end
end

end