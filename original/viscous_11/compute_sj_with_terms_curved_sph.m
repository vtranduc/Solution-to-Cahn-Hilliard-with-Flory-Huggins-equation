function sj=compute_sj_with_terms_curved_sph(...
    gbfs1,gbfs2,nnz_,dt,weights,nex,ney,nx,ny,nz,...
    determinants,diffT,wTerms,terms,...
    ne,nez,n,neight,nexney,n_eSurf,n_nSurf,nx_,...
    wTerms_surface,weights_surface,determinants_surface,normals)

sj=zeros(1,nnz_);

w=[5/18 4/9 5/18];

if isscalar(diffT)
    parfor sjth=1:1:nnz_ %Parfor must be feasible!
        sj(sjth)=compute_sjth_sj_with_terms_sph(sjth,gbfs1,gbfs2,...
            w,dt,diffT,weights,nex,ney,nx,ny,nz,determinants,wTerms,terms,...
            ne,nez,n,neight,nexney,n_eSurf,n_nSurf,nx_,...
            wTerms_surface,weights_surface,determinants_surface,normals);
    end
% else
%     parfor sjth=1:1:nnz_
%         sj(sjth)=compute_sjth_sj_gradient_with_terms_3d(sjth,gbfs1,gbfs2,...
%             w,dt,diffT,weights,nex,ney,nx,ny,nz,determinants,wTerms,terms);
%     end
end

end

% === Temperature gradient is present =============================================

% function solution=compute_sjth_sj_gradient_with_terms_3d(sjth,gbfs1,gbfs2,...
%     w,dt,diffT,weights,nex,ney,nx,ny,nz,determinants,wTerms,terms)
% 
% solution=0;
% 
% [elements,lbfs1,lbfs2]=analyze_interaction_3d(...
% 	gbfs1(sjth),gbfs2(sjth),nex,ney,nx,ny,nz);
% 
% for i_positional_element=1:1:8
% 
%     if elements(i_positional_element)~=0
%      
%         [orientation1,type1]=analyze_lbf_3d(lbfs1(i_positional_element));
%         [orientation2,type2]=analyze_lbf_3d(lbfs2(i_positional_element));
%         
%         n1xth=get_n1xth_3d(elements(i_positional_element),nex);
%         
%         for ix=1:1:3
%             for iy=1:1:3
%                 for iz=1:1:3
%                     
%                     solution=solution+w(ix)*w(iy)*w(iz)...
%                         *determinants(elements(i_positional_element),ix,iy,iz)*(...
%                         weights(elements(i_positional_element),orientation1,type1,1,ix,iy,iz)*weights(elements(i_positional_element),orientation2,type2,1,ix,iy,iz)/dt...
%                         -diffT(elements(i_positional_element))*weights(elements(i_positional_element),orientation1,type1,1,ix,iy,iz)...
%                         *(weights(elements(i_positional_element),orientation2,type2,1,ix,iy,iz)...
%                         *terms(elements(i_positional_element),4,ix,iy,iz)...
%                         +weights(elements(i_positional_element),orientation2,type2,1,ix,iy,iz)...
%                         *terms(elements(i_positional_element),5,ix,iy,iz)...
%                         +terms(elements(i_positional_element),6,ix,iy,iz)...
%                         *wTerms(elements(i_positional_element),orientation2,type2,ix,iy,iz))...
%                         +wTerms(elements(i_positional_element),orientation1,type1,ix,iy,iz)...
%                         *wTerms(elements(i_positional_element),orientation2,type2,ix,iy,iz)...
%                         ...
%                         );
%                 end
%             end
%         end
%     end
% end
% % solution=solution*dxyz;
% end

% === No temperature gradient =============================================

function solution=compute_sjth_sj_with_terms_sph(sjth,gbfs1,gbfs2,...
    w,dt,diffT,weights,nex,ney,nx,ny,nz,determinants,wTerms,terms,...
    ne,nez,n,neight,nexney,n_eSurf,n_nSurf,nx_,...
    wTerms_surface,weights_surface,determinants_surface,normals)

solution=0;

% [elements,lbfs1,lbfs2]=analyze_interaction_3d(...
% 	gbfs1(sjth),gbfs2(sjth),nex,ney,nx,ny,nz);

[elements,lbfs1,lbfs2,mx1,mx2,mx3,mx4,mx5,mx6]=...
    analyze_interaction_sph(gbfs1(sjth),gbfs2(sjth),ne,nex,ney,nez,n,nx,ny,nz,neight,...
    nexney,n_eSurf,n_nSurf,nx_);

% disp('fasdfasfasdfasfd')

if ~isnan(elements)
    for i_positional_element=1:1:8
        if elements(i_positional_element)~=0
            [orientation1,type1]=analyze_lbf_3d(lbfs1(i_positional_element));
            [orientation2,type2]=analyze_lbf_3d(lbfs2(i_positional_element));
            for ix=1:1:3
                for iy=1:1:3
                    for iz=1:1:3
                        solution=solution+w(ix)*w(iy)*w(iz)...
                            *determinants(elements(i_positional_element),ix,iy,iz)*(...
                            weights(elements(i_positional_element),orientation1,type1,1,ix,iy,iz)*weights(elements(i_positional_element),orientation2,type2,1,ix,iy,iz)/dt...
                            -diffT*weights(elements(i_positional_element),orientation1,type1,1,ix,iy,iz)...
                            *(weights(elements(i_positional_element),orientation2,type2,1,ix,iy,iz)...
                            *terms(elements(i_positional_element),4,ix,iy,iz)...
                            +weights(elements(i_positional_element),orientation2,type2,1,ix,iy,iz)...
                            *terms(elements(i_positional_element),5,ix,iy,iz)...
                            +terms(elements(i_positional_element),6,ix,iy,iz)...
                            *wTerms(elements(i_positional_element),orientation2,type2,ix,iy,iz))...
                            +wTerms(elements(i_positional_element),orientation1,type1,ix,iy,iz)...
                            *wTerms(elements(i_positional_element),orientation2,type2,ix,iy,iz)...
                            ...
                            );
                    end
                end
            end
        end
    end
end

%-------------SURFACE----------------------------------

% [node1,~]=analyze_gbs_3d(gbfs1(sjth));
% 
% if node1>n
%     if ~isnan(mx1)
%         for i_positional_element=1:1:4
%             if mx1(i_positional_element,1)~=0
%                 e_surface=mx1(i_positional_element,1)-ne;
%                 [orientation1,type1]=analyze_lbf_3d(mx1(i_positional_element,2));
%                 [orientation2,type2]=analyze_lbf_3d(mx1(i_positional_element,3));
%                 for dim1=1:1:3
%                     for dim2=1:1:3
%                         solution=solution-w(dim1)*w(dim2)...
%                             *determinants_surface(e_surface,dim1,dim2)...
%                             *wTerms_surface(e_surface,orientation2,type2,dim1,dim2)...
%                             *(normals(e_surface,1)*weights_surface(e_surface,orientation1,type1,1,dim1,dim2)...
%                             +normals(e_surface,2)*weights_surface(e_surface,orientation1,type1,2,dim1,dim2)...
%                             +normals(e_surface,3)*weights_surface(e_surface,orientation1,type1,3,dim1,dim2));
%                     end
%                 end
%             end
%         end
%     end
%     if ~isnan(mx2)
%         for i_positional_element=1:1:4
%             if mx2(i_positional_element,1)~=0
%                 e_surface=mx2(i_positional_element,1)-ne;
%                 [orientation1,type1]=analyze_lbf_3d(mx2(i_positional_element,2));
%                 [orientation2,type2]=analyze_lbf_3d(mx2(i_positional_element,3));
%                 for dim1=1:1:3
%                     for dim2=1:1:3
%                         solution=solution-w(dim1)*w(dim2)...
%                             *determinants_surface(e_surface,dim1,dim2)...
%                             *wTerms_surface(e_surface,orientation2,type2,dim1,dim2)...
%                             *(normals(e_surface,1)*weights_surface(e_surface,orientation1,type1,1,dim1,dim2)...
%                             +normals(e_surface,2)*weights_surface(e_surface,orientation1,type1,2,dim1,dim2)...
%                             +normals(e_surface,3)*weights_surface(e_surface,orientation1,type1,3,dim1,dim2));
%                     end
%                 end
%             end
%         end
%     end
%     if ~isnan(mx3)
%         for i_positional_element=1:1:4
%             if mx3(i_positional_element,1)~=0
%                 e_surface=mx3(i_positional_element,1)-ne;
%                 [orientation1,type1]=analyze_lbf_3d(mx3(i_positional_element,2));
%                 [orientation2,type2]=analyze_lbf_3d(mx3(i_positional_element,3));
%                 for dim1=1:1:3
%                     for dim2=1:1:3
%                         solution=solution-w(dim1)*w(dim2)...
%                             *determinants_surface(e_surface,dim1,dim2)...
%                             *wTerms_surface(e_surface,orientation2,type2,dim1,dim2)...
%                             *(normals(e_surface,1)*weights_surface(e_surface,orientation1,type1,1,dim1,dim2)...
%                             +normals(e_surface,2)*weights_surface(e_surface,orientation1,type1,2,dim1,dim2)...
%                             +normals(e_surface,3)*weights_surface(e_surface,orientation1,type1,3,dim1,dim2));
%                     end
%                 end
%             end
%         end
%     end
%     if ~isnan(mx4)
%         for i_positional_element=1:1:4
%             if mx4(i_positional_element,1)~=0
%                 e_surface=mx4(i_positional_element,1)-ne;
%                 [orientation1,type1]=analyze_lbf_3d(mx4(i_positional_element,2));
%                 [orientation2,type2]=analyze_lbf_3d(mx4(i_positional_element,3));
%                 for dim1=1:1:3
%                     for dim2=1:1:3
%                         solution=solution-w(dim1)*w(dim2)...
%                             *determinants_surface(e_surface,dim1,dim2)...
%                             *wTerms_surface(e_surface,orientation2,type2,dim1,dim2)...
%                             *(normals(e_surface,1)*weights_surface(e_surface,orientation1,type1,1,dim1,dim2)...
%                             +normals(e_surface,2)*weights_surface(e_surface,orientation1,type1,2,dim1,dim2)...
%                             +normals(e_surface,3)*weights_surface(e_surface,orientation1,type1,3,dim1,dim2));
%                     end
%                 end
%             end
%         end
%     end
%     if ~isnan(mx5)
%         for i_positional_element=1:1:4
%             if mx5(i_positional_element,1)~=0
%                 e_surface=mx5(i_positional_element,1)-ne;
%                 [orientation1,type1]=analyze_lbf_3d(mx5(i_positional_element,2));
%                 [orientation2,type2]=analyze_lbf_3d(mx5(i_positional_element,3));
%                 for dim1=1:1:3
%                     for dim2=1:1:3
%                         solution=solution-w(dim1)*w(dim2)...
%                             *determinants_surface(e_surface,dim1,dim2)...
%                             *wTerms_surface(e_surface,orientation2,type2,dim1,dim2)...
%                             *(normals(e_surface,1)*weights_surface(e_surface,orientation1,type1,1,dim1,dim2)...
%                             +normals(e_surface,2)*weights_surface(e_surface,orientation1,type1,2,dim1,dim2)...
%                             +normals(e_surface,3)*weights_surface(e_surface,orientation1,type1,3,dim1,dim2));
%                     end
%                 end
%             end
%         end
%     end
%     if ~isnan(mx6)
%         for i_positional_element=1:1:4
%             if mx6(i_positional_element,1)~=0
%                 e_surface=mx6(i_positional_element,1)-ne;
%                 [orientation1,type1]=analyze_lbf_3d(mx6(i_positional_element,2));
%                 [orientation2,type2]=analyze_lbf_3d(mx6(i_positional_element,3));
%                 for dim1=1:1:3
%                     for dim2=1:1:3
%                         solution=solution-w(dim1)*w(dim2)...
%                             *determinants_surface(e_surface,dim1,dim2)...
%                             *wTerms_surface(e_surface,orientation2,type2,dim1,dim2)...
%                             *(normals(e_surface,1)*weights_surface(e_surface,orientation1,type1,1,dim1,dim2)...
%                             +normals(e_surface,2)*weights_surface(e_surface,orientation1,type1,2,dim1,dim2)...
%                             +normals(e_surface,3)*weights_surface(e_surface,orientation1,type1,3,dim1,dim2));
%                     end
%                 end
%             end
%         end
%     end
% end
% 
% % if gbfs1(sjth)>neight
% %     if ~isnan(mx1)
% %         for i_positional_element=1:1:4
% %             if mx1(i_positional_element,1)~=0
% %                 [orientation1,type1]=analyze_lbf_3d(mx1(i_positional_element,2));
% %                 [orientation2,type2]=analyze_lbf_3d(mx1(i_positional_element,3));
% %                 
% % 
% %                 for iy=1:1:3
% %                     for iz=1:1:3
% %                         solution=solution+w(iy)*w(iz)...
% %                             *determinants_surface(mx1(i_positional_element,1)-ne,iy,iz)*(...
% %                             weights(mx1(i_positional_element,1),orientation1,type1,1,ix,iy,iz)...
% %                             *weights(mx1(i_positional_element,1),orientation2,type2,1,ix,iy,iz)/dt...
% %                             -diffT*weights(mx1(i_positional_element,1),orientation1,type1,1,ix,iy,iz)...
% %                             *(weights(mx1(i_positional_element,1),orientation2,type2,1,ix,iy,iz)...
% %                             *terms(mx1(i_positional_element,1),4,ix,iy,iz)...
% %                             +weights(mx1(i_positional_element,1),orientation2,type2,1,ix,iy,iz)...
% %                             *terms(mx1(i_positional_element,1),5,ix,iy,iz)...
% %                             +terms(mx1(i_positional_element,1),6,ix,iy,iz)...
% %                             *wTerms(mx1(i_positional_element,1),orientation2,type2,ix,iy,iz))...
% %                             +wTerms(mx1(i_positional_element,1),orientation1,type1,ix,iy,iz)...
% %                             *wTerms(mx1(i_positional_element,1),orientation2,type2,ix,iy,iz)...
% %                             ...
% %                             );
% %                     end
% %                 end
% % 
% %                 
% %                 
% %             end
% %         end
% %     end
% % end
% 
% 
% 
% 
% 

%----------------------------------------------------------

% return
% 
if ~isnan(mx1)
    for i_positional_element=1:1:4
        if mx1(i_positional_element,1)~=0
            [orientation1,type1]=analyze_lbf_3d(mx1(i_positional_element,2));
            [orientation2,type2]=analyze_lbf_3d(mx1(i_positional_element,3));
            for ix=1:1:3
                for iy=1:1:3
                    for iz=1:1:3
                        solution=solution+w(ix)*w(iy)*w(iz)...
                            *determinants(mx1(i_positional_element,1),ix,iy,iz)*(...
                            weights(mx1(i_positional_element,1),orientation1,type1,1,ix,iy,iz)...
                            *weights(mx1(i_positional_element,1),orientation2,type2,1,ix,iy,iz)/dt...
                            -diffT*weights(mx1(i_positional_element,1),orientation1,type1,1,ix,iy,iz)...
                            *(weights(mx1(i_positional_element,1),orientation2,type2,1,ix,iy,iz)...
                            *terms(mx1(i_positional_element,1),4,ix,iy,iz)...
                            +weights(mx1(i_positional_element,1),orientation2,type2,1,ix,iy,iz)...
                            *terms(mx1(i_positional_element,1),5,ix,iy,iz)...
                            +terms(mx1(i_positional_element,1),6,ix,iy,iz)...
                            *wTerms(mx1(i_positional_element,1),orientation2,type2,ix,iy,iz))...
                            +wTerms(mx1(i_positional_element,1),orientation1,type1,ix,iy,iz)...
                            *wTerms(mx1(i_positional_element,1),orientation2,type2,ix,iy,iz)...
                            ...
                            );
                    end
                end
            end
        end
    end
end

if ~isnan(mx2)
    for i_positional_element=1:1:4
        if mx2(i_positional_element,1)~=0
            [orientation1,type1]=analyze_lbf_3d(mx2(i_positional_element,2));
            [orientation2,type2]=analyze_lbf_3d(mx2(i_positional_element,3));
            for ix=1:1:3
                for iy=1:1:3
                    for iz=1:1:3
                        solution=solution+w(ix)*w(iy)*w(iz)...
                            *determinants(mx2(i_positional_element,1),ix,iy,iz)*(...
                            weights(mx2(i_positional_element,1),orientation1,type1,1,ix,iy,iz)...
                            *weights(mx2(i_positional_element,1),orientation2,type2,1,ix,iy,iz)/dt...
                            -diffT*weights(mx2(i_positional_element,1),orientation1,type1,1,ix,iy,iz)...
                            *(weights(mx2(i_positional_element,1),orientation2,type2,1,ix,iy,iz)...
                            *terms(mx2(i_positional_element,1),4,ix,iy,iz)...
                            +weights(mx2(i_positional_element,1),orientation2,type2,1,ix,iy,iz)...
                            *terms(mx2(i_positional_element,1),5,ix,iy,iz)...
                            +terms(mx2(i_positional_element,1),6,ix,iy,iz)...
                            *wTerms(mx2(i_positional_element,1),orientation2,type2,ix,iy,iz))...
                            +wTerms(mx2(i_positional_element,1),orientation1,type1,ix,iy,iz)...
                            *wTerms(mx2(i_positional_element,1),orientation2,type2,ix,iy,iz)...
                            ...
                            );
                    end
                end
            end
        end
    end
end

if ~isnan(mx3)
    for i_positional_element=1:1:4
        if mx3(i_positional_element,1)~=0
            [orientation1,type1]=analyze_lbf_3d(mx3(i_positional_element,2));
            [orientation2,type2]=analyze_lbf_3d(mx3(i_positional_element,3));
            for ix=1:1:3
                for iy=1:1:3
                    for iz=1:1:3
                        solution=solution+w(ix)*w(iy)*w(iz)...
                            *determinants(mx3(i_positional_element,1),ix,iy,iz)*(...
                            weights(mx3(i_positional_element,1),orientation1,type1,1,ix,iy,iz)...
                            *weights(mx3(i_positional_element,1),orientation2,type2,1,ix,iy,iz)/dt...
                            -diffT*weights(mx3(i_positional_element,1),orientation1,type1,1,ix,iy,iz)...
                            *(weights(mx3(i_positional_element,1),orientation2,type2,1,ix,iy,iz)...
                            *terms(mx3(i_positional_element,1),4,ix,iy,iz)...
                            +weights(mx3(i_positional_element,1),orientation2,type2,1,ix,iy,iz)...
                            *terms(mx3(i_positional_element,1),5,ix,iy,iz)...
                            +terms(mx3(i_positional_element,1),6,ix,iy,iz)...
                            *wTerms(mx3(i_positional_element,1),orientation2,type2,ix,iy,iz))...
                            +wTerms(mx3(i_positional_element,1),orientation1,type1,ix,iy,iz)...
                            *wTerms(mx3(i_positional_element,1),orientation2,type2,ix,iy,iz)...
                            ...
                            );
                    end
                end
            end
        end
    end
end

if ~isnan(mx4)
    for i_positional_element=1:1:4
        if mx4(i_positional_element,1)~=0
            [orientation1,type1]=analyze_lbf_3d(mx4(i_positional_element,2));
            [orientation2,type2]=analyze_lbf_3d(mx4(i_positional_element,3));
            for ix=1:1:3
                for iy=1:1:3
                    for iz=1:1:3
                        solution=solution+w(ix)*w(iy)*w(iz)...
                            *determinants(mx4(i_positional_element,1),ix,iy,iz)*(...
                            weights(mx4(i_positional_element,1),orientation1,type1,1,ix,iy,iz)...
                            *weights(mx4(i_positional_element,1),orientation2,type2,1,ix,iy,iz)/dt...
                            -diffT*weights(mx4(i_positional_element,1),orientation1,type1,1,ix,iy,iz)...
                            *(weights(mx4(i_positional_element,1),orientation2,type2,1,ix,iy,iz)...
                            *terms(mx4(i_positional_element,1),4,ix,iy,iz)...
                            +weights(mx4(i_positional_element,1),orientation2,type2,1,ix,iy,iz)...
                            *terms(mx4(i_positional_element,1),5,ix,iy,iz)...
                            +terms(mx4(i_positional_element,1),6,ix,iy,iz)...
                            *wTerms(mx4(i_positional_element,1),orientation2,type2,ix,iy,iz))...
                            +wTerms(mx4(i_positional_element,1),orientation1,type1,ix,iy,iz)...
                            *wTerms(mx4(i_positional_element,1),orientation2,type2,ix,iy,iz)...
                            ...
                            );
                    end
                end
            end
        end
    end
end

if ~isnan(mx5)
    for i_positional_element=1:1:4
        if mx5(i_positional_element,1)~=0
            [orientation1,type1]=analyze_lbf_3d(mx5(i_positional_element,2));
            [orientation2,type2]=analyze_lbf_3d(mx5(i_positional_element,3));
            for ix=1:1:3
                for iy=1:1:3
                    for iz=1:1:3
                        solution=solution+w(ix)*w(iy)*w(iz)...
                            *determinants(mx5(i_positional_element,1),ix,iy,iz)*(...
                            weights(mx5(i_positional_element,1),orientation1,type1,1,ix,iy,iz)...
                            *weights(mx5(i_positional_element,1),orientation2,type2,1,ix,iy,iz)/dt...
                            -diffT*weights(mx5(i_positional_element,1),orientation1,type1,1,ix,iy,iz)...
                            *(weights(mx5(i_positional_element,1),orientation2,type2,1,ix,iy,iz)...
                            *terms(mx5(i_positional_element,1),4,ix,iy,iz)...
                            +weights(mx5(i_positional_element,1),orientation2,type2,1,ix,iy,iz)...
                            *terms(mx5(i_positional_element,1),5,ix,iy,iz)...
                            +terms(mx5(i_positional_element,1),6,ix,iy,iz)...
                            *wTerms(mx5(i_positional_element,1),orientation2,type2,ix,iy,iz))...
                            +wTerms(mx5(i_positional_element,1),orientation1,type1,ix,iy,iz)...
                            *wTerms(mx5(i_positional_element,1),orientation2,type2,ix,iy,iz)...
                            ...
                            );
                    end
                end
            end
        end
    end
end

if ~isnan(mx6)
    for i_positional_element=1:1:4
        if mx6(i_positional_element,1)~=0
            [orientation1,type1]=analyze_lbf_3d(mx6(i_positional_element,2));
            [orientation2,type2]=analyze_lbf_3d(mx6(i_positional_element,3));
            for ix=1:1:3
                for iy=1:1:3
                    for iz=1:1:3
                        solution=solution+w(ix)*w(iy)*w(iz)...
                            *determinants(mx6(i_positional_element,1),ix,iy,iz)*(...
                            weights(mx6(i_positional_element,1),orientation1,type1,1,ix,iy,iz)...
                            *weights(mx6(i_positional_element,1),orientation2,type2,1,ix,iy,iz)/dt...
                            -diffT*weights(mx6(i_positional_element,1),orientation1,type1,1,ix,iy,iz)...
                            *(weights(mx6(i_positional_element,1),orientation2,type2,1,ix,iy,iz)...
                            *terms(mx6(i_positional_element,1),4,ix,iy,iz)...
                            +weights(mx6(i_positional_element,1),orientation2,type2,1,ix,iy,iz)...
                            *terms(mx6(i_positional_element,1),5,ix,iy,iz)...
                            +terms(mx6(i_positional_element,1),6,ix,iy,iz)...
                            *wTerms(mx6(i_positional_element,1),orientation2,type2,ix,iy,iz))...
                            +wTerms(mx6(i_positional_element,1),orientation1,type1,ix,iy,iz)...
                            *wTerms(mx6(i_positional_element,1),orientation2,type2,ix,iy,iz)...
                            ...
                            );
                    end
                end
            end
        end
    end
end
% solution=solution*dxyz;
end