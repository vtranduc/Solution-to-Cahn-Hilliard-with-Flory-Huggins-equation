function sj=compute_sj_3d(gbfs1,gbfs2,nnz_,conc_,dt,weights,nex,ney,nx,ny,nz,n1,n2,dxyz,diffT,chi,wTerms)

sj=zeros(1,nnz_);

w=[5/18 4/9 5/18];

%----
% diff=1800;
% T=0.6;
% chi=(0.5-1*(1-1/T))*ones(1,3);
% diffT=diff*T*ones(1,3);
% size(weights)
% size(conc)
%-------

if isscalar(diffT)
    parfor sjth=1:1:nnz_ %Parfor must be feasible!
        sj(sjth)=compute_sjth_sj_3d(sjth,gbfs1,gbfs2,...
            conc_,w,dt,diffT,chi,weights,nex,ney,nx,ny,nz,n1,n2,dxyz,wTerms);
    end
else
    parfor sjth=1:1:nnz_
        sj(sjth)=compute_sjth_sj_gradient_3d(sjth,gbfs1,gbfs2,...
            conc_,w,dt,diffT,chi,weights,nex,ney,nx,ny,nz,n1,n2,dxyz,wTerms);
    end
end

end

% === Temperature gradient is present=======================================

function solution=compute_sjth_sj_gradient_3d(sjth,gbfs1,gbfs2,...
    conc_,w,dt,diffT,chi,weights,nex,ney,nx,ny,nz,n1,n2,dxyz,wTerms)

solution=0;

[elements,lbfs1,lbfs2]=analyze_interaction_3d(...
	gbfs1(sjth),gbfs2(sjth),nex,ney,nx,ny,nz);

for i_positional_element=1:1:8
%     element=;
    if elements(i_positional_element)~=0
        
%         con=dim_reducer2_3d(conc_(elements(i_positional_element),1,:,:,:));
%         conx=dim_reducer2_3d(conc_(elements(i_positional_element),2,:,:,:));
%         cony=dim_reducer2_3d(conc_(elements(i_positional_element),3,:,:,:));
%         conz=dim_reducer2_3d(conc_(elements(i_positional_element),4,:,:,:));
%         conxx=dim_reducer2_3d(conc_(elements(i_positional_element),5,:,:,:));
%         conyy=dim_reducer2_3d(conc_(elements(i_positional_element),6,:,:,:));
%         conzz=dim_reducer2_3d(conc_(elements(i_positional_element),7,:,:,:));
%         
        [orientation1,type1]=analyze_lbf_3d(lbfs1(i_positional_element));
        [orientation2,type2]=analyze_lbf_3d(lbfs2(i_positional_element));
        
%         phi1=dim_reducer3_3d(weights(orientation1,type1,1,:,:,:));
%         phixx1=dim_reducer3_3d(weights(orientation1,type1,5,:,:,:));
%         phiyy1=dim_reducer3_3d(weights(orientation1,type1,6,:,:,:));
%         phizz1=dim_reducer3_3d(weights(orientation1,type1,7,:,:,:));
%         
%         phi2=dim_reducer3_3d(weights(orientation2,type2,1,:,:,:));
%         phixx2=dim_reducer3_3d(weights(orientation2,type2,5,:,:,:));
%         phiyy2=dim_reducer3_3d(weights(orientation2,type2,6,:,:,:));
%         phizz2=dim_reducer3_3d(weights(orientation2,type2,7,:,:,:));

        n1xth=get_n1xth_3d(elements(i_positional_element),nex);
        
        for ix=1:1:3
            for iy=1:1:3
                for iz=1:1:3
                    
                    solution=solution+w(ix)*w(iy)*w(iz)*(...
                        weights(orientation1,type1,1,ix,iy,iz)*weights(orientation2,type2,1,ix,iy,iz)/dt...
                        -diffT(n1xth,ix)*weights(orientation1,type1,1,ix,iy,iz)...
                        *(2*weights(orientation2,type2,1,ix,iy,iz)...
                        *((conc_(elements(i_positional_element),1,ix,iy,iz)^3*n1)^-1+((1-conc_(elements(i_positional_element),1,ix,iy,iz))^3*n2)^-1)...
                        *(conc_(elements(i_positional_element),2,ix,iy,iz)^2+conc_(elements(i_positional_element),3,ix,iy,iz)^2+conc_(elements(i_positional_element),4,ix,iy,iz)^2)...
                        +weights(orientation2,type2,1,ix,iy,iz)...
                        *(-(conc_(elements(i_positional_element),1,ix,iy,iz)^2*n1)^-1+((1-conc_(elements(i_positional_element),1,ix,iy,iz))^2*n2)^-1)...
                        *(2*(conc_(elements(i_positional_element),2,ix,iy,iz)+conc_(elements(i_positional_element),3,ix,iy,iz)+conc_(elements(i_positional_element),4,ix,iy,iz))...
                        +(conc_(elements(i_positional_element),5,ix,iy,iz)+conc_(elements(i_positional_element),6,ix,iy,iz)+conc_(elements(i_positional_element),7,ix,iy,iz)))...
                        +((conc_(elements(i_positional_element),1,ix,iy,iz)*n1)^-1+((1-conc_(elements(i_positional_element),1,ix,iy,iz))*n2)^-1-2*chi(n1xth,ix))...
                        *wTerms(orientation2,type2,ix,iy,iz))...
                        +wTerms(orientation1,type1,ix,iy,iz)...
                        *wTerms(orientation2,type2,ix,iy,iz)...
                        ...
                        );
                    
                end
            end
        end
    end
end
solution=solution*dxyz;
end

% === No temperature gradient =============================================

function solution=compute_sjth_sj_3d(sjth,gbfs1,gbfs2,...
    conc_,w,dt,diffT,chi,weights,nex,ney,nx,ny,nz,n1,n2,dxyz,wTerms)

solution=0;

[elements,lbfs1,lbfs2]=analyze_interaction_3d(...
	gbfs1(sjth),gbfs2(sjth),nex,ney,nx,ny,nz);

for i_positional_element=1:1:8
%     element=;
    if elements(i_positional_element)~=0
        
%         con=dim_reducer2_3d(conc_(elements(i_positional_element),1,:,:,:));
%         conx=dim_reducer2_3d(conc_(elements(i_positional_element),2,:,:,:));
%         cony=dim_reducer2_3d(conc_(elements(i_positional_element),3,:,:,:));
%         conz=dim_reducer2_3d(conc_(elements(i_positional_element),4,:,:,:));
%         conxx=dim_reducer2_3d(conc_(elements(i_positional_element),5,:,:,:));
%         conyy=dim_reducer2_3d(conc_(elements(i_positional_element),6,:,:,:));
%         conzz=dim_reducer2_3d(conc_(elements(i_positional_element),7,:,:,:));
%         
        [orientation1,type1]=analyze_lbf_3d(lbfs1(i_positional_element));
        [orientation2,type2]=analyze_lbf_3d(lbfs2(i_positional_element));
        
%         phi1=dim_reducer3_3d(weights(orientation1,type1,1,:,:,:));
%         phixx1=dim_reducer3_3d(weights(orientation1,type1,5,:,:,:));
%         phiyy1=dim_reducer3_3d(weights(orientation1,type1,6,:,:,:));
%         phizz1=dim_reducer3_3d(weights(orientation1,type1,7,:,:,:));
%         
%         phi2=dim_reducer3_3d(weights(orientation2,type2,1,:,:,:));
%         phixx2=dim_reducer3_3d(weights(orientation2,type2,5,:,:,:));
%         phiyy2=dim_reducer3_3d(weights(orientation2,type2,6,:,:,:));
%         phizz2=dim_reducer3_3d(weights(orientation2,type2,7,:,:,:));
        
        for ix=1:1:3
            for iy=1:1:3
                for iz=1:1:3
                    
                    solution=solution+w(ix)*w(iy)*w(iz)*(...
                        weights(orientation1,type1,1,ix,iy,iz)*weights(orientation2,type2,1,ix,iy,iz)/dt...
                        -diffT*weights(orientation1,type1,1,ix,iy,iz)...
                        *(2*weights(orientation2,type2,1,ix,iy,iz)...
                        *((conc_(elements(i_positional_element),1,ix,iy,iz)^3*n1)^-1+((1-conc_(elements(i_positional_element),1,ix,iy,iz))^3*n2)^-1)...
                        *(conc_(elements(i_positional_element),2,ix,iy,iz)^2+conc_(elements(i_positional_element),3,ix,iy,iz)^2+conc_(elements(i_positional_element),4,ix,iy,iz)^2)...
                        +weights(orientation2,type2,1,ix,iy,iz)...
                        *(-(conc_(elements(i_positional_element),1,ix,iy,iz)^2*n1)^-1+((1-conc_(elements(i_positional_element),1,ix,iy,iz))^2*n2)^-1)...
                        *(2*(conc_(elements(i_positional_element),2,ix,iy,iz)+conc_(elements(i_positional_element),3,ix,iy,iz)+conc_(elements(i_positional_element),4,ix,iy,iz))...
                        +(conc_(elements(i_positional_element),5,ix,iy,iz)+conc_(elements(i_positional_element),6,ix,iy,iz)+conc_(elements(i_positional_element),7,ix,iy,iz)))...
                        +((conc_(elements(i_positional_element),1,ix,iy,iz)*n1)^-1+((1-conc_(elements(i_positional_element),1,ix,iy,iz))*n2)^-1-2*chi)...
                        *wTerms(orientation2,type2,ix,iy,iz))...
                        +wTerms(orientation1,type1,ix,iy,iz)...
                        *wTerms(orientation2,type2,ix,iy,iz)...
                        ...
                        );
                    
                end
            end
        end
    end
end
solution=solution*dxyz;
end