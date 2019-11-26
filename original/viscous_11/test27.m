function sj=test27(gbfs1,gbfs2,nnz_,dt,weights,nex,ney,nx,ny,nz,dxyz,diffT,wTerms,terms,sj_properties)

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
        sj(sjth)=compute_sjth_sj_with_terms_3d(sjth,gbfs1,gbfs2,...
            w,dt,diffT,weights,nex,ney,nx,ny,nz,dxyz,wTerms,terms);
    end
else
    parfor sjth=1:1:nnz_
        sj(sjth)=compute_sjth_sj_gradient_with_terms_3d(sjth,gbfs1,gbfs2,...
            w,dt,diffT,weights,nex,ney,nx,ny,nz,dxyz,wTerms,terms,sj_properties);
    end
end

end

% === Temperature gradient is present =============================================

function solution=compute_sjth_sj_gradient_with_terms_3d(sjth,gbfs1,gbfs2,...
    w,dt,diffT,weights,nex,ney,nx,ny,nz,dxyz,wTerms,terms,sj_properties)

solution=0;

[elements,lbfs1,lbfs2]=analyze_interaction_3d(...
	gbfs1(sjth),gbfs2(sjth),nex,ney,nx,ny,nz);

for i_positional_element=1:1:8

    if elements(i_positional_element)~=0
     
        [orientation1,type1]=analyze_lbf_3d(lbfs1(i_positional_element));
        [orientation2,type2]=analyze_lbf_3d(lbfs2(i_positional_element));
        
        n1xth=get_n1xth_3d(elements(i_positional_element),nex);
        
        for ix=1:1:3
            for iy=1:1:3
                for iz=1:1:3
                    
                    solution=solution+w(ix)*w(iy)*w(iz)*(...
                        weights(orientation1,type1,1,ix,iy,iz)*weights(orientation2,type2,1,ix,iy,iz)/dt...
                        -diffT(n1xth,ix)*weights(orientation1,type1,1,ix,iy,iz)...
                        *(weights(orientation2,type2,1,ix,iy,iz)...
                        *terms(elements(i_positional_element),4,ix,iy,iz)...
                        +weights(orientation2,type2,1,ix,iy,iz)...
                        *terms(elements(i_positional_element),5,ix,iy,iz)...
                        +terms(elements(i_positional_element),6,ix,iy,iz)...
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

function solution=compute_sjth_sj_with_terms_3d(sjth,gbfs1,gbfs2,...
    w,dt,diffT,weights,nex,ney,nx,ny,nz,dxyz,wTerms,terms)

solution=0;

[elements,lbfs1,lbfs2]=analyze_interaction_3d(...
	gbfs1(sjth),gbfs2(sjth),nex,ney,nx,ny,nz);

for i_positional_element=1:1:8

    if elements(i_positional_element)~=0
     
        [orientation1,type1]=analyze_lbf_3d(lbfs1(i_positional_element));
        [orientation2,type2]=analyze_lbf_3d(lbfs2(i_positional_element));
        
        for ix=1:1:3
            for iy=1:1:3
                for iz=1:1:3
                    
                    solution=solution+w(ix)*w(iy)*w(iz)*(...
                        weights(orientation1,type1,1,ix,iy,iz)*weights(orientation2,type2,1,ix,iy,iz)/dt...
                        -diffT*weights(orientation1,type1,1,ix,iy,iz)...
                        *(weights(orientation2,type2,1,ix,iy,iz)...
                        *terms(elements(i_positional_element),4,ix,iy,iz)...
                        +weights(orientation2,type2,1,ix,iy,iz)...
                        *terms(elements(i_positional_element),5,ix,iy,iz)...
                        +terms(elements(i_positional_element),6,ix,iy,iz)...
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