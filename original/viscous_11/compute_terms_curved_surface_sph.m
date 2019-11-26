function terms_surface=compute_terms_curved_surface_sph(...
    ne,nex,ney,nexney,nx,ny,nz,n,nxny,n_eSurf,n_nSurf,nx_,weights_surface,c)

terms_surface=zeros(n_eSurf(6),3,3);

parfor e=1:1:n_eSurf(6)
    terms_surface(e,:,:)=get_conc_surface_sph(e,ne,nex,ney,nexney,nx,ny,nz,n,nxny,...
        n_eSurf,n_nSurf,nx_,weights_surface,c);
end

% weights_surface=zeros(n_eSurf(6),8,8,7,3,3);
end

function sol=get_conc_surface_sph(e,ne,nex,ney,nexney,nx,ny,nz,n,nxny,...
    n_eSurf,n_nSurf,nx_,weights_surface,c)
sol=zeros(3,3);
gbfs=get_gbfs_of_element_sph(e+ne,ne,nex,ney,nexney,nx,ny,nz,n,nxny,n_eSurf,n_nSurf,nx_);
for iorder=5:1:7
    igbf=0;
    for iorientation=1:1:8
        for itype=1:1:8
            igbf=igbf+1;
            for dim1=1:1:3
                for dim2=1:1:3
                    sol(dim1,dim2)=sol(dim1,dim2)...
                        +weights_surface(e,iorientation,itype,iorder,dim1,dim2)...
                        *c(gbfs(igbf));
                end
            end
        end
    end
end
end