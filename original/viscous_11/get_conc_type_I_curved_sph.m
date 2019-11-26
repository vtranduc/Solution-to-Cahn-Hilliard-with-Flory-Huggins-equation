function conc_=get_conc_type_I_curved_sph(neTotal,c,weights,ne,nexney,nex,ney,nx,ny,nz,nxny,n,n_eSurf,n_nSurf,nx_,parallel_computing)
conc_=zeros(neTotal,3,3,3);
if parallel_computing==1
    parfor e=1:1:neTotal
        conc_(e,:,:,:)=compute_elemental_weights_type_I_curved_sph(...
            e,c,weights,ne,nexney,nex,ney,nx,ny,nz,nxny,n,n_eSurf,n_nSurf,nx_);
    end
elseif parallel_computing==0
    for e=1:1:neTotal
        conc_(e,:,:,:)=compute_elemental_weights_type_I_curved_sph(...
            e,c,weights,ne,nexney,nex,ney,nx,ny,nz,nxny,n,n_eSurf,n_nSurf,nx_);
    end
end
end

function solution=compute_elemental_weights_type_I_curved_sph(e,c,weights,ne,nexney,nex,ney,nx,ny,nz,nxny,n,n_eSurf,n_nSurf,nx_,gbfs)
solution=zeros(3,3,3);
if nargin==15
    gbfs=get_gbfs_of_element_sph(e,ne,nex,ney,nexney,nx,ny,nz,n,nxny,n_eSurf,n_nSurf,nx_);
end

igbf=0;
for iorientation=1:1:8
    for itype=1:1:8
        igbf=igbf+1;
        for igpz=1:1:3
            for igpy=1:1:3
                for igpx=1:1:3
                    solution(igpx,igpy,igpz)=...
                        solution(igpx,igpy,igpz)...
                        +weights(e,iorientation,itype,1,igpx,igpy,igpz)...
                        *c(gbfs(igbf));
                end
            end
        end
    end
end

end