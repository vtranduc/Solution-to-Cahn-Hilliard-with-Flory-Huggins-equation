function solution=get_conc_curved_sph(neTotal,c,weights,ne,nexney,nex,ney,nx,ny,nz,nxny,n,n_eSurf,n_nSurf,nx_)
solution=zeros(neTotal,7,3,3,3);
parfor e=1:1:neTotal %parfor must be feasible!
    solution(e,:,:,:,:)=compute_elemental_weights_curved_sph(...
        e,c,weights,ne,nexney,nex,ney,nx,ny,nz,nxny,n,n_eSurf,n_nSurf,nx_);
end

end

function solution=compute_elemental_weights_curved_sph(e,c,weights,ne,nexney,nex,ney,nx,ny,nz,nxny,n,n_eSurf,n_nSurf,nx_,gbfs)
solution=zeros(7,3,3,3);
if nargin==15
    gbfs=get_gbfs_of_element_sph(e,ne,nex,ney,nexney,nx,ny,nz,n,nxny,n_eSurf,n_nSurf,nx_);
end

for iorder=1:1:7
    igbf=0;
    for iorientation=1:1:8
        for itype=1:1:8
            igbf=igbf+1;
            for igpz=1:1:3
                for igpy=1:1:3
                    for igpx=1:1:3
                        solution(iorder,igpx,igpy,igpz)=...
                            solution(iorder,igpx,igpy,igpz)...
                            +weights(e,iorientation,itype,iorder,igpx,igpy,igpz)...
                            *c(gbfs(igbf));
                    end
                end
            end
        end
    end
end
end